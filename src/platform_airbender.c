/*
 * The I/O half for the Airbender zkVM. Everything goes through the non-determinism CSR 0x7c0:
 *   - a read pops the next input word; the input is one word of byte length, then the bytes
 *     packed little-endian into words;
 *   - a write feeds the simulator's quasi-UART, which the platform's zkvm_platform_write speaks.
 * The public output is the eight words left in x10-x17 at the end of the run: write_output puts
 * them into airbender_output, which the platform's zkvm_platform_exit(0) loads in the exact
 * sequence the Airbender tooling looks for in the binary. The platform is the entry object of
 * the airbender platform of bflat (platform/entry.S); every guest links it.
 */
#include "u256.h"
#include "zkvm.h"

void *malloc(size_t size);

/* The platform (platform/entry.S): the output words, x10-x25, and the ends of a run. */
extern uint32_t airbender_output[16];
__attribute__((noreturn)) void zkvm_platform_exit(long code);
__attribute__((noreturn)) void airbender_finish_error(uint32_t code);
void zkvm_platform_write(int fd, const char *buf, int len);

#ifdef AIRBENDER
static inline uint32_t csr_read(void)
{
    uint32_t value;
    __asm__ volatile("csrrw %0, 0x7c0, x0" : "=r"(value));
    return value;
}
#else
uint32_t csr_read(void); /* supplied by the host input test */
#endif

/* Airbender v3's Keccak-f[1600] delegation (common_constants delegation_types/keccak_f1600.rs): per
 * round five column-parity (CSR 0x7cd), five theta/rho (0x7cc) and five chi (0x7ce) calls, then one
 * column-parity call for the delayed final iota. x10 carries the call's control word, 0 at the
 * start, and x11 the state: the 25 lanes, the five column parities and scratch, 256-byte aligned
 * in RAM. The transpiler recognizes the run by its first CSR and its length, so it is one asm
 * block. */
#define KECCAK_COLUMN_PARITY "csrrw x0, 0x7cd, x0\n"
#define KECCAK_THETA_RHO "csrrw x0, 0x7cc, x0\n"
#define KECCAK_CHI "csrrw x0, 0x7ce, x0\n"
#define KECCAK_FIVE(call) call call call call call
#define KECCAK_ROUND KECCAK_FIVE(KECCAK_COLUMN_PARITY) KECCAK_FIVE(KECCAK_THETA_RHO) KECCAK_FIVE(KECCAK_CHI)
#define KECCAK_FOUR_ROUNDS KECCAK_ROUND KECCAK_ROUND KECCAK_ROUND KECCAK_ROUND

static uint64_t KeccakState[32] __attribute__((aligned(256)));

/* 25 lanes as 50 word moves: the state is at least word aligned, and a call to memcpy costs
 * half as much again in its loop. */
static inline __attribute__((always_inline)) void copy_lanes(uint64_t *to, const uint64_t *from)
{
    uint32_t *t = (uint32_t *)to;
    const uint32_t *f = (const uint32_t *)from;
#pragma clang loop unroll(full)
    for (int i = 0; i < 50; i++)
        t[i] = f[i];
}

/* Keccak-f[1600] in place on the delegation's state at `state`. */
static inline __attribute__((always_inline)) void permute_at(uint64_t *state)
{
#ifdef AIRBENDER
    register uint32_t x11 __asm__("a1") = (uint32_t)state;
    __asm__ volatile("add x10, x0, x0\n"
                     KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS
                     KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS
                     KECCAK_COLUMN_PARITY
                     :
                     : "r"(x11)
                     : "a0", "memory");
#else
    keccak_f1600_soft(state);
#endif
}

static inline __attribute__((always_inline)) void permute(void)
{
    permute_at(KeccakState);
}

void keccak_f1600(uint64_t state[25])
{
    copy_lanes(KeccakState, state);
    permute();
    copy_lanes(state, KeccakState);
}

/* Keccak-f[1600] on a caller's state in place, with no copies: 256-byte aligned in RAM, its lanes 25-31
 * left to the delegation's column parities and scratch. */
void airbender_keccak_f1600_in_place(uint64_t state[32])
{
    permute_at(state);
}

enum
{
    RateWords = 34 /* 136 bytes */
};

/* XORs the `words` words at data into the rate. Misaligned data is read as the aligned words that
 * cover it, each output word merged from two of them; the last of those words may run up to three bytes
 * past the data, within its own aligned word. */
static inline __attribute__((always_inline)) void absorb_words(const uint8_t *data, int words, int overwrite)
{
    uint32_t *s = (uint32_t *)KeccakState;
    uint32_t misalignment = (uintptr_t)data & 3;
    /* A loop each for writing and XORing: the compiler keeps the test of `overwrite` in the loop. */
    if (misalignment == 0)
    {
        const uint32_t *w = (const uint32_t *)data;
        if (overwrite)
        {
#pragma clang loop unroll_count(8)
            for (int i = 0; i < words; i++)
                s[i] = w[i];
        }
        else
        {
#pragma clang loop unroll_count(8)
            for (int i = 0; i < words; i++)
                s[i] ^= w[i];
        }
        return;
    }
    const uint32_t *w = (const uint32_t *)(data - misalignment);
    uint32_t low = 8 * misalignment, high = 32 - low;
    uint32_t current = w[0];
    if (overwrite)
    {
#pragma clang loop unroll_count(8)
        for (int i = 0; i < words; i++)
        {
            uint32_t next = w[i + 1];
            s[i] = (current >> low) | (next << high);
            current = next;
        }
    }
    else
    {
#pragma clang loop unroll_count(8)
        for (int i = 0; i < words; i++)
        {
            uint32_t next = w[i + 1];
            s[i] ^= (current >> low) | (next << high);
            current = next;
        }
    }
}
/* absorb_words for a whole rate, with its loops unrolled in full: no counter and no branch for each group. */
static inline __attribute__((always_inline)) void absorb_rate(const uint8_t *data, int overwrite)
{
    uint32_t *s = (uint32_t *)KeccakState;
    uint32_t misalignment = (uintptr_t)data & 3;
    /* A loop each for writing and XORing: the compiler keeps the test of `overwrite` in the loop. */
    if (misalignment == 0)
    {
        const uint32_t *w = (const uint32_t *)data;
        if (overwrite)
        {
#pragma clang loop unroll(full)
            for (int i = 0; i < RateWords; i++)
                s[i] = w[i];
        }
        else
        {
#pragma clang loop unroll(full)
            for (int i = 0; i < RateWords; i++)
                s[i] ^= w[i];
        }
        return;
    }
    const uint32_t *w = (const uint32_t *)(data - misalignment);
    uint32_t low = 8 * misalignment, high = 32 - low;
    uint32_t current = w[0];
    if (overwrite)
    {
#pragma clang loop unroll(full)
        for (int i = 0; i < RateWords; i++)
        {
            uint32_t next = w[i + 1];
            s[i] = (current >> low) | (next << high);
            current = next;
        }
    }
    else
    {
#pragma clang loop unroll(full)
        for (int i = 0; i < RateWords; i++)
        {
            uint32_t next = w[i + 1];
            s[i] ^= (current >> low) | (next << high);
            current = next;
        }
    }
}

/* Keccak-256 of the len bytes at data into out, absorbed straight into the delegation's state. The first
 * block is written over the rate rather than XORed into zeros, so only the capacity, and any rate word a
 * short input leaves, starts at zero. */
__attribute__((noinline)) static void keccak256_long(const uint8_t *data, size_t len, uint8_t out[32])
{
    uint32_t *s = (uint32_t *)KeccakState;
#pragma clang loop unroll(full)
    for (int i = RateWords; i < 50; i++)
        s[i] = 0;
    int first = 1;
    while (len >= 4 * RateWords)
    {
        absorb_rate(data, first);
        first = 0;
        permute();
        data += 4 * RateWords;
        len -= 4 * RateWords;
    }

    /* The tail's whole words, then its last bytes with the 0x01 pad byte above them, and the final 0x80. */
    int words = (int)(len >> 2);
    absorb_words(data, words, first);
    if (first)
    {
#pragma clang loop unroll_count(8)
        for (int i = words; i < RateWords; i++)
            s[i] = 0;
    }
    uint32_t rest = len & 3, last = 0;
    for (uint32_t b = 0; b < rest; b++)
        last |= (uint32_t)data[4 * words + b] << (8 * b);
    s[words] ^= last ^ (1u << (8 * rest));
    s[RateWords - 1] ^= 0x80000000u;
    permute();

    if (((uintptr_t)out & 3) == 0)
    {
        for (int i = 0; i < 8; i++)
            ((uint32_t *)out)[i] = s[i];
    }
    else
    {
        memcpy(out, s, 32);
    }
}

/* The words of the rate from input laid out as EVM memory is on Airbender: each aligned group of four bytes
 * reversed, so a group loads as its bytes big-endian and a byte swap gives the lane's little-endian word. */
static inline __attribute__((always_inline)) void absorb_laid_out_words(const uint32_t *w, int words, int overwrite)
{
    uint32_t *s = (uint32_t *)KeccakState;
    if (overwrite)
    {
#pragma clang loop unroll_count(4)
        for (int i = 0; i < words; i++)
            s[i] = rev8(w[i]);
    }
    else
    {
#pragma clang loop unroll_count(4)
        for (int i = 0; i < words; i++)
            s[i] ^= rev8(w[i]);
    }
}

/* Keccak-256 of the len bytes at data, which start on an aligned group and are laid out as EVM memory: see
 * absorb_laid_out_words. Byte k of the input lies at data[k ^ 3]. */
void airbender_keccak256_laid_out(const uint8_t *data, size_t len, uint8_t out[32])
{
    uint32_t *s = (uint32_t *)KeccakState;
    const uint32_t *w = (const uint32_t *)data;
#pragma clang loop unroll(full)
    for (int i = RateWords; i < 50; i++)
        s[i] = 0;
    int first = 1;
    while (len >= 4 * RateWords)
    {
        absorb_laid_out_words(w, RateWords, first);
        first = 0;
        permute();
        w += RateWords;
        len -= 4 * RateWords;
    }

    int words = (int)(len >> 2);
    absorb_laid_out_words(w, words, first);
    if (first)
    {
#pragma clang loop unroll_count(8)
        for (int i = words; i < RateWords; i++)
            s[i] = 0;
    }
    uint32_t rest = len & 3, last = 0;
    const uint8_t *tail = (const uint8_t *)(w + words);
    for (uint32_t b = 0; b < rest; b++)
        last |= (uint32_t)tail[b ^ 3] << (8 * b);
    s[words] ^= last ^ (1u << (8 * rest));
    s[RateWords - 1] ^= 0x80000000u;
    permute();

    if (((uintptr_t)out & 3) == 0)
    {
        for (int i = 0; i < 8; i++)
            ((uint32_t *)out)[i] = s[i];
    }
    else
    {
        memcpy(out, s, 32);
    }
}

/* One block from aligned input and output, most calls, by a path small enough to need no frame; the rest by
 * keccak256_long. The last word of the input may be read past its end, within its own aligned word. */
void airbender_keccak256(const uint8_t *data, size_t len, uint8_t out[32])
{
    if (len >= 4 * RateWords || (((uintptr_t)data | (uintptr_t)out) & 3) != 0)
    {
        keccak256_long(data, len, out);
        return;
    }
    uint32_t *s = (uint32_t *)KeccakState;
    const uint32_t *w = (const uint32_t *)data;
    uint32_t words = (uint32_t)(len >> 2), rest = (uint32_t)(len & 3);
    uint32_t i = 0;
#pragma clang loop unroll_count(8)
    for (; i < words; i++)
        s[i] = w[i];
    s[i] = (rest ? w[i] & ((1u << (8 * rest)) - 1) : 0) | (1u << (8 * rest));
#pragma clang loop unroll_count(8)
    for (i++; i < 50; i++)
        s[i] = 0;
    s[RateWords - 1] ^= 0x80000000u;
    permute();
    uint32_t *o = (uint32_t *)out;
    for (i = 0; i < 8; i++)
        o[i] = s[i];
}

#ifdef AIRBENDER
/* memcpy, by --wrap (the manifest's wrap_symbols): the guest copies megabytes of witness nodes from
 * misaligned input positions. 32 bytes a round; a misaligned source is read as the aligned words that
 * cover it, each output word merged from two of them, and the last may run up to three bytes past the
 * source within its own aligned word. no_builtin keeps the byte loops from becoming memcpy calls. */
/* 2^(256 - 8m) for m = 1..31, in RAM for the delegation: a product by it shifts a word right by m bytes (MUL_HIGH)
   or left by 32 - m bytes (MUL_LOW). */
static u256d FunnelFactor[32];
/* The delegation reads only RAM, which starts at 64 MiB (the linker script's ram region); a source in ROM takes the
   word loop. */
static const uintptr_t RamStart = 64u << 20;
static u256d FunnelScratch;
static int FunnelReady;

/* Copies the 32-byte blocks of d (32-byte aligned) from s (any alignment), each merged by the delegation from the
   two aligned source words that cover it: two copies, two products by the funnel factor and an addition, where the
   word loop takes forty instructions a block. Both regions are in RAM, and the source words may extend up to 31
   bytes past either end of the source, which only reads them. Returns how many bytes it copied. */
static size_t funnel_copy(uint8_t *d, const uint8_t *s, size_t n)
{
    uint32_t m = (uintptr_t)s & 31;
    const u256d *a = (const u256d *)(s - m);
    u256d *out = (u256d *)d;
    size_t blocks = n >> 5;
    if (m == 0)
    {
        /* Four blocks a round, then one at a time: each block is the delegation and the two pointer steps, where the
           compiler would rebuild both pointers and keep a count. The delegation keeps a0 and a1. */
        register uint32_t x10 __asm__("a0") = (uint32_t)out;
        register uint32_t x11 __asm__("a1") = (uint32_t)a;
        uint32_t rounds_end = (uint32_t)(out + (blocks & ~(size_t)3)), end = (uint32_t)(out + blocks);
#define MEMCOPY_BLOCK "li a2, 0x80\n csrrw x0, 0x7ca, x0\n addi a0, a0, 32\n addi a1, a1, 32\n"
        __asm__ volatile("beq a0, %[rounds_end], 2f\n"
                         "1:\n" MEMCOPY_BLOCK MEMCOPY_BLOCK MEMCOPY_BLOCK MEMCOPY_BLOCK
                         "bne a0, %[rounds_end], 1b\n"
                         "2:\n"
                         "beq a0, %[end], 4f\n"
                         "3:\n" MEMCOPY_BLOCK
                         "bne a0, %[end], 3b\n"
                         "4:\n"
                         : "+r"(x10), "+r"(x11)
                         : [rounds_end] "r"(rounds_end), [end] "r"(end)
                         : "a2", "memory");
#undef MEMCOPY_BLOCK
        return blocks << 5;
    }
    if (!FunnelReady)
    {
        for (uint32_t k = 1; k < 32; k++)
        {
            uint32_t bit = 256 - 8 * k;
            FunnelFactor[k].w[bit >> 5] = 1u << (bit & 31);
        }
        FunnelReady = 1;
    }
    const u256d *factor = &FunnelFactor[m];
    for (size_t i = 0; i < blocks; i++)
    {
        u256_op(&out[i], &a[i], U256_MEMCOPY);
        u256_op(&out[i], factor, U256_MUL_HIGH);
        u256_op(&FunnelScratch, &a[i + 1], U256_MEMCOPY);
        u256_op(&FunnelScratch, factor, U256_MUL_LOW);
        u256_op(&out[i], &FunnelScratch, U256_ADD);
    }
    return blocks << 5;
}

__attribute__((no_builtin("memcpy"), always_inline)) static inline void copy_bytes(uint8_t *d, const uint8_t *s, size_t n)
{
    /* Stepped by the destination alone, with the source at a fixed distance: one increment a byte. */
    ptrdiff_t distance = s - d;
    uint8_t *end = d + n;
    for (; d != end; d++)
        *d = d[distance];
}

/* The word loop for a source misaligned to the word-aligned d: each output word is merged from the two aligned
   source words that cover it, and the last may run up to three bytes past the source within its own aligned word.
   Out of line: its eight live words take saved registers, which the other paths then need not save. */
__attribute__((no_builtin("memcpy"), noinline)) static void copy_merged(uint32_t *dw, const uint8_t *s, size_t n)
{
    uint32_t misalignment = (uintptr_t)s & 3;
    const uint32_t *sw = (const uint32_t *)(s - misalignment);
    uint32_t low = 8 * misalignment, high = 32 - low;
    uint32_t current = *sw++;
    for (; n >= 32; n -= 32, dw += 8, sw += 8)
    {
        uint32_t w0 = sw[0], w1 = sw[1], w2 = sw[2], w3 = sw[3], w4 = sw[4], w5 = sw[5], w6 = sw[6], w7 = sw[7];
        dw[0] = (current >> low) | (w0 << high);
        dw[1] = (w0 >> low) | (w1 << high);
        dw[2] = (w1 >> low) | (w2 << high);
        dw[3] = (w2 >> low) | (w3 << high);
        dw[4] = (w3 >> low) | (w4 << high);
        dw[5] = (w4 >> low) | (w5 << high);
        dw[6] = (w5 >> low) | (w6 << high);
        dw[7] = (w6 >> low) | (w7 << high);
        current = w7;
    }
    for (; n >= 4; n -= 4)
    {
        uint32_t next = *sw++;
        *dw++ = (current >> low) | (next << high);
        current = next;
    }
    copy_bytes((uint8_t *)dw, (const uint8_t *)sw - 4 + misalignment, n);
}

/* The word loops, for a misaligned source by copy_merged. */
__attribute__((no_builtin("memcpy"), always_inline)) static inline void copy_words(uint8_t *d, const uint8_t *s, size_t n)
{
    if (n < 16)
    {
        copy_bytes(d, s, n);
        return;
    }
    while ((uintptr_t)d & 3)
    {
        *d++ = *s++;
        n--;
    }
    uint32_t *dw = (uint32_t *)d;
    if ((uintptr_t)s & 3)
    {
        copy_merged(dw, s, n);
        return;
    }
    {
        const uint32_t *sw = (const uint32_t *)s;
        for (; n >= 32; n -= 32, dw += 8, sw += 8)
        {
            uint32_t w0 = sw[0], w1 = sw[1], w2 = sw[2], w3 = sw[3], w4 = sw[4], w5 = sw[5], w6 = sw[6], w7 = sw[7];
            dw[0] = w0;
            dw[1] = w1;
            dw[2] = w2;
            dw[3] = w3;
            dw[4] = w4;
            dw[5] = w5;
            dw[6] = w6;
            dw[7] = w7;
        }
        for (; n >= 4; n -= 4)
            *dw++ = *sw++;
        copy_bytes((uint8_t *)dw, (const uint8_t *)sw, n);
    }
}

__attribute__((no_builtin("memcpy"), noinline)) static void *copy_long(void *restrict dest, const void *restrict src, size_t n)
{
    uint8_t *d = dest;
    const uint8_t *s = src;
    /* Long enough to pay for aligning the destination to a block: the delegation copies a block in one operation when
       source and destination share their alignment within it, and merges it in five when they differ within a word;
       sharing a word's alignment only, the word loop is as quick. A copy onto itself, which callers make, would hand
       the delegation one operand twice. */
    uintptr_t apart = (uintptr_t)d ^ (uintptr_t)s;
    if (n >= 128 && (uintptr_t)s >= RamStart + 32 && d != s && ((apart & 31) == 0 || (apart & 3) != 0))
    {
        size_t head = -(uintptr_t)d & 31;
        copy_words(d, s, head);
        d += head;
        s += head;
        n -= head;
        size_t done = funnel_copy(d, s, n);
        d += done;
        s += done;
        n -= done;
    }
    copy_words(d, s, n);
    /* Hidden from the compiler, which would otherwise know the result is dest and keep dest in __wrap_memcpy across
       the call instead of jumping to it. */
    __asm__("" : "+r"(dest));
    return dest;
}

/* Most copies are a few bytes: they take the byte loop before copy_long's frame. */
__attribute__((no_builtin("memcpy"))) void *__wrap_memcpy(void *restrict dest, const void *restrict src, size_t n)
{
    if (n >= 16)
        return copy_long(dest, src, n);
    copy_bytes(dest, src, n);
    return dest;
}

void *__real_memmove(void *dest, const void *src, size_t n);

/* memmove, by --wrap: managed code copies through it, mostly a few bytes between disjoint buffers. A short copy runs
   forwards unless the destination starts inside the source; a long one goes to copy_long when the two are disjoint,
   else to the libc memmove. */
__attribute__((no_builtin("memcpy", "memmove"))) void *__wrap_memmove(void *dest, const void *src, size_t n)
{
    if (n < 16)
    {
        if ((uintptr_t)dest - (uintptr_t)src >= n)
        {
            copy_bytes(dest, src, n);
            return dest;
        }
        uint8_t *d = dest;
        const uint8_t *s = src;
        while (n--)
            d[n] = s[n];
        return dest;
    }
    if ((uintptr_t)src - (uintptr_t)dest - n <= -2 * n)
        return copy_long(dest, src, n);
    return __real_memmove(dest, src, n);
}

#endif

/* The u256 delegation for managed code. A JIT that lists it in JitZkDelegationPInvokes emits the CSR
 * access in place of a call to this. */
uint32_t airbender_u256_op(u256d *a, const u256d *b, uint32_t op)
{
    return u256_op(a, b, op);
}


/* ---- The trusted input, hashed as it is read ------------------------------------------------------
 * Nethermind's trusted input ends with a list of state nodes and a list of codes, and the guest keys each of
 * them by its keccak. Read first and hashed later, every word costs a store, then a load of the word and a load,
 * an XOR and a store of the state. Absorbed while it is read, the word is in a register already: the later load
 * goes. The digests, and for a full branch the three sponge states the guest's commit resumes from, are kept
 * for the guest to take with airbender_prehashed_nodes and airbender_prehashed_code. An input of another shape
 * is read as it always was, and the guest then hashes it itself. */
enum
{
    TrustedMagic = 0x4E49544E,
    TrustedHeaderWords = 22,
    TrustedNodesField = 16,
    TrustedCodesField = 17,
    TrustedTotalField = 18,
    TrustedItemPrefixWords = 6,
    TrustedCodePadding = 34,
    FullBranchLength = 532,
    RetainedStateWords = 64, /* a 256-byte slot: 25 lanes and the delegation's scratch */
    RetainedStates = 3
};
static struct
{
    uint8_t *node_hashes;
    uint32_t *node_states; /* in 32-bit words, as the sponge is run */
    const uint8_t *first_node, *last_node;
    const uint8_t **node_data;
    uint32_t node_count;
    uint8_t *code_hashes;
    const uint8_t **code_data;
    uint32_t code_count;
} Prehashed;

/* The capacity of `from` into `to`: lanes 16 to 27 in three delegation copies, of which the caller then writes
 * lane 16 and the rest past lane 24 is scratch. Both states are 256-byte aligned and distinct. */
static inline __attribute__((always_inline)) void copy_capacity(uint32_t *to, const uint32_t *from)
{
    u256_op((u256d *)(to + 32), (const u256d *)(from + 32), U256_MEMCOPY);
    u256_op((u256d *)(to + 40), (const u256d *)(from + 40), U256_MEMCOPY);
    u256_op((u256d *)(to + 48), (const u256d *)(from + 48), U256_MEMCOPY);
}
/* A full branch of 532 bytes, read into w: the states after each of its three whole rate blocks in the three
 * slots at st, as Nethermind's ComputeHash532Retaining leaves them, and its digest in out. */
static void read_full_branch(uint32_t *w, uint32_t *st, uint32_t *out)
{
    uint32_t *s = (uint32_t *)KeccakState;
#pragma clang loop unroll(full)
    for (int j = 0; j < RateWords; j++)
    {
        uint32_t x = csr_read();
        w[j] = x;
        st[j] = x;
    }
#pragma clang loop unroll(full)
    for (int j = RateWords; j < 50; j++)
        st[j] = 0;
    permute_at((uint64_t *)st);
    for (int block = 1; block < RetainedStates; block++)
    {
        uint32_t *next = st + RetainedStateWords;
        w += RateWords;
        copy_capacity(next, st);
#pragma clang loop unroll(full)
        for (int j = 0; j < RateWords; j++)
        {
            uint32_t x = csr_read();
            w[j] = x;
            next[j] = st[j] ^ x;
        }
        permute_at((uint64_t *)next);
        st = next;
    }
    /* The tail: 124 bytes, the 0x01 pad byte in the word after them, the 0x80 in the rate's last byte. */
    w += RateWords;
    copy_capacity(s, st);
#pragma clang loop unroll(full)
    for (int j = 0; j < 31; j++)
    {
        uint32_t x = csr_read();
        w[j] = x;
        s[j] = st[j] ^ x;
    }
    s[31] = st[31] ^ 1u;
    s[33] ^= 0x80000000u;
    permute();
#pragma clang loop unroll(full)
    for (int j = 0; j < 8; j++)
        out[j] = s[j];
}
/* Any other item of len bytes, read into w, its digest in out. Reads (len + 3) / 4 words. */
static void read_and_hash(uint32_t *w, uint32_t len, uint32_t *out)
{
    uint32_t *s = (uint32_t *)KeccakState;
#pragma clang loop unroll(full)
    for (int j = RateWords; j < 50; j++)
        s[j] = 0;
    int first = 1;
    while (len >= 4 * RateWords)
    {
        if (first)
        {
#pragma clang loop unroll(full)
            for (int j = 0; j < RateWords; j++)
            {
                uint32_t x = csr_read();
                w[j] = x;
                s[j] = x;
            }
        }
        else
        {
#pragma clang loop unroll(full)
            for (int j = 0; j < RateWords; j++)
            {
                uint32_t x = csr_read();
                w[j] = x;
                s[j] ^= x;
            }
        }
        first = 0;
        permute();
        w += RateWords;
        len -= 4 * RateWords;
    }
    uint32_t words = len >> 2, rest = len & 3, last = 0, j = 0;
    if (first)
    {
#pragma clang loop unroll_count(8)
        for (; j < words; j++)
        {
            uint32_t x = csr_read();
            w[j] = x;
            s[j] = x;
        }
    }
    else
    {
#pragma clang loop unroll_count(8)
        for (; j < words; j++)
        {
            uint32_t x = csr_read();
            w[j] = x;
            s[j] ^= x;
        }
    }
    if (rest)
    {
        uint32_t x = csr_read();
        w[words] = x;
        last = x & ((1u << (8 * rest)) - 1);
    }
    last |= 1u << (8 * rest);
    if (first)
    {
        s[words] = last;
#pragma clang loop unroll_count(8)
        for (j = words + 1; j < RateWords; j++)
            s[j] = 0;
    }
    else
    {
        s[words] ^= last;
    }
    s[RateWords - 1] ^= 0x80000000u;
    permute();
#pragma clang loop unroll(full)
    for (int k = 0; k < 8; k++)
        out[k] = s[k];
}
/* Words [i, end) of the input, with nothing else done to them. */
static uint32_t read_plain(uint32_t *words, uint32_t i, uint32_t end)
{
    for (; i + 32 <= end; i += 32)
    {
#pragma clang loop unroll(full)
        for (int k = 0; k < 32; k++)
            words[i + k] = csr_read();
    }
    for (; i < end; i++)
        words[i] = csr_read();
    return i;
}
/* Returns the first unread word. A complete list, including an empty one, sets *complete. */
static uint32_t read_list(uint32_t *words, uint32_t i, uint32_t end, uint32_t padding, int nodes, int *complete)
{
    *complete = 0;
    if (end - i < 2)
        return i;
    uint32_t n = csr_read();
    words[i] = n;
    words[i + 1] = csr_read();
    i += 2;
    if (n == 0)
    {
        *complete = i == end;
        return i;
    }
    if (n > (end - i) / TrustedItemPrefixWords)
        return i;
    size_t state_bytes = nodes ? RetainedStates * RetainedStateWords * sizeof(uint32_t) : 0;
    size_t per_node = state_bytes + 32 + sizeof(const uint8_t *);
    size_t slack = nodes ? 255 : 0;
    if (n > (SIZE_MAX - slack) / per_node)
        return i;
    uint8_t *area = malloc(n * per_node + slack);
    if (!area)
        return i;
    if (nodes)
        area = (uint8_t *)(((uintptr_t)area + 255) & ~(uintptr_t)255);
    uint32_t *states = (uint32_t *)area;
    uint32_t *hashes = (uint32_t *)(area + n * state_bytes);
    const uint8_t **data = (const uint8_t **)(hashes + (size_t)8 * n);
    const uint8_t *first = 0, *last = 0;
    for (uint32_t k = 0; k < n; k++)
    {
        if (end - i < TrustedItemPrefixWords)
            return i;
#pragma clang loop unroll(full)
        for (int p = 0; p < TrustedItemPrefixWords; p++)
            words[i + p] = csr_read();
        uint32_t len = words[i];
        i += TrustedItemPrefixWords;
        if (len > 4 * (end - i) || len > UINT32_MAX - padding - 7)
            return i;
        uint32_t stride = ((len + padding + 7) & ~7u) >> 2, body = (len + 3) >> 2;
        if (stride > end - i)
            return i;
        uint32_t *item = words + i;
        if (nodes && len == FullBranchLength)
            read_full_branch(item, states + (size_t)k * (RetainedStates * RetainedStateWords), hashes + 8 * k);
        else
            read_and_hash(item, len, hashes + 8 * k);
        for (uint32_t q = body; q < stride; q++)
            item[q] = csr_read();
        i += stride;
        if (k == 0)
            first = (const uint8_t *)item;
        last = (const uint8_t *)item;
        data[k] = (const uint8_t *)item;
    }
    if (i != end)
        return i;
    if (nodes)
    {
        Prehashed.node_hashes = (uint8_t *)hashes;
        Prehashed.node_states = states;
        Prehashed.first_node = first;
        Prehashed.last_node = last;
        Prehashed.node_data = data;
        Prehashed.node_count = n;
    }
    else
    {
        Prehashed.code_hashes = (uint8_t *)hashes;
        Prehashed.code_data = data;
        Prehashed.code_count = n;
    }
    *complete = 1;
    return i;
}
/* The digests and states of the count state nodes whose first and last lie at first and last, if those are the
 * nodes the prehashed input reader hashed. The states are 256-byte aligned, three slots a node, written for full branches only. */
int airbender_prehashed_nodes(const uint8_t *first, const uint8_t *last, uint32_t count, uint8_t **hashes, uint64_t **states)
{
    if (count == 0 || count != Prehashed.node_count || first != Prehashed.first_node || last != Prehashed.last_node)
        return 0;
    *hashes = Prehashed.node_hashes;
    *states = (uint64_t *)Prehashed.node_states;
    return 1;
}
/* ---- The tries over the nodes the input reader hashed: see mpt.c ------------------------------------------ */
static uint32_t mpt_w_count(void)
{
    return Prehashed.node_count;
}
static const uint32_t *mpt_w_hash(uint32_t k)
{
    return (const uint32_t *)(Prehashed.node_hashes + 32 * (size_t)k);
}
static const uint8_t *mpt_w_rlp(uint32_t k, uint32_t *len)
{
    const uint8_t *item = Prehashed.node_data[k];
    *len = ((const uint32_t *)item)[-TrustedItemPrefixWords];
    return item;
}
static void mpt_keccak(const uint8_t *data, uint32_t len, uint32_t out[8])
{
    keccak256_long(data, len, (uint8_t *)out);
}
/* The digest of full branch k as the 532 bytes at edited, which equal those the input reader hashed up to rate block
 * `block` and are valid from it on: the sponge resumes from the state kept after the block before. */
static void mpt_hash_edited(uint32_t k, int block, const uint8_t *original, const uint8_t *edited, uint32_t out[8])
{
    (void)original;
    uint32_t *s = (uint32_t *)KeccakState;
    const uint32_t *w = (const uint32_t *)edited + RateWords * block;
    if (block == 0)
    {
#pragma clang loop unroll(full)
        for (int j = 0; j < RateWords; j++)
            s[j] = w[j];
#pragma clang loop unroll(full)
        for (int j = RateWords; j < 50; j++)
            s[j] = 0;
    }
    else
    {
        const uint32_t *st = Prehashed.node_states
            + (size_t)k * (RetainedStates * RetainedStateWords) + (size_t)(block - 1) * RetainedStateWords;
        copy_capacity(s, st);
        if (block == RetainedStates)
        {
#pragma clang loop unroll(full)
            for (int j = 0; j < 31; j++)
                s[j] = st[j] ^ w[j];
            s[31] = st[31] ^ 1u;
            s[33] ^= 0x80000000u;
            permute();
            goto done;
        }
#pragma clang loop unroll(full)
        for (int j = 0; j < RateWords; j++)
            s[j] = st[j] ^ w[j];
    }
    permute();
    w += RateWords;
    for (block++; block < RetainedStates; block++)
    {
#pragma clang loop unroll(full)
        for (int j = 0; j < RateWords; j++)
            s[j] ^= w[j];
        permute();
        w += RateWords;
    }
#pragma clang loop unroll(full)
    for (int j = 0; j < 31; j++)
        s[j] ^= w[j];
    s[31] ^= 1u;
    s[33] ^= 0x80000000u;
    permute();
done:
#pragma clang loop unroll(full)
    for (int j = 0; j < 8; j++)
        out[j] = s[j];
}
/* len bytes from `from` to `to`, equal modulo 32, as the whole 32-byte blocks that hold them: a delegation each. */
static inline __attribute__((always_inline)) void mpt_copy(uint8_t *to, const uint8_t *from, uint32_t len)
{
    uintptr_t a = (uintptr_t)from & ~(uintptr_t)31, end = ((uintptr_t)from + len + 31) & ~(uintptr_t)31;
    uint8_t *d = to - ((uintptr_t)from - a);
    for (; a < end; a += 32, d += 32)
        u256_op((u256d *)d, (const u256d *)a, U256_MEMCOPY);
}
/* The keccak of the 20 or 32 bytes in the aligned words at w: one block, written whole. */
static void mpt_keccak_key(const uint32_t *w, uint32_t len, uint32_t out[8])
{
    uint32_t *s = (uint32_t *)KeccakState;
    if (len == 20)
    {
#pragma clang loop unroll(full)
        for (int j = 0; j < 5; j++)
            s[j] = w[j];
        s[5] = 1;
#pragma clang loop unroll(full)
        for (int j = 6; j < 50; j++)
            s[j] = 0;
    }
    else
    {
#pragma clang loop unroll(full)
        for (int j = 0; j < 8; j++)
            s[j] = w[j];
        s[8] = 1;
#pragma clang loop unroll(full)
        for (int j = 9; j < 50; j++)
            s[j] = 0;
    }
    s[33] = 0x80000000u;
    permute();
#pragma clang loop unroll(full)
    for (int j = 0; j < 8; j++)
        out[j] = s[j];
}
#define MPT_UNALIGNED_WORDS 1
#define MPT_REV8(x) rev8(x)
#include "mpt.c"
/* ---- The log bloom ------------------------------------------------------------------------------------
 * Sets in the 256-byte bloom, and in master if given, the three bits the keccak of the len bytes at data picks:
 * of each of its first three byte pairs, the low eleven bits, counted from the bloom's end. A block's logs repeat
 * a few addresses and topics, so the bits of the last inputs are kept. len is 20 or 32. */
enum
{
    BloomMemoSize = 4096
};
static struct
{
    uint32_t pre[8];
    uint32_t len;
    uint16_t bit[4];
} BloomMemo[BloomMemoSize];
void airbender_bloom_set(uint8_t *bloom, uint8_t *master, const uint8_t *data, uint32_t len)
{
    uint32_t pre[8];
    if (len == 20)
        mpt_load_address(pre, data);
    else
        mpt_load32(pre, data);
    uint32_t x = pre[0] ^ pre[1] ^ pre[4] ^ pre[7];
    x ^= x >> 13;
    x ^= x >> 7;
    __typeof__(BloomMemo[0]) *e = &BloomMemo[x & (BloomMemoSize - 1)];
    if (e->len != len
        || ((e->pre[0] ^ pre[0]) | (e->pre[1] ^ pre[1]) | (e->pre[2] ^ pre[2]) | (e->pre[3] ^ pre[3])
            | (e->pre[4] ^ pre[4]) | (e->pre[5] ^ pre[5]) | (e->pre[6] ^ pre[6]) | (e->pre[7] ^ pre[7])) != 0)
    {
        uint32_t h[8];
#pragma clang loop unroll(full)
        for (int i = 0; i < 8; i++)
            e->pre[i] = pre[i];
        e->len = len;
        mpt_keccak_key(e->pre, len, h);
        /* Byte pairs (0,1), (2,3), (4,5) of the digest, big-endian. */
        e->bit[0] = (uint16_t)(((h[0] & 7) << 8) | ((h[0] >> 8) & 0xff));
        e->bit[1] = (uint16_t)((((h[0] >> 16) & 7) << 8) | (h[0] >> 24));
        e->bit[2] = (uint16_t)(((h[1] & 7) << 8) | ((h[1] >> 8) & 0xff));
    }
#pragma clang loop unroll(full)
    for (int i = 0; i < 3; i++)
    {
        uint32_t bit = e->bit[i], at = 255 - (bit >> 3), mask = 1u << (bit & 7);
        bloom[at] |= (uint8_t)mask;
        if (master)
            master[at] |= (uint8_t)mask;
    }
}
/* The nodes the input reader hashed, which the tries are over: the data of the first and of the last, and their count. */
int airbender_mpt_nodes(const uint8_t **first, const uint8_t **last, uint32_t *count)
{
    *first = Prehashed.first_node;
    *last = Prehashed.last_node;
    *count = Prehashed.node_count;
    return Prehashed.node_count != 0;
}
/* Whether the tries can serve these nodes: the count nodes from first to last are the ones the input reader hashed. */
int airbender_mpt_attach(const uint8_t *first, const uint8_t *last, uint32_t count)
{
    return count != 0 && count == Prehashed.node_count && first == Prehashed.first_node
        && last == Prehashed.last_node;
}
/* The digest of the code at data, if it is the index-th code the input reader hashed. */
int airbender_prehashed_code(const uint8_t *data, uint32_t index, uint8_t *out)
{
    if (index >= Prehashed.code_count || Prehashed.code_data[index] != data)
        return 0;
    const uint32_t *h = (const uint32_t *)(Prehashed.code_hashes + 32 * (size_t)index);
    if (((uintptr_t)out & 3) == 0)
    {
#pragma clang loop unroll(full)
        for (int i = 0; i < 8; i++)
            ((uint32_t *)out)[i] = h[i];
    }
    else
    {
        memcpy(out, h, 32);
    }
    return 1;
}

static struct
{
    uint8_t *start;
    uint32_t count;
} DivisionHints;
/* The division hints the input reader found after the trusted input: count pairs of a quotient and a remainder. */
void airbender_division_hints(uint8_t **start, uint32_t *count)
{
    *start = DivisionHints.start;
    *count = DivisionHints.count;
}
/* Return aligned operands for the delegation while preserving the complete raw input. */
static uint32_t *hint_payload(uint32_t *words, uint32_t count)
{
    if (((uintptr_t)words & 31) == 0 || count == 0)
        return words;
    size_t bytes = (size_t)count * sizeof(*words);
    if (bytes > SIZE_MAX - 31)
        return 0;
    uint8_t *raw = malloc(bytes + 31);
    if (!raw)
        return 0;
    uint32_t *aligned = (uint32_t *)(((uintptr_t)raw + 31) & ~(uintptr_t)31);
    memcpy(aligned, words, bytes);
    return aligned;
}

static void read_hints(uint32_t *words, uint32_t count)
{
    if (count < 2 || words[0] > (count - 2) / 16)
        return;
    uint32_t n = words[0], size = 16 * n;
    DivisionHints.start = (uint8_t *)hint_payload(words + 2, size);
    if (DivisionHints.start)
        DivisionHints.count = n;
    words += 2 + size;
    count -= 2 + size;
    if (count < 2 || words[0] > (count - 2) / 8)
        return;
    uint32_t m = words[0];
    if (!m)
    {
        FieldHints.active = 1; /* An empty tape requests recording. */
        return;
    }
    const u256d *tape = (const u256d *)hint_payload(words + 2, 8 * m);
    if (tape)
        FieldHints = (field_hints){tape, tape + m, 1};
}

static inline __attribute__((always_inline)) void read_guest_input(const uint8_t **buf_ptr, size_t *buf_size, int prehash)
{
    uint32_t len = csr_read();
    uint32_t count = len / 4 + (len % 4 != 0);
    if (count > UINT32_MAX / 4)
        airbender_fail(0x202);
    uint32_t *words = malloc((size_t)count * 4);
    if (!words && count)
        airbender_fail(0x202);
    uint32_t i = 0;
    if (prehash && count >= TrustedHeaderWords)
    {
        i = read_plain(words, 0, TrustedHeaderWords);
        uint32_t nodes = words[TrustedNodesField], codes = words[TrustedCodesField], total = words[TrustedTotalField];
        if (words[0] == TrustedMagic && words[1] == 1 && total <= len && ((len | total) & 7) == 0
            && ((nodes | codes) & 7) == 0 && nodes >= 4 * TrustedHeaderWords && nodes < codes && codes < total)
        {
            int complete;
            i = read_plain(words, i, nodes >> 2);
            i = read_list(words, i, codes >> 2, 0, 1, &complete);
            if (complete)
                i = read_list(words, i, total >> 2, TrustedCodePadding, 0, &complete);
            if (complete && i < count)
            {
                read_plain(words, i, count);
                read_hints(words + i, count - i);
                i = count;
            }
        }
    }
    /* Thirty-two words a round: the input is megabytes, and the loop's own instructions would be a third of
     * each word's cost. Nothing past the input is read: the source has no more words to give. */
    read_plain(words, i, count);
    *buf_ptr = (const uint8_t *)words;
    *buf_size = len;
}

void read_input(const uint8_t **buf_ptr, size_t *buf_size)
{
    read_guest_input(buf_ptr, buf_size, 0);
}

void airbender_read_prehashed_input(const uint8_t **buf_ptr, size_t *buf_size)
{
    read_guest_input(buf_ptr, buf_size, 1);
}

void write_output(const uint8_t *output, size_t size)
{
    if (size > 32)
        airbender_fail(0x201);
    for (size_t i = 0; i < 8; i++)
    {
        uint32_t w = 0;
        for (size_t b = 0; b < 4; b++)
            if (4 * i + b < size)
                w |= (uint32_t)output[4 * i + b] << (8 * b);
        airbender_output[i] = w;
    }
}

void sys_write(uint32_t fd, const uint8_t *ptr, size_t len)
{
    zkvm_platform_write((int)fd, (const char *)ptr, (int)len);
}

/* The end of a run: success publishes airbender_output, anything else proves nothing. */
void airbender_exit(int status)
{
    zkvm_platform_exit(status);
}

void airbender_fail(uint32_t code)
{
    airbender_finish_error(code);
}

/* 64-bit unsigned division for RV32IM: the compiler-rt routine shifts and subtracts a bit at a time, some
 * 29 instructions for each bit of the quotient. These divide by 16-bit digits with divu (Hacker's Delight,
 * divlu and its 64-by-64 form). The manifest wraps __udivdi3 and __umoddi3 onto them. */
static inline __attribute__((always_inline)) uint32_t leading_zeros(uint32_t x)
{
    uint32_t n = 0;
    if (x <= 0x0000ffffu) { n += 16; x <<= 16; }
    if (x <= 0x00ffffffu) { n += 8; x <<= 8; }
    if (x <= 0x0fffffffu) { n += 4; x <<= 4; }
    if (x <= 0x3fffffffu) { n += 2; x <<= 2; }
    if (x <= 0x7fffffffu) { n += 1; }
    return n;
}
/* (u1:u0) / v for u1 < v; the remainder to *r. */
static inline __attribute__((always_inline)) uint32_t divide_64_by_32(uint32_t u1, uint32_t u0, uint32_t v, uint32_t *r)
{
    const uint32_t b = 65536;
    uint32_t s = leading_zeros(v);
    v <<= s;
    uint32_t vn1 = v >> 16, vn0 = v & 0xffff;
    uint32_t un32 = s ? (u1 << s) | (u0 >> (32 - s)) : u1;
    uint32_t un10 = u0 << s;
    uint32_t un1 = un10 >> 16, un0 = un10 & 0xffff;
    uint32_t q1 = un32 / vn1, rhat = un32 - q1 * vn1;
    while (q1 >= b || q1 * vn0 > b * rhat + un1)
    {
        q1--;
        rhat += vn1;
        if (rhat >= b)
            break;
    }
    uint32_t un21 = un32 * b + un1 - q1 * v;
    uint32_t q0 = un21 / vn1;
    rhat = un21 - q0 * vn1;
    while (q0 >= b || q0 * vn0 > b * rhat + un0)
    {
        q0--;
        rhat += vn1;
        if (rhat >= b)
            break;
    }
    *r = (un21 * b + un0 - q0 * v) >> s;
    return q1 * b + q0;
}
static inline __attribute__((always_inline)) uint64_t divide_64(uint64_t u, uint64_t v, uint64_t *rem)
{
    uint32_t u1 = (uint32_t)(u >> 32), u0 = (uint32_t)u, v1 = (uint32_t)(v >> 32), v0 = (uint32_t)v, r;
    if (v1 == 0)
    {
        if (u1 == 0)
        {
            uint32_t q = u0 / v0;
            *rem = u0 - q * v0;
            return q;
        }
        uint32_t q1 = 0;
        if (u1 >= v0)
        {
            q1 = u1 / v0;
            u1 -= q1 * v0;
        }
        uint32_t q0 = divide_64_by_32(u1, u0, v0, &r);
        *rem = r;
        return ((uint64_t)q1 << 32) | q0;
    }
    uint32_t n = leading_zeros(v1);
    uint32_t vn = (uint32_t)((v << n) >> 32);
    uint64_t half = u >> 1;
    uint32_t q1 = divide_64_by_32((uint32_t)(half >> 32), (uint32_t)half, vn, &r);
    uint64_t q = ((uint64_t)q1 << n) >> 31;
    if (q != 0)
        q--;
    uint64_t left = u - q * v;
    if (left >= v)
    {
        q++;
        left -= v;
    }
    *rem = left;
    return q;
}
uint64_t __wrap___udivdi3(uint64_t u, uint64_t v)
{
    uint64_t r;
    return divide_64(u, v, &r);
}
uint64_t __wrap___umoddi3(uint64_t u, uint64_t v)
{
    uint64_t r;
    divide_64(u, v, &r);
    return r;
}
