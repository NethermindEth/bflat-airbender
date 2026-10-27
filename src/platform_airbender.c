/*
 * The I/O half for the Airbender zkVM. Everything goes through the non-determinism CSR 0x7c0:
 *   - a read pops the next input word; the input is one word of byte length, then the bytes
 *     packed little-endian into words;
 *   - a write feeds the simulator's quasi-UART: 0xffffffff, the number of words that follow,
 *     the byte length, then the bytes packed little-endian into words. Debug output only; the
 *     prover ignores it.
 * The public output is the eight words left in x10-x17 at the end of the run, which
 * airbender_finish_success (platform_airbender.S) loads in the exact sequence the Airbender
 * tooling looks for in the binary.
 */
#include "u256.h"
#include "zkvm.h"

void *malloc(size_t size);
__attribute__((noreturn)) void airbender_finish_success(const uint32_t words[16]);
__attribute__((noreturn)) void airbender_finish_error(uint32_t code);

/* x10-x17 the output, x18-x25 zero, as the recursion convention expects. */
static uint32_t Output[16];

static inline uint32_t csr_read(void)
{
    uint32_t value;
    __asm__ volatile("csrrw %0, 0x7c0, x0" : "=r"(value));
    return value;
}

static inline void csr_write(uint32_t value)
{
    __asm__ volatile("csrrw x0, 0x7c0, %0" : : "r"(value));
}

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

/* Keccak-f[1600] on KeccakState in place. */
static inline __attribute__((always_inline)) void permute(void)
{
    register uint32_t x11 __asm__("a1") = (uint32_t)KeccakState;
    __asm__ volatile("add x10, x0, x0\n"
                     KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS
                     KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS KECCAK_FOUR_ROUNDS
                     KECCAK_COLUMN_PARITY
                     :
                     : "r"(x11)
                     : "a0", "memory");
}

void keccak_f1600(uint64_t state[25])
{
    copy_lanes(KeccakState, state);
    permute();
    copy_lanes(state, KeccakState);
}

enum
{
    RateWords = 34 /* 136 bytes */
};

/* XORs the `words` words at data into the rate. Misaligned data is read as the aligned words that
 * cover it, each output word merged from two of them; the last of those words may run up to three bytes
 * past the data, within its own aligned word. */
static inline __attribute__((always_inline)) void absorb_words(const uint8_t *data, int words)
{
    uint32_t *s = (uint32_t *)KeccakState;
    uint32_t misalignment = (uintptr_t)data & 3;
    if (misalignment == 0)
    {
        const uint32_t *w = (const uint32_t *)data;
        for (int i = 0; i < words; i++)
            s[i] ^= w[i];
        return;
    }
    const uint32_t *w = (const uint32_t *)(data - misalignment);
    uint32_t low = 8 * misalignment, high = 32 - low;
    uint32_t current = w[0];
    for (int i = 0; i < words; i++)
    {
        uint32_t next = w[i + 1];
        s[i] ^= (current >> low) | (next << high);
        current = next;
    }
}

/* Keccak-256 of the len bytes at data into out, absorbed straight into the delegation's state. */
void airbender_keccak256(const uint8_t *data, size_t len, uint8_t out[32])
{
    uint32_t *s = (uint32_t *)KeccakState;
    for (int i = 0; i < 50; i++)
        s[i] = 0;
    while (len >= 4 * RateWords)
    {
        absorb_words(data, RateWords);
        permute();
        data += 4 * RateWords;
        len -= 4 * RateWords;
    }

    /* The tail's whole words, then its last bytes with the 0x01 pad byte above them, and the final 0x80. */
    int words = (int)(len >> 2);
    absorb_words(data, words);
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

/* The u256 delegation for managed code. A JIT that lists it in JitZkDelegationPInvokes emits the CSR
 * access in place of a call to this. */
uint32_t airbender_u256_op(u256d *a, const u256d *b, uint32_t op)
{
    return u256_op(a, b, op);
}

void read_input(const uint8_t **buf_ptr, size_t *buf_size)
{
    uint32_t len = csr_read();
    uint32_t *words = malloc(((size_t)len + 3) & ~(size_t)3);
    for (uint32_t i = 0; i < (len + 3) / 4; i++)
        words[i] = csr_read();
    *buf_ptr = (const uint8_t *)words;
    *buf_size = len;
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
        Output[i] = w;
    }
}

void sys_write(uint32_t fd, const uint8_t *ptr, size_t len)
{
    (void)fd;
    csr_write(0xffffffff);
    csr_write((uint32_t)(1 + (len + 3) / 4));
    csr_write((uint32_t)len);
    for (size_t i = 0; i < len; i += 4)
    {
        uint32_t w = 0;
        for (size_t b = 0; b < 4 && i + b < len; b++)
            w |= (uint32_t)ptr[i + b] << (8 * b);
        csr_write(w);
    }
}

/* The end of a run: success publishes Output, anything else proves nothing. */
void airbender_exit(int status)
{
    if (status == 0)
        airbender_finish_success(Output);
    airbender_finish_error((uint32_t)status);
}

void airbender_fail(uint32_t code)
{
    airbender_finish_error(code);
}
