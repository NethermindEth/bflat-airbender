/* Host harness for mpt.c: a portable keccak and a witness held in arrays. */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL, 0x8000000080008000ULL,
    0x000000000000808bULL, 0x0000000080000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008aULL, 0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL, 0x000000000000800aULL, 0x800000008000000aULL,
    0x8000000080008081ULL, 0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL};
static const int ROT[24] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44};
static const int PIL[24] = {10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1};

static void keccakf(uint64_t s[25])
{
    for (int r = 0; r < 24; r++)
    {
        uint64_t bc[5], t;
        for (int i = 0; i < 5; i++)
            bc[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
        for (int i = 0; i < 5; i++)
        {
            t = bc[(i + 4) % 5] ^ ((bc[(i + 1) % 5] << 1) | (bc[(i + 1) % 5] >> 63));
            for (int j = 0; j < 25; j += 5)
                s[j + i] ^= t;
        }
        t = s[1];
        for (int i = 0; i < 24; i++)
        {
            int j = PIL[i];
            uint64_t b = s[j];
            s[j] = (t << ROT[i]) | (t >> (64 - ROT[i]));
            t = b;
        }
        for (int j = 0; j < 25; j += 5)
        {
            for (int i = 0; i < 5; i++)
                bc[i] = s[j + i];
            for (int i = 0; i < 5; i++)
                s[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }
        s[0] ^= RC[r];
    }
}

void t_keccak(const uint8_t *data, uint32_t len, uint8_t *out)
{
    uint64_t s[25];
    uint8_t block[136];
    memset(s, 0, sizeof s);
    while (len >= 136)
    {
        for (int i = 0; i < 17; i++)
        {
            uint64_t w;
            memcpy(&w, data + 8 * i, 8);
            s[i] ^= w;
        }
        keccakf(s);
        data += 136;
        len -= 136;
    }
    memset(block, 0, 136);
    memcpy(block, data, len);
    block[len] ^= 0x01;
    block[135] ^= 0x80;
    for (int i = 0; i < 17; i++)
    {
        uint64_t w;
        memcpy(&w, block + 8 * i, 8);
        s[i] ^= w;
    }
    keccakf(s);
    memcpy(out, s, 32);
}

static uint32_t WCount, WCap;
static uint32_t *WHash;
static const uint8_t **WRlp;
static uint32_t *WLen;
static int EditedCalls;

static uint32_t mpt_w_count(void) { return WCount; }
static const uint32_t *mpt_w_hash(uint32_t k) { return WHash + 8 * k; }
static const uint8_t *mpt_w_rlp(uint32_t k, uint32_t *len) { *len = WLen[k]; return WRlp[k]; }
static void mpt_keccak(const uint8_t *data, uint32_t len, uint32_t out[8]) { t_keccak(data, len, (uint8_t *)out); }
static void mpt_hash_edited(uint32_t k, int block, const uint8_t *original, const uint8_t *edited, uint32_t out[8])
{
    uint8_t whole[532];
    (void)k;
    EditedCalls++;
    memcpy(whole, original, (size_t)136 * block);
    memcpy(whole + 136 * block, edited + 136 * block, 532 - (size_t)136 * block);
    t_keccak(whole, 532, (uint8_t *)out);
}
static void mpt_keccak_key(const uint32_t *words, uint32_t len, uint32_t out[8]) { t_keccak((const uint8_t *)words, len, (uint8_t *)out); }
static void mpt_copy(uint8_t *to, const uint8_t *from, uint32_t len)
{
    /* As the delegation copies: whole 32-byte blocks around the range. */
    uintptr_t a = (uintptr_t)from & ~(uintptr_t)31, end = ((uintptr_t)from + len + 31) & ~(uintptr_t)31;
    memcpy(to - ((uintptr_t)from - a), (const void *)a, end - a);
}

#include "mpt.c"

void t_reset(void)
{
    WCount = 0;
    MptSlots = 0;
    MptInfo = 0;
    MptNext = 0;
    memset(MptMemo, 0, sizeof MptMemo);
}

void t_add_node(const uint8_t *rlp, uint32_t len)
{
    if (WCount == WCap)
    {
        WCap = WCap ? 2 * WCap : 1024;
        WHash = realloc(WHash, (size_t)WCap * 32);
        WRlp = realloc(WRlp, (size_t)WCap * sizeof *WRlp);
        WLen = realloc(WLen, (size_t)WCap * 4);
    }
    /* As in the input: word aligned, with bytes before and after it. */
    uint8_t *copy = (uint8_t *)malloc(len + 128) + 40;
    memcpy(copy, rlp, len);
    WRlp[WCount] = copy;
    WLen[WCount] = len;
    t_keccak(rlp, len, (uint8_t *)(WHash + 8 * WCount));
    WCount++;
}

int t_edited_calls(void) { return EditedCalls; }
