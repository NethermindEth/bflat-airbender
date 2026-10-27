/* The BLAKE2b compression function F (EIP-152), all values little-endian. */
#include "zkvm.h"

static const uint64_t IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL,
};

static const uint8_t Sigma[10][16] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3},
    {11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4}, {7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8},
    {9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13}, {2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9},
    {12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11}, {13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10},
    {6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5}, {10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0},
};

static inline uint64_t rotr64(uint64_t x, unsigned n)
{
    return (x >> n) | (x << (64 - n));
}

static inline uint64_t load64(const uint8_t *p)
{
    uint64_t v = 0;
    for (int b = 7; b >= 0; b--)
        v = (v << 8) | p[b];
    return v;
}

static inline void store64(uint8_t *p, uint64_t v)
{
    for (int b = 0; b < 8; b++)
        p[b] = (uint8_t)(v >> (8 * b));
}

#define G(a, b, c, d, x, y)                                                                                            \
    do                                                                                                                 \
    {                                                                                                                  \
        v[a] = v[a] + v[b] + (x);                                                                                      \
        v[d] = rotr64(v[d] ^ v[a], 32);                                                                                \
        v[c] = v[c] + v[d];                                                                                            \
        v[b] = rotr64(v[b] ^ v[c], 24);                                                                                \
        v[a] = v[a] + v[b] + (y);                                                                                      \
        v[d] = rotr64(v[d] ^ v[a], 16);                                                                                \
        v[c] = v[c] + v[d];                                                                                            \
        v[b] = rotr64(v[b] ^ v[c], 63);                                                                                \
    } while (0)

zkvm_status zkvm_blake2f(uint32_t rounds, uint8_t h[64], const uint8_t m[128], const uint8_t t[16], uint8_t f)
{
    uint64_t hw[8], mw[16], v[16];
    for (int i = 0; i < 8; i++)
        hw[i] = load64(h + 8 * i);
    for (int i = 0; i < 16; i++)
        mw[i] = load64(m + 8 * i);

    for (int i = 0; i < 8; i++)
    {
        v[i] = hw[i];
        v[i + 8] = IV[i];
    }
    v[12] ^= load64(t);
    v[13] ^= load64(t + 8);
    if (f)
        v[14] = ~v[14];

    for (uint32_t r = 0; r < rounds; r++)
    {
        const uint8_t *s = Sigma[r % 10];
        G(0, 4, 8, 12, mw[s[0]], mw[s[1]]);
        G(1, 5, 9, 13, mw[s[2]], mw[s[3]]);
        G(2, 6, 10, 14, mw[s[4]], mw[s[5]]);
        G(3, 7, 11, 15, mw[s[6]], mw[s[7]]);
        G(0, 5, 10, 15, mw[s[8]], mw[s[9]]);
        G(1, 6, 11, 12, mw[s[10]], mw[s[11]]);
        G(2, 7, 8, 13, mw[s[12]], mw[s[13]]);
        G(3, 4, 9, 14, mw[s[14]], mw[s[15]]);
    }

    for (int i = 0; i < 8; i++)
        store64(h + 8 * i, hw[i] ^ v[i] ^ v[i + 8]);
    return ZKVM_EOK;
}
