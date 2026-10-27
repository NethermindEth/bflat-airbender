/* SHA-256 (FIPS 180-4). */

#include "zkvm.h"

static const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
};

static inline uint32_t rotr32(uint32_t x, unsigned n)
{
    return (x >> n) | (x << (32 - n));
}

#define BSIG0(x) (rotr32(x, 2) ^ rotr32(x, 13) ^ rotr32(x, 22))
#define BSIG1(x) (rotr32(x, 6) ^ rotr32(x, 11) ^ rotr32(x, 25))
#define SSIG0(x) (rotr32(x, 7) ^ rotr32(x, 18) ^ ((x) >> 3))
#define SSIG1(x) (rotr32(x, 17) ^ rotr32(x, 19) ^ ((x) >> 10))

/* One round with the working variables renamed rather than moved: eight rounds bring them back
 * to their places. kw is K[i] + W[i]. Maj(a, b, c) = b ^ ((a ^ b) & (b ^ c)), and this round's a ^ b
 * is the next one's b ^ c: bc carries it. */
#define ROUND(a, b, c, d, e, f, g, h, kw)                                                          \
    do                                                                                             \
    {                                                                                              \
        uint32_t t1 = (h) + BSIG1(e) + ((g) ^ ((e) & ((f) ^ (g)))) + (kw);                         \
        uint32_t ab = (a) ^ (b);                                                                   \
        (d) += t1;                                                                                 \
        (h) = t1 + BSIG0(a) + ((b) ^ (ab & bc));                                                   \
        bc = ab;                                                                                   \
    } while (0)

/* The 64 rounds over kw[i] = K[i] + W[i]. */
static void rounds(uint32_t h[8], const uint32_t kw[64])
{
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    uint32_t bc = b ^ c;
    for (int i = 0; i < 64; i += 8)
    {
        ROUND(a, b, c, d, e, f, g, hh, kw[i]);
        ROUND(hh, a, b, c, d, e, f, g, kw[i + 1]);
        ROUND(g, hh, a, b, c, d, e, f, kw[i + 2]);
        ROUND(f, g, hh, a, b, c, d, e, kw[i + 3]);
        ROUND(e, f, g, hh, a, b, c, d, kw[i + 4]);
        ROUND(d, e, f, g, hh, a, b, c, kw[i + 5]);
        ROUND(c, d, e, f, g, hh, a, b, kw[i + 6]);
        ROUND(b, c, d, e, f, g, hh, a, kw[i + 7]);
    }
    h[0] += a;
    h[1] += b;
    h[2] += c;
    h[3] += d;
    h[4] += e;
    h[5] += f;
    h[6] += g;
    h[7] += hh;
}

/* kw[i] = K[i] + W[i] for the message schedule W of a block whose first sixteen words are w. */
static void schedule(uint32_t kw[64], uint32_t w[16])
{
    for (int i = 0; i < 16; i++)
        kw[i] = K[i] + w[i];
    /* Sixteen steps a round, so every index into the window is a constant. */
    for (int i = 16; i < 64; i += 16)
    {
#pragma clang loop unroll(full)
        for (int j = 0; j < 16; j++)
        {
            uint32_t next = w[j] + SSIG0(w[(j + 1) & 15]) + w[(j + 9) & 15] + SSIG1(w[(j + 14) & 15]);
            w[j] = next;
            kw[i + j] = K[i + j] + next;
        }
    }
}

static void compress(uint32_t h[8], const uint8_t block[64])
{
    uint32_t w[16], kw[64];
    for (int i = 0; i < 16; i++)
        w[i] = ((uint32_t)block[4 * i] << 24) | ((uint32_t)block[4 * i + 1] << 16) |
               ((uint32_t)block[4 * i + 2] << 8) | block[4 * i + 3];
    schedule(kw, w);
    rounds(h, kw);
}

/* The padding block of a 64-byte message, the SSZ pair hash, is a constant: its schedule is
 * computed once. */
static uint32_t PadKw64[64];
static int PadKw64Ready;

static void compress_pad64(uint32_t h[8])
{
    if (!PadKw64Ready)
    {
        uint32_t w[16] = {0x80000000};
        w[15] = 512;
        schedule(PadKw64, w);
        PadKw64Ready = 1;
    }
    rounds(h, PadKw64);
}

zkvm_status zkvm_sha256(const uint8_t *data, size_t len, uint8_t output[32])
{
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    uint64_t bits = (uint64_t)len * 8;

    while (len >= 64)
    {
        compress(h, data);
        data += 64;
        len -= 64;
    }

    if (len == 0 && bits == 512)
        compress_pad64(h);
    else
    {
        uint8_t tail[128] = {0};
        memcpy(tail, data, len);
        tail[len] = 0x80;
        size_t tailLen = len < 56 ? 64 : 128;
        for (int i = 0; i < 8; i++)
            tail[tailLen - 1 - i] = (uint8_t)(bits >> (8 * i));
        compress(h, tail);
        if (tailLen == 128)
            compress(h, tail + 64);
    }

    for (int i = 0; i < 8; i++)
    {
        output[4 * i] = (uint8_t)(h[i] >> 24);
        output[4 * i + 1] = (uint8_t)(h[i] >> 16);
        output[4 * i + 2] = (uint8_t)(h[i] >> 8);
        output[4 * i + 3] = (uint8_t)h[i];
    }
    return ZKVM_EOK;
}
