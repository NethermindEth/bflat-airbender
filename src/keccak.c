/* Keccak-f[1600] and Keccak-256 (the original padding, 0x01, as Ethereum uses). */

#include "zkvm.h"

static const uint64_t RoundConstants[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL, 0x8000000080008000ULL,
    0x000000000000808bULL, 0x0000000080000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008aULL, 0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL, 0x000000000000800aULL, 0x800000008000000aULL,
    0x8000000080008081ULL, 0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL,
};

static const uint8_t Rotations[24] = {
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44,
};

static const uint8_t PiLanes[24] = {
    10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1,
};

static inline uint64_t rotl64(uint64_t x, unsigned n)
{
    return (x << n) | (x >> (64 - n));
}

void keccak_f1600(uint64_t a[25])
{
    for (int round = 0; round < 24; round++)
    {
        uint64_t c[5];
        for (int x = 0; x < 5; x++)
            c[x] = a[x] ^ a[x + 5] ^ a[x + 10] ^ a[x + 15] ^ a[x + 20];
        for (int x = 0; x < 5; x++)
        {
            uint64_t d = c[(x + 4) % 5] ^ rotl64(c[(x + 1) % 5], 1);
            for (int y = 0; y < 25; y += 5)
                a[y + x] ^= d;
        }

        uint64_t carried = a[1];
        for (int i = 0; i < 24; i++)
        {
            int lane = PiLanes[i];
            uint64_t next = a[lane];
            a[lane] = rotl64(carried, Rotations[i]);
            carried = next;
        }

        for (int y = 0; y < 25; y += 5)
        {
            uint64_t row[5];
            for (int x = 0; x < 5; x++)
                row[x] = a[y + x];
            for (int x = 0; x < 5; x++)
                a[y + x] = row[x] ^ (~row[(x + 1) % 5] & row[(x + 2) % 5]);
        }

        a[0] ^= RoundConstants[round];
    }
}

/* The ZisK-named permutation entry the managed Keccak reaches in zkVM builds. */
void syscall_keccak_f(uint64_t state[25])
{
    keccak_f1600(state);
}

zkvm_status zkvm_keccak256(const uint8_t *data, size_t len, uint8_t output[32])
{
    enum { Rate = 136 };
    uint64_t state[25] = {0};
    uint8_t block[Rate];

    for (;;)
    {
        size_t take = len < Rate ? len : Rate;
        memcpy(block, data, take);
        if (take < Rate)
        {
            memset(block + take, 0, Rate - take);
            block[take] ^= 0x01;
            block[Rate - 1] ^= 0x80;
        }
        for (int i = 0; i < Rate / 8; i++)
        {
            uint64_t lane = 0;
            for (int b = 7; b >= 0; b--)
                lane = (lane << 8) | block[i * 8 + b];
            state[i] ^= lane;
        }
        keccak_f1600(state);
        if (take < Rate)
            break;
        data += Rate;
        len -= Rate;
    }

    for (int i = 0; i < 4; i++)
        for (int b = 0; b < 8; b++)
            output[i * 8 + b] = (uint8_t)(state[i] >> (8 * b));
    return ZKVM_EOK;
}
