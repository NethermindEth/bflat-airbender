/* The u256 delegation's results computed in C, for the platforms without it (see u256.h). */
#include "u256.h"
#include "zkvm.h"

uint32_t u256_op(u256d *a, const u256d *b, uint32_t op)
{
    uint32_t carry = (op & U256_CARRY) ? 1 : 0, flag = 0;
    uint32_t r[8];

    if (op & (U256_ADD | U256_SUB | U256_SUB_AND_NEGATE))
    {
        const uint32_t *l = (op & U256_SUB_AND_NEGATE) ? b->w : a->w;
        const uint32_t *s = (op & U256_SUB_AND_NEGATE) ? a->w : b->w;
        uint64_t c = carry;
        for (int i = 0; i < 8; i++)
        {
            if (op & U256_ADD)
            {
                uint64_t t = (uint64_t)l[i] + s[i] + c;
                r[i] = (uint32_t)t;
                c = t >> 32;
            }
            else
            {
                uint64_t t = (uint64_t)l[i] - s[i] - c;
                r[i] = (uint32_t)t;
                c = (t >> 63) & 1;
            }
        }
        flag = (uint32_t)c;
    }
    else if (op & (U256_MUL_LOW | U256_MUL_HIGH))
    {
        uint32_t t[16] = {0};
        for (int i = 0; i < 8; i++)
        {
            uint64_t c = 0;
            for (int j = 0; j < 8; j++)
            {
                c += (uint64_t)a->w[i] * b->w[j] + t[i + j];
                t[i + j] = (uint32_t)c;
                c >>= 32;
            }
            t[i + 8] = (uint32_t)c;
        }
        for (int i = 0; i < 8; i++)
            r[i] = (op & U256_MUL_LOW) ? t[i] : t[i + 8];
        if (op & U256_MUL_LOW)
            for (int i = 8; i < 16; i++)
                flag |= t[i] != 0;
    }
    else if (op & U256_EQ)
    {
        flag = 1;
        for (int i = 0; i < 8; i++)
        {
            r[i] = a->w[i];
            flag &= a->w[i] == b->w[i];
        }
    }
    else
    {
        uint64_t c = carry;
        for (int i = 0; i < 8; i++)
        {
            uint64_t t = (uint64_t)b->w[i] + c;
            r[i] = (uint32_t)t;
            c = t >> 32;
        }
        flag = carry ? (uint32_t)c : 0;
    }

    for (int i = 0; i < 8; i++)
        a->w[i] = r[i];
    return flag;
}
