/*
 * MODEXP (EIP-198) on 32-bit limbs: left-to-right square-and-multiply, each product reduced
 * with Knuth's algorithm D.
 *
 * Operands are capped at 1024 bytes, the EIP-7823 bound from Osaka on, which keeps every
 * buffer static: the runtime has no allocator to lend. Larger operands fail.
 */

#include "zkvm.h"

enum
{
    MaxBytes = 1024,
    MaxLimbs = MaxBytes / 4
};

static uint32_t ModLimbs[MaxLimbs];
static uint32_t Acc[MaxLimbs];
static uint32_t Base[MaxLimbs];
static uint32_t Product[2 * MaxLimbs + 1];
static uint32_t NormMod[MaxLimbs];

/* Big-endian bytes to little-endian limbs; returns the limb count without leading zero limbs. */
static size_t load(uint32_t *r, size_t cap, const uint8_t *in, size_t len)
{
    memset(r, 0, cap * sizeof(uint32_t));
    for (size_t i = 0; i < len; i++)
    {
        size_t bit = (len - 1 - i) * 8;
        r[bit / 32] |= (uint32_t)in[i] << (bit % 32);
    }
    size_t n = (len + 3) / 4;
    while (n > 0 && r[n - 1] == 0)
        n--;
    return n;
}

/* r = num mod m, num has numLen limbs and is destroyed; m has mLen > 0 limbs, top limb non-zero. */
static void reduce(uint32_t *r, uint32_t *num, size_t numLen, const uint32_t *m, size_t mLen)
{
    if (numLen < mLen)
    {
        memset(r, 0, mLen * sizeof(uint32_t));
        memcpy(r, num, numLen * sizeof(uint32_t));
        return;
    }

    if (mLen == 1)
    {
        uint64_t rem = 0;
        for (size_t i = numLen; i-- > 0;)
            rem = ((rem << 32) | num[i]) % m[0];
        r[0] = (uint32_t)rem;
        return;
    }

    /* D1: normalise so the divisor's top bit is set. */
    unsigned shift = (unsigned)__builtin_clz(m[mLen - 1]);
    for (size_t i = mLen; i-- > 0;)
        NormMod[i] = (m[i] << shift) | (shift && i > 0 ? m[i - 1] >> (32 - shift) : 0);
    num[numLen] = shift ? num[numLen - 1] >> (32 - shift) : 0;
    for (size_t i = numLen; i-- > 0;)
        num[i] = (num[i] << shift) | (shift && i > 0 ? num[i - 1] >> (32 - shift) : 0);

    for (size_t j = numLen - mLen + 1; j-- > 0;)
    {
        /* D3: estimate the quotient digit from the top two limbs. */
        uint64_t top = ((uint64_t)num[j + mLen] << 32) | num[j + mLen - 1];
        uint64_t qhat = top / NormMod[mLen - 1];
        uint64_t rhat = top % NormMod[mLen - 1];
        while (qhat >> 32 || qhat * NormMod[mLen - 2] > ((rhat << 32) | num[j + mLen - 2]))
        {
            qhat--;
            rhat += NormMod[mLen - 1];
            if (rhat >> 32)
                break;
        }

        /* D4: multiply and subtract. */
        int64_t borrow = 0;
        uint64_t carry = 0;
        for (size_t i = 0; i < mLen; i++)
        {
            uint64_t p = qhat * NormMod[i] + carry;
            carry = p >> 32;
            int64_t t = (int64_t)num[i + j] - (int64_t)(uint32_t)p + borrow;
            num[i + j] = (uint32_t)t;
            borrow = t >> 32;
        }
        int64_t t = (int64_t)num[j + mLen] - (int64_t)carry + borrow;
        num[j + mLen] = (uint32_t)t;

        /* D6: add back when the estimate was one too large. */
        if (t < 0)
        {
            uint64_t c = 0;
            for (size_t i = 0; i < mLen; i++)
            {
                c += (uint64_t)num[i + j] + NormMod[i];
                num[i + j] = (uint32_t)c;
                c >>= 32;
            }
            num[j + mLen] += (uint32_t)c;
        }
    }

    /* D8: unnormalise the remainder. */
    for (size_t i = 0; i < mLen; i++)
        r[i] = (num[i] >> shift) | (shift && i + 1 < mLen + 1 ? (uint32_t)((uint64_t)num[i + 1] << (32 - shift)) : 0);
}

static void mul_mod(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *m, size_t mLen)
{
    memset(Product, 0, (2 * mLen + 1) * sizeof(uint32_t));
    for (size_t i = 0; i < mLen; i++)
    {
        uint64_t carry = 0;
        for (size_t j = 0; j < mLen; j++)
        {
            carry += (uint64_t)a[i] * b[j] + Product[i + j];
            Product[i + j] = (uint32_t)carry;
            carry >>= 32;
        }
        Product[i + mLen] = (uint32_t)carry;
    }
    reduce(r, Product, 2 * mLen, m, mLen);
}

zkvm_status zkvm_modexp(const uint8_t *base, size_t base_len, const uint8_t *exp, size_t exp_len,
                        const uint8_t *modulus, size_t mod_len, uint8_t *output)
{
    if (base_len > MaxBytes || exp_len > MaxBytes || mod_len > MaxBytes)
        return ZKVM_EFAIL;

    memset(output, 0, mod_len);
    size_t mLen = load(ModLimbs, MaxLimbs, modulus, mod_len);
    if (mLen == 0)
        return ZKVM_EOK;

    size_t baseLen = load(Product, 2 * MaxLimbs + 1, base, base_len);
    reduce(Base, Product, baseLen, ModLimbs, mLen);

    /* 1 mod m: zero when m is 1. */
    memset(Acc, 0, mLen * sizeof(uint32_t));
    Acc[0] = (mLen == 1 && ModLimbs[0] == 1) ? 0 : 1;

    for (size_t i = 0; i < exp_len; i++)
    {
        for (int bit = 7; bit >= 0; bit--)
        {
            mul_mod(Acc, Acc, Acc, ModLimbs, mLen);
            if ((exp[i] >> bit) & 1)
                mul_mod(Acc, Acc, Base, ModLimbs, mLen);
        }
    }

    for (size_t i = 0; i < mod_len; i++)
    {
        size_t bit = (mod_len - 1 - i) * 8;
        output[i] = bit / 32 < mLen ? (uint8_t)(Acc[bit / 32] >> (bit % 32)) : 0;
    }
    return ZKVM_EOK;
}
