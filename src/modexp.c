/*
 * MODEXP (EIP-198) on 32-bit limbs: left-to-right square-and-multiply, each product reduced
 * with Knuth's algorithm D.
 *
 * Operands are capped at 1024 bytes, the EIP-7823 bound from Osaka on, which keeps every
 * buffer static: the runtime has no allocator to lend. Larger operands fail.
 */

#include "u256.h"
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

/*
 * An odd modulus of at most 1024 bits: Montgomery exponentiation on the u256 delegation, R = 2^(256 k)
 * for k 256-bit digits, by 4-bit windows. A product accumulates digit by digit in T, 2k + 1 digits, and
 * every row adds the multiple of the modulus that clears its lowest digit (CIOS on 256-bit words). R mod m
 * and R^2 mod m come from the division above.
 */
enum
{
    MaxDigits = 4
};

static u256d MontN[MaxDigits], MontNPrime, MontZero2, MontX, MontM2, MontT[2 * MaxDigits + 1], MontD[MaxDigits + 1];
static u256d MontTable[16][MaxDigits], MontValue[MaxDigits], MontR2[MaxDigits], MontOne[MaxDigits];

/* T[d..] += *v. */
static void accumulate_digit(int d, const u256d *v)
{
    uint32_t carry = u256_op(&MontT[d], v, U256_ADD);
    while (carry)
        carry = u256_op(&MontT[++d], &MontZero2, U256_ADD | U256_CARRY);
}

/* r = a b / R mod N, for a, b < N; r may alias a or b. */
static void mont_mul_k(u256d *r, const u256d *a, const u256d *b, int k)
{
    for (int d = 0; d <= 2 * k; d++)
        u256_op(&MontT[d], &MontZero2, U256_MEMCOPY);
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            u256_op(&MontX, &a[j], U256_MEMCOPY);
            u256_op(&MontX, &b[i], U256_MUL_LOW);
            accumulate_digit(i + j, &MontX);
            u256_op(&MontX, &a[j], U256_MEMCOPY);
            u256_op(&MontX, &b[i], U256_MUL_HIGH);
            accumulate_digit(i + j + 1, &MontX);
        }
        u256_op(&MontM2, &MontT[i], U256_MEMCOPY);
        u256_op(&MontM2, &MontNPrime, U256_MUL_LOW);
        for (int j = 0; j < k; j++)
        {
            u256_op(&MontX, &MontM2, U256_MEMCOPY);
            u256_op(&MontX, &MontN[j], U256_MUL_LOW);
            accumulate_digit(i + j, &MontX);
            u256_op(&MontX, &MontM2, U256_MEMCOPY);
            u256_op(&MontX, &MontN[j], U256_MUL_HIGH);
            accumulate_digit(i + j + 1, &MontX);
        }
    }

    /* The result, T[k..2k], is below 2N: subtract N once when that does not borrow. */
    uint32_t borrow = 0;
    for (int d = 0; d < k; d++)
    {
        u256_op(&MontD[d], &MontT[k + d], U256_MEMCOPY);
        borrow = u256_op(&MontD[d], &MontN[d], U256_SUB | (borrow ? U256_CARRY : 0));
    }
    const u256d *src = (borrow && u256_op(&MontT[2 * k], &MontZero2, U256_EQ)) ? &MontT[k] : MontD;
    for (int d = 0; d < k; d++)
        u256_op(&r[d], &src[d], U256_MEMCOPY);
}

static void load_digits(u256d *r, const uint32_t *limbs, size_t len, int k)
{
    for (int d = 0; d < k; d++)
        for (int w = 0; w < 8; w++)
        {
            size_t i = (size_t)(8 * d + w);
            r[d].w[w] = i < len ? limbs[i] : 0;
        }
}

static void modexp_mont(uint32_t *acc, const uint32_t *base, const uint8_t *exp, size_t exp_len, size_t mLen)
{
    int k = (int)((mLen + 7) / 8);
    load_digits(MontN, ModLimbs, mLen, k);

    /* -N^-1 mod 2^256 by Newton on the low digit: it is its own inverse mod 8, and each step doubles the bits. */
    static u256d x, t, two;
    x = MontN[0];
    for (int step = 0; step < 7; step++)
    {
        t = MontN[0];
        u256_op(&t, &x, U256_MUL_LOW);
        two = (u256d){{2}};
        u256_op(&two, &t, U256_SUB);
        u256_op(&x, &two, U256_MUL_LOW);
    }
    MontNPrime = (u256d){{0}};
    u256_op(&MontNPrime, &x, U256_SUB);

    /* R and R^2 mod N. */
    size_t rLimbs = (size_t)(8 * k);
    memset(Product, 0, (2 * rLimbs + 1) * sizeof(uint32_t));
    Product[rLimbs] = 1;
    reduce(Acc, Product, rLimbs + 1, ModLimbs, mLen);
    load_digits(MontOne, Acc, mLen, k);
    memset(Product, 0, (2 * rLimbs + 1) * sizeof(uint32_t));
    Product[2 * rLimbs] = 1;
    reduce(Acc, Product, 2 * rLimbs + 1, ModLimbs, mLen);
    load_digits(MontR2, Acc, mLen, k);

    load_digits(MontValue, base, mLen, k);
    for (int d = 0; d < k; d++)
        MontTable[0][d] = MontOne[d];
    mont_mul_k(MontTable[1], MontValue, MontR2, k);
    for (int i = 2; i < 16; i++)
        mont_mul_k(MontTable[i], MontTable[i - 1], MontTable[1], k);

    for (int d = 0; d < k; d++)
        MontValue[d] = MontOne[d];
    int started = 0;
    for (size_t i = 0; i < exp_len; i++)
    {
        for (int half = 1; half >= 0; half--)
        {
            uint32_t bits = (exp[i] >> (4 * half)) & 15;
            if (started)
            {
                for (int s = 0; s < 4; s++)
                    mont_mul_k(MontValue, MontValue, MontValue, k);
            }
            if (bits)
            {
                mont_mul_k(MontValue, MontValue, MontTable[bits], k);
                started = 1;
            }
        }
    }

    /* Out of the domain: times 1. */
    static u256d one[MaxDigits];
    for (int d = 0; d < k; d++)
        one[d] = (u256d){{d == 0 ? 1u : 0u}};
    mont_mul_k(MontValue, MontValue, one, k);
    for (size_t i = 0; i < mLen; i++)
        acc[i] = MontValue[i / 8].w[i % 8];
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

    if (mLen <= 8 * MaxDigits && (ModLimbs[0] & 1) && !(mLen == 1 && ModLimbs[0] == 1))
    {
        modexp_mont(Acc, Base, exp, exp_len, mLen);
    }
    else
    {
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
    }

    for (size_t i = 0; i < mod_len; i++)
    {
        size_t bit = (mod_len - 1 - i) * 8;
        output[i] = bit / 32 < mLen ? (uint8_t)(Acc[bit / 32] >> (bit % 32)) : 0;
    }
    return ZKVM_EOK;
}
