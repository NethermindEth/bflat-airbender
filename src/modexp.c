/*
 * MODEXP (EIP-198) on 32-bit limbs: left-to-right square-and-multiply, each product reduced
 * with Knuth's algorithm D.
 *
 * Operands are capped at 1024 bytes, the EIP-7823 bound from Osaka on, which keeps every
 * buffer static: the runtime has no allocator to lend. Larger operands fail.
 */

#include "mont.h"
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

/* high:low / divisor, for a divisor with its top bit set and above high; the remainder in *remainder.
   Hacker's Delight divlu: two 32-bit divisions by the divisor's top half, each corrected at most twice. */
static uint32_t divide_normalized(uint32_t high, uint32_t low, uint32_t divisor, uint32_t *remainder)
{
    const uint32_t base = 1u << 16;
    uint32_t vn1 = divisor >> 16, vn0 = divisor & 0xFFFF;
    uint32_t un1 = low >> 16, un0 = low & 0xFFFF;

    uint32_t q1 = high / vn1;
    uint32_t rhat = high - q1 * vn1;
    while (q1 >= base || q1 * vn0 > ((rhat << 16) | un1))
    {
        q1--;
        rhat += vn1;
        if (rhat >= base)
            break;
    }

    uint32_t un21 = (high << 16) + un1 - q1 * divisor;
    uint32_t q0 = un21 / vn1;
    rhat = un21 - q0 * vn1;
    while (q0 >= base || q0 * vn0 > ((rhat << 16) | un0))
    {
        q0--;
        rhat += vn1;
        if (rhat >= base)
            break;
    }

    *remainder = (un21 << 16) + un0 - q0 * divisor;
    return (q1 << 16) | q0;
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
        /* D3: estimate the quotient digit from the top two limbs. The top limb is at most the divisor's, so the
           estimate fits a limb once capped, and comes from two 32-bit divisions rather than a 64-bit one. */
        uint32_t qhat;
        uint64_t rhat;
        if (num[j + mLen] >= NormMod[mLen - 1])
        {
            qhat = UINT32_MAX;
            rhat = (uint64_t)num[j + mLen - 1] + NormMod[mLen - 1];
        }
        else
        {
            uint32_t r;
            qhat = divide_normalized(num[j + mLen], num[j + mLen - 1], NormMod[mLen - 1], &r);
            rhat = r;
        }
        while (!(rhat >> 32) && (uint64_t)qhat * NormMod[mLen - 2] > ((rhat << 32) | num[j + mLen - 2]))
        {
            qhat--;
            rhat += NormMod[mLen - 1];
        }

        /* D4: multiply and subtract. */
        int64_t borrow = 0;
        uint64_t carry = 0;
        for (size_t i = 0; i < mLen; i++)
        {
            uint64_t p = (uint64_t)qhat * NormMod[i] + carry;
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

static mont_ctx MontSingle;
static u256d MontN[MaxDigits], MontNPrime, MontZero2, MontX, MontM2, MontT[2 * MaxDigits + 1], MontD[MaxDigits + 1];
static u256d MontTable[16][MaxDigits], MontValue[MaxDigits], MontR2[MaxDigits], MontOne[MaxDigits];

/* T[d..] += *v. */
static void accumulate_digit(int d, const u256d *v)
{
    uint32_t carry = u256_op(&MontT[d], v, U256_ADD);
    while (carry)
        carry = u256_op(&MontT[++d], &MontZero2, U256_ADD | U256_CARRY);
}

static void mont_mul_multi(u256d *r, const u256d *a, const u256d *b, int k)
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

/* Inline so mont_mul can resolve operand aliasing at each square or multiply. */
MONT_INLINE void mont_mul_k(u256d *r, const u256d *a, const u256d *b, int k)
{
    if (k == 1)
        mont_mul(r, a, b, &MontSingle);
    else
        mont_mul_multi(r, a, b, k);
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

MONT_INLINE void mont_pow_windows(const uint8_t *exp, size_t exp_len, int k)
{
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
}

/* MontSingle.ninv = n^-1 mod 2^256 for the low digit n of the modulus, by Newton: n is its own inverse mod 8, and each
   step doubles the bits. */
static void mont_inverse(const u256d *n)
{
    static u256d t, two;
    MontSingle.ninv = *n;
    for (int step = 0; step < 7; step++)
    {
        t = *n;
        u256_op(&t, &MontSingle.ninv, U256_MUL_LOW);
        two = (u256d){{2}};
        u256_op(&two, &t, U256_SUB);
        u256_op(&MontSingle.ninv, &two, U256_MUL_LOW);
    }
}

static void modexp_mont(uint32_t *acc, const uint32_t *base, const uint8_t *exp, size_t exp_len, size_t mLen)
{
    int k = (int)((mLen + 7) / 8);
    u256d *n = k == 1 ? &MontSingle.n : MontN;
    load_digits(n, ModLimbs, mLen, k);

    mont_inverse(&n[0]);
    if (k > 1)
    {
        MontNPrime = (u256d){{0}};
        u256_op(&MontNPrime, &MontSingle.ninv, U256_SUB);
    }

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
    /* Specialize the window loop once so each product uses its known digit path. */
    if (k == 1)
        mont_pow_windows(exp, exp_len, 1);
    else
        mont_pow_windows(exp, exp_len, k);
    for (size_t i = 0; i < mLen; i++)
        acc[i] = MontValue[i / 8].w[i % 8];
}

/* Primes of at most 256 bits, little-endian limbs, for which x^(m - 2) mod m is the inverse of x (Fermat): BN254's
   scalar and base fields, secp256k1's and P-256's base and scalar fields, BLS12-381's scalar field. SNARK verifiers
   invert modulo the BN254 scalar field this way. */
static const uint32_t InversePrimes[][8] = {
    {0xf0000001, 0x43e1f593, 0x79b97091, 0x2833e848, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72},
    {0xd87cfd47, 0x3c208c16, 0x6871ca8d, 0x97816a91, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72},
    {0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff},
    {0xd0364141, 0xbfd25e8c, 0xaf48a03b, 0xbaaedce6, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff},
    {0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0xffffffff},
    {0xfc632551, 0xf3b9cac2, 0xa7179e84, 0xbce6faad, 0xffffffff, 0xffffffff, 0x00000000, 0xffffffff},
    {0x00000001, 0xffffffff, 0xfffe5bfe, 0x53bda402, 0x09a1d805, 0x3339d808, 0x299d7d48, 0x73eda753},
};

enum
{
    InverseNotHinted = -1, /* not an inversion modulo one of InversePrimes, or no hint tape */
    InverseComputed = 0,   /* an inversion whose hint is missing or wrong: the caller computes it and gives it out */
    InverseDone = 1        /* acc holds the result */
};

/*
 * x^(m - 2) mod m, with m in ModLimbs (mLen limbs) one of InversePrimes, is the inverse of x. With a hint tape, the next
 * value y is the result when it is below m and x y = 1 mod m, checked as mont(x, y) = mont(1, 1): two Montgomery
 * products where the exponentiation takes some 340. x = 0 gives 0 and takes no hint.
 */
static int inverse_hinted(uint32_t *acc, const uint32_t *base, const uint8_t *exp, size_t exp_len, size_t mLen)
{
    if (!FieldHints.active || mLen > 8)
        return InverseNotHinted;
    uint32_t m[8] = {0};
    for (size_t i = 0; i < mLen; i++)
        m[i] = ModLimbs[i];
    int known = 0;
    for (size_t p = 0; p < sizeof(InversePrimes) / sizeof(InversePrimes[0]) && !known; p++)
    {
        uint32_t diff = 0;
        for (int i = 0; i < 8; i++)
            diff |= m[i] ^ InversePrimes[p][i];
        known = diff == 0;
    }
    if (!known)
        return InverseNotHinted;

    /* The exponent, big-endian with any number of leading zero bytes, must be m - 2. */
    uint32_t e[8], borrow = 2;
    for (int i = 0; i < 8; i++)
    {
        e[i] = m[i] - borrow;
        borrow = m[i] < borrow;
    }
    for (size_t i = 0; i < exp_len; i++)
    {
        size_t bit = (exp_len - 1 - i) * 8;
        if (exp[i] != (bit < 256 ? (uint8_t)(e[bit / 32] >> (bit % 32)) : 0))
            return InverseNotHinted;
    }
    for (size_t bit = exp_len * 8; bit < 256; bit += 8)
        if ((uint8_t)(e[bit / 32] >> (bit % 32)) != 0)
            return InverseNotHinted;

    uint32_t nonzero = 0;
    for (size_t i = 0; i < mLen; i++)
        nonzero |= base[i];
    if (!nonzero)
    {
        memset(acc, 0, mLen * sizeof(uint32_t));
        return InverseDone;
    }

    if (FieldHints.next == FieldHints.end)
        return InverseComputed;
    const u256d *y = FieldHints.next++;
    int below = 0;
    for (int i = 7; i >= 0; i--)
        if (y->w[i] != m[i])
        {
            below = y->w[i] < m[i];
            break;
        }
    if (!below)
        return InverseComputed;

    static u256d value, hint, product, unit;
    load_digits(&MontSingle.n, ModLimbs, mLen, 1);
    mont_inverse(&MontSingle.n);
    load_digits(&value, base, mLen, 1);
    u256_op(&hint, y, U256_MEMCOPY);
    mont_mul(&product, &value, &hint, &MontSingle);
    unit = (u256d){{1}};
    mont_mul(&unit, &unit, &unit, &MontSingle);
    if (!u256_op(&product, &unit, U256_EQ))
        return InverseComputed;
    for (size_t i = 0; i < mLen; i++)
        acc[i] = hint.w[i];
    return InverseDone;
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
        int inverse = inverse_hinted(Acc, Base, exp, exp_len, mLen);
        if (inverse != InverseDone)
            modexp_mont(Acc, Base, exp, exp_len, mLen);
        if (inverse == InverseComputed)
        {
            static u256d computed;
            load_digits(&computed, Acc, mLen, 1);
            field_hint_emit(&computed);
        }
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
