/*
 * Field, point and pairing arithmetic shared by the BN254 and BLS12-381 accelerators; see ec.h.
 * Every routine works on the curve selected by ec_use, so one copy of the code serves both.
 */

#include "ec.h"
#include "u256.h"

/* Composite operations stay out of line: inlined, the arithmetic would double the code. */
#define NOINLINE __attribute__((noinline))

const curve *C;

/* ---- Fp ------------------------------------------------------------------------------------- */

static int limbs_cmp(const uint32_t *a, const uint32_t *b, uint32_t n)
{
    for (uint32_t i = n; i-- > 0;)
        if (a[i] != b[i])
            return a[i] < b[i] ? -1 : 1;
    return 0;
}

/*
 * The field operations run on the u256 delegation, which takes its operands from RAM, 32-byte aligned
 * and distinct, and writes its result over the first. An element is two 256-bit digits (ec.h), only
 * the first in use for BN254, and the curve tables are in RAM, so the operands are read in place.
 */
static u256d M, X, T[4], Zero, Pinv, PInverse, P0, P1;
static fp Scratch;
static const curve *Prepared;

static inline u256d *digit(const fp *a, int i)
{
    return (u256d *)(uintptr_t)&a->l[8 * i];
}

/* Pinv = -p^-1 mod 2^256, PInverse = p^-1 mod 2^256 and p as digits, once per curve. */
static __attribute__((noinline)) void prepare(void)
{
    static u256d x, t, two;
    for (int i = 0; i < 8; i++)
    {
        P0.w[i] = C->p.l[i];
        P1.w[i] = C->p.l[i + 8];
    }
    /* p^-1 mod 2^256 by Newton: p is its own inverse mod 8, and each step doubles the bits. */
    x = P0;
    for (int step = 0; step < 7; step++)
    {
        u256_op(&t, &P0, U256_MEMCOPY);
        u256_op(&t, &x, U256_MUL_LOW);
        two = (u256d){{2}};
        u256_op(&two, &t, U256_SUB);
        u256_op(&x, &two, U256_MUL_LOW);
    }
    PInverse = x;
    Pinv = (u256d){{0}};
    u256_op(&Pinv, &x, U256_SUB);
    Prepared = C;
}

void ec_prepare(void)
{
    if (Prepared != C)
        prepare();
}

/* r = a, digit by digit; r and a may be the same. */
static inline void fp_copy(fp *r, const fp *a, uint32_t n)
{
    if (r == a)
        return;
    u256_op(digit(r, 0), digit(a, 0), U256_MEMCOPY);
    if (n > 8)
        u256_op(digit(r, 1), digit(a, 1), U256_MEMCOPY);
}

/* The top limb of BN254's p, the one prime of eight limbs here. */
enum
{
    Bn254PTop = 0x30644e72
};

/* r -= p when r >= p, for r below 2p. */
static inline __attribute__((always_inline)) void fp_reduce(fp *r, uint32_t n)
{
    /* A top limb below p's settles it without the trial subtraction, and one above makes r > p, so p comes off in
       place; only an equal top limb needs the trial. */
    uint32_t top = n == 8 ? Bn254PTop : C->p.l[n - 1];
    if (r->l[n - 1] < top)
        return;
    if (r->l[n - 1] > top)
    {
        uint32_t below = u256_op(digit(r, 0), &P0, U256_SUB);
        if (n > 8)
            u256_op(digit(r, 1), &P1, U256_SUB | (below ? U256_CARRY : 0));
        return;
    }
    u256_op(&T[0], digit(r, 0), U256_MEMCOPY);
    uint32_t borrow = u256_op(&T[0], &P0, U256_SUB);
    if (n > 8)
    {
        u256_op(&T[1], digit(r, 1), U256_MEMCOPY);
        borrow = u256_op(&T[1], &P1, U256_SUB | (borrow ? U256_CARRY : 0));
    }
    if (borrow)
        return;
    u256_op(digit(r, 0), &T[0], U256_MEMCOPY);
    if (n > 8)
        u256_op(digit(r, 1), &T[1], U256_MEMCOPY);
}

/* Both primes leave the top bit of their last limb clear, so a sum of two reduced elements
   never carries out of n limbs. Inline for the Fp2 operations, which pass a literal n. */
static inline __attribute__((always_inline)) void fp_add_n(fp *r, const fp *a, const fp *b, uint32_t n)
{
    if (r == b)
    {
        /* r += a; when a is r too, r doubles a copy of itself. */
        b = a;
        a = r;
        if (b == r)
        {
            fp_copy(&Scratch, r, n);
            b = &Scratch;
        }
    }
    fp_copy(r, a, n);
    uint32_t carry = u256_op(digit(r, 0), digit(b, 0), U256_ADD);
    if (n > 8)
        u256_op(digit(r, 1), digit(b, 1), U256_ADD | (carry ? U256_CARRY : 0));
    fp_reduce(r, n);
}

void fp_add(fp *r, const fp *a, const fp *b)
{
    fp_add_n(r, a, b, C->n);
}

static inline __attribute__((always_inline)) void fp_sub_n(fp *r, const fp *a, const fp *b, uint32_t n)
{
    if (r == b)
    {
        fp_copy(&Scratch, b, n);
        b = &Scratch;
    }
    fp_copy(r, a, n);
    uint32_t borrow = u256_op(digit(r, 0), digit(b, 0), U256_SUB);
    if (n > 8)
        borrow = u256_op(digit(r, 1), digit(b, 1), U256_SUB | (borrow ? U256_CARRY : 0));
    if (!borrow)
        return;
    /* Below zero: add p back, the wrap past 2^(256 * digits) undoing the borrow. The borrow is 1, U256_ADD. */
    uint32_t carry = u256_op(digit(r, 0), &P0, borrow);
    if (n > 8)
        u256_op(digit(r, 1), &P1, U256_ADD | (carry ? U256_CARRY : 0));
}

void fp_sub(fp *r, const fp *a, const fp *b)
{
    fp_sub_n(r, a, b, C->n);
}

/*
 * Montgomery product a*b/R mod p. BN254 (R = 2^256) takes the 512-bit product and one 256-bit
 * reduction step, all in place. BLS12-381 (R = 2^384) splits an element into a 256-bit digit and one
 * below 2^128: the product accumulates in T, three digits, and the reduction adds multiples of p
 * that clear T from the bottom, one 256-bit step and one of 128 bits.
 */
static void fp_mul_wide(fp *r, const fp *a, const fp *b);

static inline __attribute__((always_inline)) void fp_mul_n(fp *r, const fp *a, const fp *b, uint32_t n)
{
    const u256d *a0 = digit(a, 0), *b0 = digit(b, 0);
    if (n == 8)
    {
        /* With m = lo(ab) p^-1, m p has the low half of ab, so (ab - m p) / R = hi(ab) - hi(m p): a difference in
           (-p, p) that needs no carry from the low halves. hi(ab) is formed in r itself, unless r is b alone. */
        u256d *rd = digit(r, 0);
        if (rd == b0 && rd != a0)
        {
            const u256d *t = a0;
            a0 = b0;
            b0 = t;
        }
        if (b0 == rd)
        {
            u256_op(&T[1], b0, U256_MEMCOPY);
            b0 = &T[1];
        }
        u256_op(&T[0], a0, U256_MEMCOPY);
        u256_op(&T[0], b0, U256_MUL_LOW);
        u256_op(&T[0], &PInverse, U256_MUL_LOW);
        u256_op(&T[0], &P0, U256_MUL_HIGH);
        if (rd != a0)
            u256_op(rd, a0, U256_MEMCOPY);
        u256_op(rd, b0, U256_MUL_HIGH);
        /* A borrow leaves 1 in x12, which is U256_ADD. */
        uint32_t borrow = u256_op(rd, &T[0], U256_SUB);
        if (borrow)
            u256_op(rd, &P0, borrow);
        return;
    }
    fp_mul_wide(r, a, b);
}

void fp_mul(fp *r, const fp *a, const fp *b)
{
    fp_mul_n(r, a, b, C->n);
}

/* The BLS12-381 product of fp_mul. */
static NOINLINE void fp_mul_wide(fp *r, const fp *a, const fp *b)
{
    const u256d *a0 = digit(a, 0), *b0 = digit(b, 0);
    /* The product lands in three digits, T0 + T1 W + T2 W^2 (W = 2^256), each carry passed on as the next
       addition's carry in. a1 and b1 are below 2^128, so a1 b1 fits a digit and T2 never carries out. */
    const u256d *a1 = digit(a, 1), *b1 = digit(b, 1);
    u256d *t0 = &T[0], *t1 = &T[1], *t2 = &T[2];
    u256_op(t0, a0, U256_MEMCOPY);
    u256_op(t0, b0, U256_MUL_LOW);
    u256_op(t1, a0, U256_MEMCOPY);
    u256_op(t1, b0, U256_MUL_HIGH);
    u256_op(t2, a0, U256_MEMCOPY);
    u256_op(t2, b1, U256_MUL_HIGH);
    u256_op(&X, a0, U256_MEMCOPY);
    u256_op(&X, b1, U256_MUL_LOW);
    uint32_t c1 = u256_op(t1, &X, U256_ADD);
    u256_op(&X, a1, U256_MEMCOPY);
    u256_op(&X, b0, U256_MUL_LOW);
    uint32_t c2 = u256_op(t1, &X, U256_ADD);
    u256_op(&X, a1, U256_MEMCOPY);
    u256_op(&X, b0, U256_MUL_HIGH);
    u256_op(t2, &X, U256_ADD | (c1 ? U256_CARRY : 0));
    u256_op(&X, a1, U256_MEMCOPY);
    u256_op(&X, b1, U256_MUL_LOW);
    u256_op(t2, &X, U256_ADD | (c2 ? U256_CARRY : 0));

    /* A 256-bit step: m = T0 Pinv, made in T0 itself, and T0 + lo(m P0) is W exactly when T0 != 0. p = P0 + P1 W
       with P1 below 2^128. */
    uint32_t low = u256_op(t0, &Zero, U256_EQ) ? 0 : U256_CARRY;
    u256_op(t0, &Pinv, U256_MUL_LOW);
    u256_op(&X, t0, U256_MEMCOPY);
    u256_op(&X, &P0, U256_MUL_HIGH);
    c1 = u256_op(t1, &X, U256_ADD | low);
    u256_op(&X, t0, U256_MEMCOPY);
    u256_op(&X, &P1, U256_MUL_LOW);
    c2 = u256_op(t1, &X, U256_ADD);
    u256_op(t0, &P1, U256_MUL_HIGH);
    u256_op(t2, t0, U256_ADD | (c1 ? U256_CARRY : 0));
    if (c2)
        u256_op(t2, &Zero, U256_ADD | U256_CARRY);

    /* A 128-bit step: m = T1 Pinv mod 2^128 clears the low half of T1; m P1 is below 2^253. */
    u256_op(&M, t1, U256_MEMCOPY);
    u256_op(&M, &Pinv, U256_MUL_LOW);
    for (int i = 4; i < 8; i++)
        M.w[i] = 0;
    u256_op(&X, &M, U256_MEMCOPY);
    u256_op(&X, &P0, U256_MUL_LOW);
    c1 = u256_op(t1, &X, U256_ADD);
    u256_op(&X, &M, U256_MEMCOPY);
    u256_op(&X, &P0, U256_MUL_HIGH);
    u256_op(&M, &P1, U256_MUL_LOW);
    u256_op(t2, &X, U256_ADD | (c1 ? U256_CARRY : 0));
    u256_op(t2, &M, U256_ADD);

    /* The product over R = 2^384 is T from bit 384 on: words 4-7 of T1, then T2; below 2p. */
    for (int i = 0; i < 4; i++)
        r->l[i] = t1->w[4 + i];
    for (int i = 0; i < 8; i++)
        r->l[4 + i] = t2->w[i];
    for (int i = 12; i < FpMaxLimbs; i++)
        r->l[i] = 0;
    fp_reduce(r, 12);
}

bool fp_is_zero(const fp *a)
{
    uint32_t acc = 0;
    for (uint32_t i = 0; i < C->n; i++)
        acc |= a->l[i];
    return acc == 0;
}

/* Reads a big-endian element of 4n bytes; false when it is not below p. */
bool fp_from_be(fp *r, const uint8_t *in)
{
    uint32_t n = C->n;
    fp raw = {{0}};
    for (uint32_t i = 0; i < n; i++)
    {
        const uint8_t *b = in + 4 * (n - 1 - i);
        raw.l[i] = ((uint32_t)b[0] << 24) | ((uint32_t)b[1] << 16) | ((uint32_t)b[2] << 8) | b[3];
    }
    if (limbs_cmp(raw.l, C->p.l, n) >= 0)
        return false;
    fp_mul(r, &raw, &C->r2);
    return true;
}

void fp_from_mont(fp *r, const fp *a)
{
    fp one = {{1}};
    fp_mul(r, a, &one);
}

void fp_to_be(uint8_t *out, const fp *a)
{
    uint32_t n = C->n;
    fp t;
    fp_from_mont(&t, a);
    for (uint32_t i = 0; i < n; i++)
    {
        uint8_t *b = out + 4 * (n - 1 - i);
        b[0] = (uint8_t)(t.l[i] >> 24);
        b[1] = (uint8_t)(t.l[i] >> 16);
        b[2] = (uint8_t)(t.l[i] >> 8);
        b[3] = (uint8_t)t.l[i];
    }
}

/* (p + add) >> shift, the exponents of inversion and square roots. */
static void p_exp(fp *e, int32_t add, unsigned shift)
{
    uint32_t n = C->n;
    int64_t carry = add;
    for (uint32_t i = 0; i < n; i++)
    {
        carry += C->p.l[i];
        e->l[i] = (uint32_t)carry;
        carry >>= 32;
    }
    for (uint32_t i = 0; shift && i < n; i++)
        e->l[i] = (e->l[i] >> shift) | (i + 1 < n ? e->l[i + 1] << (32 - shift) : 0);
}

/* ---- Fp and Fp2 as one degree-d field ------------------------------------------------------- */

NOINLINE void fq_add(fp2 *r, const fp2 *a, const fp2 *b, int d)
{
    if (C->n == 8)
    {
        fp_add_n(&r->c[0], &a->c[0], &b->c[0], 8);
        if (d == 2)
            fp_add_n(&r->c[1], &a->c[1], &b->c[1], 8);
        return;
    }
    fp_add_n(&r->c[0], &a->c[0], &b->c[0], 12);
    if (d == 2)
        fp_add_n(&r->c[1], &a->c[1], &b->c[1], 12);
}

NOINLINE void fq_sub(fp2 *r, const fp2 *a, const fp2 *b, int d)
{
    if (C->n == 8)
    {
        fp_sub_n(&r->c[0], &a->c[0], &b->c[0], 8);
        if (d == 2)
            fp_sub_n(&r->c[1], &a->c[1], &b->c[1], 8);
        return;
    }
    fp_sub_n(&r->c[0], &a->c[0], &b->c[0], 12);
    if (d == 2)
        fp_sub_n(&r->c[1], &a->c[1], &b->c[1], 12);
}

NOINLINE void fq_neg(fp2 *r, const fp2 *a, int d)
{
    fp2 zero;
    memset(&zero, 0, sizeof(zero));
    fq_sub(r, &zero, a, d);
}

/*
 * The BN254 sum r = a + b for a product's factor, r distinct from a and b, left below 2p unreduced: with 4p < 2^256,
 * factors below 2p keep hi(ab) below p, all fp_mul_n needs.
 */
static inline __attribute__((always_inline)) void fp_add_factor8(fp *r, const fp *a, const fp *b)
{
    u256_op(digit(r, 0), digit(a, 0), U256_MEMCOPY);
    u256_op(digit(r, 0), digit(b, 0), U256_ADD);
}

/* The BN254 product of fp_mul_n for r distinct from a and b, which needs none of its operand checks. */
static inline __attribute__((always_inline)) void fp_mul_new8(fp *r, const fp *a, const fp *b)
{
    u256d *rd = digit(r, 0);
    u256_op(&T[0], digit(a, 0), U256_MEMCOPY);
    u256_op(&T[0], digit(b, 0), U256_MUL_LOW);
    u256_op(&T[0], &PInverse, U256_MUL_LOW);
    u256_op(&T[0], &P0, U256_MUL_HIGH);
    u256_op(rd, digit(a, 0), U256_MEMCOPY);
    u256_op(rd, digit(b, 0), U256_MUL_HIGH);
    uint32_t borrow = u256_op(rd, &T[0], U256_SUB);
    if (borrow)
        u256_op(rd, &P0, borrow);
}

/* The BN254 difference of fp_sub_n for r distinct from a and b. */
static inline __attribute__((always_inline)) void fp_sub_new8(fp *r, const fp *a, const fp *b)
{
    u256d *rd = digit(r, 0);
    u256_op(rd, digit(a, 0), U256_MEMCOPY);
    uint32_t borrow = u256_op(rd, digit(b, 0), U256_SUB);
    if (borrow)
        u256_op(rd, &P0, borrow);
}

/* Fp2 by Karatsuba: (a0 + a1 u)(b0 + b1 u) = a0 b0 - a1 b1 + ((a0 + a1)(b0 + b1) - a0 b0 - a1 b1) u. */
static inline __attribute__((always_inline)) void fq_mul_n(fp2 *r, const fp2 *a, const fp2 *b, uint32_t n)
{
    fp t0, t1, t2, t3;
    if (n == 8)
    {
        /* The temporaries are distinct from a, b and r, and r is written once nothing reads a or b. */
        if (a == b)
        {
            fp_add_factor8(&t0, &a->c[0], &a->c[1]);
            fp_sub_new8(&t1, &a->c[0], &a->c[1]);
            fp_mul_new8(&t2, &a->c[0], &a->c[1]);
            fp_mul_new8(&r->c[0], &t0, &t1);
            fp_add_n(&r->c[1], &t2, &t2, 8);
            return;
        }
        fp_mul_new8(&t0, &a->c[0], &b->c[0]);
        fp_mul_new8(&t1, &a->c[1], &b->c[1]);
        fp_add_factor8(&t2, &a->c[0], &a->c[1]);
        fp_add_factor8(&t3, &b->c[0], &b->c[1]);
        fp_mul_n(&t2, &t2, &t3, 8);
        fp_sub_new8(&r->c[0], &t0, &t1);
        fp_sub_n(&t2, &t2, &t0, 8);
        fp_sub_new8(&r->c[1], &t2, &t1);
        return;
    }
    if (a == b)
    {
        /* A square by the complex method, u^2 = -1: (a0 + a1)(a0 - a1) + 2 a0 a1 u, two products for three. */
        fp_add_n(&t0, &a->c[0], &a->c[1], n);
        fp_sub_n(&t1, &a->c[0], &a->c[1], n);
        fp_mul_n(&t2, &a->c[0], &a->c[1], n);
        fp_mul_n(&r->c[0], &t0, &t1, n);
        fp_add_n(&r->c[1], &t2, &t2, n);
        return;
    }
    fp_mul_n(&t0, &a->c[0], &b->c[0], n);
    fp_mul_n(&t1, &a->c[1], &b->c[1], n);
    fp_add_n(&t2, &a->c[0], &a->c[1], n);
    fp_add_n(&t3, &b->c[0], &b->c[1], n);
    fp_mul_n(&t2, &t2, &t3, n);
    fp_sub_n(&r->c[0], &t0, &t1, n);
    fp_sub_n(&t2, &t2, &t0, n);
    fp_sub_n(&r->c[1], &t2, &t1, n);
}

/* Each limb count in a function of its own, which saves only the registers that count's product needs. */
static NOINLINE void fq_mul8(fp2 *r, const fp2 *a, const fp2 *b)
{
    fq_mul_n(r, a, b, 8);
}

static NOINLINE void fq_mul12(fp2 *r, const fp2 *a, const fp2 *b)
{
    fq_mul_n(r, a, b, 12);
}

/* The field operations inline, each with the limb count as a constant. */
NOINLINE void fq_mul(fp2 *r, const fp2 *a, const fp2 *b, int d)
{
    if (d == 1)
        fp_mul(&r->c[0], &a->c[0], &b->c[0]);
    else if (C->n == 8)
        fq_mul8(r, a, b);
    else
        fq_mul12(r, a, b);
}

NOINLINE bool fq_is_zero(const fp2 *a, int d)
{
    for (int i = 0; i < d; i++)
        if (!fp_is_zero(&a->c[i]))
            return false;
    return true;
}

NOINLINE bool fq_eq(const fp2 *a, const fp2 *b, int d)
{
    fp2 t;
    fq_sub(&t, a, b, d);
    return fq_is_zero(&t, d);
}

static void fq_set_one(fp2 *r)
{
    memset(r, 0, sizeof(*r));
    r->c[0] = C->one;
}

/*
 * The degree-d field operations of the point formulas and of fq_pow. With n = 8, BN254's, they run inline, as a call
 * through fq_* and its degree and limb-count dispatch costs about three times the product it makes; with n = 0 they
 * take the out-of-line fq_*, which keeps the code small. Elements are reduced, so equality compares limbs.
 */
#define PT_OP static inline __attribute__((always_inline))

PT_OP void pt_fmul(fp2 *r, const fp2 *a, const fp2 *b, int d, uint32_t n)
{
    if (n != 8)
        fq_mul(r, a, b, d);
    else if (d == 1)
        fp_mul_n(&r->c[0], &a->c[0], &b->c[0], 8);
    else
        fq_mul8(r, a, b);
}

PT_OP void pt_fadd(fp2 *r, const fp2 *a, const fp2 *b, int d, uint32_t n)
{
    if (n != 8)
    {
        fq_add(r, a, b, d);
        return;
    }
    for (int i = 0; i < d; i++)
        fp_add_n(&r->c[i], &a->c[i], &b->c[i], 8);
}

PT_OP void pt_fsub(fp2 *r, const fp2 *a, const fp2 *b, int d, uint32_t n)
{
    if (n != 8)
    {
        fq_sub(r, a, b, d);
        return;
    }
    for (int i = 0; i < d; i++)
        fp_sub_n(&r->c[i], &a->c[i], &b->c[i], 8);
}

PT_OP bool pt_fzero(const fp2 *a, int d, uint32_t n)
{
    if (n != 8)
        return fq_is_zero(a, d);
    uint32_t acc = 0;
    for (int i = 0; i < d; i++)
        for (int j = 0; j < 8; j++)
            acc |= a->c[i].l[j];
    return acc == 0;
}

PT_OP bool pt_feq(const fp2 *a, const fp2 *b, int d, uint32_t n)
{
    if (n != 8)
        return fq_eq(a, b, d);
    uint32_t acc = 0;
    for (int i = 0; i < d; i++)
        for (int j = 0; j < 8; j++)
            acc |= a->c[i].l[j] ^ b->c[i].l[j];
    return acc == 0;
}

PT_OP bool pt_fone(const fp2 *a, int d, uint32_t n)
{
    uint32_t limbs = n ? n : C->n, acc = 0;
    for (uint32_t i = 0; i < limbs; i++)
        acc |= (a->c[0].l[i] ^ C->one.l[i]) | (d == 2 ? a->c[1].l[i] : 0);
    return acc == 0;
}

/* a^e by 4-bit windows from the top non-zero one: 14 products for the table, then one per window. */
PT_OP void fq_pow_in(fp2 *r, const fp2 *a, const fp *e, int d, uint32_t n)
{
    fp2 table[16], acc;
    fq_set_one(&table[0]);
    table[1] = *a;
    for (int i = 2; i < 16; i++)
        pt_fmul(&table[i], &table[i - 1], a, d, n);

    int window = 8 * (int)C->n - 1;
    while (window > 0 && ((e->l[window / 8] >> (4 * (window % 8))) & 15) == 0)
        window--;
    acc = table[(e->l[window / 8] >> (4 * (window % 8))) & 15];
    for (window--; window >= 0; window--)
    {
        for (int s = 0; s < 4; s++)
            pt_fmul(&acc, &acc, &acc, d, n);
        uint32_t bits = (e->l[window / 8] >> (4 * (window % 8))) & 15;
        if (bits)
            pt_fmul(&acc, &acc, &table[bits], d, n);
    }
    *r = acc;
}

NOINLINE static void fq_pow(fp2 *r, const fp2 *a, const fp *e, int d)
{
    if (d == 1 && C->n == 8)
        fq_pow_in(r, a, e, 1, 8);
    else
        fq_pow_in(r, a, e, d, 0);
}

/* Inverse by Fermat in Fp; in Fp2 through the norm, 1/(a0 + a1 u) = (a0 - a1 u)/(a0^2 + a1^2). */
NOINLINE void fq_inv(fp2 *r, const fp2 *a, int d)
{
    fp e;
    p_exp(&e, -2, 0);
    if (d == 1)
    {
        fq_pow(r, a, &e, 1);
        return;
    }
    fp2 t;
    fp s;
    fp_mul(&t.c[0], &a->c[0], &a->c[0]);
    fp_mul(&s, &a->c[1], &a->c[1]);
    fp_add(&t.c[0], &t.c[0], &s);
    fq_pow(&t, &t, &e, 1);
    fp_mul(&r->c[0], &a->c[0], &t.c[0]);
    fp_mul(&s, &a->c[1], &t.c[0]);
    memset(&t, 0, sizeof(t));
    fp_sub(&r->c[1], &t.c[0], &s);
}

NOINLINE void fp2_conj(fp2 *r, const fp2 *a)
{
    fp zero = {{0}};
    r->c[0] = a->c[0];
    fp_sub(&r->c[1], &zero, &a->c[1]);
}

/* 1/2 in Montgomery form, formed once per curve: half of R mod p, which is C->one. */
static const fp *fp_half(void)
{
    static fp half;
    static const curve *formed_for;
    if (formed_for != C)
    {
        uint32_t n = C->n;
        fp t = C->one;
        uint64_t carry = 0;
        if (t.l[0] & 1)
        {
            for (uint32_t i = 0; i < n; i++)
            {
                carry += (uint64_t)t.l[i] + C->p.l[i];
                t.l[i] = (uint32_t)carry;
                carry >>= 32;
            }
        }
        for (uint32_t i = 0; i < n; i++)
            t.l[i] = (t.l[i] >> 1) | (i + 1 < n ? t.l[i + 1] << 31 : (uint32_t)carry << 31);
        half = t;
        formed_for = C;
    }
    return &half;
}

/* r = a^((p+1)/4) in Fp, held in an Fp2's first coefficient; true when it squares back to a. */
static bool fp_sqrt_in(fp2 *r, const fp2 *a)
{
    fp e;
    fp2 t;
    p_exp(&e, 1, 2);
    fq_pow(r, a, &e, 1);
    fq_mul(&t, r, r, 1);
    return fq_eq(&t, a, 1);
}

/*
 * A square root, false when there is none. Both primes are 3 mod 4: in Fp the root is
 * a^((p+1)/4). In Fp2 = Fp(u), u^2 = -1, a = a0 + a1 u with a1 != 0 is a square exactly when its
 * norm a0^2 + a1^2 is one in Fp, with root alpha; then one of delta = (a0 +- alpha) / 2 is a square
 * in Fp, and with t = delta^((p-3)/4), gamma = t delta is its root and 1/gamma is t, so the root of a
 * is gamma + (a1 t / 2) u. Two or three exponentiations in Fp, and one for a non-square, where the
 * powers in Fp2 took two. Either candidate is checked by squaring.
 */
NOINLINE bool fq_sqrt(fp2 *r, const fp2 *a, int d)
{
    fp2 x, t;
    if (d == 1)
    {
        fp_sqrt_in(&x, a);
    }
    else if (fp_is_zero(&a->c[1]))
    {
        /* a0 or -a0 is a square in Fp, as -1 is not: the root is sqrt(a0) or sqrt(-a0) u. */
        fp2 b;
        memset(&b, 0, sizeof(b));
        b.c[0] = a->c[0];
        memset(&x, 0, sizeof(x));
        if (!fp_sqrt_in(&t, &b))
        {
            fq_neg(&b, &b, 1);
            fp_sqrt_in(&t, &b);
            x.c[1] = t.c[0];
        }
        else
            x.c[0] = t.c[0];
    }
    else
    {
        fp2 n, alpha, delta, power, gamma;
        fp s, e;
        memset(&n, 0, sizeof(n));
        fp_mul(&n.c[0], &a->c[0], &a->c[0]);
        fp_mul(&s, &a->c[1], &a->c[1]);
        fp_add(&n.c[0], &n.c[0], &s);
        if (!fp_sqrt_in(&alpha, &n))
        {
            *r = *a;
            return false;
        }
        p_exp(&e, -3, 2);
        memset(&delta, 0, sizeof(delta));
        fp_add(&delta.c[0], &a->c[0], &alpha.c[0]);
        fp_mul(&delta.c[0], &delta.c[0], fp_half());
        fq_pow(&power, &delta, &e, 1);
        fq_mul(&gamma, &power, &delta, 1);
        fq_mul(&t, &gamma, &gamma, 1);
        if (!fq_eq(&t, &delta, 1))
        {
            fp_sub(&delta.c[0], &a->c[0], &alpha.c[0]);
            fp_mul(&delta.c[0], &delta.c[0], fp_half());
            fq_pow(&power, &delta, &e, 1);
            fq_mul(&gamma, &power, &delta, 1);
        }
        x.c[0] = gamma.c[0];
        fp_mul(&x.c[1], &a->c[1], &power.c[0]);
        fp_mul(&x.c[1], &x.c[1], fp_half());
    }
    fq_mul(&t, &x, &x, d);
    bool ok = fq_eq(&t, a, d);
    *r = x;
    return ok;
}

/* (a0 k - a1) + (a1 k + a0) u, with ka the products by k, or a itself for k = 1. r may be a. */
static inline __attribute__((always_inline)) void fp2_mul_xi_n(fp2 *r, const fp2 *a, const fp2 *ka, uint32_t n)
{
    fp t;
    fp_sub_n(&t, &ka->c[0], &a->c[1], n);
    fp_add_n(&r->c[1], &ka->c[1], &a->c[0], n);
    fp_copy(&r->c[0], &t, n);
}

/* 9 in Montgomery form, for BN254's xi = 9 + u; in RAM, where the delegation reads its operands. */
static fp Nine = {{0x410d7ff7, 0xf60647ce, 0xd31bd011, 0x2f3d6f4d, 0x3940c6d1, 0x2943337e, 0xa7e39857, 0x1d9598e8}};

/* a * (9 + u) on BN254, with 9 a as two products by 9: an addition costs a third of a product, and 9 a by additions
   takes four per coefficient. r may be a. */
static inline __attribute__((always_inline)) void fp2_mul_xi8(fp2 *r, const fp2 *a)
{
    fp2 k;
    fp_mul_new8(&k.c[0], &a->c[0], &Nine);
    fp_mul_new8(&k.c[1], &a->c[1], &Nine);
    fp2_mul_xi_n(r, a, &k, 8);
}

/* a * xi with xi = k + u: (a0 k - a1) + (a0 + a1 k) u. r may be a. */
NOINLINE static void fp2_mul_xi(fp2 *r, const fp2 *a)
{
    if (C->xi == 1)
        fp2_mul_xi_n(r, a, a, 12);
    else
        fp2_mul_xi8(r, a);
}

/* ---- points --------------------------------------------------------------------------------- */

static void pt_set_inf(point *p)
{
    memset(p, 0, sizeof(*p));
}

bool pt_is_inf(const point *p, int d)
{
    return fq_is_zero(&p->z, d);
}

/* Doubling for a = 0 (dbl-2009-l). */
PT_OP void pt_dbl_in(point *r, const point *p, int d, uint32_t n)
{
    if (pt_fzero(&p->z, d, n) || pt_fzero(&p->y, d, n))
    {
        pt_set_inf(r);
        return;
    }
    fp2 xx, yy, yyyy, s, m, t;
    pt_fmul(&xx, &p->x, &p->x, d, n);
    pt_fmul(&yy, &p->y, &p->y, d, n);
    pt_fmul(&yyyy, &yy, &yy, d, n);

    pt_fmul(&s, &p->x, &yy, d, n);
    pt_fadd(&s, &s, &s, d, n);
    pt_fadd(&s, &s, &s, d, n);

    pt_fadd(&m, &xx, &xx, d, n);
    pt_fadd(&m, &m, &xx, d, n);

    /* Each coordinate is written once nothing reads its input any more, so r may be p with no copies. */
    pt_fmul(&r->z, &p->y, &p->z, d, n);
    pt_fadd(&r->z, &r->z, &r->z, d, n);

    pt_fmul(&r->x, &m, &m, d, n);
    pt_fsub(&r->x, &r->x, &s, d, n);
    pt_fsub(&r->x, &r->x, &s, d, n);

    pt_fsub(&t, &s, &r->x, d, n);
    pt_fmul(&r->y, &m, &t, d, n);
    pt_fadd(&t, &yyyy, &yyyy, d, n);
    pt_fadd(&t, &t, &t, d, n);
    pt_fadd(&t, &t, &t, d, n);
    pt_fsub(&r->y, &r->y, &t, d, n);
}

NOINLINE void pt_dbl(point *r, const point *p, int d)
{
    if (C->n != 8)
        pt_dbl_in(r, p, d, 0);
    else if (d == 1)
        pt_dbl_in(r, p, 1, 8);
    else
        pt_dbl_in(r, p, 2, 8);
}

/* General addition (add-2007-bl without the doubling shortcut), falling back to doubling. */
PT_OP void pt_add_in(point *r, const point *p, const point *q, int d, uint32_t n)
{
    if (pt_fzero(&p->z, d, n))
    {
        *r = *q;
        return;
    }
    if (pt_fzero(&q->z, d, n))
    {
        *r = *p;
        return;
    }
    fp2 z1z1, z2z2, u1v, u2, s1v, s2, h, rr, hh, hhh, v, x3, y3, z3, t;
    /* An affine q (z = 1, as loaded) leaves p's coordinates as they are: five products fewer. */
    bool q_affine = pt_fone(&q->z, d, n);
    const fp2 *u1 = &p->x, *s1 = &p->y;
    pt_fmul(&z1z1, &p->z, &p->z, d, n);
    if (!q_affine)
    {
        pt_fmul(&z2z2, &q->z, &q->z, d, n);
        pt_fmul(&u1v, &p->x, &z2z2, d, n);
        pt_fmul(&t, &q->z, &z2z2, d, n);
        pt_fmul(&s1v, &p->y, &t, d, n);
        u1 = &u1v;
        s1 = &s1v;
    }
    pt_fmul(&u2, &q->x, &z1z1, d, n);
    pt_fmul(&t, &p->z, &z1z1, d, n);
    pt_fmul(&s2, &q->y, &t, d, n);

    if (pt_feq(u1, &u2, d, n))
    {
        if (pt_feq(s1, &s2, d, n))
            pt_dbl(r, p, d);
        else
            pt_set_inf(r);
        return;
    }

    pt_fsub(&h, &u2, u1, d, n);
    pt_fsub(&rr, &s2, s1, d, n);
    pt_fmul(&hh, &h, &h, d, n);
    pt_fmul(&hhh, &hh, &h, d, n);
    pt_fmul(&v, u1, &hh, d, n);

    pt_fmul(&x3, &rr, &rr, d, n);
    pt_fsub(&x3, &x3, &hhh, d, n);
    pt_fsub(&x3, &x3, &v, d, n);
    pt_fsub(&x3, &x3, &v, d, n);

    pt_fsub(&t, &v, &x3, d, n);
    pt_fmul(&y3, &rr, &t, d, n);
    pt_fmul(&t, s1, &hhh, d, n);
    pt_fsub(&y3, &y3, &t, d, n);

    if (q_affine)
        pt_fmul(&z3, &p->z, &h, d, n);
    else
    {
        pt_fmul(&z3, &p->z, &q->z, d, n);
        pt_fmul(&z3, &z3, &h, d, n);
    }

    r->x = x3;
    r->y = y3;
    r->z = z3;
}

NOINLINE void pt_add(point *r, const point *p, const point *q, int d)
{
    if (C->n != 8)
        pt_add_in(r, p, q, d, 0);
    else if (d == 1)
        pt_add_in(r, p, q, 1, 8);
    else
        pt_add_in(r, p, q, 2, 8);
}

void pt_neg(point *r, const point *p, int d)
{
    *r = *p;
    fq_neg(&r->y, &p->y, d);
}

static inline uint32_t scalar_bit(const uint8_t *k, size_t klen, size_t i)
{
    return i < 8 * klen ? (k[klen - 1 - i / 8] >> (i % 8)) & 1 : 0;
}

/*
 * The digits of the big-endian k, least significant first, for a double-and-add loop: its non-adjacent form (-1, 0
 * or 1, no two adjacent non-zero) when that saves additions net of the extra step a longer form takes, else its
 * bits. Returns the digit count, at most 8 klen + 1, the top digit 1; none for k = 0.
 */
static int scalar_digits(int8_t *dg, const uint8_t *k, size_t klen)
{
    int naf_len = 0, naf_weight = 0, bin_len = 0, bin_weight = 0;
    uint32_t carry = 0;
    for (size_t i = 0; i < 8 * klen || carry; i++)
    {
        uint32_t bit = scalar_bit(k, klen, i), v = bit + carry;
        int8_t z = 0;
        if (v == 1)
        {
            carry = scalar_bit(k, klen, i + 1);
            z = carry ? -1 : 1;
            naf_len = (int)i + 1;
            naf_weight++;
        }
        else
            carry = v >> 1;
        dg[i] = z;
        if (bit)
        {
            bin_len = (int)i + 1;
            bin_weight++;
        }
    }
    if (naf_weight + (naf_len - bin_len) < bin_weight)
        return naf_len;
    for (int i = 0; i < bin_len; i++)
        dg[i] = (int8_t)scalar_bit(k, klen, (size_t)i);
    return bin_len;
}

/* |x|, the curve family parameter, big-endian. */
static void x_bytes(uint8_t k[8])
{
    for (int i = 0; i < 8; i++)
        k[i] = (uint8_t)(C->x >> (56 - 8 * i));
}

enum
{
    MaxScalarBytes = 32
};

/* k*p for a big-endian scalar of at most MaxScalarBytes, by double-and-add on scalar_digits. */
void pt_mul(point *r, const point *p, const uint8_t *k, size_t klen, int d)
{
    int8_t dg[8 * MaxScalarBytes + 1];
    int n = scalar_digits(dg, k, klen);
    if (n == 0)
    {
        pt_set_inf(r);
        return;
    }
    point base = *p, neg, acc = *p;
    pt_neg(&neg, p, d);
    for (int i = n - 2; i >= 0; i--)
    {
        pt_dbl(&acc, &acc, d);
        if (dg[i])
            pt_add(&acc, &acc, dg[i] > 0 ? &base : &neg, d);
    }
    *r = acc;
}

/* kp*p + kq*q for big-endian scalars of klen bytes each, at most MaxScalarBytes, by one double-and-add over both. */
void pt_mul2(point *r, const point *p, const uint8_t *kp, const point *q, const uint8_t *kq, size_t klen, int d)
{
    int8_t dp[8 * MaxScalarBytes + 1], dq[8 * MaxScalarBytes + 1];
    int np = scalar_digits(dp, kp, klen), nq = scalar_digits(dq, kq, klen);
    point negp, negq, acc;
    pt_neg(&negp, p, d);
    pt_neg(&negq, q, d);
    pt_set_inf(&acc);
    for (int i = (np > nq ? np : nq) - 1; i >= 0; i--)
    {
        pt_dbl(&acc, &acc, d);
        if (i < np && dp[i])
            pt_add(&acc, &acc, dp[i] > 0 ? p : &negp, d);
        if (i < nq && dq[i])
            pt_add(&acc, &acc, dq[i] > 0 ? q : &negq, d);
    }
    *r = acc;
}

/* x*p, x the signed curve parameter. */
void pt_mul_x(point *r, const point *p, int d)
{
    uint8_t k[8];
    x_bytes(k);
    pt_mul(r, p, k, sizeof(k), d);
    if (C->xneg)
        pt_neg(r, r, d);
}

bool pt_eq(const point *p, const point *q, int d)
{
    if (pt_is_inf(p, d) || pt_is_inf(q, d))
        return pt_is_inf(p, d) && pt_is_inf(q, d);
    fp2 z1z1, z2z2, a, b;
    fq_mul(&z1z1, &p->z, &p->z, d);
    fq_mul(&z2z2, &q->z, &q->z, d);
    fq_mul(&a, &p->x, &z2z2, d);
    fq_mul(&b, &q->x, &z1z1, d);
    if (!fq_eq(&a, &b, d))
        return false;
    fq_mul(&z1z1, &z1z1, &p->z, d);
    fq_mul(&z2z2, &z2z2, &q->z, d);
    fq_mul(&a, &p->y, &z2z2, d);
    fq_mul(&b, &q->y, &z1z1, d);
    return fq_eq(&a, &b, d);
}

void pt_set_affine(point *p, int d)
{
    if (pt_is_inf(p, d))
        return;
    fp2 one;
    fq_set_one(&one);
    if (fq_eq(&p->z, &one, d))
        return;
    fp2 zi, zi2;
    fq_inv(&zi, &p->z, d);
    fq_mul(&zi2, &zi, &zi, d);
    fq_mul(&p->x, &p->x, &zi2, d);
    fq_mul(&zi2, &zi2, &zi, d);
    fq_mul(&p->y, &p->y, &zi2, d);
    fq_set_one(&p->z);
}

/* For an affine point: y^2 == x^3 + b, on G1's curve (d = 1) or the twist (d = 2). */
bool pt_on_curve(const point *p, int d)
{
    if (pt_is_inf(p, d))
        return true;
    fp2 lhs, rhs;
    fq_mul(&lhs, &p->y, &p->y, d);
    fq_mul(&rhs, &p->x, &p->x, d);
    fq_mul(&rhs, &rhs, &p->x, d);
    fq_add(&rhs, &rhs, d == 1 ? &C->b1 : &C->b2, d);
    return fq_eq(&lhs, &rhs, d);
}

/*
 * Reads an affine point x || y of 2d big-endian field elements, all zero meaning infinity. An
 * Fp2 element is c0 || c1, or c1 || c0 with swap (the EIP-197 order). False when an element is
 * not below p; the curve equation is left to pt_on_curve.
 */
bool pt_load(point *p, const uint8_t *in, int d, bool swap)
{
    uint32_t fb = 4 * C->n;
    memset(p, 0, sizeof(*p));
    fp2 *coord[2] = {&p->x, &p->y};
    for (int i = 0; i < 2 * d; i++)
    {
        int c = d == 2 && swap ? 1 - i % 2 : i % 2;
        if (!fp_from_be(&coord[i / d]->c[d == 1 ? 0 : c], in + i * fb))
            return false;
    }
    if (!fq_is_zero(&p->x, d) || !fq_is_zero(&p->y, d))
        fq_set_one(&p->z);
    return true;
}

void pt_store(uint8_t *out, const point *p, int d, bool swap)
{
    uint32_t fb = 4 * C->n;
    point a = *p;
    pt_set_affine(&a, d);
    memset(out, 0, 2 * d * fb);
    if (pt_is_inf(&a, d))
        return;
    const fp2 *coord[2] = {&a.x, &a.y};
    for (int i = 0; i < 2 * d; i++)
    {
        int c = d == 2 && swap ? 1 - i % 2 : i % 2;
        fp_to_be(out + i * fb, &coord[i / d]->c[d == 1 ? 0 : c]);
    }
}

/*
 * The p-power Frobenius carried to the twist: untwisting, raising coordinates to the p-th power
 * and twisting back gives (conj(x) g^2, conj(y) g^3) with g = gamma on a D-type twist and
 * 1/gamma on an M-type one. Jacobian coordinates take the same map with z conjugated.
 */
void pt_frob(point *r, const point *p, bool inverse)
{
    /* g^2 and g^3 per curve and direction, formed once: 1/gamma is an exponentiation. */
    static struct
    {
        const curve *cv;
        fp2 g2, g3;
    } factors[2][2];
    __typeof__(factors[0][0]) *k = &factors[C->mtwist][inverse];
    if (k->cv != C)
    {
        fp2 g = C->gamma;
        if (inverse)
            fq_inv(&g, &g, 2);
        fq_mul(&k->g2, &g, &g, 2);
        fq_mul(&k->g3, &k->g2, &g, 2);
        k->cv = C;
    }
    fp2_conj(&r->x, &p->x);
    fp2_conj(&r->y, &p->y);
    fp2_conj(&r->z, &p->z);
    fq_mul(&r->x, &r->x, &k->g2, 2);
    fq_mul(&r->y, &r->y, &k->g3, 2);
}

/* ---- Fp12 ----------------------------------------------------------------------------------- */

/*
 * The tower and Miller loop functions on the hot path are written once for n limbs, and BY_LIMBS makes F8 and F12, an
 * instance for each curve with n a constant, and F, which calls the instance of the curve in use. BN254's (eight limbs)
 * Fp2 sums and products by xi then run inline, as a call costs about as much as the sum, and BLS12-381 calls.
 */
#define TOWER static inline __attribute__((always_inline))
#define BY_LIMBS(F, PARAMS, ...)                                                                                       \
    static NOINLINE void F##8 PARAMS                                                                                   \
    {                                                                                                                  \
        F##_n(__VA_ARGS__, 8);                                                                                         \
    }                                                                                                                  \
    static NOINLINE void F##12 PARAMS                                                                                  \
    {                                                                                                                  \
        F##_n(__VA_ARGS__, 12);                                                                                        \
    }                                                                                                                  \
    static inline void F PARAMS                                                                                        \
    {                                                                                                                  \
        if (C->n == 8)                                                                                                 \
            F##8(__VA_ARGS__);                                                                                         \
        else                                                                                                           \
            F##12(__VA_ARGS__);                                                                                        \
    }

TOWER void fp2_add_n(fp2 *r, const fp2 *a, const fp2 *b, uint32_t n)
{
    if (n != 8)
    {
        fq_add(r, a, b, 2);
        return;
    }
    fp_add_n(&r->c[0], &a->c[0], &b->c[0], 8);
    fp_add_n(&r->c[1], &a->c[1], &b->c[1], 8);
}

TOWER void fp2_sub_n(fp2 *r, const fp2 *a, const fp2 *b, uint32_t n)
{
    if (n != 8)
    {
        fq_sub(r, a, b, 2);
        return;
    }
    fp_sub_n(&r->c[0], &a->c[0], &b->c[0], 8);
    fp_sub_n(&r->c[1], &a->c[1], &b->c[1], 8);
}

TOWER void fp2_mul_n(fp2 *r, const fp2 *a, const fp2 *b, uint32_t n)
{
    if (n == 8)
        fq_mul8(r, a, b);
    else
        fq_mul12(r, a, b);
}

/* a * xi; r may be a. */
TOWER void fp2_xi_n(fp2 *r, const fp2 *a, uint32_t n)
{
    if (n == 8)
        fp2_mul_xi8(r, a);
    else
        fp2_mul_xi(r, a);
}

/* r = a b for r distinct from a and b. */
TOWER void fp_mul_new_n(fp *r, const fp *a, const fp *b, uint32_t n)
{
    if (n == 8)
        fp_mul_new8(r, a, b);
    else
        fp_mul(r, a, b);
}

void fp12_one(fp12 *r)
{
    memset(r, 0, sizeof(*r));
    r->c[0].c[0] = C->one;
}

/*
 * Fp12 as Fp6[w]/(w^2 - v) over Fp6 = Fp2[v]/(v^3 - xi), v = w^2: the even coefficients of the flat
 * form are the Fp6 element A and the odd ones B, for A + B w. Karatsuba at both levels takes 18 Fp2
 * products for a product and 12 for a square, against 36 and 21 for the flat schoolbook. An fp12
 * stores A then B, so the coefficient of w^i is at W(i), and the halves are read in place.
 */
typedef struct
{
    fp2 c[3]; /* c[0] + c[1] v + c[2] v^2 */
} fp6;

#define W(i) ((i) % 2 ? 3 + (i) / 2 : (i) / 2)

static inline fp6 *fp12_even(const fp12 *a)
{
    return (fp6 *)(uintptr_t)&a->c[0];
}

static inline fp6 *fp12_odd(const fp12 *a)
{
    return (fp6 *)(uintptr_t)&a->c[3];
}

TOWER void fp6_add_n(fp6 *r, const fp6 *a, const fp6 *b, uint32_t n)
{
    for (int i = 0; i < 3; i++)
        fp2_add_n(&r->c[i], &a->c[i], &b->c[i], n);
}

TOWER void fp6_sub_n(fp6 *r, const fp6 *a, const fp6 *b, uint32_t n)
{
    for (int i = 0; i < 3; i++)
        fp2_sub_n(&r->c[i], &a->c[i], &b->c[i], n);
}

/* r = a v: (a0 + a1 v + a2 v^2) v = xi a2 + a0 v + a1 v^2. */
TOWER void fp6_mul_v_n(fp6 *r, const fp6 *a, uint32_t n)
{
    fp2 t;
    fp2_xi_n(&t, &a->c[2], n);
    r->c[2] = a->c[1];
    r->c[1] = a->c[0];
    r->c[0] = t;
}

TOWER void fp6_mul_n(fp6 *r, const fp6 *a, const fp6 *b, uint32_t n)
{
    fp2 t0, t1, t2, s, u;
    fp6 out;
    fp2_mul_n(&t0, &a->c[0], &b->c[0], n);
    fp2_mul_n(&t1, &a->c[1], &b->c[1], n);
    fp2_mul_n(&t2, &a->c[2], &b->c[2], n);

    /* c0 = t0 + xi ((a1 + a2)(b1 + b2) - t1 - t2) */
    fp2_add_n(&s, &a->c[1], &a->c[2], n);
    fp2_add_n(&u, &b->c[1], &b->c[2], n);
    fp2_mul_n(&s, &s, &u, n);
    fp2_sub_n(&s, &s, &t1, n);
    fp2_sub_n(&s, &s, &t2, n);
    fp2_xi_n(&s, &s, n);
    fp2_add_n(&out.c[0], &t0, &s, n);

    /* c1 = (a0 + a1)(b0 + b1) - t0 - t1 + xi t2 */
    fp2_add_n(&s, &a->c[0], &a->c[1], n);
    fp2_add_n(&u, &b->c[0], &b->c[1], n);
    fp2_mul_n(&s, &s, &u, n);
    fp2_sub_n(&s, &s, &t0, n);
    fp2_sub_n(&s, &s, &t1, n);
    fp2_xi_n(&u, &t2, n);
    fp2_add_n(&out.c[1], &s, &u, n);

    /* c2 = (a0 + a2)(b0 + b2) - t0 - t2 + t1 */
    fp2_add_n(&s, &a->c[0], &a->c[2], n);
    fp2_add_n(&u, &b->c[0], &b->c[2], n);
    fp2_mul_n(&s, &s, &u, n);
    fp2_sub_n(&s, &s, &t0, n);
    fp2_sub_n(&s, &s, &t2, n);
    fp2_add_n(&out.c[2], &s, &t1, n);
    *r = out;
}
BY_LIMBS(fp6_mul, (fp6 *r, const fp6 *a, const fp6 *b), r, a, b)

TOWER void fp12_mul_n(fp12 *r, const fp12 *a, const fp12 *b, uint32_t n)
{
    const fp6 *x = fp12_even(a), *y = fp12_odd(a);
    fp6 s, t, u;
    if (a == b)
    {
        /* (A + B w)^2 = (A^2 + v B^2) + 2 A B w, where (A + B)(A + v B) = A^2 + v B^2 + (1 + v) A B. */
        fp6_mul(&u, x, y);
        fp6_mul_v_n(&t, y, n);
        fp6_add_n(&t, &t, x, n);
        fp6_add_n(&s, x, y, n);
        fp6_mul(&s, &s, &t);
        fp6_sub_n(&s, &s, &u, n);
        fp6_mul_v_n(&t, &u, n);
        fp6_sub_n(fp12_even(r), &s, &t, n);
        fp6_add_n(fp12_odd(r), &u, &u, n);
        return;
    }

    /* (A + B w)(C + D w) = (A C + v B D) + ((A + B)(C + D) - A C - B D) w */
    const fp6 *c = fp12_even(b), *d = fp12_odd(b);
    fp6 ac, bd;
    fp6_mul(&ac, x, c);
    fp6_mul(&bd, y, d);
    fp6_add_n(&s, x, y, n);
    fp6_add_n(&t, c, d, n);
    fp6_mul(&s, &s, &t);
    fp6_sub_n(&s, &s, &ac, n);
    fp6_sub_n(fp12_odd(r), &s, &bd, n);
    fp6_mul_v_n(&t, &bd, n);
    fp6_add_n(fp12_even(r), &t, &ac, n);
}
BY_LIMBS(fp12_mul, (fp12 *r, const fp12 *a, const fp12 *b), r, a, b)

/* The p^6-power Frobenius, w -> -w; it inverts elements of the cyclotomic subgroup. */
/* r = a (b0 + b1 v): 5 Fp2 products by Karatsuba. r may be a. */
TOWER void fp6_mul_01_n(fp6 *r, const fp6 *a, const fp2 *b0, const fp2 *b1, uint32_t n)
{
    fp2 t0, t1, c0, c1, c2;
    fp2_mul_n(&t0, &a->c[0], b0, n);
    fp2_mul_n(&t1, &a->c[1], b1, n);
    fp2_add_n(&c1, &a->c[0], &a->c[1], n);
    fp2_add_n(&c2, b0, b1, n);
    fp2_mul_n(&c1, &c1, &c2, n);
    fp2_sub_n(&c1, &c1, &t0, n);
    fp2_sub_n(&c1, &c1, &t1, n);
    fp2_mul_n(&c0, &a->c[2], b1, n);
    fp2_xi_n(&c0, &c0, n);
    fp2_add_n(&c0, &c0, &t0, n);
    fp2_mul_n(&c2, &a->c[2], b0, n);
    fp2_add_n(&c2, &c2, &t1, n);
    r->c[0] = c0;
    r->c[1] = c1;
    r->c[2] = c2;
}
BY_LIMBS(fp6_mul_01, (fp6 *r, const fp6 *a, const fp2 *b0, const fp2 *b1), r, a, b0, b1)

/* r = a b1 v: 3 Fp2 products. r may be a. */
TOWER void fp6_mul_1_n(fp6 *r, const fp6 *a, const fp2 *b1, uint32_t n)
{
    fp2 c0, c1, c2;
    fp2_mul_n(&c0, &a->c[2], b1, n);
    fp2_xi_n(&c0, &c0, n);
    fp2_mul_n(&c1, &a->c[0], b1, n);
    fp2_mul_n(&c2, &a->c[1], b1, n);
    r->c[0] = c0;
    r->c[1] = c1;
    r->c[2] = c2;
}
BY_LIMBS(fp6_mul_1, (fp6 *r, const fp6 *a, const fp2 *b1), r, a, b1)

/*
 * A line through points of the twist, evaluated at the G1 point p: cy y_P + cx x_P w + c0 w^3 on a D-type twist and
 * c0 + cx x_P w^2 + cy y_P w^3 on an M-type one, up to factors in Fp2 and powers of w, which the final exponentiation
 * sends to 1. y holds cy y_P and x holds cx x_P.
 */
typedef struct
{
    fp2 y, x, c0;
} line;

/*
 * f *= l. As L0 + L1 w over Fp6, one of L0 and L1 is a single Fp2 term and the other has two, so Karatsuba takes 13
 * Fp2 products where a dense product takes 18.
 */
TOWER void fp12_mul_line_n(fp12 *f, const line *l, uint32_t n)
{
    fp6 *x = fp12_even(f), *y = fp12_odd(f);
    fp6 a, b, s;
    fp2 sum;
    fp6_add_n(&s, x, y, n);
    if (C->mtwist)
    {
        /* L0 = c0 + (cx x_P) v, L1 = (cy y_P) v */
        fp6_mul_01(&a, x, &l->c0, &l->x);
        fp6_mul_1(&b, y, &l->y);
        fp2_add_n(&sum, &l->x, &l->y, n);
        fp6_mul_01(&s, &s, &l->c0, &sum);
    }
    else
    {
        /* L0 = cy y_P, L1 = cx x_P + c0 v */
        for (int i = 0; i < 3; i++)
            fp2_mul_n(&a.c[i], &x->c[i], &l->y, n);
        fp6_mul_01(&b, y, &l->x, &l->c0);
        fp2_add_n(&sum, &l->y, &l->x, n);
        fp6_mul_01(&s, &s, &sum, &l->c0);
    }
    /* (F0 + F1 w)(L0 + L1 w) = (F0 L0 + v F1 L1) + ((F0 + F1)(L0 + L1) - F0 L0 - F1 L1) w */
    fp6_sub_n(&s, &s, &a, n);
    fp6_sub_n(y, &s, &b, n);
    fp6_mul_v_n(&b, &b, n);
    fp6_add_n(x, &a, &b, n);
}
BY_LIMBS(fp12_mul_line, (fp12 *f, const line *l), f, l)

NOINLINE static void fp12_conj(fp12 *r, const fp12 *a)
{
    *r = *a;
    for (int i = 1; i < 6; i += 2)
        fq_neg(&r->c[W(i)], &a->c[W(i)], 2);
}

/* The p-power Frobenius: c_i w^i -> conj(c_i) gamma^i w^i. */
NOINLINE static void fp12_frob(fp12 *r, const fp12 *a)
{
    fp2 g;
    fq_set_one(&g);
    for (int i = 0; i < 6; i++)
    {
        fp2_conj(&r->c[W(i)], &a->c[W(i)]);
        fq_mul(&r->c[W(i)], &r->c[W(i)], &g, 2);
        fq_mul(&g, &g, &C->gamma, 2);
    }
}

/*
 * 1/a = conj(a) / (a conj(a)), where a conj(a) lies in Fp6 = Fp2[v], v = w^2, v^3 = xi; an Fp6
 * element b0 + b1 v + b2 v^2 inverts through its norm to Fp2.
 */
NOINLINE static void fp12_inv(fp12 *r, const fp12 *a)
{
    fp12 c, t;
    fp12_conj(&c, a);
    fp12_mul(&t, a, &c);
    fp2 b0 = t.c[W(0)], b1 = t.c[W(2)], b2 = t.c[W(4)], s0, s1, s2, n, m;

    fq_mul(&s0, &b0, &b0, 2);
    fq_mul(&m, &b1, &b2, 2);
    fp2_mul_xi(&m, &m);
    fq_sub(&s0, &s0, &m, 2);

    fq_mul(&s1, &b2, &b2, 2);
    fp2_mul_xi(&s1, &s1);
    fq_mul(&m, &b0, &b1, 2);
    fq_sub(&s1, &s1, &m, 2);

    fq_mul(&s2, &b1, &b1, 2);
    fq_mul(&m, &b0, &b2, 2);
    fq_sub(&s2, &s2, &m, 2);

    fq_mul(&n, &b2, &s1, 2);
    fq_mul(&m, &b1, &s2, 2);
    fq_add(&n, &n, &m, 2);
    fp2_mul_xi(&n, &n);
    fq_mul(&m, &b0, &s0, 2);
    fq_add(&n, &n, &m, 2);
    fq_inv(&n, &n, 2);

    memset(&t, 0, sizeof(t));
    fq_mul(&t.c[W(0)], &s0, &n, 2);
    fq_mul(&t.c[W(2)], &s1, &n, 2);
    fq_mul(&t.c[W(4)], &s2, &n, 2);
    fp12_mul(r, &c, &t);
}

/* (z + y v)^2 in Fp4 = Fp2[v]/(v^2 - xi), as t + s v: t = z^2 + xi y^2 and s = (z + y)^2 - z^2 - y^2, by squares. */
TOWER void fp4_sqr_n(fp2 *t, fp2 *s, const fp2 *z, const fp2 *y, uint32_t n)
{
    fp2 zz, yy;
    fp2_mul_n(&zz, z, z, n);
    fp2_mul_n(&yy, y, y, n);
    fp2_add_n(s, z, y, n);
    fp2_mul_n(s, s, s, n);
    fp2_sub_n(s, s, &zz, n);
    fp2_sub_n(s, s, &yy, n);
    fp2_xi_n(t, &yy, n);
    fp2_add_n(t, t, &zz, n);
}

/* *r = 3 t + 2 z when plus, else 3 t - 2 z; r may be z. */
TOWER void cyclotomic_part_n(fp2 *r, const fp2 *t, const fp2 *z, bool plus, uint32_t n)
{
    if (plus)
        fp2_add_n(r, t, z, n);
    else
        fp2_sub_n(r, t, z, n);
    fp2_add_n(r, r, r, n);
    fp2_add_n(r, r, t, n);
}

/*
 * r = a^2 for a in the cyclotomic subgroup, by Granger and Scott (ePrint 2009/565): three Fp4 squarings, 18 Fp
 * products where a general square takes 36. The arrangement follows arkworks' cyclotomic_square, whose Fp6 parts
 * are the even and the odd coefficients of w here. The squarings read every coefficient before any is written, and
 * each result takes only its own coefficient after them, so r may be a with no copies of either.
 */
TOWER void fp12_cyclotomic_sqr_n(fp12 *r, const fp12 *a, uint32_t n)
{
    const fp2 *z0 = &a->c[W(0)], *z4 = &a->c[W(2)], *z3 = &a->c[W(4)], *z2 = &a->c[W(1)], *z1 = &a->c[W(3)],
              *z5 = &a->c[W(5)];
    fp2 t0, t1, t2, t3, t4, t5, x;
    fp4_sqr_n(&t0, &t1, z0, z1, n);
    fp4_sqr_n(&t2, &t3, z2, z3, n);
    fp4_sqr_n(&t4, &t5, z4, z5, n);
    cyclotomic_part_n(&r->c[W(0)], &t0, z0, false, n);
    cyclotomic_part_n(&r->c[W(3)], &t1, z1, true, n);
    fp2_xi_n(&x, &t5, n);
    cyclotomic_part_n(&r->c[W(1)], &x, z2, true, n);
    cyclotomic_part_n(&r->c[W(4)], &t4, z3, false, n);
    cyclotomic_part_n(&r->c[W(2)], &t2, z4, false, n);
    cyclotomic_part_n(&r->c[W(5)], &t3, z5, true, n);
}
BY_LIMBS(fp12_cyclotomic_sqr, (fp12 *r, const fp12 *a), r, a)

/*
 * a^x for the signed curve parameter x; a must be in the cyclotomic subgroup, where the squares are cheaper and the
 * conjugate is the inverse that a negative digit multiplies by.
 */
NOINLINE static void fp12_pow_x(fp12 *r, const fp12 *a)
{
    uint8_t k[8];
    int8_t dg[8 * sizeof(k) + 1];
    x_bytes(k);
    fp12 acc = *a, base = *a, inv;
    fp12_conj(&inv, a);
    for (int i = scalar_digits(dg, k, sizeof(k)) - 2; i >= 0; i--)
    {
        fp12_cyclotomic_sqr(&acc, &acc);
        if (dg[i])
            fp12_mul(&acc, &acc, dg[i] > 0 ? &base : &inv);
    }
    if (C->xneg)
        fp12_conj(&acc, &acc);
    *r = acc;
}

/* ---- pairing -------------------------------------------------------------------------------- */

/* Evaluates the line with coefficients cy, cx and c0 at the G1 point p. */
TOWER void line_eval_n(line *l, const fp2 *cy, const fp2 *cx, const fp2 *c0, const point *p, uint32_t n)
{
    for (int i = 0; i < 2; i++)
    {
        fp_mul_new_n(&l->y.c[i], &cy->c[i], &p->y.c[0], n);
        fp_mul_new_n(&l->x.c[i], &cx->c[i], &p->x.c[0], n);
    }
    l->c0 = *c0;
}

static void fq_half(fp2 *r, const fp2 *a)
{
    fp_mul(&r->c[0], &a->c[0], fp_half());
    fp_mul(&r->c[1], &a->c[1], fp_half());
}

/*
 * t = 2t, and the tangent at t evaluated at p, for t in homogeneous projective coordinates (x = X/Z, y = Y/Z) on the
 * twist y^2 = x^3 + b' (Costello, Lange and Naehrig, ePrint 2009/615, as arkworks doubles): with B = Y^2, C = Z^2,
 * E = 3 b' C, F = 3E, H = 2 Y Z, I = E - B and J = X^2, 2t = (X Y (B - F) / 2, ((B + F) / 2)^2 - 3 E^2, B H), and the
 * tangent, scaled by -Z^2, has cy = -H, cx = 3J and c0 = I (y^2 - 3b' = 3x^3 - 2y^2 on the curve).
 */
TOWER void miller_dbl_n(line *l, point *t, const point *p, uint32_t n)
{
    fp2 a, b, c, e, f, g, h, i, j, s;
    fp2_mul_n(&a, &t->x, &t->y, n);
    fq_half(&a, &a);
    fp2_mul_n(&b, &t->y, &t->y, n);
    fp2_mul_n(&c, &t->z, &t->z, n);
    fp2_add_n(&e, &c, &c, n);
    fp2_add_n(&e, &e, &c, n);
    fp2_mul_n(&e, &e, &C->b2, n);
    fp2_add_n(&f, &e, &e, n);
    fp2_add_n(&f, &f, &e, n);
    fp2_add_n(&g, &b, &f, n);
    fq_half(&g, &g);
    fp2_add_n(&h, &t->y, &t->z, n);
    fp2_mul_n(&h, &h, &h, n);
    fp2_sub_n(&h, &h, &b, n);
    fp2_sub_n(&h, &h, &c, n);
    fp2_sub_n(&i, &e, &b, n);
    fp2_mul_n(&j, &t->x, &t->x, n);

    fp2_sub_n(&s, &b, &f, n);
    fp2_mul_n(&t->x, &a, &s, n);
    fp2_mul_n(&s, &e, &e, n);
    fp2_mul_n(&t->y, &g, &g, n);
    fp2_sub_n(&t->y, &t->y, &s, n);
    fp2_sub_n(&t->y, &t->y, &s, n);
    fp2_sub_n(&t->y, &t->y, &s, n);
    fp2_mul_n(&t->z, &b, &h, n);

    fq_neg(&h, &h, 2);
    fp2_add_n(&s, &j, &j, n);
    fp2_add_n(&s, &s, &j, n);
    line_eval_n(l, &h, &s, &i, p, n);
}
BY_LIMBS(miller_dbl, (line *l, point *t, const point *p), l, t, p)

/*
 * f *= the chord through t and the affine q evaluated at p, and t += q, t projective as in miller_dbl (arkworks'
 * mixed addition): with O = Y - y_Q Z and L = X - x_Q Z (the slope is O/L), E = L^3, F = Z O^2, G = X L^2 and
 * H = E + F - 2G, t + q = (L H, O (G - H) - E Y, Z E), and the chord, scaled by L Z, has cy = L, cx = -O and
 * c0 = O x_Q - L y_Q. t is neither q nor -q: the loop reaches no such multiple of a point of order r.
 */
TOWER void miller_add_n(fp12 *f, point *t, const point *q, const point *p, uint32_t n)
{
    fp2 o, el, c, d, e, g, h, s;
    line l;
    fp2_mul_n(&s, &q->y, &t->z, n);
    fp2_sub_n(&o, &t->y, &s, n);
    fp2_mul_n(&s, &q->x, &t->z, n);
    fp2_sub_n(&el, &t->x, &s, n);
    fp2_mul_n(&c, &o, &o, n);
    fp2_mul_n(&d, &el, &el, n);
    fp2_mul_n(&e, &el, &d, n);
    fp2_mul_n(&c, &t->z, &c, n);
    fp2_mul_n(&g, &t->x, &d, n);
    fp2_add_n(&h, &e, &c, n);
    fp2_sub_n(&h, &h, &g, n);
    fp2_sub_n(&h, &h, &g, n);

    fp2_mul_n(&t->x, &el, &h, n);
    fp2_sub_n(&s, &g, &h, n);
    fp2_mul_n(&s, &s, &o, n);
    fp2_mul_n(&g, &e, &t->y, n);
    fp2_sub_n(&t->y, &s, &g, n);
    fp2_mul_n(&t->z, &t->z, &e, n);

    fp2_mul_n(&s, &o, &q->x, n);
    fp2_mul_n(&g, &el, &q->y, n);
    fp2_sub_n(&s, &s, &g, n);
    fq_neg(&o, &o, 2);
    line_eval_n(&l, &el, &o, &s, p, n);
    fp12_mul_line(f, &l);
}
BY_LIMBS(miller_add, (fp12 *f, point *t, const point *q, const point *p), f, t, q, p)

/*
 * Multiplies f by the Miller functions of the optimal ate pairing of each pair g1[i] (G1) and
 * g2[i] (G2), at most PairingBatch of them, in one loop that squares its accumulator once a step
 * for all; a pair with a point at infinity pairs to 1 and leaves f alone. BN254 loops over 6x + 2
 * and closes with the lines to pi(q) and -pi^2(q). BLS12-381 loops over |x| = -x; the function of
 * x itself is the inverse of that one up to a vertical line, so every pair's value is inverted
 * alike and a product of pairings is 1 either both ways or neither.
 */
void pairing_add_batch(fp12 *f, const point *g1, const point *g2, int count)
{
    point pa[PairingBatch], qa[PairingBatch], qn[PairingBatch], t[PairingBatch];
    int k = 0;
    for (int i = 0; i < count; i++)
    {
        if (pt_is_inf(&g1[i], 1) || pt_is_inf(&g2[i], 2))
            continue;
        pa[k] = g1[i];
        qa[k] = g2[i];
        pt_set_affine(&pa[k], 1);
        pt_set_affine(&qa[k], 2);
        pt_neg(&qn[k], &qa[k], 2);
        t[k] = qa[k];
        k++;
    }
    if (k == 0)
        return;
    fp12 acc;
    line l;
    fp12_one(&acc);
    /* A digit -1 adds -q: the Miller function for -1 is a vertical line, which the final exponentiation removes. */
    int8_t dg[8 * sizeof(C->loop) + 1];
    int top = scalar_digits(dg, C->loop, sizeof(C->loop)) - 1;
    for (int i = top - 1; i >= 0; i--)
    {
        /* The first square is of 1. */
        if (i != top - 1)
            fp12_mul(&acc, &acc, &acc);
        for (int j = 0; j < k; j++)
        {
            miller_dbl(&l, &t[j], &pa[j]);
            fp12_mul_line(&acc, &l);
        }
        if (dg[i])
            for (int j = 0; j < k; j++)
                miller_add(&acc, &t[j], dg[i] > 0 ? &qa[j] : &qn[j], &pa[j]);
    }
    if (!C->mtwist)
    {
        for (int j = 0; j < k; j++)
        {
            point q1, q2;
            pt_frob(&q1, &qa[j], false);
            pt_frob(&q2, &q1, false);
            pt_neg(&q2, &q2, 2);
            miller_add(&acc, &t[j], &q1, &pa[j]);
            miller_add(&acc, &t[j], &q2, &pa[j]);
        }
    }
    fp12_mul(f, f, &acc);
}

/*
 * Hard part of the BN254 final exponentiation, Algorithm 6 of Duquesne and Ghammam (ePrint
 * 2015/192) as arkworks has it: a power of f coprime to r times (p^4 - p^2 + 1)/r.
 */
NOINLINE static void bn_hard(fp12 *r, const fp12 *f)
{
    fp12 y1, y3, y4, y8, y9, t;
    fp12_pow_x(&t, f);
    fp12_conj(&t, &t);
    fp12_cyclotomic_sqr(&y1, &t);
    fp12_cyclotomic_sqr(&y3, &y1);
    fp12_mul(&y3, &y3, &y1);
    fp12_pow_x(&y4, &y3);
    fp12_conj(&y4, &y4);
    fp12_cyclotomic_sqr(&t, &y4);
    fp12_pow_x(&t, &t);
    fp12_mul(&y8, &t, &y4);
    fp12_conj(&t, &y3);
    fp12_mul(&y8, &y8, &t);
    fp12_mul(&y9, &y8, &y1);
    fp12_mul(&y4, &y8, &y4);
    fp12_mul(&y4, &y4, f);
    fp12_frob(&t, &y9);
    fp12_mul(r, &t, &y4);
    fp12_frob(&t, &y8);
    fp12_frob(&t, &t);
    fp12_mul(r, r, &t);
    fp12_conj(&t, f);
    fp12_mul(&t, &t, &y9);
    fp12_frob(&t, &t);
    fp12_frob(&t, &t);
    fp12_frob(&t, &t);
    fp12_mul(r, r, &t);
}

/*
 * Hard part of the BLS12-381 final exponentiation by 3(p^4 - p^2 + 1)/r =
 * (x - 1)^2 (x + p)(x^2 + p^2 - 1) + 3 (Hayashida, Hayasaka and Teruya, ePrint 2020/875).
 */
NOINLINE static void bls_hard(fp12 *r, const fp12 *f)
{
    fp12 a, b, t;
    fp12_pow_x(&t, f);
    fp12_conj(&a, f);
    fp12_mul(&a, &a, &t);
    fp12_pow_x(&t, &a);
    fp12_conj(&a, &a);
    fp12_mul(&a, &a, &t);
    fp12_pow_x(&t, &a);
    fp12_frob(&b, &a);
    fp12_mul(&b, &b, &t);
    fp12_pow_x(&t, &b);
    fp12_pow_x(&t, &t);
    fp12_frob(&a, &b);
    fp12_frob(&a, &a);
    fp12_mul(&t, &t, &a);
    fp12_conj(&b, &b);
    fp12_mul(&t, &t, &b);
    fp12_cyclotomic_sqr(&a, f);
    fp12_mul(&a, &a, f);
    fp12_mul(r, &t, &a);
}

static bool fp12_is_one(const fp12 *a)
{
    fp12 one;
    fp12_one(&one);
    for (int i = 0; i < 6; i++)
        if (!fq_eq(&a->c[i], &one.c[i], 2))
            return false;
    return true;
}

/* Whether f^((p^12 - 1)/r) == 1; f is destroyed. */
bool pairing_is_one(fp12 *f)
{
    fp12 t;
    /* Easy part, f^((p^6 - 1)(p^2 + 1)), which lands f in the cyclotomic subgroup. */
    fp12_inv(&t, f);
    fp12_conj(f, f);
    fp12_mul(f, f, &t);
    fp12_frob(&t, f);
    fp12_frob(&t, &t);
    fp12_mul(f, f, &t);
    if (C->mtwist)
        bls_hard(&t, f);
    else
        bn_hard(&t, f);
    return fp12_is_one(&t);
}
