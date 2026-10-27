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

static uint32_t limbs_sub(uint32_t *r, const uint32_t *a, const uint32_t *b, uint32_t n)
{
    int64_t borrow = 0;
    for (uint32_t i = 0; i < n; i++)
    {
        borrow += (int64_t)a[i] - b[i];
        r[i] = (uint32_t)borrow;
        borrow >>= 32;
    }
    return borrow != 0;
}

/* Both primes leave the top bit of their last limb clear, so a sum of two reduced elements
   never carries out of n limbs. */
void fp_add(fp *r, const fp *a, const fp *b)
{
    uint32_t n = C->n;
    uint64_t carry = 0;
    for (uint32_t i = 0; i < n; i++)
    {
        carry += (uint64_t)a->l[i] + b->l[i];
        r->l[i] = (uint32_t)carry;
        carry >>= 32;
    }
    if (limbs_cmp(r->l, C->p.l, n) >= 0)
        limbs_sub(r->l, r->l, C->p.l, n);
}

void fp_sub(fp *r, const fp *a, const fp *b)
{
    uint32_t n = C->n;
    if (limbs_sub(r->l, a->l, b->l, n))
    {
        uint64_t carry = 0;
        for (uint32_t i = 0; i < n; i++)
        {
            carry += (uint64_t)r->l[i] + C->p.l[i];
            r->l[i] = (uint32_t)carry;
            carry >>= 32;
        }
    }
}

/* Montgomery multiplication on the u256 delegation. An element is split into 256-bit digits, a
 * low one and, for BLS12-381, a high one below 2^128; the product accumulates in T, four digits.
 * The reduction adds multiples of p that clear T from the bottom: one 256-bit step, which is the
 * whole of it for BN254 (R = 2^256), and for BLS12-381 (R = 2^384) a second step that clears only
 * 128 bits. The delegation takes its operands from RAM, 32-byte aligned and distinct, so they all
 * live in these statics. */
static u256d A0, A1, B0, B1, M, X, T[4], Zero, Pinv, P0, P1;
static const curve *Prepared;

/* Pinv = -p^-1 mod 2^256 and p as digits, once per curve. */
static void prepare(void)
{
    static u256d x, t, two;
    for (int i = 0; i < 8; i++)
    {
        P0.w[i] = C->p.l[i];
        P1.w[i] = i + 8 < (int)C->n ? C->p.l[i + 8] : 0;
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
    Pinv = (u256d){{0}};
    u256_op(&Pinv, &x, U256_SUB);
    Prepared = C;
}

/* T[k..] += *v. */
static void accumulate(int k, const u256d *v)
{
    uint32_t carry = u256_op(&T[k], v, U256_ADD);
    for (int i = k + 1; carry && i < 4; i++)
        carry = u256_op(&T[i], &Zero, U256_ADD | U256_CARRY);
}

/* T += a * b * 2^(256 k); the high half only when the product can reach it. */
static void accumulate_product(int k, const u256d *a, const u256d *b, bool high)
{
    u256_op(&X, a, U256_MEMCOPY);
    u256_op(&X, b, U256_MUL_LOW);
    accumulate(k, &X);
    if (high)
    {
        u256_op(&X, a, U256_MEMCOPY);
        u256_op(&X, b, U256_MUL_HIGH);
        accumulate(k + 1, &X);
    }
}

/* Montgomery product a*b/R mod p. */
void fp_mul(fp *r, const fp *a, const fp *b)
{
    uint32_t n = C->n;
    if (Prepared != C)
        prepare();
    for (int i = 0; i < 8; i++)
    {
        A0.w[i] = a->l[i];
        B0.w[i] = b->l[i];
        A1.w[i] = i + 8 < (int)n ? a->l[i + 8] : 0;
        B1.w[i] = i + 8 < (int)n ? b->l[i + 8] : 0;
    }
    for (int i = 0; i < 4; i++)
        u256_op(&T[i], &Zero, U256_MEMCOPY);

    fp res = {{0}};
    if (n == 8)
    {
        accumulate_product(0, &A0, &B0, true);
        u256_op(&M, &T[0], U256_MEMCOPY);
        u256_op(&M, &Pinv, U256_MUL_LOW);
        accumulate_product(0, &M, &P0, true);
        memcpy(res.l, T[1].w, 8 * sizeof(uint32_t));
    }
    else
    {
        accumulate_product(0, &A0, &B0, true);
        accumulate_product(1, &A0, &B1, true);
        accumulate_product(1, &A1, &B0, true);
        accumulate_product(2, &A1, &B1, false);

        u256_op(&M, &T[0], U256_MEMCOPY);
        u256_op(&M, &Pinv, U256_MUL_LOW);
        accumulate_product(0, &M, &P0, true);
        accumulate_product(1, &M, &P1, true);

        u256_op(&M, &T[1], U256_MEMCOPY);
        u256_op(&M, &Pinv, U256_MUL_LOW);
        for (int i = 4; i < 8; i++)
            M.w[i] = 0;
        accumulate_product(1, &M, &P0, true);
        accumulate_product(2, &M, &P1, false);

        memcpy(res.l, T[1].w + 4, 4 * sizeof(uint32_t));
        memcpy(res.l + 4, T[2].w, 8 * sizeof(uint32_t));
    }
    if (limbs_cmp(res.l, C->p.l, n) >= 0)
        limbs_sub(res.l, res.l, C->p.l, n);
    *r = res;
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
    for (int i = 0; i < d; i++)
        fp_add(&r->c[i], &a->c[i], &b->c[i]);
}

NOINLINE void fq_sub(fp2 *r, const fp2 *a, const fp2 *b, int d)
{
    for (int i = 0; i < d; i++)
        fp_sub(&r->c[i], &a->c[i], &b->c[i]);
}

NOINLINE void fq_neg(fp2 *r, const fp2 *a, int d)
{
    fp2 zero;
    memset(&zero, 0, sizeof(zero));
    fq_sub(r, &zero, a, d);
}

/* Fp2 by Karatsuba: (a0 + a1 u)(b0 + b1 u) = a0 b0 - a1 b1 + ((a0 + a1)(b0 + b1) - a0 b0 - a1 b1) u. */
NOINLINE void fq_mul(fp2 *r, const fp2 *a, const fp2 *b, int d)
{
    if (d == 1)
    {
        fp_mul(&r->c[0], &a->c[0], &b->c[0]);
        return;
    }
    fp t0, t1, t2, t3;
    fp_mul(&t0, &a->c[0], &b->c[0]);
    fp_mul(&t1, &a->c[1], &b->c[1]);
    fp_add(&t2, &a->c[0], &a->c[1]);
    fp_add(&t3, &b->c[0], &b->c[1]);
    fp_mul(&t2, &t2, &t3);
    fp_sub(&r->c[0], &t0, &t1);
    fp_sub(&t2, &t2, &t0);
    fp_sub(&r->c[1], &t2, &t1);
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

NOINLINE static void fq_pow(fp2 *r, const fp2 *a, const fp *e, int d)
{
    fp2 acc, base = *a;
    fq_set_one(&acc);
    for (int i = 32 * (int)C->n - 1; i >= 0; i--)
    {
        fq_mul(&acc, &acc, &acc, d);
        if ((e->l[i / 32] >> (i % 32)) & 1)
            fq_mul(&acc, &acc, &base, d);
    }
    *r = acc;
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

/*
 * A square root, false when there is none. Both primes are 3 mod 4: in Fp the root is
 * a^((p+1)/4); in Fp2 it is Algorithm 9 of Adj and Rodriguez-Henriquez, "Square root computation
 * over even extension fields" (2012). Either candidate is checked by squaring.
 */
NOINLINE bool fq_sqrt(fp2 *r, const fp2 *a, int d)
{
    fp e;
    fp2 x, t;
    if (d == 1)
    {
        p_exp(&e, 1, 2);
        fq_pow(&x, a, &e, 1);
    }
    else
    {
        fp2 a1, alpha, minus_one;
        p_exp(&e, -3, 2);
        fq_pow(&a1, a, &e, 2);
        fq_mul(&alpha, &a1, &a1, 2);
        fq_mul(&alpha, &alpha, a, 2);
        fq_mul(&x, &a1, a, 2);
        fq_set_one(&t);
        fq_neg(&minus_one, &t, 2);
        if (fq_eq(&alpha, &minus_one, 2))
        {
            /* x = u x0 = -x0.c1 + x0.c0 u */
            fq_neg(&t, &x, 2);
            t.c[0] = t.c[1];
            t.c[1] = x.c[0];
            x = t;
        }
        else
        {
            fq_add(&alpha, &alpha, &t, 2);
            p_exp(&e, -1, 1);
            fq_pow(&alpha, &alpha, &e, 2);
            fq_mul(&x, &alpha, &x, 2);
        }
    }
    fq_mul(&t, &x, &x, d);
    bool ok = fq_eq(&t, a, d);
    *r = x;
    return ok;
}

/* a * xi with xi = k + u: (a0 k - a1) + (a0 + a1 k) u. */
NOINLINE static void fp2_mul_xi(fp2 *r, const fp2 *a)
{
    fp2 k = *a, t;
    for (uint32_t i = 1; i < C->xi; i++)
        fq_add(&k, &k, a, 2);
    fp_sub(&t.c[0], &k.c[0], &a->c[1]);
    fp_add(&t.c[1], &k.c[1], &a->c[0]);
    *r = t;
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
NOINLINE void pt_dbl(point *r, const point *p, int d)
{
    if (pt_is_inf(p, d) || fq_is_zero(&p->y, d))
    {
        pt_set_inf(r);
        return;
    }
    fp2 xx, yy, yyyy, s, m, t, x3, y3, z3;
    fq_mul(&xx, &p->x, &p->x, d);
    fq_mul(&yy, &p->y, &p->y, d);
    fq_mul(&yyyy, &yy, &yy, d);

    fq_mul(&s, &p->x, &yy, d);
    fq_add(&s, &s, &s, d);
    fq_add(&s, &s, &s, d);

    fq_add(&m, &xx, &xx, d);
    fq_add(&m, &m, &xx, d);

    fq_mul(&x3, &m, &m, d);
    fq_sub(&x3, &x3, &s, d);
    fq_sub(&x3, &x3, &s, d);

    fq_sub(&t, &s, &x3, d);
    fq_mul(&y3, &m, &t, d);
    fq_add(&t, &yyyy, &yyyy, d);
    fq_add(&t, &t, &t, d);
    fq_add(&t, &t, &t, d);
    fq_sub(&y3, &y3, &t, d);

    fq_mul(&z3, &p->y, &p->z, d);
    fq_add(&z3, &z3, &z3, d);

    r->x = x3;
    r->y = y3;
    r->z = z3;
}

/* General addition (add-2007-bl without the doubling shortcut), falling back to doubling. */
NOINLINE void pt_add(point *r, const point *p, const point *q, int d)
{
    if (pt_is_inf(p, d))
    {
        *r = *q;
        return;
    }
    if (pt_is_inf(q, d))
    {
        *r = *p;
        return;
    }
    fp2 z1z1, z2z2, u1, u2, s1, s2, h, rr, hh, hhh, v, x3, y3, z3, t;
    fq_mul(&z1z1, &p->z, &p->z, d);
    fq_mul(&z2z2, &q->z, &q->z, d);
    fq_mul(&u1, &p->x, &z2z2, d);
    fq_mul(&u2, &q->x, &z1z1, d);
    fq_mul(&t, &q->z, &z2z2, d);
    fq_mul(&s1, &p->y, &t, d);
    fq_mul(&t, &p->z, &z1z1, d);
    fq_mul(&s2, &q->y, &t, d);

    if (fq_eq(&u1, &u2, d))
    {
        if (fq_eq(&s1, &s2, d))
            pt_dbl(r, p, d);
        else
            pt_set_inf(r);
        return;
    }

    fq_sub(&h, &u2, &u1, d);
    fq_sub(&rr, &s2, &s1, d);
    fq_mul(&hh, &h, &h, d);
    fq_mul(&hhh, &hh, &h, d);
    fq_mul(&v, &u1, &hh, d);

    fq_mul(&x3, &rr, &rr, d);
    fq_sub(&x3, &x3, &hhh, d);
    fq_sub(&x3, &x3, &v, d);
    fq_sub(&x3, &x3, &v, d);

    fq_sub(&t, &v, &x3, d);
    fq_mul(&y3, &rr, &t, d);
    fq_mul(&t, &s1, &hhh, d);
    fq_sub(&y3, &y3, &t, d);

    fq_mul(&z3, &p->z, &q->z, d);
    fq_mul(&z3, &z3, &h, d);

    r->x = x3;
    r->y = y3;
    r->z = z3;
}

void pt_neg(point *r, const point *p, int d)
{
    *r = *p;
    fq_neg(&r->y, &p->y, d);
}

/* k*p for a big-endian scalar of any length, by double-and-add. */
void pt_mul(point *r, const point *p, const uint8_t *k, size_t klen, int d)
{
    point base = *p, acc;
    pt_set_inf(&acc);
    for (size_t i = 0; i < klen; i++)
        for (int b = 7; b >= 0; b--)
        {
            pt_dbl(&acc, &acc, d);
            if ((k[i] >> b) & 1)
                pt_add(&acc, &acc, &base, d);
        }
    *r = acc;
}

/* x*p, x the signed curve parameter. */
void pt_mul_x(point *r, const point *p, int d)
{
    uint8_t k[8];
    for (int i = 0; i < 8; i++)
        k[i] = (uint8_t)(C->x >> (56 - 8 * i));
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
    fp2 g, g2, g3;
    g = C->gamma;
    if (inverse)
        fq_inv(&g, &g, 2);
    fq_mul(&g2, &g, &g, 2);
    fq_mul(&g3, &g2, &g, 2);
    fp2_conj(&r->x, &p->x);
    fp2_conj(&r->y, &p->y);
    fp2_conj(&r->z, &p->z);
    fq_mul(&r->x, &r->x, &g2, 2);
    fq_mul(&r->y, &r->y, &g3, 2);
}

/* ---- Fp12 ----------------------------------------------------------------------------------- */

void fp12_one(fp12 *r)
{
    memset(r, 0, sizeof(*r));
    r->c[0].c[0] = C->one;
}

/* Schoolbook product modulo w^6 - xi. Zero coefficients of b are skipped, which makes the
   product by a line (three non-zero coefficients) half the cost; a square (a == b) takes each
   cross product once, doubled, for 21 Fp2 products instead of 36. */
NOINLINE static void fp12_mul(fp12 *r, const fp12 *a, const fp12 *b)
{
    fp2 t[11], m;
    memset(t, 0, sizeof(t));
    for (int j = 0; j < 6; j++)
    {
        if (fq_is_zero(&b->c[j], 2))
            continue;
        for (int i = a == b ? j : 0; i < 6; i++)
        {
            fq_mul(&m, &a->c[i], &b->c[j], 2);
            if (a == b && i != j)
                fq_add(&m, &m, &m, 2);
            fq_add(&t[i + j], &t[i + j], &m, 2);
        }
    }
    for (int k = 10; k >= 6; k--)
    {
        fp2_mul_xi(&m, &t[k]);
        fq_add(&t[k - 6], &t[k - 6], &m, 2);
    }
    memcpy(r->c, t, sizeof(r->c));
}

/* The p^6-power Frobenius, w -> -w; it inverts elements of the cyclotomic subgroup. */
NOINLINE static void fp12_conj(fp12 *r, const fp12 *a)
{
    *r = *a;
    for (int i = 1; i < 6; i += 2)
        fq_neg(&r->c[i], &a->c[i], 2);
}

/* The p-power Frobenius: c_i w^i -> conj(c_i) gamma^i w^i. */
NOINLINE static void fp12_frob(fp12 *r, const fp12 *a)
{
    fp2 g;
    fq_set_one(&g);
    for (int i = 0; i < 6; i++)
    {
        fp2_conj(&r->c[i], &a->c[i]);
        fq_mul(&r->c[i], &r->c[i], &g, 2);
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
    fp2 b0 = t.c[0], b1 = t.c[2], b2 = t.c[4], s0, s1, s2, n, m;

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
    fq_mul(&t.c[0], &s0, &n, 2);
    fq_mul(&t.c[2], &s1, &n, 2);
    fq_mul(&t.c[4], &s2, &n, 2);
    fp12_mul(r, &c, &t);
}

/* a^x for the signed curve parameter x; a must be in the cyclotomic subgroup when x < 0. */
NOINLINE static void fp12_pow_x(fp12 *r, const fp12 *a)
{
    fp12 acc, base = *a;
    fp12_one(&acc);
    for (int i = 63; i >= 0; i--)
    {
        fp12_mul(&acc, &acc, &acc);
        if ((C->x >> i) & 1)
            fp12_mul(&acc, &acc, &base);
    }
    if (C->xneg)
        fp12_conj(&acc, &acc);
    *r = acc;
}

/* ---- pairing -------------------------------------------------------------------------------- */

/*
 * A line through points of the twist, evaluated at the G1 point p, is cy*y_P + cx*x_P*w + c0*w^3
 * on a D-type twist and c0 + cx*x_P*w^2 + cy*y_P*w^3 on an M-type one, up to factors in Fp2 and
 * powers of w, which the final exponentiation sends to 1.
 */
NOINLINE static void line_place(fp12 *l, const fp2 *cy, const fp2 *cx, const fp2 *c0, const point *p)
{
    fp2 a, b;
    for (int i = 0; i < 2; i++)
    {
        fp_mul(&a.c[i], &cy->c[i], &p->y.c[0]);
        fp_mul(&b.c[i], &cx->c[i], &p->x.c[0]);
    }
    memset(l, 0, sizeof(*l));
    if (C->mtwist)
    {
        l->c[0] = *c0;
        l->c[2] = b;
        l->c[3] = a;
    }
    else
    {
        l->c[0] = a;
        l->c[1] = b;
        l->c[3] = *c0;
    }
}

/* The tangent at t = (X, Y, Z): with slope 3x^2/2y, scaled by 2y Z^6, cy = 2 Y Z^3,
   cx = -3 X^2 Z^2 and c0 = 3 X^3 - 2 Y^2. */
NOINLINE static void line_dbl(fp12 *l, const point *t, const point *p)
{
    fp2 x2, z2, cy, cx, c0, s;
    fq_mul(&x2, &t->x, &t->x, 2);
    fq_mul(&z2, &t->z, &t->z, 2);

    fq_mul(&cy, &t->y, &z2, 2);
    fq_mul(&cy, &cy, &t->z, 2);
    fq_add(&cy, &cy, &cy, 2);

    fq_mul(&s, &x2, &z2, 2);
    fq_add(&cx, &s, &s, 2);
    fq_add(&cx, &cx, &s, 2);
    fq_neg(&cx, &cx, 2);

    fq_mul(&s, &x2, &t->x, 2);
    fq_add(&c0, &s, &s, 2);
    fq_add(&c0, &c0, &s, 2);
    fq_mul(&s, &t->y, &t->y, 2);
    fq_add(&s, &s, &s, 2);
    fq_sub(&c0, &c0, &s, 2);

    line_place(l, &cy, &cx, &c0, p);
}

/* The chord through t and the affine q: with A = X Z - x_Q Z^3 and B = Y - y_Q Z^3 (the slope
   is B/A), cy = A, cx = -B and c0 = B x_Q - A y_Q. */
NOINLINE static void line_add(fp12 *l, const point *t, const point *q, const point *p)
{
    fp2 z3, a, b, c0, s;
    fq_mul(&z3, &t->z, &t->z, 2);
    fq_mul(&z3, &z3, &t->z, 2);

    fq_mul(&a, &t->x, &t->z, 2);
    fq_mul(&s, &q->x, &z3, 2);
    fq_sub(&a, &a, &s, 2);

    fq_mul(&s, &q->y, &z3, 2);
    fq_sub(&b, &t->y, &s, 2);

    fq_mul(&c0, &b, &q->x, 2);
    fq_mul(&s, &q->y, &a, 2);
    fq_sub(&c0, &c0, &s, 2);

    fq_neg(&b, &b, 2);
    line_place(l, &a, &b, &c0, p);
}

NOINLINE static void miller_add(fp12 *f, point *t, const point *q, const point *p)
{
    fp12 l;
    line_add(&l, t, q, p);
    fp12_mul(f, f, &l);
    pt_add(t, t, q, 2);
}

/*
 * Multiplies f by the Miller function of the optimal ate pairing of p (G1) and q (G2); a pair
 * with a point at infinity pairs to 1 and leaves f alone. BN254 loops over 6x + 2 and closes
 * with the lines to pi(q) and -pi^2(q). BLS12-381 loops over |x| = -x; the function of x itself
 * is the inverse of that one up to a vertical line, so every pair's value is inverted alike and
 * a product of pairings is 1 either both ways or neither.
 */
void pairing_add(fp12 *f, const point *g1, const point *g2)
{
    if (pt_is_inf(g1, 1) || pt_is_inf(g2, 2))
        return;
    fp12 acc, l;
    point pa = *g1, qa = *g2;
    pt_set_affine(&pa, 1);
    pt_set_affine(&qa, 2);
    const point *p = &pa, *q = &qa;
    point t = qa;
    fp12_one(&acc);
    int top = 8 * (int)sizeof(C->loop) - 1;
    while (!((C->loop[sizeof(C->loop) - 1 - top / 8] >> (top % 8)) & 1))
        top--;
    for (int i = top - 1; i >= 0; i--)
    {
        fp12_mul(&acc, &acc, &acc);
        line_dbl(&l, &t, p);
        fp12_mul(&acc, &acc, &l);
        pt_dbl(&t, &t, 2);
        if ((C->loop[sizeof(C->loop) - 1 - i / 8] >> (i % 8)) & 1)
            miller_add(&acc, &t, q, p);
    }
    if (!C->mtwist)
    {
        point q1, q2;
        pt_frob(&q1, q, false);
        pt_frob(&q2, &q1, false);
        pt_neg(&q2, &q2, 2);
        miller_add(&acc, &t, &q1, p);
        miller_add(&acc, &t, &q2, p);
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
    fp12_mul(&y1, &t, &t);
    fp12_mul(&y3, &y1, &y1);
    fp12_mul(&y3, &y3, &y1);
    fp12_pow_x(&y4, &y3);
    fp12_conj(&y4, &y4);
    fp12_mul(&t, &y4, &y4);
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
    fp12_mul(&a, f, f);
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
