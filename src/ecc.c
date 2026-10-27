/*
 * secp256k1 public-key recovery and secp256r1 (P-256) verification on 32-bit limbs.
 *
 * Every modulus used here (the field primes and group orders of both curves) has the form
 * 2^256 - c with c < 2^225, so a 512-bit product reduces by folding: hi * 2^256 + lo is
 * congruent to hi * c + lo, which shrinks hi each round. One routine serves all four moduli.
 * Points are Jacobian; the double-scalar product u1*G + u2*Q is computed with Shamir's trick.
 * Nothing here is constant time: a zkVM guest has no side channel.
 */

#include "zkvm.h"

typedef uint32_t u256[8]; /* little-endian limbs */

typedef struct
{
    u256 m;
    u256 c; /* 2^256 - m */
} modulus;

typedef struct
{
    modulus p;
    modulus n;
    u256 a;
    u256 b;
    u256 gx;
    u256 gy;
} curve;

typedef struct
{
    u256 x, y, z; /* z == 0 is the point at infinity */
} jpoint;

static void u256_from_be(u256 r, const uint8_t in[32])
{
    for (int i = 0; i < 8; i++)
        r[i] = ((uint32_t)in[31 - 4 * i - 3] << 24) | ((uint32_t)in[31 - 4 * i - 2] << 16) |
               ((uint32_t)in[31 - 4 * i - 1] << 8) | in[31 - 4 * i];
}

static void u256_to_be(uint8_t out[32], const u256 a)
{
    for (int i = 0; i < 8; i++)
    {
        out[31 - 4 * i] = (uint8_t)a[i];
        out[31 - 4 * i - 1] = (uint8_t)(a[i] >> 8);
        out[31 - 4 * i - 2] = (uint8_t)(a[i] >> 16);
        out[31 - 4 * i - 3] = (uint8_t)(a[i] >> 24);
    }
}

static bool u256_is_zero(const u256 a)
{
    uint32_t acc = 0;
    for (int i = 0; i < 8; i++)
        acc |= a[i];
    return acc == 0;
}

static int u256_cmp(const u256 a, const u256 b)
{
    for (int i = 7; i >= 0; i--)
        if (a[i] != b[i])
            return a[i] < b[i] ? -1 : 1;
    return 0;
}

static uint32_t u256_add(u256 r, const u256 a, const u256 b)
{
    uint64_t carry = 0;
    for (int i = 0; i < 8; i++)
    {
        carry += (uint64_t)a[i] + b[i];
        r[i] = (uint32_t)carry;
        carry >>= 32;
    }
    return (uint32_t)carry;
}

static uint32_t u256_sub(u256 r, const u256 a, const u256 b)
{
    int64_t borrow = 0;
    for (int i = 0; i < 8; i++)
    {
        borrow += (int64_t)a[i] - b[i];
        r[i] = (uint32_t)borrow;
        borrow >>= 32;
    }
    return (uint32_t)(borrow != 0);
}

static void mul_wide(uint32_t t[16], const u256 a, const u256 b)
{
    memset(t, 0, 16 * sizeof(uint32_t));
    for (int i = 0; i < 8; i++)
    {
        uint64_t carry = 0;
        for (int j = 0; j < 8; j++)
        {
            carry += (uint64_t)a[i] * b[j] + t[i + j];
            t[i + j] = (uint32_t)carry;
            carry >>= 32;
        }
        t[i + 8] = (uint32_t)carry;
    }
}

static void mod_reduce(u256 r, uint32_t t[16], const modulus *m)
{
    for (;;)
    {
        u256 hi;
        memcpy(hi, t + 8, sizeof(hi));
        if (u256_is_zero(hi))
            break;
        uint32_t prod[16];
        mul_wide(prod, hi, m->c);
        uint64_t carry = 0;
        for (int i = 0; i < 16; i++)
        {
            carry += (uint64_t)prod[i] + (i < 8 ? t[i] : 0);
            t[i] = (uint32_t)carry;
            carry >>= 32;
        }
    }
    memcpy(r, t, sizeof(u256));
    while (u256_cmp(r, m->m) >= 0)
        u256_sub(r, r, m->m);
}

static void mod_mul(u256 r, const u256 a, const u256 b, const modulus *m)
{
    uint32_t t[16];
    mul_wide(t, a, b);
    mod_reduce(r, t, m);
}

static void mod_add(u256 r, const u256 a, const u256 b, const modulus *m)
{
    uint32_t carry = u256_add(r, a, b);
    if (carry || u256_cmp(r, m->m) >= 0)
        u256_sub(r, r, m->m);
}

static void mod_sub(u256 r, const u256 a, const u256 b, const modulus *m)
{
    if (u256_sub(r, a, b))
        u256_add(r, r, m->m);
}

static void mod_pow(u256 r, const u256 base, const u256 e, const modulus *m)
{
    u256 acc = {1};
    for (int i = 255; i >= 0; i--)
    {
        mod_mul(acc, acc, acc, m);
        if ((e[i / 32] >> (i % 32)) & 1)
            mod_mul(acc, acc, base, m);
    }
    memcpy(r, acc, sizeof(u256));
}

/* Inverse by Fermat; every modulus here is prime. */
static void mod_inv(u256 r, const u256 a, const modulus *m)
{
    u256 e, two = {2};
    u256_sub(e, m->m, two);
    mod_pow(r, a, e, m);
}

static void jdouble(jpoint *r, const jpoint *p, const curve *cv)
{
    const modulus *f = &cv->p;
    if (u256_is_zero(p->z) || u256_is_zero(p->y))
    {
        memset(r, 0, sizeof(*r));
        return;
    }
    u256 xx, yy, yyyy, zz, s, m, t, x3, y3, z3;
    mod_mul(xx, p->x, p->x, f);
    mod_mul(yy, p->y, p->y, f);
    mod_mul(yyyy, yy, yy, f);
    mod_mul(zz, p->z, p->z, f);

    mod_mul(s, p->x, yy, f);
    mod_add(s, s, s, f);
    mod_add(s, s, s, f);

    mod_add(m, xx, xx, f);
    mod_add(m, m, xx, f);
    if (!u256_is_zero(cv->a))
    {
        mod_mul(t, zz, zz, f);
        mod_mul(t, t, cv->a, f);
        mod_add(m, m, t, f);
    }

    mod_mul(x3, m, m, f);
    mod_sub(x3, x3, s, f);
    mod_sub(x3, x3, s, f);

    mod_sub(t, s, x3, f);
    mod_mul(y3, m, t, f);
    mod_add(t, yyyy, yyyy, f);
    mod_add(t, t, t, f);
    mod_add(t, t, t, f);
    mod_sub(y3, y3, t, f);

    mod_mul(z3, p->y, p->z, f);
    mod_add(z3, z3, z3, f);

    memcpy(r->x, x3, sizeof(u256));
    memcpy(r->y, y3, sizeof(u256));
    memcpy(r->z, z3, sizeof(u256));
}

static void jadd(jpoint *r, const jpoint *p, const jpoint *q, const curve *cv)
{
    const modulus *f = &cv->p;
    if (u256_is_zero(p->z))
    {
        *r = *q;
        return;
    }
    if (u256_is_zero(q->z))
    {
        *r = *p;
        return;
    }
    u256 z1z1, z2z2, u1, u2, s1, s2, h, rr, hh, hhh, v, x3, y3, z3, t;
    mod_mul(z1z1, p->z, p->z, f);
    mod_mul(z2z2, q->z, q->z, f);
    mod_mul(u1, p->x, z2z2, f);
    mod_mul(u2, q->x, z1z1, f);
    mod_mul(t, q->z, z2z2, f);
    mod_mul(s1, p->y, t, f);
    mod_mul(t, p->z, z1z1, f);
    mod_mul(s2, q->y, t, f);

    if (u256_cmp(u1, u2) == 0)
    {
        if (u256_cmp(s1, s2) == 0)
            jdouble(r, p, cv);
        else
            memset(r, 0, sizeof(*r));
        return;
    }

    mod_sub(h, u2, u1, f);
    mod_sub(rr, s2, s1, f);
    mod_mul(hh, h, h, f);
    mod_mul(hhh, hh, h, f);
    mod_mul(v, u1, hh, f);

    mod_mul(x3, rr, rr, f);
    mod_sub(x3, x3, hhh, f);
    mod_sub(x3, x3, v, f);
    mod_sub(x3, x3, v, f);

    mod_sub(t, v, x3, f);
    mod_mul(y3, rr, t, f);
    mod_mul(t, s1, hhh, f);
    mod_sub(y3, y3, t, f);

    mod_mul(z3, p->z, q->z, f);
    mod_mul(z3, z3, h, f);

    memcpy(r->x, x3, sizeof(u256));
    memcpy(r->y, y3, sizeof(u256));
    memcpy(r->z, z3, sizeof(u256));
}

/* u1*G + u2*Q. */
static void double_mul(jpoint *r, const u256 u1, const u256 u2, const jpoint *q, const curve *cv)
{
    jpoint g, gq, acc;
    memcpy(g.x, cv->gx, sizeof(u256));
    memcpy(g.y, cv->gy, sizeof(u256));
    memset(g.z, 0, sizeof(u256));
    g.z[0] = 1;
    jadd(&gq, &g, q, cv);

    memset(&acc, 0, sizeof(acc));
    for (int i = 255; i >= 0; i--)
    {
        jdouble(&acc, &acc, cv);
        int b1 = (u1[i / 32] >> (i % 32)) & 1;
        int b2 = (u2[i / 32] >> (i % 32)) & 1;
        if (b1 && b2)
            jadd(&acc, &acc, &gq, cv);
        else if (b1)
            jadd(&acc, &acc, &g, cv);
        else if (b2)
            jadd(&acc, &acc, q, cv);
    }
    *r = acc;
}

static void to_affine(u256 x, u256 y, const jpoint *p, const curve *cv)
{
    const modulus *f = &cv->p;
    u256 zinv, zinv2, zinv3;
    mod_inv(zinv, p->z, f);
    mod_mul(zinv2, zinv, zinv, f);
    mod_mul(zinv3, zinv2, zinv, f);
    mod_mul(x, p->x, zinv2, f);
    mod_mul(y, p->y, zinv3, f);
}

static void set_modulus(modulus *m, const uint8_t be[32])
{
    u256 zero = {0};
    u256_from_be(m->m, be);
    u256_sub(m->c, zero, m->m);
}

static void load_curve(curve *cv, const uint8_t p[32], const uint8_t n[32], const uint8_t a[32], const uint8_t b[32],
                       const uint8_t gx[32], const uint8_t gy[32])
{
    set_modulus(&cv->p, p);
    set_modulus(&cv->n, n);
    u256_from_be(cv->a, a);
    u256_from_be(cv->b, b);
    u256_from_be(cv->gx, gx);
    u256_from_be(cv->gy, gy);
}

static const uint8_t K1_P[32] = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xfc, 0x2f};
static const uint8_t K1_N[32] = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xba, 0xae, 0xdc, 0xe6, 0xaf, 0x48, 0xa0, 0x3b, 0xbf, 0xd2, 0x5e, 0x8c, 0xd0, 0x36, 0x41, 0x41};
static const uint8_t K1_A[32] = {0};
static const uint8_t K1_B[32] = {[31] = 7};
static const uint8_t K1_GX[32] = {
    0x79, 0xbe, 0x66, 0x7e, 0xf9, 0xdc, 0xbb, 0xac, 0x55, 0xa0, 0x62, 0x95, 0xce, 0x87, 0x0b, 0x07,
    0x02, 0x9b, 0xfc, 0xdb, 0x2d, 0xce, 0x28, 0xd9, 0x59, 0xf2, 0x81, 0x5b, 0x16, 0xf8, 0x17, 0x98};
static const uint8_t K1_GY[32] = {
    0x48, 0x3a, 0xda, 0x77, 0x26, 0xa3, 0xc4, 0x65, 0x5d, 0xa4, 0xfb, 0xfc, 0x0e, 0x11, 0x08, 0xa8,
    0xfd, 0x17, 0xb4, 0x48, 0xa6, 0x85, 0x54, 0x19, 0x9c, 0x47, 0xd0, 0x8f, 0xfb, 0x10, 0xd4, 0xb8};

static const uint8_t R1_P[32] = {
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
static const uint8_t R1_N[32] = {
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xbc, 0xe6, 0xfa, 0xad, 0xa7, 0x17, 0x9e, 0x84, 0xf3, 0xb9, 0xca, 0xc2, 0xfc, 0x63, 0x25, 0x51};
static const uint8_t R1_A[32] = {
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfc};
static const uint8_t R1_B[32] = {
    0x5a, 0xc6, 0x35, 0xd8, 0xaa, 0x3a, 0x93, 0xe7, 0xb3, 0xeb, 0xbd, 0x55, 0x76, 0x98, 0x86, 0xbc,
    0x65, 0x1d, 0x06, 0xb0, 0xcc, 0x53, 0xb0, 0xf6, 0x3b, 0xce, 0x3c, 0x3e, 0x27, 0xd2, 0x60, 0x4b};
static const uint8_t R1_GX[32] = {
    0x6b, 0x17, 0xd1, 0xf2, 0xe1, 0x2c, 0x42, 0x47, 0xf8, 0xbc, 0xe6, 0xe5, 0x63, 0xa4, 0x40, 0xf2,
    0x77, 0x03, 0x7d, 0x81, 0x2d, 0xeb, 0x33, 0xa0, 0xf4, 0xa1, 0x39, 0x45, 0xd8, 0x98, 0xc2, 0x96};
static const uint8_t R1_GY[32] = {
    0x4f, 0xe3, 0x42, 0xe2, 0xfe, 0x1a, 0x7f, 0x9b, 0x8e, 0xe7, 0xeb, 0x4a, 0x7c, 0x0f, 0x9e, 0x16,
    0x2b, 0xce, 0x33, 0x57, 0x6b, 0x31, 0x5e, 0xce, 0xcb, 0xb6, 0x40, 0x68, 0x37, 0xbf, 0x51, 0xf5};

static bool in_scalar_range(const u256 v, const curve *cv)
{
    return !u256_is_zero(v) && u256_cmp(v, cv->n.m) < 0;
}

/* The message hash as a scalar: its 256 bits taken mod n. */
static void hash_to_scalar(u256 e, const uint8_t msg[32], const curve *cv)
{
    u256_from_be(e, msg);
    if (u256_cmp(e, cv->n.m) >= 0)
        u256_sub(e, e, cv->n.m);
}

static bool on_curve(const u256 x, const u256 y, const curve *cv)
{
    const modulus *f = &cv->p;
    u256 lhs, rhs, t;
    mod_mul(lhs, y, y, f);
    mod_mul(rhs, x, x, f);
    mod_mul(rhs, rhs, x, f);
    mod_mul(t, cv->a, x, f);
    mod_add(rhs, rhs, t, f);
    mod_add(rhs, rhs, cv->b, f);
    return u256_cmp(lhs, rhs) == 0;
}

zkvm_status zkvm_secp256k1_ecrecover(const uint8_t msg[32], const uint8_t sig[64], uint8_t recid,
                                     uint8_t output[64])
{
    curve cv;
    load_curve(&cv, K1_P, K1_N, K1_A, K1_B, K1_GX, K1_GY);
    const modulus *f = &cv.p;
    const modulus *n = &cv.n;

    u256 r, s;
    u256_from_be(r, sig);
    u256_from_be(s, sig + 32);
    if (recid > 3 || !in_scalar_range(r, &cv) || !in_scalar_range(s, &cv))
        return ZKVM_EFAIL;

    u256 x;
    memcpy(x, r, sizeof(u256));
    if (recid & 2)
    {
        if (u256_add(x, x, n->m))
            return ZKVM_EFAIL;
    }
    if (u256_cmp(x, f->m) >= 0)
        return ZKVM_EFAIL;

    /* y^2 = x^3 + 7; p = 3 mod 4, so a root is alpha^((p+1)/4). */
    u256 alpha, beta, check, e;
    mod_mul(alpha, x, x, f);
    mod_mul(alpha, alpha, x, f);
    mod_add(alpha, alpha, cv.b, f);
    u256 one = {1};
    u256_add(e, f->m, one);
    for (int i = 0; i < 8; i++)
        e[i] = (e[i] >> 2) | (i < 7 ? e[i + 1] << 30 : 0);
    mod_pow(beta, alpha, e, f);
    mod_mul(check, beta, beta, f);
    if (u256_cmp(check, alpha) != 0)
        return ZKVM_EFAIL;
    if ((beta[0] & 1) != (recid & 1))
    {
        u256 zero = {0};
        mod_sub(beta, zero, beta, f);
    }

    jpoint rp;
    memcpy(rp.x, x, sizeof(u256));
    memcpy(rp.y, beta, sizeof(u256));
    memset(rp.z, 0, sizeof(u256));
    rp.z[0] = 1;

    /* Q = r^-1 (s R - e G) */
    u256 rinv, u1, u2, zero = {0};
    hash_to_scalar(e, msg, &cv);
    mod_inv(rinv, r, n);
    mod_sub(u1, zero, e, n);
    mod_mul(u1, u1, rinv, n);
    mod_mul(u2, s, rinv, n);

    jpoint q;
    double_mul(&q, u1, u2, &rp, &cv);
    if (u256_is_zero(q.z))
        return ZKVM_EFAIL;

    u256 qx, qy;
    to_affine(qx, qy, &q, &cv);
    u256_to_be(output, qx);
    u256_to_be(output + 32, qy);
    return ZKVM_EOK;
}

zkvm_status zkvm_secp256r1_verify(const uint8_t msg[32], const uint8_t sig[64], const uint8_t pubkey[64],
                                  bool *verified)
{
    curve cv;
    load_curve(&cv, R1_P, R1_N, R1_A, R1_B, R1_GX, R1_GY);
    const modulus *n = &cv.n;
    *verified = false;

    u256 r, s, qx, qy;
    u256_from_be(r, sig);
    u256_from_be(s, sig + 32);
    u256_from_be(qx, pubkey);
    u256_from_be(qy, pubkey + 32);
    if (!in_scalar_range(r, &cv) || !in_scalar_range(s, &cv))
        return ZKVM_EOK;
    if (u256_cmp(qx, cv.p.m) >= 0 || u256_cmp(qy, cv.p.m) >= 0 || (u256_is_zero(qx) && u256_is_zero(qy)) ||
        !on_curve(qx, qy, &cv))
        return ZKVM_EOK;

    u256 e, w, u1, u2;
    hash_to_scalar(e, msg, &cv);
    mod_inv(w, s, n);
    mod_mul(u1, e, w, n);
    mod_mul(u2, r, w, n);

    jpoint q, rp;
    memcpy(q.x, qx, sizeof(u256));
    memcpy(q.y, qy, sizeof(u256));
    memset(q.z, 0, sizeof(u256));
    q.z[0] = 1;
    double_mul(&rp, u1, u2, &q, &cv);
    if (u256_is_zero(rp.z))
        return ZKVM_EOK;

    u256 x, y;
    to_affine(x, y, &rp, &cv);
    if (u256_cmp(x, n->m) >= 0)
        u256_sub(x, x, n->m);
    *verified = u256_cmp(x, r) == 0;
    return ZKVM_EOK;
}
