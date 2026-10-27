/*
 * secp256k1 public-key recovery and secp256r1 (P-256) verification.
 *
 * Field and scalar arithmetic is Montgomery arithmetic on the u256 delegation (mont.c); field
 * elements live in the Montgomery domain from input decoding to output encoding. Points are
 * Jacobian; the double-scalar product u1*G + u2*Q is computed with Shamir's trick. Nothing here is
 * constant time: a zkVM guest has no side channel.
 */
#include "mont.h"
#include "zkvm.h"

typedef struct
{
    mont_ctx p;
    mont_ctx n;
    u256d a;  /* Montgomery form */
    u256d b;  /* Montgomery form */
    u256d gx; /* Montgomery form */
    u256d gy; /* Montgomery form */
    int a_is_zero;
} curve;

typedef struct
{
    u256d x, y, z; /* Montgomery form; z == 0 is the point at infinity */
} jpoint;

static void from_be(u256d *r, const uint8_t in[32])
{
    for (int i = 0; i < 8; i++)
        r->w[i] = ((uint32_t)in[31 - 4 * i - 3] << 24) | ((uint32_t)in[31 - 4 * i - 2] << 16) |
                  ((uint32_t)in[31 - 4 * i - 1] << 8) | in[31 - 4 * i];
}

static void to_be(uint8_t out[32], const u256d *a)
{
    for (int i = 0; i < 8; i++)
    {
        out[31 - 4 * i] = (uint8_t)a->w[i];
        out[31 - 4 * i - 1] = (uint8_t)(a->w[i] >> 8);
        out[31 - 4 * i - 2] = (uint8_t)(a->w[i] >> 16);
        out[31 - 4 * i - 3] = (uint8_t)(a->w[i] >> 24);
    }
}

static int is_zero(const u256d *a)
{
    uint32_t acc = 0;
    for (int i = 0; i < 8; i++)
        acc |= a->w[i];
    return acc == 0;
}

static int equal(const u256d *a, const u256d *b)
{
    uint32_t acc = 0;
    for (int i = 0; i < 8; i++)
        acc |= a->w[i] ^ b->w[i];
    return acc == 0;
}

/* a < b as plain integers. */
static int less(const u256d *a, const u256d *b)
{
    for (int i = 7; i >= 0; i--)
        if (a->w[i] != b->w[i])
            return a->w[i] < b->w[i];
    return 0;
}

#define MUL(r, x, y) mont_mul((r), (x), (y), f)
#define ADD(r, x, y) mont_add((r), (x), (y), f)
#define SUB(r, x, y) mont_sub((r), (x), (y), f)

static void jdouble(jpoint *r, const jpoint *p, const curve *cv)
{
    const mont_ctx *f = &cv->p;
    if (is_zero(&p->z) || is_zero(&p->y))
    {
        r->z = (u256d){{0}};
        return;
    }
    u256d xx, yy, yyyy, zz, s, m, t, x3, y3, z3;
    MUL(&xx, &p->x, &p->x);
    MUL(&yy, &p->y, &p->y);
    MUL(&yyyy, &yy, &yy);
    MUL(&zz, &p->z, &p->z);

    MUL(&s, &p->x, &yy);
    ADD(&s, &s, &s);
    ADD(&s, &s, &s);

    ADD(&m, &xx, &xx);
    ADD(&m, &m, &xx);
    if (!cv->a_is_zero)
    {
        MUL(&t, &zz, &zz);
        MUL(&t, &t, &cv->a);
        ADD(&m, &m, &t);
    }

    MUL(&x3, &m, &m);
    SUB(&x3, &x3, &s);
    SUB(&x3, &x3, &s);

    SUB(&t, &s, &x3);
    MUL(&y3, &m, &t);
    ADD(&t, &yyyy, &yyyy);
    ADD(&t, &t, &t);
    ADD(&t, &t, &t);
    SUB(&y3, &y3, &t);

    MUL(&z3, &p->y, &p->z);
    ADD(&z3, &z3, &z3);

    r->x = x3;
    r->y = y3;
    r->z = z3;
}

static void jadd(jpoint *r, const jpoint *p, const jpoint *q, const curve *cv)
{
    const mont_ctx *f = &cv->p;
    if (is_zero(&p->z))
    {
        *r = *q;
        return;
    }
    if (is_zero(&q->z))
    {
        *r = *p;
        return;
    }
    u256d z1z1, z2z2, u1, u2, s1, s2, h, rr, hh, hhh, v, x3, y3, z3, t;
    MUL(&z1z1, &p->z, &p->z);
    MUL(&z2z2, &q->z, &q->z);
    MUL(&u1, &p->x, &z2z2);
    MUL(&u2, &q->x, &z1z1);
    MUL(&t, &q->z, &z2z2);
    MUL(&s1, &p->y, &t);
    MUL(&t, &p->z, &z1z1);
    MUL(&s2, &q->y, &t);

    if (equal(&u1, &u2))
    {
        if (equal(&s1, &s2))
            jdouble(r, p, cv);
        else
            r->z = (u256d){{0}};
        return;
    }

    SUB(&h, &u2, &u1);
    SUB(&rr, &s2, &s1);
    MUL(&hh, &h, &h);
    MUL(&hhh, &hh, &h);
    MUL(&v, &u1, &hh);

    MUL(&x3, &rr, &rr);
    SUB(&x3, &x3, &hhh);
    SUB(&x3, &x3, &v);
    SUB(&x3, &x3, &v);

    SUB(&t, &v, &x3);
    MUL(&y3, &rr, &t);
    MUL(&t, &s1, &hhh);
    SUB(&y3, &y3, &t);

    MUL(&z3, &p->z, &q->z);
    MUL(&z3, &z3, &h);

    r->x = x3;
    r->y = y3;
    r->z = z3;
}

/* u1*G + u2*Q, for plain (non-Montgomery) scalars. */
static void double_mul(jpoint *r, const u256d *u1, const u256d *u2, const jpoint *q, const curve *cv)
{
    jpoint g, gq, acc;
    g.x = cv->gx;
    g.y = cv->gy;
    g.z = cv->p.one;
    jadd(&gq, &g, q, cv);

    acc.z = (u256d){{0}};
    for (int i = 255; i >= 0; i--)
    {
        jdouble(&acc, &acc, cv);
        int b1 = (u1->w[i / 32] >> (i % 32)) & 1;
        int b2 = (u2->w[i / 32] >> (i % 32)) & 1;
        if (b1 && b2)
            jadd(&acc, &acc, &gq, cv);
        else if (b1)
            jadd(&acc, &acc, &g, cv);
        else if (b2)
            jadd(&acc, &acc, q, cv);
    }
    *r = acc;
}

/* Inverse by Fermat; every modulus here is prime. */
static void inverse(u256d *r, const u256d *a, const mont_ctx *m)
{
    static u256d e, two;
    e = m->n;
    two = (u256d){{2}};
    u256_op(&e, &two, U256_SUB);
    mont_pow(r, a, &e, m);
}

/* Affine coordinates, plain. */
static void to_affine(u256d *x, u256d *y, const jpoint *p, const curve *cv)
{
    const mont_ctx *f = &cv->p;
    u256d zinv, zinv2, zinv3;
    inverse(&zinv, &p->z, f);
    MUL(&zinv2, &zinv, &zinv);
    MUL(&zinv3, &zinv2, &zinv);
    MUL(x, &p->x, &zinv2);
    MUL(y, &p->y, &zinv3);
    mont_from(x, x, f);
    mont_from(y, y, f);
}

typedef struct
{
    uint8_t p[32], n[32], a[32], b[32], gx[32], gy[32];
} curve_params;

static void load_curve(curve *cv, const curve_params *cp)
{
    u256d v;
    uint32_t words[8];
    from_be(&v, cp->p);
    for (int i = 0; i < 8; i++)
        words[i] = v.w[i];
    mont_init(&cv->p, words);
    from_be(&v, cp->n);
    for (int i = 0; i < 8; i++)
        words[i] = v.w[i];
    mont_init(&cv->n, words);
    from_be(&v, cp->a);
    cv->a_is_zero = is_zero(&v);
    mont_to(&cv->a, &v, &cv->p);
    from_be(&v, cp->b);
    mont_to(&cv->b, &v, &cv->p);
    from_be(&v, cp->gx);
    mont_to(&cv->gx, &v, &cv->p);
    from_be(&v, cp->gy);
    mont_to(&cv->gy, &v, &cv->p);
}

static const curve_params K1 = {
    .p = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
          0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xfc, 0x2f},
    .n = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
          0xba, 0xae, 0xdc, 0xe6, 0xaf, 0x48, 0xa0, 0x3b, 0xbf, 0xd2, 0x5e, 0x8c, 0xd0, 0x36, 0x41, 0x41},
    .a = {0},
    .b = {[31] = 7},
    .gx = {0x79, 0xbe, 0x66, 0x7e, 0xf9, 0xdc, 0xbb, 0xac, 0x55, 0xa0, 0x62, 0x95, 0xce, 0x87, 0x0b, 0x07,
           0x02, 0x9b, 0xfc, 0xdb, 0x2d, 0xce, 0x28, 0xd9, 0x59, 0xf2, 0x81, 0x5b, 0x16, 0xf8, 0x17, 0x98},
    .gy = {0x48, 0x3a, 0xda, 0x77, 0x26, 0xa3, 0xc4, 0x65, 0x5d, 0xa4, 0xfb, 0xfc, 0x0e, 0x11, 0x08, 0xa8,
           0xfd, 0x17, 0xb4, 0x48, 0xa6, 0x85, 0x54, 0x19, 0x9c, 0x47, 0xd0, 0x8f, 0xfb, 0x10, 0xd4, 0xb8},
};

static const curve_params R1 = {
    .p = {0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
          0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
    .n = {0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
          0xbc, 0xe6, 0xfa, 0xad, 0xa7, 0x17, 0x9e, 0x84, 0xf3, 0xb9, 0xca, 0xc2, 0xfc, 0x63, 0x25, 0x51},
    .a = {0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
          0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfc},
    .b = {0x5a, 0xc6, 0x35, 0xd8, 0xaa, 0x3a, 0x93, 0xe7, 0xb3, 0xeb, 0xbd, 0x55, 0x76, 0x98, 0x86, 0xbc,
          0x65, 0x1d, 0x06, 0xb0, 0xcc, 0x53, 0xb0, 0xf6, 0x3b, 0xce, 0x3c, 0x3e, 0x27, 0xd2, 0x60, 0x4b},
    .gx = {0x6b, 0x17, 0xd1, 0xf2, 0xe1, 0x2c, 0x42, 0x47, 0xf8, 0xbc, 0xe6, 0xe5, 0x63, 0xa4, 0x40, 0xf2,
           0x77, 0x03, 0x7d, 0x81, 0x2d, 0xeb, 0x33, 0xa0, 0xf4, 0xa1, 0x39, 0x45, 0xd8, 0x98, 0xc2, 0x96},
    .gy = {0x4f, 0xe3, 0x42, 0xe2, 0xfe, 0x1a, 0x7f, 0x9b, 0x8e, 0xe7, 0xeb, 0x4a, 0x7c, 0x0f, 0x9e, 0x16,
           0x2b, 0xce, 0x33, 0x57, 0x6b, 0x31, 0x5e, 0xce, 0xcb, 0xb6, 0x40, 0x68, 0x37, 0xbf, 0x51, 0xf5},
};

/* The curves, set up once: a context costs a few hundred delegations. */
static curve CurveK1, CurveR1;
static int CurveK1Ready, CurveR1Ready;

static const curve *get_curve(curve *cv, int *ready, const curve_params *cp)
{
    if (!*ready)
    {
        load_curve(cv, cp);
        *ready = 1;
    }
    return cv;
}

static int in_scalar_range(const u256d *v, const curve *cv)
{
    return !is_zero(v) && less(v, &cv->n.n);
}

/* The message hash as a scalar, mod n, in the Montgomery domain of n. */
static void hash_to_scalar(u256d *e, const uint8_t msg[32], const curve *cv)
{
    from_be(e, msg);
    if (!less(e, &cv->n.n))
        u256_op(e, &cv->n.n, U256_SUB);
    mont_to(e, e, &cv->n);
}

static int on_curve(const u256d *x, const u256d *y, const curve *cv)
{
    const mont_ctx *f = &cv->p;
    u256d lhs, rhs, t;
    MUL(&lhs, y, y);
    MUL(&rhs, x, x);
    MUL(&rhs, &rhs, x);
    MUL(&t, &cv->a, x);
    ADD(&rhs, &rhs, &t);
    ADD(&rhs, &rhs, &cv->b);
    return equal(&lhs, &rhs);
}

zkvm_status zkvm_secp256k1_ecrecover(const uint8_t msg[32], const uint8_t sig[64], uint8_t recid,
                                     uint8_t output[64])
{
    const curve *cv = get_curve(&CurveK1, &CurveK1Ready, &K1);
    const mont_ctx *f = &cv->p;
    const mont_ctx *n = &cv->n;

    u256d r, s;
    from_be(&r, sig);
    from_be(&s, sig + 32);
    if (recid > 3 || !in_scalar_range(&r, cv) || !in_scalar_range(&s, cv))
        return ZKVM_EFAIL;

    u256d x = r;
    if (recid & 2)
    {
        u256d nn = n->n;
        if (u256_op(&x, &nn, U256_ADD))
            return ZKVM_EFAIL;
    }
    if (!less(&x, &f->n))
        return ZKVM_EFAIL;

    /* y^2 = x^3 + 7; p = 3 mod 4, so a root is alpha^((p+1)/4). */
    u256d xm, alpha, beta, check, e, one = {{1}};
    mont_to(&xm, &x, f);
    MUL(&alpha, &xm, &xm);
    MUL(&alpha, &alpha, &xm);
    ADD(&alpha, &alpha, &cv->b);
    e = f->n;
    u256_op(&e, &one, U256_ADD);
    for (int i = 0; i < 8; i++)
        e.w[i] = (e.w[i] >> 2) | (i < 7 ? e.w[i + 1] << 30 : 0);
    mont_pow(&beta, &alpha, &e, f);
    MUL(&check, &beta, &beta);
    if (!equal(&check, &alpha))
        return ZKVM_EFAIL;
    u256d plain;
    mont_from(&plain, &beta, f);
    if ((plain.w[0] & 1) != (recid & 1))
    {
        u256d zero = {{0}};
        SUB(&beta, &zero, &beta);
    }

    jpoint rp;
    rp.x = xm;
    rp.y = beta;
    rp.z = f->one;

    /* Q = r^-1 (s R - e G) */
    u256d rm, sm, rinv, u1, u2, zero = {{0}};
    hash_to_scalar(&e, msg, cv);
    mont_to(&rm, &r, n);
    mont_to(&sm, &s, n);
    inverse(&rinv, &rm, n);
    mont_sub(&u1, &zero, &e, n);
    mont_mul(&u1, &u1, &rinv, n);
    mont_mul(&u2, &sm, &rinv, n);
    mont_from(&u1, &u1, n);
    mont_from(&u2, &u2, n);

    jpoint q;
    double_mul(&q, &u1, &u2, &rp, cv);
    if (is_zero(&q.z))
        return ZKVM_EFAIL;

    u256d qx, qy;
    to_affine(&qx, &qy, &q, cv);
    to_be(output, &qx);
    to_be(output + 32, &qy);
    return ZKVM_EOK;
}

zkvm_status zkvm_secp256r1_verify(const uint8_t msg[32], const uint8_t sig[64], const uint8_t pubkey[64],
                                  bool *verified)
{
    const curve *cv = get_curve(&CurveR1, &CurveR1Ready, &R1);
    const mont_ctx *n = &cv->n;
    *verified = false;

    u256d r, s, qx, qy;
    from_be(&r, sig);
    from_be(&s, sig + 32);
    from_be(&qx, pubkey);
    from_be(&qy, pubkey + 32);
    if (!in_scalar_range(&r, cv) || !in_scalar_range(&s, cv))
        return ZKVM_EOK;
    if (!less(&qx, &cv->p.n) || !less(&qy, &cv->p.n) || (is_zero(&qx) && is_zero(&qy)))
        return ZKVM_EOK;

    jpoint q;
    mont_to(&q.x, &qx, &cv->p);
    mont_to(&q.y, &qy, &cv->p);
    q.z = cv->p.one;
    if (!on_curve(&q.x, &q.y, cv))
        return ZKVM_EOK;

    u256d e, rm, sm, w, u1, u2;
    hash_to_scalar(&e, msg, cv);
    mont_to(&rm, &r, n);
    mont_to(&sm, &s, n);
    inverse(&w, &sm, n);
    mont_mul(&u1, &e, &w, n);
    mont_mul(&u2, &rm, &w, n);
    mont_from(&u1, &u1, n);
    mont_from(&u2, &u2, n);

    jpoint rp;
    double_mul(&rp, &u1, &u2, &q, cv);
    if (is_zero(&rp.z))
        return ZKVM_EOK;

    u256d x, y;
    to_affine(&x, &y, &rp, cv);
    if (!less(&x, &n->n))
        u256_op(&x, &n->n, U256_SUB);
    *verified = equal(&x, &r);
    return ZKVM_EOK;
}
