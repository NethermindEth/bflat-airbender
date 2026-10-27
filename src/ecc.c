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

/* ---- secp256k1: GLV, wNAF and mixed additions ---------------------------------------------- */
/*
 * Recovery computes u1*G + u2*R. With the endomorphism lambda*(x, y) = (beta*x, y), each scalar
 * splits into two of about 128 bits (libsecp256k1's split_lambda), so one doubling chain of ~129
 * steps serves four terms. Their wNAF digits add odd multiples from tables: window 8 against G and
 * lambda*G, built once and kept affine so each addition is a mixed one; window 5 against R and
 * lambda*R, built per call. The square root and the field inverse use libsecp256k1's addition
 * chains. Every operand the delegation reads lives in RAM, so the constants are copied there.
 */

typedef struct
{
    u256d x, y; /* Montgomery form */
} apoint;

enum
{
    WindowG = 8,
    EntriesG = 1 << (WindowG - 2),
    WindowR = 5,
    EntriesR = 1 << (WindowR - 2),
    MaxDigits = 131 /* the parts are below 2^129; one more position takes the last carry */
};

static const uint32_t GlvG1[8] = {0x45dbb031, 0xe893209a, 0x71e8ca7f, 0x3daa8a14, 0x9284eb15, 0xe86c90e4, 0xa7d46bcd, 0x3086d221};
static const uint32_t GlvG2[8] = {0x8ac47f71, 0x1571b4ae, 0x9df506c6, 0x221208ac, 0x0abfe4c4, 0x6f547fa9, 0x010e8828, 0xe4437ed6};
static const uint32_t MinusB1[8] = {0x0abfe4c3, 0x6f547fa9, 0x010e8828, 0xe4437ed6, 0, 0, 0, 0};
static const uint32_t MinusB2[8] = {0x3db1562c, 0xd765cda8, 0x0774346d, 0x8a280ac5, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff};
static const uint32_t MinusLambda[8] = {0xb51283cf, 0xe0cfc810, 0x8ec739c2, 0xa880b9fc, 0x77ed9ba4, 0x5ad9e3fd, 0x3fa3cf1f, 0xac9c52b3};
static const uint32_t Beta[8] = {0x719501ee, 0xc1396c28, 0x12f58995, 0x9cf04975, 0xac3434e9, 0x6e64479e, 0x657c0710, 0x7ae96a2b};
static const uint32_t HalfN[8] = {0x681b20a0, 0xdfe92f46, 0x57a4501d, 0x5d576e73, 0xffffffff, 0xffffffff, 0xffffffff, 0x7fffffff};

static struct
{
    int ready;
    u256d g1, g2, half_n;                       /* plain */
    u256d minus_b1, minus_b2, minus_lambda;     /* Montgomery form mod n */
    u256d beta;                                 /* Montgomery form mod p */
    apoint g[EntriesG], lambda_g[EntriesG];     /* odd multiples of G and lambda*G */
} Glv;

static void load_words(u256d *r, const uint32_t w[8])
{
    for (int i = 0; i < 8; i++)
        r->w[i] = w[i];
}

/* r = 2p for a = 0 (dbl-2009-l): 2M + 5S. r may be p. */
static void jdouble0(jpoint *r, const jpoint *p, const mont_ctx *f)
{
    if (is_zero(&p->z) || is_zero(&p->y))
    {
        r->z = (u256d){{0}};
        return;
    }
    u256d a, b, c, d, e, t, x3;
    MUL(&a, &p->x, &p->x);
    MUL(&b, &p->y, &p->y);
    MUL(&c, &b, &b);
    ADD(&d, &p->x, &b);
    MUL(&d, &d, &d);
    SUB(&d, &d, &a);
    SUB(&d, &d, &c);
    ADD(&d, &d, &d);
    ADD(&e, &a, &a);
    ADD(&e, &e, &a);
    MUL(&t, &e, &e);
    SUB(&x3, &t, &d);
    SUB(&x3, &x3, &d);
    MUL(&r->z, &p->y, &p->z);
    ADD(&r->z, &r->z, &r->z);
    SUB(&t, &d, &x3);
    MUL(&t, &e, &t);
    ADD(&c, &c, &c);
    ADD(&c, &c, &c);
    ADD(&c, &c, &c);
    SUB(&r->y, &t, &c);
    r->x = x3;
}

/* r = p + q for p not at infinity and q affine (madd-2007-bl): 7M + 4S. r may be p. */
static void jadd_affine(jpoint *r, const jpoint *p, const apoint *q, const mont_ctx *f)
{
    u256d z1z1, u2, s2, h, hh, i, j, rr, v, x3;
    MUL(&z1z1, &p->z, &p->z);
    MUL(&u2, &q->x, &z1z1);
    MUL(&s2, &q->y, &p->z);
    MUL(&s2, &s2, &z1z1);
    SUB(&h, &u2, &p->x);
    SUB(&rr, &s2, &p->y);
    if (is_zero(&h))
    {
        if (is_zero(&rr))
            jdouble0(r, p, f);
        else
            r->z = (u256d){{0}};
        return;
    }
    MUL(&hh, &h, &h);
    ADD(&i, &hh, &hh);
    ADD(&i, &i, &i);
    MUL(&j, &h, &i);
    ADD(&rr, &rr, &rr);
    MUL(&v, &p->x, &i);
    MUL(&x3, &rr, &rr);
    SUB(&x3, &x3, &j);
    SUB(&x3, &x3, &v);
    SUB(&x3, &x3, &v);
    ADD(&u2, &p->z, &h);
    MUL(&u2, &u2, &u2);
    SUB(&u2, &u2, &z1z1);
    SUB(&u2, &u2, &hh);
    SUB(&v, &v, &x3);
    MUL(&v, &rr, &v);
    MUL(&s2, &p->y, &j);
    ADD(&s2, &s2, &s2);
    SUB(&r->y, &v, &s2);
    r->x = x3;
    r->z = u2;
}

/* r = a^(2^k). */
static void sqr_n(u256d *r, const u256d *a, int k, const mont_ctx *f)
{
    if (r != a)
        *r = *a;
    while (k-- > 0)
        MUL(r, r, r);
}

/* x223 = a^(2^223 - 1), x22 and x2 alongside: the common prefix of libsecp256k1's field chains. */
static void chain_x223(u256d *x223, u256d *x22, u256d *x2, const u256d *a, const mont_ctx *f)
{
    u256d x3, x6, x9, x11, x44, x88, x176, x220;
    MUL(x2, a, a);
    MUL(x2, x2, a);
    MUL(&x3, x2, x2);
    MUL(&x3, &x3, a);
    sqr_n(&x6, &x3, 3, f);
    MUL(&x6, &x6, &x3);
    sqr_n(&x9, &x6, 3, f);
    MUL(&x9, &x9, &x3);
    sqr_n(&x11, &x9, 2, f);
    MUL(&x11, &x11, x2);
    sqr_n(x22, &x11, 11, f);
    MUL(x22, x22, &x11);
    sqr_n(&x44, x22, 22, f);
    MUL(&x44, &x44, x22);
    sqr_n(&x88, &x44, 44, f);
    MUL(&x88, &x88, &x44);
    sqr_n(&x176, &x88, 88, f);
    MUL(&x176, &x176, &x88);
    sqr_n(&x220, &x176, 44, f);
    MUL(&x220, &x220, &x44);
    sqr_n(x223, &x220, 3, f);
    MUL(x223, x223, &x3);
}

/* r = a^((p+1)/4), a square root of a when one exists. */
static void field_sqrt(u256d *r, const u256d *a, const mont_ctx *f)
{
    u256d x223, x22, x2;
    chain_x223(&x223, &x22, &x2, a, f);
    sqr_n(r, &x223, 23, f);
    MUL(r, r, &x22);
    sqr_n(r, r, 6, f);
    MUL(r, r, &x2);
    sqr_n(r, r, 2, f);
}

/* r = a^(p-2) = a^-1. */
static void field_inverse(u256d *r, const u256d *a, const mont_ctx *f)
{
    u256d x223, x22, x2;
    chain_x223(&x223, &x22, &x2, a, f);
    sqr_n(r, &x223, 23, f);
    MUL(r, r, &x22);
    sqr_n(r, r, 5, f);
    MUL(r, r, a);
    sqr_n(r, r, 3, f);
    MUL(r, r, &x2);
    sqr_n(r, r, 2, f);
    MUL(r, r, a);
}

/* The odd multiples of G (1G, 3G, ..., (2 EntriesG - 1)G) in affine form, by one inversion. */
static void glv_init(const curve *cv)
{
    const mont_ctx *f = &cv->p;
    const mont_ctx *n = &cv->n;
    static jpoint pts[EntriesG];
    static u256d prefix[EntriesG];
    u256d u, inv, zinv, zinv2, zinv3;
    jpoint twice;

    load_words(&Glv.g1, GlvG1);
    load_words(&Glv.g2, GlvG2);
    load_words(&Glv.half_n, HalfN);
    load_words(&u, MinusB1);
    mont_to(&Glv.minus_b1, &u, n);
    load_words(&u, MinusB2);
    mont_to(&Glv.minus_b2, &u, n);
    load_words(&u, MinusLambda);
    mont_to(&Glv.minus_lambda, &u, n);
    load_words(&u, Beta);
    mont_to(&Glv.beta, &u, f);

    pts[0].x = cv->gx;
    pts[0].y = cv->gy;
    pts[0].z = f->one;
    jdouble0(&twice, &pts[0], f);
    for (int i = 1; i < EntriesG; i++)
        jadd(&pts[i], &pts[i - 1], &twice, cv);

    prefix[0] = pts[0].z;
    for (int i = 1; i < EntriesG; i++)
        MUL(&prefix[i], &prefix[i - 1], &pts[i].z);
    field_inverse(&inv, &prefix[EntriesG - 1], f);
    for (int i = EntriesG - 1; i >= 0; i--)
    {
        if (i > 0)
        {
            MUL(&zinv, &inv, &prefix[i - 1]);
            MUL(&inv, &inv, &pts[i].z);
        }
        else
            zinv = inv;
        MUL(&zinv2, &zinv, &zinv);
        MUL(&zinv3, &zinv2, &zinv);
        MUL(&Glv.g[i].x, &pts[i].x, &zinv2);
        MUL(&Glv.g[i].y, &pts[i].y, &zinv3);
        MUL(&Glv.lambda_g[i].x, &Glv.g[i].x, &Glv.beta);
        Glv.lambda_g[i].y = Glv.g[i].y;
    }
    Glv.ready = 1;
}

/* r = round(a * b / 2^384), the GLV rounding step. */
static void mul_shift_384(u256d *r, const u256d *a, const u256d *b)
{
    static u256d hi;
    u256_op(&hi, a, U256_MEMCOPY);
    u256_op(&hi, b, U256_MUL_HIGH);
    uint64_t carry = hi.w[3] >> 31;
    for (int i = 0; i < 8; i++)
    {
        carry += i < 4 ? hi.w[i + 4] : 0;
        r->w[i] = (uint32_t)carry;
        carry >>= 32;
    }
}

/* k = k1 + k2 * lambda (mod n), each part returned as a magnitude below 2^129 and a sign. */
static void split_lambda(u256d *k1, int *neg1, u256d *k2, int *neg2, const u256d *k, const mont_ctx *n)
{
    u256d c1, c2, t;
    mul_shift_384(&c1, k, &Glv.g1);
    mul_shift_384(&c2, k, &Glv.g2);
    mont_mul(&c1, &c1, &Glv.minus_b1, n);
    mont_mul(&c2, &c2, &Glv.minus_b2, n);
    mont_add(k2, &c1, &c2, n);
    mont_mul(k1, k2, &Glv.minus_lambda, n);
    mont_add(k1, k1, k, n);

    *neg1 = less(&Glv.half_n, k1);
    if (*neg1)
    {
        t = n->n;
        u256_op(&t, k1, U256_SUB);
        *k1 = t;
    }
    *neg2 = less(&Glv.half_n, k2);
    if (*neg2)
    {
        t = n->n;
        u256_op(&t, k2, U256_SUB);
        *k2 = t;
    }
}

/* The width-w NAF of s: odd digits in (-2^(w-1), 2^(w-1)), at least w - 1 zeros between two. */
static void wnaf(int8_t digits[MaxDigits], const u256d *s, int w)
{
    memset(digits, 0, MaxDigits);
    int bit = 0;
    uint32_t carry = 0;
    while (bit < MaxDigits)
    {
        if (((s->w[bit >> 5] >> (bit & 31)) & 1) == carry)
        {
            bit++;
            continue;
        }
        int now = MaxDigits - bit < w ? MaxDigits - bit : w;
        int limb = bit >> 5;
        int offset = bit & 31;
        uint32_t bits = s->w[limb] >> offset;
        if (offset + now > 32)
            bits |= s->w[limb + 1] << (32 - offset);
        int32_t word = (int32_t)(bits & ((1u << now) - 1)) + (int32_t)carry;
        carry = ((uint32_t)word >> (w - 1)) & 1;
        word -= (int32_t)(carry << w);
        digits[bit] = (int8_t)word;
        bit += now;
    }
}

/* acc += the table entry of digit (negated when the signs disagree); *inf tracks acc at infinity. */
static void add_affine_term(jpoint *acc, int *inf, const apoint *table, int digit, int negate, const mont_ctx *f)
{
    if (digit < 0)
    {
        digit = -digit;
        negate = !negate;
    }
    apoint q = table[digit >> 1];
    if (negate)
    {
        u256d zero = {{0}};
        SUB(&q.y, &zero, &q.y);
    }
    if (*inf)
    {
        acc->x = q.x;
        acc->y = q.y;
        acc->z = f->one;
        *inf = 0;
        return;
    }
    jadd_affine(acc, acc, &q, f);
    *inf = is_zero(&acc->z);
}

static void add_jacobian_term(jpoint *acc, int *inf, const jpoint *table, int digit, int negate, const curve *cv)
{
    const mont_ctx *f = &cv->p;
    if (digit < 0)
    {
        digit = -digit;
        negate = !negate;
    }
    jpoint q = table[digit >> 1];
    if (negate)
    {
        u256d zero = {{0}};
        SUB(&q.y, &zero, &q.y);
    }
    if (*inf)
    {
        *acc = q;
        *inf = is_zero(&acc->z);
        return;
    }
    jadd(acc, acc, &q, cv);
    *inf = is_zero(&acc->z);
}

/* u1*G + u2*R for plain scalars below n and R affine (z = 1). */
static void glv_double_mul(jpoint *r, const u256d *u1, const u256d *u2, const jpoint *rp, const curve *cv)
{
    const mont_ctx *f = &cv->p;
    if (!Glv.ready)
        glv_init(cv);

    u256d a1, a2, b1, b2;
    int na1, na2, nb1, nb2;
    split_lambda(&a1, &na1, &a2, &na2, u1, &cv->n);
    split_lambda(&b1, &nb1, &b2, &nb2, u2, &cv->n);
    int8_t d[4][MaxDigits];
    wnaf(d[0], &a1, WindowG);
    wnaf(d[1], &a2, WindowG);
    wnaf(d[2], &b1, WindowR);
    wnaf(d[3], &b2, WindowR);

    jpoint tr[EntriesR], tlr[EntriesR], twice;
    tr[0] = *rp;
    jdouble0(&twice, rp, f);
    for (int i = 1; i < EntriesR; i++)
        jadd(&tr[i], &tr[i - 1], &twice, cv);
    for (int i = 0; i < EntriesR; i++)
    {
        tlr[i] = tr[i];
        MUL(&tlr[i].x, &tr[i].x, &Glv.beta);
    }

    int top = MaxDigits - 1;
    while (top >= 0 && !d[0][top] && !d[1][top] && !d[2][top] && !d[3][top])
        top--;
    jpoint acc;
    int inf = 1;
    for (int i = top; i >= 0; i--)
    {
        if (!inf)
        {
            jdouble0(&acc, &acc, f);
            inf = is_zero(&acc.z);
        }
        if (d[0][i])
            add_affine_term(&acc, &inf, Glv.g, d[0][i], na1, f);
        if (d[1][i])
            add_affine_term(&acc, &inf, Glv.lambda_g, d[1][i], na2, f);
        if (d[2][i])
            add_jacobian_term(&acc, &inf, tr, d[2][i], nb1, cv);
        if (d[3][i])
            add_jacobian_term(&acc, &inf, tlr, d[3][i], nb2, cv);
    }
    if (inf)
        acc.z = (u256d){{0}};
    *r = acc;
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
    u256d xm, alpha, beta, check, e;
    mont_to(&xm, &x, f);
    MUL(&alpha, &xm, &xm);
    MUL(&alpha, &alpha, &xm);
    ADD(&alpha, &alpha, &cv->b);
    field_sqrt(&beta, &alpha, f);
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
    glv_double_mul(&q, &u1, &u2, &rp, cv);
    if (is_zero(&q.z))
        return ZKVM_EFAIL;

    u256d qx, qy, zinv, zinv2, zinv3;
    field_inverse(&zinv, &q.z, f);
    MUL(&zinv2, &zinv, &zinv);
    MUL(&zinv3, &zinv2, &zinv);
    MUL(&qx, &q.x, &zinv2);
    MUL(&qy, &q.y, &zinv3);
    mont_from(&qx, &qx, f);
    mont_from(&qy, &qy, f);
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
