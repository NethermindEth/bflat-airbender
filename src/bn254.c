/*
 * BN254 (alt_bn128) accelerators: EIP-196 point addition and scalar multiplication on G1 and the
 * EIP-197 pairing check, on the shared arithmetic of ec.c.
 *
 * Encodings follow the precompiles, which hand their input through unchanged: a G1 point is
 * x || y, a G2 point x || y with each Fp2 element imaginary part first (c1 || c0), every field
 * element 32 bytes big-endian, and all zeros the point at infinity. ZKVM_EFAIL means invalid
 * input, which fails the precompile: a coordinate not below p, a point off its curve, or a G2
 * point outside the order-r subgroup (G1 is the whole curve). A scalar is any 256-bit value.
 */

#include "ec.h"

/* In RAM, where the delegation reads its operands. */
static curve Bn254 = {
    .n = 8,
    .inv = 0xe4866389,
    .p = {{0xd87cfd47, 0x3c208c16, 0x6871ca8d, 0x97816a91, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72}},
    .one = {{0xc58f0d9d, 0xd35d438d, 0xf5c70b3d, 0x0a78eb28, 0x7879462c, 0x666ea36f, 0x9a07df2f, 0x0e0a77c1}},
    .r2 = {{0x538afa89, 0xf32cfc5b, 0xd44501fb, 0xb5e71911, 0x0a417ff6, 0x47ab1eff, 0xcab8351f, 0x06d89f71}},
    .b1 = {{{{0x50ad28d7, 0x7a17caa9, 0xe15521b9, 0x1f6ac17a, 0x696bd284, 0x334bea4e, 0xce179d8e, 0x2a1f6744}}, {{0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000}}}},
    .b2 = {{{{0x77b802a8, 0x3bf938e3, 0x3633535d, 0x020b1b27, 0x49755260, 0x26b7edf0, 0x4384a86d, 0x2514c632}}, {{0xd1dcff67, 0x38e7eccc, 0x93ce0d3e, 0x65f0b37d, 0x22ac00aa, 0xd749d0dd, 0x4a688d4d, 0x0141b9ce}}}},
    .gamma = {{{{0x33144907, 0xaf9ba696, 0x87afb78a, 0xca6b1d73, 0xf08a2087, 0x11bded5e, 0x1a1f3a7c, 0x02f34d75}}, {{0x4c492d72, 0xa222ae23, 0x565de15b, 0xd00f02a4, 0x53dfc926, 0xdc2ff3a2, 0xb3899551, 0x10a75716}}}},
    .xi = 9,
    .mtwist = false,
    .x = 4965661367192848881u,
    .xneg = false,
    .loop = {0x01, 0x9d, 0x79, 0x70, 0x39, 0xbe, 0x76, 0x3b, 0xa8},
};

/*
 * GLV on G1 (Gallant, Lambert and Vanstone): (x, y) -> (beta x, y), beta a cube root of 1 mod p, is multiplication by
 * lambda, a cube root of 1 mod r, so k P = k1 P + k2 (beta x, y). The short basis (A1, -B), (B, B2) of the lattice of
 * (a, b) with a + b lambda = 0 mod r gives k1 = k - c1 A1 - c2 B and k2 = c1 B - c2 B2 with c1 = floor(k G1 / 2^256)
 * and c2 = floor(k G2 / 2^256), G1 and G2 the basis scaled by 2^256 / r: any c1 and c2 keep k1 + k2 lambda = k mod r,
 * and these leave both halves at most 128 bits for any 256-bit k. Little-endian 32-bit limbs; Beta in Montgomery form.
 */
static fp Beta = {{0x13e80b9c, 0x3350c88e, 0xdb5e56b9, 0x7dce557c, 0xb615564a, 0x6001b4b8, 0x020217e0, 0x2682e617}};
static const uint32_t GlvA1[4] = {0x7d4f1128, 0x8211bbeb, 0xeeb859fc, 0x6f4d8248};
static const uint32_t GlvB[2] = {0x94d213e3, 0x89d32568};
static const uint32_t GlvB2[4] = {0x1221250b, 0x0be4e154, 0xeeb859fd, 0x6f4d8248};
static const uint32_t GlvG1[5] = {0x00ff6565, 0x5398fd03, 0xa773d2d2, 0x4ccef014, 0x00000002};
static const uint32_t GlvG2[3] = {0xc7e0b3d7, 0xd91d232e, 0x00000002};

enum
{
    GlvLimbs = 9 /* k and the halves, signed, mod 2^288 */
};

/* r = a b, na + nb limbs. */
static void mp_mul(uint32_t *r, const uint32_t *a, int na, const uint32_t *b, int nb)
{
    memset(r, 0, 4 * (size_t)(na + nb));
    for (int i = 0; i < na; i++)
    {
        uint64_t carry = 0;
        for (int j = 0; j < nb; j++)
        {
            carry += (uint64_t)a[i] * b[j] + r[i + j];
            r[i + j] = (uint32_t)carry;
            carry >>= 32;
        }
        r[i + nb] = (uint32_t)carry;
    }
}

/* r -= a mod 2^(32 GlvLimbs), a of na limbs. */
static void mp_sub(uint32_t *r, const uint32_t *a, int na)
{
    uint32_t borrow = 0;
    for (int i = 0; i < GlvLimbs; i++)
    {
        uint32_t ai = i < na ? a[i] : 0, d = r[i] - ai - borrow;
        borrow = r[i] < ai || (r[i] == ai && borrow);
        r[i] = d;
    }
}

/* |v| as 20 big-endian bytes; true when v is negative. */
static bool mp_abs_be(uint8_t out[20], uint32_t *v)
{
    bool neg = v[GlvLimbs - 1] >> 31;
    if (neg)
    {
        uint32_t m[GlvLimbs];
        memcpy(m, v, sizeof(m));
        memset(v, 0, sizeof(m));
        mp_sub(v, m, GlvLimbs);
    }
    for (int i = 0; i < 20; i++)
        out[19 - i] = (uint8_t)(v[i / 4] >> (8 * (i % 4)));
    return neg;
}

static void g1_mul_glv(point *r, const point *p, const uint8_t scalar[32])
{
    uint32_t k[8], prod[13], c1[5], c2[3], k1[GlvLimbs] = {0}, k2[GlvLimbs] = {0}, t[9];
    for (int i = 0; i < 8; i++)
    {
        const uint8_t *b = scalar + 4 * (7 - i);
        k[i] = ((uint32_t)b[0] << 24) | ((uint32_t)b[1] << 16) | ((uint32_t)b[2] << 8) | b[3];
    }
    mp_mul(prod, k, 8, GlvG1, 5);
    memcpy(c1, prod + 8, sizeof(c1));
    mp_mul(prod, k, 8, GlvG2, 3);
    memcpy(c2, prod + 8, sizeof(c2));

    memcpy(k1, k, sizeof(k));
    mp_mul(t, c1, 5, GlvA1, 4);
    mp_sub(k1, t, 9);
    mp_mul(t, c2, 3, GlvB, 2);
    mp_sub(k1, t, 5);
    mp_mul(t, c1, 5, GlvB, 2);
    memcpy(k2, t, 4 * 7);
    mp_mul(t, c2, 3, GlvB2, 4);
    mp_sub(k2, t, 7);

    uint8_t e1[20], e2[20];
    point p1 = *p, p2 = *p;
    fp_mul(&p2.x.c[0], &p->x.c[0], &Beta);
    if (mp_abs_be(e1, k1))
        pt_neg(&p1, &p1, 1);
    if (mp_abs_be(e2, k2))
        pt_neg(&p2, &p2, 1);
    pt_mul2(r, &p1, e1, &p2, e2, sizeof(e1), 1);
}

static bool load_g1(point *p, const uint8_t in[64])
{
    return pt_load(p, in, 1, false) && pt_on_curve(p, 1);
}

/*
 * G2 membership without multiplying by r: q is in G2 iff [x+1]q + psi([x]q) + psi^2([x]q) ==
 * psi^3([2x]q), psi the Frobenius on the twist (Scott; El Housni, Guillevic and Piellard, ePrint
 * 2022/348 section 5.1, as gnark checks it).
 */
static bool load_g2(point *q, const uint8_t in[128])
{
    if (!pt_load(q, in, 2, true) || !pt_on_curve(q, 2))
        return false;
    point a, t, lhs, rhs;
    pt_mul_x(&a, q, 2);
    pt_add(&lhs, &a, q, 2);
    pt_frob(&t, &a, false);
    pt_add(&lhs, &lhs, &t, 2);
    pt_frob(&t, &t, false);
    pt_add(&lhs, &lhs, &t, 2);
    pt_dbl(&rhs, &a, 2);
    for (int i = 0; i < 3; i++)
        pt_frob(&rhs, &rhs, false);
    return pt_eq(&lhs, &rhs, 2);
}

zkvm_status zkvm_bn254_g1_add(const uint8_t p1[64], const uint8_t p2[64], uint8_t result[64])
{
    ec_use(&Bn254);
    point a, b;
    if (!load_g1(&a, p1) || !load_g1(&b, p2))
        return ZKVM_EFAIL;
    pt_add(&a, &a, &b, 1);
    pt_store(result, &a, 1, false);
    return ZKVM_EOK;
}

zkvm_status zkvm_bn254_g1_mul(const uint8_t point_in[64], const uint8_t scalar[32], uint8_t result[64])
{
    ec_use(&Bn254);
    point a;
    if (!load_g1(&a, point_in))
        return ZKVM_EFAIL;
    g1_mul_glv(&a, &a, scalar);
    pt_store(result, &a, 1, false);
    return ZKVM_EOK;
}

/* Each pair is a G1 point (64 bytes) then a G2 point (128 bytes). No pairs is a passing check. */
zkvm_status zkvm_bn254_pairing(const uint8_t *pairs, size_t num_pairs, bool *verified)
{
    ec_use(&Bn254);
    fp12 f;
    fp12_one(&f);
    *verified = false;
    for (size_t i = 0; i < num_pairs; i += PairingBatch)
    {
        point p[PairingBatch], q[PairingBatch];
        int count = num_pairs - i < PairingBatch ? (int)(num_pairs - i) : PairingBatch;
        for (int j = 0; j < count; j++)
        {
            const uint8_t *in = pairs + 192 * (i + j);
            if (!load_g1(&p[j], in) || !load_g2(&q[j], in + 64))
                return ZKVM_EFAIL;
        }
        pairing_add_batch(&f, p, q, count);
    }
    *verified = pairing_is_one(&f);
    return ZKVM_EOK;
}
