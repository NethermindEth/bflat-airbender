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
    pt_mul(&a, &a, scalar, 32, 1);
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
