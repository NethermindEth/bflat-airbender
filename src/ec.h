/*
 * Shared elliptic-curve arithmetic of the BN254 and BLS12-381 accelerators.
 *
 * Both curves are y^2 = x^3 + b over a prime field with a pairing-friendly sextic twist, so one
 * implementation serves both, parameterised by the curve table in use (ec_use):
 *   - Fp: Montgomery arithmetic on 32-bit limbs, 8 limbs for BN254 and 12 for BLS12-381;
 *   - Fp2 = Fp[u]/(u^2 + 1), the same for both;
 *   - Fp12 = Fp2[w]/(w^6 - xi), flat rather than towered, with xi = 9 + u (BN254) or 1 + u
 *     (BLS12-381), so one schoolbook product covers both;
 *   - Jacobian points over Fp (G1, degree 1) or Fp2 (G2, degree 2) by the same routines;
 *   - the optimal ate pairing check: a Miller loop per pair and one final exponentiation.
 * Everything favours code size over speed, and nothing is constant time: a zkVM guest has no
 * side channel. The tables hold Montgomery-form constants generated from the curve parameters.
 */
#ifndef AIRBENDER_EC_H
#define AIRBENDER_EC_H

#include "zkvm.h"

enum
{
    FpMaxLimbs = 12
};

typedef struct
{
    uint32_t l[FpMaxLimbs]; /* little-endian limbs, Montgomery form; limbs past n are zero */
} fp;

typedef struct
{
    fp c[2]; /* c[0] + c[1]*u; a degree-1 element uses c[0] only */
} fp2;

typedef struct
{
    fp2 c[6]; /* sum of c[i] * w^i */
} fp12;

typedef struct
{
    fp2 x, y, z; /* Jacobian; z == 0 is the point at infinity */
} point;

typedef struct
{
    uint32_t n;        /* limbs in use */
    uint32_t inv;      /* -p^-1 mod 2^32 */
    fp p;              /* the field prime, plain */
    fp one;            /* R mod p, R = 2^(32n): 1 in Montgomery form */
    fp r2;             /* R^2 mod p, which converts into Montgomery form */
    fp2 b1;            /* b of G1: y^2 = x^3 + b */
    fp2 b2;            /* b' of the twist carrying G2 */
    fp2 gamma;         /* xi^((p-1)/6), so that w^p = gamma * w */
    uint32_t xi;       /* xi = xi + u */
    bool mtwist;       /* M-type twist (BLS12-381); BN254 is D-type */
    uint64_t x;        /* |x|, the curve family parameter */
    bool xneg;         /* x < 0 (BLS12-381) */
    uint8_t loop[9];   /* the Miller loop scalar, big-endian: 6x + 2 (BN254) or |x| */
} curve;

extern const curve *C;

static inline void ec_use(const curve *cv)
{
    C = cv;
}

void fp_add(fp *r, const fp *a, const fp *b);
void fp_sub(fp *r, const fp *a, const fp *b);
void fp_mul(fp *r, const fp *a, const fp *b);
bool fp_is_zero(const fp *a);
bool fp_from_be(fp *r, const uint8_t *in);
void fp_to_be(uint8_t *out, const fp *a);
void fp_from_mont(fp *r, const fp *a);

/* Degree-d (1 or 2) field operations: Fp on c[0] alone, or Fp2. */
void fq_add(fp2 *r, const fp2 *a, const fp2 *b, int d);
void fq_sub(fp2 *r, const fp2 *a, const fp2 *b, int d);
void fq_neg(fp2 *r, const fp2 *a, int d);
void fq_mul(fp2 *r, const fp2 *a, const fp2 *b, int d);
void fq_inv(fp2 *r, const fp2 *a, int d);
bool fq_sqrt(fp2 *r, const fp2 *a, int d);
bool fq_is_zero(const fp2 *a, int d);
bool fq_eq(const fp2 *a, const fp2 *b, int d);
void fp2_conj(fp2 *r, const fp2 *a);

void pt_dbl(point *r, const point *p, int d);
void pt_add(point *r, const point *p, const point *q, int d);
void pt_neg(point *r, const point *p, int d);
void pt_mul(point *r, const point *p, const uint8_t *k, size_t klen, int d);
void pt_mul_x(point *r, const point *p, int d);
bool pt_eq(const point *p, const point *q, int d);
bool pt_is_inf(const point *p, int d);
void pt_set_affine(point *p, int d);
bool pt_on_curve(const point *p, int d);
bool pt_load(point *p, const uint8_t *in, int d, bool swap);
void pt_store(uint8_t *out, const point *p, int d, bool swap);
void pt_frob(point *r, const point *p, bool inverse);

void pairing_add(fp12 *f, const point *p, const point *q);
bool pairing_is_one(fp12 *f);
void fp12_one(fp12 *r);

#endif
