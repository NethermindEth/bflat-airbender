/*
 * Montgomery arithmetic modulo an odd 256-bit modulus on the u256 delegation, R = 2^256.
 * Values are little-endian u256d; everything but mont_init is a handful of delegations.
 */
#ifndef AIRBENDER_MONT_H
#define AIRBENDER_MONT_H

#include "u256.h"

typedef struct
{
    u256d n;      /* the modulus, odd */
    u256d nprime; /* -n^-1 mod R */
    u256d r2;     /* R^2 mod n, to enter the domain */
    u256d one;    /* R mod n, the domain's 1 */
} mont_ctx;

/* Fills ctx for the modulus n (little-endian words). */
void mont_init(mont_ctx *ctx, const uint32_t n[8]);

/* Inline: a call and its return cost as much as the operation's own glue. Copies go through the
   delegation too (MEMCOPY), one instruction instead of sixteen. */
extern u256d MontLo, MontHi, MontM, MontT, MontZero;
#define MONT_INLINE static inline __attribute__((always_inline))

/* *v -= n when that does not borrow, or unconditionally when `wrapped` says *v overflowed. */
MONT_INLINE void mont_reduce_once(u256d *v, uint32_t wrapped, const mont_ctx *ctx)
{
    if (wrapped)
    {
        u256_op(v, &ctx->n, U256_SUB);
        return;
    }
    u256_op(&MontT, v, U256_MEMCOPY);
    if (!u256_op(&MontT, &ctx->n, U256_SUB))
        u256_op(v, &MontT, U256_MEMCOPY);
}

/* r = a * b * R^-1 mod n, for a, b < n. r may alias a or b. */
MONT_INLINE void mont_mul(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx)
{
    u256_op(&MontLo, a, U256_MEMCOPY);
    u256_op(&MontLo, b, U256_MUL_LOW);
    u256_op(&MontHi, a, U256_MEMCOPY);
    u256_op(&MontHi, b, U256_MUL_HIGH);

    /* m = lo * n' mod R makes lo + m*n a multiple of R: its low half is R exactly when lo != 0. */
    u256_op(&MontM, &MontLo, U256_MEMCOPY);
    u256_op(&MontM, &ctx->nprime, U256_MUL_LOW);
    u256_op(&MontM, &ctx->n, U256_MUL_HIGH);

    uint32_t carry = u256_op(&MontLo, &MontZero, U256_EQ) ? 0 : U256_CARRY;
    uint32_t wrapped = u256_op(&MontHi, &MontM, U256_ADD | carry);
    mont_reduce_once(&MontHi, wrapped, ctx);
    u256_op(r, &MontHi, U256_MEMCOPY);
}

/* r = a + b mod n and r = a - b mod n, for a, b < n. r may alias a or b. */
MONT_INLINE void mont_add(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx)
{
    u256_op(&MontHi, a, U256_MEMCOPY);
    uint32_t wrapped = u256_op(&MontHi, b, U256_ADD);
    mont_reduce_once(&MontHi, wrapped, ctx);
    u256_op(r, &MontHi, U256_MEMCOPY);
}

MONT_INLINE void mont_sub(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx)
{
    u256_op(&MontHi, a, U256_MEMCOPY);
    if (u256_op(&MontHi, b, U256_SUB))
        u256_op(&MontHi, &ctx->n, U256_ADD);
    u256_op(r, &MontHi, U256_MEMCOPY);
}

/* Into the domain (a < n) and out of it. */
void mont_to(u256d *r, const u256d *a, const mont_ctx *ctx);
void mont_from(u256d *r, const u256d *a, const mont_ctx *ctx);

/* r = a^e in the domain; e is a plain 256-bit exponent. r may alias a. */
void mont_pow(u256d *r, const u256d *a, const u256d *e, const mont_ctx *ctx);

#endif
