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

/* r = a * b * R^-1 mod n, for a, b < n. r may alias a or b. */
void mont_mul(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx);

/* r = a + b mod n and r = a - b mod n, for a, b < n. r may alias a or b. */
void mont_add(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx);
void mont_sub(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx);

/* Into the domain (a < n) and out of it. */
void mont_to(u256d *r, const u256d *a, const mont_ctx *ctx);
void mont_from(u256d *r, const u256d *a, const mont_ctx *ctx);

/* r = a^e in the domain; e is a plain 256-bit exponent. r may alias a. */
void mont_pow(u256d *r, const u256d *a, const u256d *e, const mont_ctx *ctx);

#endif
