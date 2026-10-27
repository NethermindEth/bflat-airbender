/*
 * Montgomery arithmetic on the u256 delegation (see mont.h).
 *
 * The delegation reads both operands from RAM, so nothing here takes an operand from a constant:
 * scratch values are static (RAM) and contexts are filled at run time.
 */
#include "mont.h"

/* The scratch of the inline operations in mont.h, in RAM where the delegation reads its operands. */
u256d MontLo, MontHi, MontM, MontT, MontZero;

void mont_to(u256d *r, const u256d *a, const mont_ctx *ctx)
{
    mont_mul(r, a, &ctx->r2, ctx);
}

void mont_from(u256d *r, const u256d *a, const mont_ctx *ctx)
{
    static u256d one_plain;
    one_plain = (u256d){{1}};
    mont_mul(r, a, &one_plain, ctx);
}

void mont_pow(u256d *r, const u256d *a, const u256d *e, const mont_ctx *ctx)
{
    static u256d base, acc;
    base = *a;
    acc = ctx->one;
    int started = 0;
    for (int i = 255; i >= 0; i--)
    {
        if (started)
            mont_mul(&acc, &acc, &acc, ctx);
        if ((e->w[i / 32] >> (i % 32)) & 1)
        {
            mont_mul(&acc, &acc, &base, ctx);
            started = 1;
        }
    }
    *r = acc;
}

void mont_init(mont_ctx *ctx, const uint32_t n[8])
{
    static u256d x, t, u, two;
    for (int i = 0; i < 8; i++)
        ctx->n.w[i] = n[i];

    /* n^-1 mod R by Newton: n is its own inverse mod 8, and each step doubles the bits. */
    x = ctx->n;
    for (int step = 0; step < 7; step++)
    {
        t = ctx->n;
        u256_op(&t, &x, U256_MUL_LOW);
        two = (u256d){{2}};
        u256_op(&two, &t, U256_SUB);
        u256_op(&x, &two, U256_MUL_LOW);
    }
    ctx->nprime = (u256d){{0}};
    u256_op(&ctx->nprime, &x, U256_SUB);

    /* R mod n = (R - n) mod n. */
    u = (u256d){{0}};
    u256_op(&u, &ctx->n, U256_SUB);
    for (;;)
    {
        t = u;
        if (u256_op(&t, &ctx->n, U256_SUB))
            break;
        u = t;
    }
    ctx->one = u;

    /* R^2 mod n: double R mod n 256 times. */
    for (int i = 0; i < 256; i++)
    {
        t = u;
        uint32_t wrapped = u256_op(&u, &t, U256_ADD);
        mont_reduce_once(&u, wrapped, ctx);
    }
    ctx->r2 = u;
}
