/*
 * Montgomery arithmetic on the u256 delegation (see mont.h).
 *
 * The delegation reads both operands from RAM, so nothing here takes an operand from a constant:
 * scratch values are static (RAM) and contexts are filled at run time.
 */
#include "mont.h"

/* Copies go through the delegation too (MEMCOPY), one instruction instead of sixteen. */
static u256d Lo, Hi, M, T, Zero;

/* *v -= n when that does not borrow, or unconditionally when `wrapped` says *v overflowed. */
static void reduce_once(u256d *v, uint32_t wrapped, const mont_ctx *ctx)
{
    if (wrapped)
    {
        u256_op(v, &ctx->n, U256_SUB);
        return;
    }
    u256_op(&T, v, U256_MEMCOPY);
    if (!u256_op(&T, &ctx->n, U256_SUB))
        u256_op(v, &T, U256_MEMCOPY);
}

void mont_mul(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx)
{
    u256_op(&Lo, a, U256_MEMCOPY);
    u256_op(&Lo, b, U256_MUL_LOW);
    u256_op(&Hi, a, U256_MEMCOPY);
    u256_op(&Hi, b, U256_MUL_HIGH);

    /* m = lo * n' mod R makes lo + m*n a multiple of R: its low half is R exactly when lo != 0. */
    u256_op(&M, &Lo, U256_MEMCOPY);
    u256_op(&M, &ctx->nprime, U256_MUL_LOW);
    u256_op(&M, &ctx->n, U256_MUL_HIGH);

    uint32_t carry = u256_op(&Lo, &Zero, U256_EQ) ? 0 : U256_CARRY;
    uint32_t wrapped = u256_op(&Hi, &M, U256_ADD | carry);
    reduce_once(&Hi, wrapped, ctx);
    u256_op(r, &Hi, U256_MEMCOPY);
}

void mont_add(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx)
{
    u256_op(&Hi, a, U256_MEMCOPY);
    uint32_t wrapped = u256_op(&Hi, b, U256_ADD);
    reduce_once(&Hi, wrapped, ctx);
    u256_op(r, &Hi, U256_MEMCOPY);
}

void mont_sub(u256d *r, const u256d *a, const u256d *b, const mont_ctx *ctx)
{
    u256_op(&Hi, a, U256_MEMCOPY);
    if (u256_op(&Hi, b, U256_SUB))
        u256_op(&Hi, &ctx->n, U256_ADD);
    u256_op(r, &Hi, U256_MEMCOPY);
}

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
        reduce_once(&u, wrapped, ctx);
    }
    ctx->r2 = u;
}
