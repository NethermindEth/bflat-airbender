/*
 * Airbender's u256 delegation (CSR 0x7ca): one 256-bit operation per call, on two operands that
 * must sit in RAM (not ROM), 32-byte aligned and distinct. The result replaces `a`; the call
 * returns the operation's flag: the carry or borrow out of ADD/SUB, whether MUL_LOW overflowed,
 * whether EQ found the operands equal.
 *
 * The airbender platform issues the delegation; the linux platform computes the same results in C
 * (u256_soft.c), so the code built on it runs and is tested unchanged under qemu.
 */
#ifndef AIRBENDER_U256_H
#define AIRBENDER_U256_H

#include <stdint.h>

/* Little-endian 32-bit words, aligned as the delegation requires. */
typedef struct
{
    uint32_t w[8];
} __attribute__((aligned(32))) u256d;

enum
{
    U256_ADD = 1u << 0,
    U256_SUB = 1u << 1,
    U256_SUB_AND_NEGATE = 1u << 2, /* a = b - a */
    U256_MUL_LOW = 1u << 3,
    U256_MUL_HIGH = 1u << 4,
    U256_EQ = 1u << 5,
    U256_CARRY = 1u << 6, /* with ADD/SUB: carry or borrow in */
    U256_MEMCOPY = 1u << 7 /* a = b */
};

uint32_t u256_op(u256d *a, const u256d *b, uint32_t op);

#endif
