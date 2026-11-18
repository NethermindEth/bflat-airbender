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

#ifdef AIRBENDER
/* Inline, as the delegation itself is one instruction: a call would cost several times it. */
static inline uint32_t u256_op(u256d *a, const u256d *b, uint32_t op)
{
    register uint32_t x10 __asm__("a0") = (uint32_t)a;
    register uint32_t x11 __asm__("a1") = (uint32_t)b;
    register uint32_t x12 __asm__("a2") = op;
    __asm__ volatile("csrrw x0, 0x7ca, x0" : "+r"(x12) : "r"(x10), "r"(x11) : "memory");
    return x12;
}
#else
uint32_t u256_op(u256d *a, const u256d *b, uint32_t op);
#endif

/* The bytes of x reversed: Zbb's rev8, which Airbender's machines decode on their own, without the rest of Zbb.
 * One instruction where the shifts and masks are eight or more. */
static inline __attribute__((always_inline)) uint32_t rev8(uint32_t x)
{
#ifdef AIRBENDER
    uint32_t r;
    __asm__(".insn i 0x13, 5, %0, %1, 0x698" : "=r"(r) : "r"(x));
    return r;
#else
    return __builtin_bswap32(x);
#endif
}

/* The field hint tape: values the prover gives and the code checks, one after another, such as the slope of an
 * affine addition, checked by one product where computing it takes an inversion. `active` says the input has a tape
 * section at all; an exhausted or wrong tape only costs the computation. */
typedef struct
{
    const u256d *next, *end;
    uint32_t active;
} field_hints;
extern field_hints FieldHints;
/* Called with each value that was computed for want of a hint, in tape order: the point a tape is recorded at. */
void field_hint_emit(const u256d *v);

#endif
