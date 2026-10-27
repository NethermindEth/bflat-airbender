/*
 * Shared declarations of the Airbender zkVM runtime: the eth-act accelerator status codes and the
 * helpers the accelerator sources share. The public entry points follow
 * eth-act/zkvm-standards c-interface-accelerators/zkvm_accelerators.h.
 */
#ifndef AIRBENDER_ZKVM_H
#define AIRBENDER_ZKVM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* The bare-metal build has no libc headers; the runtime provides both. */
void *memcpy(void *dst, const void *src, size_t n);
void *memset(void *dst, int c, size_t n);

typedef enum
{
    ZKVM_EOK = 0,
    ZKVM_EFAIL = -1
} zkvm_status;

/* Ends the program as failed. The accelerators the guest never reaches on the blocks it proves
   call it rather than return a wrong answer. */
__attribute__((noreturn)) void airbender_fail(uint32_t code);

void keccak_f1600(uint64_t state[25]);

#endif
