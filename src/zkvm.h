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

/* The Keccak-f[1600] permutation: the platform's, a delegation on Airbender. */
void keccak_f1600(uint64_t state[25]);
#ifdef AIRBENDER
void airbender_keccak256(const uint8_t *data, size_t len, uint8_t out[32]);
#endif
void keccak_f1600_soft(uint64_t state[25]);

/* Recovers the uncompressed public keys (x || y, big-endian) of count signatures: msgs holds 32-byte hashes,
   sigs 64-byte r || s, recids one byte each. ok[i] is 1 when outputs[i] holds signature i's key, as
   zkvm_secp256k1_ecrecover would give it, and 0 when that call would fail. The inverses of every r and of every
   recovered z are taken together, an exponentiation each for the whole batch rather than per signature. */
void airbender_secp256k1_ecrecover_batch(const uint8_t *msgs, const uint8_t *sigs, const uint8_t *recids,
                                         uint8_t *outputs, uint8_t *ok, uint32_t count);

#endif
