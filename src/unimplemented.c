/*
 * Accelerators not implemented yet. Each ends the run as failed with its own code rather than
 * returning a status: a status of failure would be taken as an invalid input, and the guest
 * would prove a wrong result instead of stopping.
 */
#include "zkvm.h"

#define UNIMPLEMENTED(name, code)                                                                                      \
    zkvm_status name(void)                                                                                             \
    {                                                                                                                  \
        airbender_fail(code);                                                                                          \
    }

UNIMPLEMENTED(zkvm_kzg_point_eval, 0x104)
UNIMPLEMENTED(zkvm_bls12_g1_add, 0x105)
UNIMPLEMENTED(zkvm_bls12_g1_msm, 0x106)
UNIMPLEMENTED(zkvm_bls12_g2_add, 0x107)
UNIMPLEMENTED(zkvm_bls12_g2_msm, 0x108)
UNIMPLEMENTED(zkvm_bls12_pairing, 0x109)
UNIMPLEMENTED(zkvm_bls12_map_fp_to_g1, 0x10a)
UNIMPLEMENTED(zkvm_bls12_map_fp2_to_g2, 0x10b)
