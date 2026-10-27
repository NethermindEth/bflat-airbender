/*
 * The two ways an Airbender run ends.
 *
 * airbender_finish_success(a0 = 16 words): loads them into x10-x25 and spins. The loads from s10
 * and the final `j .` are the EXIT_SEQUENCE (execution_utils) that the Airbender tooling finds in
 * the binary to know the final pc, so the encoding must stay exactly this and appear only here.
 *
 * airbender_finish_error(a0 = code): writes the cycle CSR, which no circuit satisfies, so the run
 * cannot be proven. The code stays in a0 for whoever reads the simulator's registers.
 */
    .option norvc
    .section .text.airbender_finish, "ax", @progbits
    .balign 4

    .global airbender_finish_success
    .type airbender_finish_success, @function
airbender_finish_success:
    mv      s10, a0
    lw      a0, 0x0(s10)
    lw      a1, 0x4(s10)
    lw      a2, 0x8(s10)
    lw      a3, 0xc(s10)
    lw      a4, 0x10(s10)
    lw      a5, 0x14(s10)
    lw      a6, 0x18(s10)
    lw      a7, 0x1c(s10)
    lw      s2, 0x20(s10)
    lw      s3, 0x24(s10)
    lw      s4, 0x28(s10)
    lw      s5, 0x2c(s10)
    lw      s6, 0x30(s10)
    lw      s7, 0x34(s10)
    lw      s8, 0x38(s10)
    lw      s9, 0x3c(s10)
1:  j       1b
    .size airbender_finish_success, .-airbender_finish_success

    .global airbender_finish_error
    .type airbender_finish_error, @function
airbender_finish_error:
    csrrw   x0, cycle, x0
2:  j       2b
    .size airbender_finish_error, .-airbender_finish_error
