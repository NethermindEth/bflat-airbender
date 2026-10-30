/*
 * The entry point of bflat guests for ZKsync Airbender, and the platform
 * interface that bflat's pal calls (docs/manifest.md in bflat-riscv64).
 *
 * Airbender (matter-labs/zksync-airbender) proves RV32IM: ilp32, no A, no C,
 * and every load and store naturally aligned. A program is a flat image whose
 * instructions the circuit reads from a ROM at address 0; RAM starts zeroed,
 * so nothing initialised can live there until _start copies it in.
 *
 * The platform owns how a run ends and the console:
 *   - zkvm_platform_exit(0) publishes airbender_output: it loads the 16 words
 *     into x10-x25 and spins (airbender_finish_success). libairbender's
 *     write_output fills the words; without it they stay zero.
 *   - zkvm_platform_exit(code != 0) writes the cycle CSR, which no circuit
 *     satisfies, so the run cannot be proven (airbender_finish_error).
 *   - zkvm_platform_write prints through the quasi-UART of the
 *     non-determinism CSR (0x7c0): 0xffffffff, the number of words that
 *     follow, the byte length, then the bytes packed little-endian into
 *     words. Debug output only; the prover ignores it.
 *
 * Copyright (C) 2026 Demerzel Solutions Limited (Nethermind)
 */
    .option norvc
    .option norelax

.section .rodata
    .balign 4
app_str:
    .asciz "app"

    .section .data
    .balign 4
argv_vec:
    .word app_str
    .word 0

/* Where a jump to address 0 goes once the program runs; see _start. */
    .balign 4
zero_jump_target:
    .word zero_jump_unset

/* The published output: x10-x17 the output, x18-x25 zero, as the recursion
 * convention expects. */
    .section .bss.airbender_output, "aw", @nobits
    .balign 4
    .global airbender_output
    .type   airbender_output, @object
airbender_output:
    .zero   64
    .size   airbender_output, 64

.section .text.init, "ax", @progbits
.global _start
    .cfi_sections .eh_frame
_start:
    .cfi_startproc
    .cfi_def_cfa sp, 0
    .cfi_undefined ra
    /* Every register starts zeroed, so sp is zero only at reset. Running code
     * reaches address 0 only through a jump to an empty code pointer, which
     * goes on to airbender_set_zero_jump_target's target with every register
     * as the jump left it: a dispatch table that starts zeroed, as bump
     * allocated memory does, needs no fill and no test for an empty entry. */
    bnez    sp, zero_jump
    la      gp, _global_pointer
    la      sp, _init_stack_top

    /* Initialised RAM: one block linked to run in RAM and loaded in ROM, which
     * the delegations cannot read. Eight words a round, then the rest one by
     * one: without .NET's data rehydration the block is hundreds of KiB. */
    la      t0, __ram_init_load
    la      t1, __ram_init_start
    la      t2, __ram_init_end
    addi    t4, t2, -32
    bgtu    t1, t4, 3f
1:  lw      t3, 0(t0)
    lw      t5, 4(t0)
    lw      t6, 8(t0)
    lw      a2, 12(t0)
    lw      a3, 16(t0)
    lw      a4, 20(t0)
    lw      a5, 24(t0)
    lw      a6, 28(t0)
    sw      t3, 0(t1)
    sw      t5, 4(t1)
    sw      t6, 8(t1)
    sw      a2, 12(t1)
    sw      a3, 16(t1)
    sw      a4, 20(t1)
    sw      a5, 24(t1)
    sw      a6, 28(t1)
    addi    t0, t0, 32
    addi    t1, t1, 32
    bleu    t1, t4, 1b
3:  bgeu    t1, t2, 2f
    lw      t3, 0(t0)
    sw      t3, 0(t1)
    addi    t0, t0, 4
    addi    t1, t1, 4
    j       3b
2:
    /* .bss needs no clearing: RAM starts zeroed, in the circuit (its memory
     * argument fixes the initial value at 0) and in the simulator alike, and
     * nothing runs before this entry point. It is megabytes. */
    /* noos_start_main (bflat's noos module) replaces musl's
     * __libc_start_main: thread pointer, .init_array, then the .NET start-up.
     * It ends through musl's exit(), which pal wraps into zkvm_platform_exit. */
    li      a0, 1            # argc = 1
    la      a1, argv_vec     # argv = &argv_vec[0]
    call    noos_start_main
5:  j       5b               # unreachable: noos_start_main never returns
    .cfi_endproc
    .size   _start, .-_start

zero_jump:
    lui     t0, %hi(zero_jump_target)
    lw      t0, %lo(zero_jump_target)(t0)
    jr      t0

/* A jump to address 0 before any target is set is a fault. */
zero_jump_unset:
    unimp

/* void airbender_set_zero_jump_target(void *target) */
.section .text.airbender_set_zero_jump_target, "ax", @progbits
.global airbender_set_zero_jump_target
.type   airbender_set_zero_jump_target, @function
airbender_set_zero_jump_target:
    lui     t0, %hi(zero_jump_target)
    sw      a0, %lo(zero_jump_target)(t0)
    ret
.size   airbender_set_zero_jump_target, .-airbender_set_zero_jump_target

/* void zkvm_platform_exit(long code): never returns. */
.section .text.zkvm_platform_exit, "ax", @progbits
.global zkvm_platform_exit
.type   zkvm_platform_exit, @function
zkvm_platform_exit:
    bnez    a0, airbender_finish_error
    la      a0, airbender_output
    j       airbender_finish_success
.size   zkvm_platform_exit, .-zkvm_platform_exit

/* void zkvm_platform_write(int fd, const char *buf, int len): fd is not
 * used; Airbender has one console. */
.section .text.zkvm_platform_write, "ax", @progbits
.global zkvm_platform_write
.type   zkvm_platform_write, @function
zkvm_platform_write:
    blez    a2, 4f
    li      t0, -1
    csrrw   x0, 0x7c0, t0
    addi    t0, a2, 3
    srli    t0, t0, 2
    addi    t0, t0, 1         # the words that follow: the length and the data
    csrrw   x0, 0x7c0, t0
    csrrw   x0, 0x7c0, a2
1:  li      t0, 0             # the word
    li      t1, 0             # the shift of the next byte
    li      t2, 32
2:  lbu     t3, 0(a1)
    sll     t3, t3, t1
    or      t0, t0, t3
    addi    a1, a1, 1
    addi    a2, a2, -1
    addi    t1, t1, 8
    beqz    a2, 3f
    bne     t1, t2, 2b
3:  csrrw   x0, 0x7c0, t0
    bnez    a2, 1b
4:  ret
.size   zkvm_platform_write, .-zkvm_platform_write

/*
 * The two ways an Airbender run ends.
 *
 * airbender_finish_success(a0 = 16 words): loads them into x10-x25 and
 * spins. The loads from s10 and the final `j .` are the EXIT_SEQUENCE
 * (execution_utils) that the Airbender tooling finds in the binary to know
 * the final pc, so the encoding must stay exactly this and appear only here.
 *
 * airbender_finish_error(a0 = code): writes the cycle CSR, which no circuit
 * satisfies, so the run cannot be proven. The code stays in a0 for whoever
 * reads the simulator's registers.
 */
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
