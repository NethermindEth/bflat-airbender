# bflat-airbender

The zkVM runtime library for guests that bflat builds for ZKsync Airbender
(`bflat build --arch riscv32 --libc airbender`): `libairbender.a`, for RV32IM
with the ilp32 ABI.

It supplies what the guest reaches through `__Internal` P/Invokes
(Nethermind.Zkvm.Abstractions):

- Input and output over Airbender's non-determinism CSR (0x7c0): the input
  words, the UART, and the finish and error sequences Airbender recognizes.
- The precompiles in C: Keccak-256, SHA-256, RIPEMD-160, BLAKE2 F,
  secp256k1 and P-256 signatures, modexp, BN254, BLS12-381 and KZG point
  evaluation. The field arithmetic runs on Airbender's u256 delegation
  (CSR 0x7ca).

## Build

Needs clang with the riscv32 target and llvm-ar.

```bash
build/build.sh airbender   # output/airbender/libairbender.a and libairbender.bflat.manifest
build/build.sh linux       # output/linux/libairbender.a, for tests under qemu-riscv32
```

`SYSROOT` names the rv32 musl sysroot that the linux platform compiles
against. The airbender platform is freestanding. bflat takes the library with
`--extlib output/airbender/libairbender.bflat.manifest`.

## Test

`tests/run.py <harness-command...>` checks the precompiles against
pure-Python references, through `tests/harness.c` built against the linux
library and run under qemu-riscv32. That build has no u256 delegation, so
`u256_soft.c` computes its results.
