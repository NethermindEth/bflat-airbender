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

```bash
build/ci.sh
```

builds both platforms in the image of the `Dockerfile`, runs the gates, and
then the tests under `qemu-riscv32` (`QEMU` names another). The image is the
x64 cross-build container of the .NET runtime that dotnet-riscv and uGC build
in, so clang, lld and the LLVM tools are those of the riscv32 runtime; the rv32im
musl sysroot and the gates come from dotnet-riscv at `DOTNET_RISCV_REF`.
In the image, or with clang, ld.lld and the LLVM tools at hand:

```bash
build/build.sh airbender   # output/airbender/libairbender.a and libairbender.bflat.manifest
build/build.sh linux       # output/linux/libairbender.a and harness, for tests under qemu-riscv32
build/check.sh airbender linux
```

`SYSROOT` names the rv32 musl sysroot, with the compiler-rt builtins, that the
linux platform compiles and links against. The airbender platform is
freestanding: bflat links it with musl's memory functions and the builtins,
and the build fails if it needs anything else. bflat takes the library with
`--extlib output/airbender/libairbender.bflat.manifest`.

The archive holds one object, in which only the entry points are global: the
`zkvm_` accelerators, the I/O of Nethermind.Zkvm.Abstractions (`read_input`,
`write_output`, `sys_write`), the `airbender_` calls of bflat and of the guest,
and the `__wrap_` targets of the manifest. The helpers stay local, so that
their names cannot clash with those of the guest's other native code.

`build/check.sh` fails if an object is not RV32 ilp32 without the C extension
(`check_ilp32.sh` of dotnet-riscv) or carries a compressed, atomic or fence
instruction (`tools/isa_scan.py`), as Airbender decodes RV32IM only.

## Test

`tests/run.py <harness-command...>` checks the precompiles against
pure-Python references, through `tests/harness.c` built against the linux
library and run under qemu-riscv32 (`python3 tests/run.py qemu-riscv32
output/linux/harness`). That build has no u256 delegation, so `u256_soft.c`
computes its results. It fails if the harness answers fewer cases than it was
given.

## Release

The Release workflow (run by hand, with the major.minor version) builds and
tests as CI does and releases `v<version>.<run>`: `libairbender.a`,
`libairbender.bflat.manifest`, and `libairbender-riscv32.tar.gz` with both.
Two more assets are `bflat-signature.sha256`, the SHA-256 of each other
asset, and `bflat-signature.sha256.sig`, its signature. The sign job, the
reusable workflow `sign.yml` of bflat-riscv64, signs the list with the secret
`BFLAT_SIGNING_KEY`; it is the only job that has the key, and it runs no
build. The release job checks the assets against the list and the signature
(`tools/sign_package.py check-files`, with the variable
`BFLAT_SIGNING_PUBLIC_KEY`) before it makes the release.
`tools/sign_package.py` is a copy of the script in bflat-riscv64; keep it the
same.
