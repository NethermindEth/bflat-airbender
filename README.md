# bflat-airbender

The ZKsync Airbender platform of bflat, and the zkVM runtime library of its
guests, for RV32IM with the ilp32 ABI
(`bflat build --arch riscv32 --libc airbender`).

The platform (`platform/`, see the manifest format in bflat-riscv64's
`docs/manifest.md`) is what bflat links every Airbender guest with:

- `airbender_entry.o` (`entry.S`): `_start`, which copies the initialized
  data from ROM to RAM and enters the runtime, and the platform interface of
  bflat's pal: `zkvm_platform_exit` (success publishes the 16 output words in
  x10-x25 through the EXIT_SEQUENCE the Airbender tooling looks for; failure
  writes the cycle CSR, which no circuit satisfies) and `zkvm_platform_write`
  (the quasi-UART of the non-determinism CSR).
- `airbender.ld`: the memory map. The ROM at address 0 holds the image, 4 MiB
  in v3; the link fails if it does not fit.
- `airbender.bflat.manifest`: the platform manifest. The modules of bflat it
  requires, the ISA (no C, A, F or D), the RyuJIT options (unsigned RV32M only,
  no unaligned access, the zkVM tuning measured on the Nethermind guest), no
  unwind tables, and the post-link steps: the FENCE and signed RV32M rewrites
  of the native code, and the flat images `<output>.bin` and `<output>.text`
  that Airbender loads.

`libairbender.a`, with its manifest `libairbender.bflat.manifest`, supplies
what the guest reaches through `__Internal` P/Invokes
(Nethermind.Zkvm.Abstractions, and the `airbender_` calls of the Nethermind
guest):

- Input and output over Airbender's non-determinism CSR (0x7c0): the input
  words and the output words (`read_input`, `write_output`), and the console
  through the platform.
- The precompiles in C: Keccak-256, SHA-256, RIPEMD-160, BLAKE2 F,
  secp256k1 and P-256 signatures, modexp, BN254, BLS12-381 and KZG point
  evaluation. The field arithmetic runs on Airbender's u256 delegation
  (CSR 0x7ca), Keccak on its Keccak-f delegation.

A guest that uses no accelerator links with the platform alone.

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
build/build.sh airbender   # output/airbender: the platform, libairbender.a and libairbender.bflat.manifest
build/build.sh linux       # output/linux/libairbender.a and harness, for tests under qemu-riscv32
build/check.sh airbender linux
```

`SYSROOT` names the rv32 musl sysroot, with the compiler-rt builtins, that the
linux platform compiles and links against. The airbender platform is
freestanding: bflat links it with musl's memory functions and the builtins,
and the build fails if it needs anything else. bflat takes the library with
`--extlib output/airbender/libairbender.bflat.manifest`: it reads every
`.bflat.manifest` of the directory, the platform among them.
`platform/build.sh <dir>` builds the platform alone, with clang only.

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

`tests/platform/run-tls.sh` checks that the memory map keeps the initialized
TLS in the ROM copy of the RAM block.

## Release

The Release workflow (run by hand, with the major.minor version) builds and
tests as CI does and releases `v<version>.<run>`: the platform
(`airbender_entry.o`, `airbender.ld`, `airbender.bflat.manifest`),
`libairbender.a` and `libairbender.bflat.manifest`, and
`libairbender-riscv32.tar.gz` with all of them.
Two more assets are `bflat-signature.sha256`, the SHA-256 of each other
asset, and `bflat-signature.sha256.sig`, its signature. The sign job, the
reusable workflow `sign.yml` of bflat-riscv64, signs the list with the secret
`BFLAT_SIGNING_KEY`; it is the only job that has the key, and it runs no
build. The release job checks the assets against the list and the signature
(`tools/sign_package.py check-files`, with the variable
`BFLAT_SIGNING_PUBLIC_KEY`) before it makes the release.
`tools/sign_package.py` is a copy of the script in bflat-riscv64; keep it the
same.
