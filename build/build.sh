#!/bin/bash
# Builds libairbender.a for rv32im/ilp32 into output/<platform>/, where <platform> is
#   airbender - the bare-metal zkVM guest (CSR input, register output);
#   linux     - a Linux process under qemu-riscv32, for testing the guest off the zkVM.
# The linux platform also builds output/linux/harness, the test harness of tests/run.py.
# Needs clang with the riscv32 target, ld.lld and the LLVM tools; SYSROOT is the rv32 musl sysroot
# (with the compiler-rt builtins) the linux platform compiles and links against. The image of
# the Dockerfile has them all: build/ci.sh.
set -euo pipefail

TOP="$(cd "$(dirname "$0")/.." && pwd -P)"
platform="${1:-airbender}"
SYSROOT="${SYSROOT:-/rv32}"

common=(--target=riscv32-unknown-elf -march=rv32im -mabi=ilp32 -mno-relax -O2 -ffunction-sections -fdata-sections
        -fno-stack-protector -Wall -Wextra -Werror)
case "$platform" in
    airbender) flags=("${common[@]}" -DAIRBENDER -ffreestanding -nostdinc -isystem "$(clang -print-resource-dir)/include") ;;
    linux) flags=("${common[@]/--target=riscv32-unknown-elf/--target=riscv32-unknown-linux-musl}" --sysroot="$SYSROOT") ;;
    *) echo "unknown platform $platform" >&2; exit 1 ;;
esac

sources=(keccak sha256 ripemd160 blake2f ecc modexp ec bn254 bls12 mont "platform_$platform")
# The u256 delegation exists only on Airbender; elsewhere u256_soft computes its results.
[ "$platform" = airbender ] || sources+=(u256_soft)
out="$TOP/output/$platform"
obj="$out/obj"
mkdir -p "$obj"
objects=()
for s in "${sources[@]}"; do
    clang "${flags[@]}" -c "$TOP/src/$s.c" -o "$obj/$s.o"
    objects+=("$obj/$s.o")
done
# What bflat --extlib reads: the platform (its entry object, linker script and manifest) and the manifest of the
# bindings, beside the library. The entry object has the finish and error sequences Airbender recognizes, and the
# console.
if [ "$platform" = airbender ]; then
    bash "$TOP/platform/build.sh" "$out"
    cp "$TOP/bflat-manifest.json" "$out/libairbender.bflat.manifest"
fi
# The archive holds one object, in which only the entry points stay global: the guest links the library with musl,
# uGC and its own native code, and the names of the helpers (fp_mul, pt_add, C, ...) would clash with theirs. The
# entry points are the zkvm_ accelerators, the I/O of Nethermind.Zkvm.Abstractions, the airbender_ calls of bflat
# and of the guest, and the __wrap_ targets of the manifest.
exports='^(zkvm_|airbender_|syscall_|__wrap_)|^(read_input|write_output|sys_write)$'
ld.lld -r -o "$obj/libairbender.o" "${objects[@]}"
llvm-nm --defined-only -g --format=just-symbols "$obj/libairbender.o" | grep -E "$exports" > "$obj/exports.txt" || true
llvm-objcopy --keep-global-symbols="$obj/exports.txt" "$obj/libairbender.o"
rm -f "$out/libairbender.a"
llvm-ar rcs "$out/libairbender.a" "$obj/libairbender.o"

# What bflat and the guest call must be there; on Airbender, what the library calls must be what the link of
# bflat provides: the memory functions of musl, the compiler-rt builtins and the platform interface of the entry.
required=(read_input write_output sys_write zkvm_keccak256 zkvm_sha256 zkvm_secp256k1_ecrecover)
[ "$platform" = airbender ] && required+=(airbender_exit airbender_fail airbender_u256_op airbender_keccak256
    airbender_secp256k1_ecrecover_batch airbender_keccak_f1600_in_place)
for s in "${required[@]}"; do
    grep -qx "$s" "$obj/exports.txt" || { echo "libairbender.a does not define $s" >&2; exit 1; }
done
if [ "$platform" = airbender ]; then
    unexpected=$(llvm-nm -u --format=just-symbols "$obj/libairbender.o" \
        | grep -vxE '(malloc|calloc|realloc|free|memcpy|memmove|memset|memcmp|__real_memcpy)|__[a-z]+[dst]i[0-9]|(zkvm_platform_exit|zkvm_platform_write|airbender_finish_error|airbender_output)' || true)
    [ -z "$unexpected" ] || { echo "libairbender.a needs what the guest link does not provide:" $unexpected >&2; exit 1; }
fi
echo "built $out/libairbender.a: $(wc -l < "$obj/exports.txt") entry points"

# The test harness of tests/run.py, a static Linux executable for qemu-riscv32.
if [ "$platform" = linux ]; then
    # stdio_noldbl replaces musl's formatted I/O, which needs quad-precision helpers riscv32 has none of.
    clang "${flags[@]}" -c "$TOP/tests/harness.c" -o "$obj/harness.o"
    clang "${flags[@]}" -c "$TOP/tests/stdio_noldbl.c" -o "$obj/stdio_noldbl.o"
    lib="$SYSROOT/usr/lib"
    ld.lld -static -o "$out/harness" "$lib/crt1.o" "$lib/crti.o" "$obj/harness.o" "$obj/stdio_noldbl.o" "$out/libairbender.a" \
        "$lib/libc.a" "$lib/libclang_rt.builtins-riscv32.a" "$lib/crtn.o"
    echo "built $out/harness"
fi
