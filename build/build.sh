#!/bin/bash
# Builds libairbender.a for rv32im/ilp32 into output/<platform>/, where <platform> is
#   airbender - the bare-metal zkVM guest (CSR input, register output);
#   linux     - a Linux process under qemu-riscv32, for testing the guest off the zkVM.
# Needs clang with the riscv32 target and llvm-ar; SYSROOT is the rv32 musl sysroot the
# linux platform compiles against.
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
# The finish and error sequences Airbender recognizes, which only assembly can lay down.
if [ "$platform" = airbender ]; then
    clang --target=riscv32-unknown-elf -march=rv32im -mabi=ilp32 -c "$TOP/src/platform_airbender_asm.S" -o "$obj/platform_airbender_asm.o"
    objects+=("$obj/platform_airbender_asm.o")
    # What bflat --extlib reads: the static library and the target it is for.
    cat > "$out/libairbender.bflat.manifest" <<'EOF'
{
  "name": "libairbender",
  "builds": [
    { "arch": "riscv32", "os": "linux", "libc": "airbender", "static_lib": "libairbender.a" }
  ]
}
EOF
fi
rm -f "$out/libairbender.a"
llvm-ar rcs "$out/libairbender.a" "${objects[@]}"
echo "built $out/libairbender.a"
