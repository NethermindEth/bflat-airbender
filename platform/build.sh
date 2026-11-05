#!/bin/bash
# Builds the Airbender platform for bflat guests into the output directory
# (the first argument; default: output/airbender in the repository):
#
#   - airbender_entry.o: entry.S assembled; the _start of the guest, the end
#     of a run and the console (the platform interface of bflat's pal)
#   - airbender.ld: the memory map (the linker script)
#   - airbender.bflat.manifest: the platform manifest (docs/manifest.md in
#     bflat-riscv64)
#
# It needs only clang with the riscv32 target, not the bindings library, so a
# guest that uses no accelerator links with the platform alone. build/build.sh
# airbender calls it for the package.
set -euo pipefail

PLATFORM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
OUTPUT_DIR="${1:-$(cd "${PLATFORM_DIR}/.." && pwd -P)/output/airbender}"
mkdir -p "${OUTPUT_DIR}"

# rv32im and the ilp32 ABI: what the guest links, without the C extension.
clang --target=riscv32-unknown-elf -march=rv32im -mabi=ilp32 -mno-relax -c \
    "${PLATFORM_DIR}/entry.S" -o "${OUTPUT_DIR}/airbender_entry.o"
cp "${PLATFORM_DIR}/airbender.ld" "${PLATFORM_DIR}/airbender.bflat.manifest" "${OUTPUT_DIR}/"

echo "Platform: ${OUTPUT_DIR}/airbender_entry.o airbender.ld airbender.bflat.manifest"
