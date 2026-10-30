#!/bin/bash
# The gates on what build/build.sh made: every object is RV32 ilp32 without the C extension
# (check_ilp32.sh of dotnet-riscv), and carries no compressed, atomic or fence instruction
# (tools/isa_scan.py of dotnet-riscv), as Airbender decodes RV32IM only. DOTNET_RISCV names the
# dotnet-riscv checkout; the image of the Dockerfile has one.
# Usage: build/check.sh <platform>...
set -euo pipefail

TOP="$(cd "$(dirname "$0")/.." && pwd -P)"
DOTNET_RISCV="${DOTNET_RISCV:?DOTNET_RISCV must name a dotnet-riscv checkout}"
[ "$#" -gt 0 ] || { echo "usage: $0 <platform>..." >&2; exit 2; }

dirs=() files=()
for platform in "$@"; do
    out="$TOP/output/$platform"
    [ -f "$out/libairbender.a" ] || { echo "no $out/libairbender.a: run build/build.sh $platform" >&2; exit 1; }
    dirs+=("$out")
    files+=("$out/libairbender.a")
    [ ! -f "$out/harness" ] || files+=("$out/harness")
done
"$DOTNET_RISCV/check_ilp32.sh" "${dirs[@]}"
OBJDUMP=llvm-objdump ISA_SCAN_FENCE=1 python3 "$DOTNET_RISCV/tools/isa_scan.py" "${files[@]}"
