#!/bin/bash
# The inversion hint of src/modexp.c on the host, with the address and undefined-behaviour sanitizers: the u256
# delegation's C model stands in for the delegation, and the test owns the hint tape. Usage: tests/modexp/run.sh
set -euo pipefail

TOP="$(cd "$(dirname "$0")/../.." && pwd -P)"
out="$TOP/output/host"
mkdir -p "$out"
${CC:-cc} -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -I "$TOP/src" \
    -o "$out/inverse_hint" "$TOP/tests/modexp/inverse_hint.c"
ASAN_OPTIONS=detect_leaks=0 "$out/inverse_hint"
