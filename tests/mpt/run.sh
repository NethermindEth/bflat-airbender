#!/bin/bash
# The trie engine of src/mpt.c on the host, with the address and undefined-behaviour sanitizers:
# mpt_host.c gives it a portable Keccak and a witness held in arrays, mpt_test.py compares its
# reads, writes and roots with a trie built from scratch. Usage: tests/mpt/run.sh [rounds] [seed]
set -euo pipefail

TOP="$(cd "$(dirname "$0")/../.." && pwd -P)"
out="$TOP/output/host"
mkdir -p "$out"
${CC:-cc} -O1 -g -shared -fPIC -fsanitize=address,undefined -fno-sanitize-recover=undefined \
    -I "$TOP/src" -o "$out/libmpt.so" "$TOP/tests/mpt/mpt_host.c"
asan="$(${CC:-cc} -print-file-name=libasan.so)"
LD_PRELOAD="$asan" ASAN_OPTIONS=detect_leaks=0 python3 "$TOP/tests/mpt/mpt_test.py" "$out/libmpt.so" "${1:-200}" "${2:-1}"
