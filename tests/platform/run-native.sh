#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../.."
out=output/host
mkdir -p "$out"
cc=${CC:-clang}
flags=(-O2 -ffunction-sections -fdata-sections -Wl,--gc-sections -I src)
"$cc" "${flags[@]}" -fsanitize=address,undefined -fno-sanitize-recover=undefined \
    tests/platform/input.c src/keccak.c src/u256_soft.c -o "$out/input"
ASAN_OPTIONS=detect_leaks=0 "$out/input"
"$cc" "${flags[@]}" tests/pairing-witness.c src/{ecc,bn254,bls12,mont,u256_soft}.c -o "$out/pairing-witness"
# One known non-identity pair on each curve, from the precompile corpus.
awk '$3 == "0" && (($1 == "bnpair" && length($2) == 384) || ($1 == "blspair" && length($2) == 576)) && !seen[$1]++ { print $1, $2 }' \
    tests/vectors/ec.txt | "$out/pairing-witness"
