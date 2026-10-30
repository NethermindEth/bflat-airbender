#!/bin/bash
# The build of CI and of the release: both platforms and the gates in the image of the Dockerfile,
# then tests/run.py with the harness under qemu-riscv32 (QEMU names it) on the host.
set -euo pipefail

TOP="$(cd "$(dirname "$0")/.." && pwd -P)"
image=bflat-airbender-build

docker build --platform linux/amd64 -t "$image" "$TOP"
docker run --platform linux/amd64 --rm -v "$TOP:$TOP" -w "$TOP" "$image" bash -euc '
    bash build/build.sh airbender
    bash build/build.sh linux
    bash build/check.sh airbender linux
    bash tests/platform/run-tls.sh
    chmod -R a+rwX output'
python3 "$TOP/tests/run.py" "${QEMU:-qemu-riscv32}" "$TOP/output/linux/harness"
