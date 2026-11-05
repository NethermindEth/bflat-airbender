# The Airbender CLI (the cli of matter-labs/zksync-airbender), which runs an
# Airbender guest on the transpiler VM:
#
#   docker run --rm -v "$PWD:/n" nethermindeth/airbender-runner:<commit> run \
#       --bin /n/nethermind.bin --text /n/nethermind.text \
#       --input-file /n/input.hex --cycles 12000000000 [--expected-results ...]
#
# AIRBENDER_COMMIT is the revision with the 4 MiB ROM, the unsigned RV32M
# machine and the Keccak-f delegation that libairbender calls (see README.md).
# cycles.patch prints the cycles of a run on stderr ("Cycles: <n>"), which
# the CLI does not: a cycle takes 4 steps of the timestamp, which starts at 4.
# The repository pins its nightly toolchain (rust-toolchain.toml), which
# rustup installs by itself.
FROM rust:1.91-bookworm AS build
ARG AIRBENDER_COMMIT
RUN apt-get update \
 && apt-get install -y --no-install-recommends pkg-config libssl-dev \
 && rm -rf /var/lib/apt/lists/*
COPY cycles.patch /cycles.patch
RUN git init -q /src && cd /src \
 && git fetch -q --depth 1 https://github.com/matter-labs/zksync-airbender "${AIRBENDER_COMMIT}" \
 && git checkout -q FETCH_HEAD \
 && git apply /cycles.patch \
 && cargo build --release -p cli

FROM debian:bookworm-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates libssl3 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/cli /usr/local/bin/cli
ENTRYPOINT ["/usr/local/bin/cli"]
