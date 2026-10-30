# The build image of libairbender: build/ci.sh.
#
# The x64 cross-build container of the .NET 11 runtime pipelines, which
# dotnet-riscv builds the riscv32 runtime and uGC builds ugc.o in: clang, lld
# and the LLVM tools at one pinned version. At DOTNET_RISCV_REF, the
# dotnet-riscv scripts make the rv32im musl sysroot the linux platform and the
# test harness compile against (SYSROOT), and provide the gates of
# build/check.sh (DOTNET_RISCV).
FROM mcr.microsoft.com/dotnet-buildtools/prereqs:azurelinux-3.0-net11.0-cross-amd64@sha256:ca817af24ff68c801a3e31f16a0fbf6bb7ba57764b0aeae5be260dbc70ae25cb

ARG DOTNET_RISCV_REF=57af2904b5dfa86900065b12f5d0ece0d7664024

RUN git init -q /dotnet-riscv && cd /dotnet-riscv \
    && git fetch -q --depth 1 https://github.com/NethermindEth/dotnet-riscv "${DOTNET_RISCV_REF}" \
    && git checkout -q FETCH_HEAD \
    && ./build_musl_rv32im.sh /rv32 \
    && ./build_builtins_rv32im.sh /rv32

ENV SYSROOT=/rv32
ENV DOTNET_RISCV=/dotnet-riscv
