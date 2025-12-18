# syntax=docker/dockerfile:1

# ============ BUILD STAGE ============
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git pkg-config clang libclang-dev \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    curl && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install maturin with patchelf support
RUN pip install --no-cache-dir "maturin[patchelf]" uv

WORKDIR /build

# Copy project files
COPY pyproject.toml README.md LICENSE ./
COPY rust/ rust/
COPY src/ src/

# Build the wheel
RUN maturin build --release --out /wheels --manifest-path rust/Cargo.toml

# ============ RUNTIME STAGE ============
FROM python:3.10-slim AS runtime

# OCI Labels for GitHub Container Registry
LABEL org.opencontainers.image.title="krewlyzer"
LABEL org.opencontainers.image.description="A comprehensive toolkit for ctDNA fragmentomics analysis from GRCh37 aligned BAM files"
LABEL org.opencontainers.image.url="https://github.com/msk-access/krewlyzer"
LABEL org.opencontainers.image.source="https://github.com/msk-access/krewlyzer"
LABEL org.opencontainers.image.vendor="MSK-ACCESS"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.authors="Ronak Shah <shahr2@mskcc.org>"

# Runtime dependencies only (amd64 has pre-built wheels)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 libssl3 zlib1g libbz2-1.0 liblzma5 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm /tmp/*.whl

ENTRYPOINT ["krewlyzer"]
CMD ["--help"]
