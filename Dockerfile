# syntax=docker/dockerfile:1
FROM python:3.10-slim

# Install OS-level dependencies for pybedtools, pysam, skmisc, numpy, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    pkg-config \
    gfortran \
    clang \
    libclang-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libreadline-dev \
    libffi-dev \
    libxml2-dev \
    libxslt1-dev \
    libopenblas-dev \
    liblapack-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Rust toolchain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install uv (fast Python package/dependency manager)
RUN pip install --no-cache-dir uv

# Set workdir
WORKDIR /app

# Copy project files
COPY . /app

# Install dependencies and package in a uv venv
RUN uv venv .venv && \
    . .venv/bin/activate && \
    uv pip install .

# Set environment variables for uv venv activation
ENV PATH="/app/.venv/bin:$PATH"
ENV PYTHONPATH="/app"

# Default command: show CLI help
CMD ["krewlyzer", "--help"]
