# Multi-stage build to optimize final image size
FROM rust:1.95 AS builder

# Install system dependencies required for compilation
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    cmake \
    libclang-dev \
    clang \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /app

# Copy project files
COPY Cargo.toml Cargo.lock ./
COPY src/ ./src/

# Compile in release mode. The image runs on machines other than the one that
# builds it, so no target-cpu=native (an image built on a newer CPU would stop
# with "illegal instruction" on an older one): x86-64-v2 (SSE4.2) runs on every
# x86-64 server of the last 15 years; parasail picks its SIMD kernels at run
# time. For the fastest local build use RUSTFLAGS="-C target-cpu=native".
ARG RUST_TARGET_CPU=x86-64-v2
RUN RUSTFLAGS="-C target-cpu=${RUST_TARGET_CPU}" cargo build --release

# Final stage with smaller image
FROM debian:bookworm-slim

# Install only necessary runtime dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN useradd -r -s /bin/false cgdist

# Copy compiled binaries (cgdist; cgdist-cache for cache stores; cgdist-diff
# to inspect one allele pair: use --entrypoint to run the helpers)
COPY --from=builder /app/target/release/cgdist /usr/local/bin/cgdist
COPY --from=builder /app/target/release/cgdist-cache /usr/local/bin/cgdist-cache
COPY --from=builder /app/target/release/cgdist-diff /usr/local/bin/cgdist-diff

# Ensure binaries are executable
RUN chmod +x /usr/local/bin/cgdist /usr/local/bin/cgdist-cache /usr/local/bin/cgdist-diff

# Create directory for data
RUN mkdir -p /data && chown cgdist:cgdist /data

# Switch to non-root user
USER cgdist

# Working directory for data
WORKDIR /data

# Entry point
ENTRYPOINT ["/usr/local/bin/cgdist"]

# Default help if no arguments are passed
CMD ["--help"]

# Metadata
LABEL maintainer="andrea.deruvo@gssi.it"
LABEL description="cgDist: Ultra-fast SNP/indel-level distance calculator for core genome MLST analysis"
LABEL version="0.1.0"
LABEL org.opencontainers.image.title="cgDist"
LABEL org.opencontainers.image.description="High-performance tool for calculating pairwise SNP and indel distances between bacterial isolates using core genome MLST allelic profiles with sequence alignment"
LABEL org.opencontainers.image.vendor="GenPat-IT Bioinformatics"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.source="https://github.com/genpat-it/cgDist"
LABEL org.opencontainers.image.documentation="https://github.com/genpat-it/cgDist"