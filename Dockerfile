# Multi-stage build to optimize final image size
FROM rust:1.82 as builder

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

# Compile in release mode with native optimizations
RUN RUSTFLAGS="-C target-cpu=native" cargo build --release

# Final stage with smaller image
FROM debian:bookworm-slim

# Install only necessary runtime dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN useradd -r -s /bin/false cgdist

# Copy compiled binary
COPY --from=builder /app/target/release/cgdist /usr/local/bin/cgdist

# Ensure binary is executable
RUN chmod +x /usr/local/bin/cgdist

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
LABEL description="cgdist: Ultra-fast SNP/indel-level distance calculator for core genome MLST analysis"
LABEL version="1.0.0"
LABEL org.opencontainers.image.title="cgdist"
LABEL org.opencontainers.image.description="High-performance tool for calculating pairwise SNP and indel distances between bacterial isolates using core genome MLST allelic profiles with sequence alignment"
LABEL org.opencontainers.image.vendor="Bioinformatics Laboratory"
LABEL org.opencontainers.image.licenses="MIT"