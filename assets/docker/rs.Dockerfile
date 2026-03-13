# syntax=docker/dockerfile:1.7
# version 0.0.1

# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm as builder
WORKDIR /apps

# INFO: moving codebase
COPY ../rust/splicing/Cargo.toml ../rust/splicing/Cargo.lock ./splicing/
COPY ../rust/splicing/src ./splicing/src

# INFO: build
RUN cargo build --release --manifest-path ./splicing/Cargo.toml

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*


# INFO: copy binaries to bin section
COPY --from=builder /apps/splicing/target/release/splicing /usr/local/bin/splicing

# INFO: run help
RUN splicing --help
