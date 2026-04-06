# syntax=docker/dockerfile:1.7
# version 0.0.3

# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm as builder
WORKDIR /apps

# INFO: copy sources from the repo-root build context
COPY assets/rust/splicing/Cargo.toml assets/rust/splicing/Cargo.lock ./splicing/
COPY assets/rust/splicing/src ./splicing/src

COPY assets/rust/aparent/Cargo.toml assets/rust/aparent/Cargo.lock ./aparent/
COPY assets/rust/aparent/src ./aparent/src


# INFO: build
RUN cargo build --release --manifest-path ./splicing/Cargo.toml
RUN cargo build --release --manifest-path ./aparent/Cargo.toml

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    openssh-client \
    rsync \
    && rm -rf /var/lib/apt/lists/*


# INFO: copy binaries to bin section
COPY --from=builder /apps/splicing/target/release/splicing /usr/local/bin/splicing
COPY --from=builder /apps/aparent/target/release/aparent /usr/local/bin/aparent

# INFO: run help
RUN splicing --help
RUN aparent --help
