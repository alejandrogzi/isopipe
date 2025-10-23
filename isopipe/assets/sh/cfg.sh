#!/usr/bin/env bash

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.0.1

set -euo pipefail

# INFO: check if a command exists
check_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "Error: '$1' is not installed or not in PATH." >&2
        exit 1
    fi

    echo "INFO: '$1' is installed."
}

echo "INFO: [1/9] Checking dependencies..."
check_command cargo
check_command rustc
check_command uv

echo "INFO: [2/9] Building isopipe and installing it locally..."
cd isopipe
cargo build --release && cargo install --path .

echo "INFO: [3/9] Building isotools..."
cd ../isotools/isotools
cargo build --release

echo "INFO: [4/9] Building orf + setting up Python env..."
cd ../../isopipe/assets/rust/orf
cargo build --release

cd tai
uv venv
# shellcheck disable=SC1091
source .venv/bin/activate
uv pip install ".[cpu]"

echo "INFO: [5/9] Building extract..."
cd ../../extract
cargo build --release

echo "INFO: [6/9] Building collapse..."
cd ../collapse
cargo build --release

echo "INFO: [7/9] Building orfipy + venv"
cd ../py/orfipy
uv venv
source .venv/bin/activate
uv pip install "."

echo "INFO: [8/9] Building predict venv"
cd ../predict
uv venv
source .venv/bin/activate
uv pip sync

echo "INFO: [9/10] Returning to repo root..."
cd ../../../../

echo "INFO: [10/10] Running test..."
echo "INFO: Configuration completed successfully!"
