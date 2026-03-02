#!/usr/bin/env bash
# run_swam-g_local.sh
# Run SWAM-g end-to-end on a local workstation.
# Usage: bash run_swam-g_local.sh [input_dir] [output_dir]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INPUT="${1:-$(pwd)/input}"
OUTPUT="${2:-$(pwd)/output}"

echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo ""

snakemake \
    --profile config/local/ \
    --config in_dir="$INPUT" out_dir="$OUTPUT"
