#!/usr/bin/env bash
# run_ohm-g_local.sh
# Run OHM-g end-to-end on a local workstation.
# Usage: bash run_ohm-g_local.sh [input_dir] [output_dir] [extra snakemake args...]

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash run_ohm-g_local.sh [input_dir] [output_dir] [extra snakemake args...]

Edit config/local/config.yaml to choose optional analyses and PD lookup settings.
Any extra arguments are forwarded to snakemake unchanged, so one-off overrides can
still be passed with --config key=value when needed.
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INPUT="${1:-$(pwd)/input}"
OUTPUT="${2:-$(pwd)/output}"
shift $(( $# >= 1 ? 1 : 0 ))
shift $(( $# >= 1 ? 1 : 0 ))

SNAKEMAKE_ARGS=("$@")

echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo ""

snakemake \
    --profile config/local/ \
    --config in_dir="$INPUT" out_dir="$OUTPUT" \
    "${SNAKEMAKE_ARGS[@]}"
