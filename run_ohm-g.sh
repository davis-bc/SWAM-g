#!/usr/bin/env bash
# Run OHM-g on Slurm with the default resilient-style profile and emit a
# concise post-run status summary.
# Usage: bash run_ohm-g.sh [input_dir] [output_dir] [extra snakemake args...]

set -uo pipefail

usage() {
    cat <<'EOF'
Usage: bash run_ohm-g.sh [input_dir] [output_dir] [extra snakemake args...]

Edit config/slurm/config.yaml to choose optional analyses and PD lookup settings.
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

mkdir -p "$OUTPUT/logs/run_status"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
DRIVER_LOG="$OUTPUT/logs/run_status/snakemake.$TIMESTAMP.log"

echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo "Profile: config/slurm/"
echo "Driver log: $DRIVER_LOG"
echo ""

snakemake \
    --profile config/slurm/ \
    --config in_dir="$INPUT" out_dir="$OUTPUT" \
    --local-cores 1 \
    --quiet rules \
    "${SNAKEMAKE_ARGS[@]}" 2>&1 | tee "$DRIVER_LOG"

SNAKEMAKE_EXIT=${PIPESTATUS[0]}

python3 workflow/scripts/build_resilient_report.py \
    --input-dir "$INPUT" \
    --output-dir "$OUTPUT" \
    --driver-log "$DRIVER_LOG"

REPORT_EXIT=$?

echo ""
echo "Run report:        $OUTPUT/logs/run_status/run_report.txt"
echo "Rule status table: $OUTPUT/logs/run_status/sample_rule_status.tsv"

if [ "$SNAKEMAKE_EXIT" -ne 0 ]; then
    exit "$SNAKEMAKE_EXIT"
fi

exit "$REPORT_EXIT"
