#!/usr/bin/env bash
# Run SWAM-g in resilient HPC mode and always emit a post-run status summary.
# Usage: bash run_swam-g_resilient.sh [input_dir] [output_dir] [extra snakemake args...]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INPUT="${1:-$(pwd)/input}"
OUTPUT="${2:-$(pwd)/output}"
shift $(( $# >= 1 ? 1 : 0 ))
shift $(( $# >= 1 ? 1 : 0 ))

mkdir -p "$OUTPUT/run_status"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
DRIVER_LOG="$OUTPUT/run_status/snakemake-resilient.$TIMESTAMP.log"

echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo "Profile: config/slurm/resilient/"
echo "Driver log: $DRIVER_LOG"
echo ""

snakemake \
    --profile config/slurm/resilient/ \
    --config in_dir="$INPUT" out_dir="$OUTPUT" \
    --local-cores 1 \
    "$@" 2>&1 | tee "$DRIVER_LOG"

SNAKEMAKE_EXIT=${PIPESTATUS[0]}

python3 workflow/scripts/build_resilient_report.py \
    --input-dir "$INPUT" \
    --output-dir "$OUTPUT" \
    --driver-log "$DRIVER_LOG"

REPORT_EXIT=$?

echo ""
echo "Run report:         $OUTPUT/run_status/run_report.txt"
echo "Sample summary:     $OUTPUT/run_status/sample_status.tsv"
echo "Rule status table:  $OUTPUT/run_status/sample_rule_status.tsv"
echo "Retry samples:      $OUTPUT/run_status/retry_samples.txt"
echo "Review samples:     $OUTPUT/run_status/review_samples.txt"
echo "Retry manifest:     $OUTPUT/run_status/retry_manifest.tsv"
echo "Review manifest:    $OUTPUT/run_status/review_manifest.tsv"
echo "Completed manifest: $OUTPUT/run_status/completed_manifest.tsv"

if [ "$SNAKEMAKE_EXIT" -ne 0 ]; then
    exit "$SNAKEMAKE_EXIT"
fi

exit "$REPORT_EXIT"
