#!/usr/bin/env bash
# Run SWAM-g on Slurm with the default resilient-style profile and emit a
# concise post-run status summary.
# Usage: bash run_swam-g.sh [input_dir] [output_dir] [SWAM-g options] [extra snakemake args...]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

INPUT="${1:-$(pwd)/input}"
OUTPUT="${2:-$(pwd)/output}"
shift $(( $# >= 1 ? 1 : 0 ))
shift $(( $# >= 1 ? 1 : 0 ))

usage() {
    cat <<'EOF'
Usage: bash run_swam-g.sh [input_dir] [output_dir] [SWAM-g options] [extra snakemake args...]

SWAM-g options:
  --pd-lookup=true|false
  --pd-backend=ftp|table
  --pd-comparator-limit=N
  --pd-sample-metadata-tsv=PATH
  --pd-isolates-tsv=PATH
  --pd-exceptions-tsv=PATH
  -h, --help

Any other arguments are forwarded to snakemake unchanged.
EOF
}

SNAKEMAKE_ARGS=()
CONFIG_OVERRIDES=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --pd-lookup=*)
            CONFIG_OVERRIDES+=("pd_lookup=${1#*=}")
            ;;
        --pd-lookup)
            shift
            CONFIG_OVERRIDES+=("pd_lookup=${1:?Missing value for --pd-lookup}")
            ;;
        --pd-backend=*)
            CONFIG_OVERRIDES+=("pd_backend=${1#*=}")
            ;;
        --pd-backend)
            shift
            CONFIG_OVERRIDES+=("pd_backend=${1:?Missing value for --pd-backend}")
            ;;
        --pd-comparator-limit=*)
            CONFIG_OVERRIDES+=("pd_comparator_limit=${1#*=}")
            ;;
        --pd-comparator-limit)
            shift
            CONFIG_OVERRIDES+=("pd_comparator_limit=${1:?Missing value for --pd-comparator-limit}")
            ;;
        --pd-sample-metadata-tsv=*)
            CONFIG_OVERRIDES+=("pd_sample_metadata_tsv=${1#*=}")
            ;;
        --pd-sample-metadata-tsv)
            shift
            CONFIG_OVERRIDES+=("pd_sample_metadata_tsv=${1:?Missing value for --pd-sample-metadata-tsv}")
            ;;
        --pd-isolates-tsv=*)
            CONFIG_OVERRIDES+=("pd_isolates_tsv=${1#*=}")
            ;;
        --pd-isolates-tsv)
            shift
            CONFIG_OVERRIDES+=("pd_isolates_tsv=${1:?Missing value for --pd-isolates-tsv}")
            ;;
        --pd-exceptions-tsv=*)
            CONFIG_OVERRIDES+=("pd_exceptions_tsv=${1#*=}")
            ;;
        --pd-exceptions-tsv)
            shift
            CONFIG_OVERRIDES+=("pd_exceptions_tsv=${1:?Missing value for --pd-exceptions-tsv}")
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            ;;
    esac
    shift
done

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
    --config in_dir="$INPUT" out_dir="$OUTPUT" "${CONFIG_OVERRIDES[@]}" \
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
