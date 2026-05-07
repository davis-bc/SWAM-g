#!/usr/bin/env bash
# Launch the OHM-g Slurm profile configuration TUI.
# Usage: bash run_slurm_config_tui.sh [path/to/config.yaml]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON_BIN="${PYTHON:-python3}"
CONFIG_PATH="${1:-config/slurm/config.yaml}"

exec "$PYTHON_BIN" workflow/scripts/slurm_config_tui.py "$CONFIG_PATH"
