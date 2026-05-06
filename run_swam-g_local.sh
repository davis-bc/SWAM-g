#!/usr/bin/env bash
# Backward-compatible wrapper for the old SWAM-g local entrypoint name.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "run_swam-g_local.sh is deprecated; use run_ohm-g_local.sh instead." >&2
exec bash "$SCRIPT_DIR/run_ohm-g_local.sh" "$@"
