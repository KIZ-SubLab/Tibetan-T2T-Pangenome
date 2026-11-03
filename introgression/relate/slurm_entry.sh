#!/usr/bin/env bash
set -euo pipefail

CONFIG=$(readlink -f "$1")
# Slurm executes a copied entry script from its spool directory. Resolve the
# main script through BASE from the submitted config.
# shellcheck disable=SC1090
source "$CONFIG"
: "${BASE:?BASE is required in config}"
exec bash "$BASE/scripts/run_relate_clues.sh" "$CONFIG"

