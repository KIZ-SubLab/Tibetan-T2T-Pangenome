#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash submit_relate_clues.sh CONFIG.env" >&2
  exit 2
fi
CONFIG=$(readlink -f "$1")
# shellcheck disable=SC1090
source "$CONFIG"
: "${BASE:?BASE is required in config}"
: "${OUTPUT_DIR:?OUTPUT_DIR is required in config}"
: "${ANALYSIS_TAG:?ANALYSIS_TAG is required in config}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
bash "$SCRIPT_DIR/validate_config.sh" "$CONFIG"
mkdir -p "$OUTPUT_DIR/logs"
job_name=$(printf 'relcl_%s' "$ANALYSIS_TAG" | tr -c 'A-Za-z0-9_.-' '_' | cut -c1-100)
job_id=$(sbatch --parsable \
  --job-name="$job_name" \
  --account="${SLURM_ACCOUNT:-acc2159rgq}" \
  --partition="${SLURM_PARTITION:-wzacnormal02}" \
  --nodes=1 --ntasks=1 \
  --cpus-per-task="${THREADS:-20}" \
  --mem="${SLURM_MEM:-128G}" \
  --time="${SLURM_TIME:-72:00:00}" \
  --output="$OUTPUT_DIR/logs/%x.%j.out" \
  --error="$OUTPUT_DIR/logs/%x.%j.err" \
  "$SCRIPT_DIR/slurm_entry.sh" "$CONFIG")
printf '%s\n' "$job_id"
printf 'squeue -j %s\n' "$job_id"
