#!/usr/bin/env bash
# Genotype every sample listed in the sample sheet.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

[[ -s "${SAMPLE_SHEET}" ]] || { echo "missing sample sheet: ${SAMPLE_SHEET}"; exit 1; }

while IFS=$'\t' read -r SAMPLE R1 R2; do
    [[ -z "${SAMPLE:-}" ]] && continue
    [[ "${SAMPLE}" == \#* ]] && continue
    [[ "${SAMPLE}" == "sample" ]] && continue
    bash "${HERE}/04_genotype_sample.sh" "${SAMPLE}" "${R1}" "${R2}"
done < "${SAMPLE_SHEET}"
