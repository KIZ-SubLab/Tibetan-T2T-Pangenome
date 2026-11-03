#!/usr/bin/env bash

set -eo pipefail

ANALYSIS_ROOT="${ANALYSIS_ROOT:?ANALYSIS_ROOT is required}"
ARCHAIC_LABEL="${ARCHAIC_LABEL:?ARCHAIC_LABEL is required}"
ARCHAIC_CODE="${ARCHAIC_CODE:?ARCHAIC_CODE is required}"
INCLUDE_CEU="${INCLUDE_CEU:-0}"
MANIFEST="${ANALYSIS_ROOT}/sv_manifest.tsv"

source /work/home/songjiaxu/anaconda3/etc/profile.d/conda.sh
conda activate plink
set -u

Rscript "${ANALYSIS_ROOT}/scripts/summarize_batch.R" \
  "${ANALYSIS_ROOT}" "${MANIFEST}" "${ARCHAIC_LABEL}" "${ARCHAIC_CODE}"
Rscript "${ANALYSIS_ROOT}/scripts/summarize_hetero200.R" \
  "${ANALYSIS_ROOT}" "${MANIFEST}" "${ARCHAIC_LABEL}" "${ARCHAIC_CODE}"
Rscript "${ANALYSIS_ROOT}/scripts/summarize_direct_phased200.R" \
  "${ANALYSIS_ROOT}" "${MANIFEST}" "${ARCHAIC_LABEL}" "${ARCHAIC_CODE}"

if [[ "${INCLUDE_CEU}" == "1" ]]; then
  Rscript "${ANALYSIS_ROOT}/scripts_ceu35/summarize_ceu35.R" \
    "${ANALYSIS_ROOT}" "${MANIFEST}"
fi

Rscript "${ANALYSIS_ROOT}/scripts/summarize_integrated.R" \
  "${ANALYSIS_ROOT}" "${MANIFEST}" "${ARCHAIC_LABEL}" \
  "${ARCHAIC_CODE}" "${INCLUDE_CEU}"

