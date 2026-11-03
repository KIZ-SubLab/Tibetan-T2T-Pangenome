#!/usr/bin/env bash
set -euo pipefail

: "${BASE:?Set BASE to the TMRCA analysis directory}"
ENV_PREFIX=${ENV_PREFIX:-/work/home/songjiaxu/anaconda3/envs/tmrca_fig5b}
export PATH="$ENV_PREFIX/bin:$PATH"

ANALYSIS=${ANALYSIS:-$BASE/msmc2_TBN2_hap2_vs_hap1}
POSTERIOR="$ANALYSIS/decode/TBN2.target.hap2_vs_hap1.posterior.tsv"
PARAMETERS="$ANALYSIS/analysis_parameters.tsv"
OUTDIR="$ANALYSIS/figure"
mkdir -p "$OUTDIR"

python "$BASE/scripts/05_plot_tmrca.py" \
  --posterior "$POSTERIOR" \
  --parameters "$PARAMETERS" \
  --physical-mutation-rate 1.25e-8 \
  --generation-years 29 \
  --y-max-mya 2.0 \
  --title "" \
  --output-prefix "$OUTDIR/TBN2_hap2_vs_hap1_TMRCA_Fig5B"
