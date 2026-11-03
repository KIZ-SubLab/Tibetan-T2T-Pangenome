#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash run_ghost_fd.sh CONFIG.env" >&2
  exit 2
fi
CONFIG=$(readlink -f "$1")
# shellcheck disable=SC1090
source "$CONFIG"

REGION_ONLY=${REGION_ONLY:-NO}
case "${REGION_ONLY^^}" in
  YES|TRUE|1) REGION_ONLY=YES ;;
  NO|FALSE|0) REGION_ONLY=NO ;;
  *) echo "ERROR: REGION_ONLY must be YES or NO" >&2; exit 2 ;;
esac

required=(BASE ENV_DIR EPAS_FD_BASE REGION_LABEL ANALYSIS_TAG BED KG_VCF FAM TBN_VCF DEN_VCF NEA_VCF ANCESTRAL_FASTA OUTPUT_DIR WINDOW_SNPS STEP_SNPS WINDOW_JACKKNIFE_BLOCK_SNPS REGION_JACKKNIFE_BP)
for var in "${required[@]}"; do
  [[ -n "${!var:-}" ]] || { echo "ERROR: missing config variable $var" >&2; exit 2; }
done
for path in "$BED" "$KG_VCF" "$FAM" "$TBN_VCF" "$DEN_VCF" "$NEA_VCF" "$ANCESTRAL_FASTA"; do
  [[ -s "$path" ]] || { echo "ERROR: missing or empty input: $path" >&2; exit 2; }
done
for vcf in "$KG_VCF" "$TBN_VCF" "$DEN_VCF" "$NEA_VCF"; do
  [[ -s "$vcf.tbi" || -s "$vcf.csi" ]] || { echo "ERROR: indexed VCF required: $vcf" >&2; exit 2; }
done

RESULTS="$OUTPUT_DIR/results"
LOGS="$OUTPUT_DIR/logs"
mkdir -p "$RESULTS" "$LOGS"
cp "$CONFIG" "$OUTPUT_DIR/run_config_snapshot.env"
WINDOW_FILE="$RESULTS/windows_${WINDOW_SNPS}SNP_step${STEP_SNPS}.tsv"
rm -f "$RESULTS/validation.txt" "$RESULTS/analysis_mode.txt"
if [[ "$REGION_ONLY" == YES ]]; then
  rm -f "$WINDOW_FILE" "$RESULTS/genomics_general_B.csv" \
    "$RESULTS/genomics_general_C.csv" \
    "$RESULTS/${ANALYSIS_TAG}_ABC_D_fd.png" \
    "$RESULTS/${ANALYSIS_TAG}_ABC_D_fd.pdf"
fi

build_cmd=(
  "$ENV_DIR/bin/python" "$BASE/scripts/build_ghost_fd.py"
  --kg-vcf "$KG_VCF" --fam "$FAM" --tbn-vcf "$TBN_VCF"
  --den-vcf "$DEN_VCF" --nea-vcf "$NEA_VCF" --bed "$BED"
  --ancestral-fasta "$ANCESTRAL_FASTA" --outdir "$RESULTS"
  --window-snps "$WINDOW_SNPS" --step-snps "$STEP_SNPS"
  --window-jackknife-block-snps "$WINDOW_JACKKNIFE_BLOCK_SNPS"
  --region-jackknife-bp "$REGION_JACKKNIFE_BP"
)
if [[ "$REGION_ONLY" == YES ]]; then
  build_cmd+=(--region-only)
fi

"${build_cmd[@]}" \
  > "$LOGS/build.stdout.log" 2> "$LOGS/build.stderr.log"

validate_extra=()
if [[ "$REGION_ONLY" == YES ]]; then
  validate_extra+=(--region-only)
else
  "$ENV_DIR/bin/python" "$BASE/scripts/make_gg_input.py" \
    --sites "$RESULTS/sites_common_ancestral_polarized.tsv.gz" \
    --geno "$RESULTS/genomics_general_common.haplo.geno.gz" \
    --pops "$RESULTS/genomics_general_common.pops.tsv" \
    > "$LOGS/make_gg.stdout.log" 2> "$LOGS/make_gg.stderr.log"

  GG="$EPAS_FD_BASE/software/genomics_general/ABBABABAwindows.py"
  for test in B C; do
    if [[ "$test" == B ]]; then P3=DEN; else P3=NEA; fi
    "$ENV_DIR/bin/python" "$GG" \
      --windType sites --windSize "$WINDOW_SNPS" \
      --overlap "$((WINDOW_SNPS-STEP_SNPS))" --minSites "$WINDOW_SNPS" \
      -P1 YRI -P2 TBN -P3 "$P3" -O ANC \
      --popsFile "$RESULTS/genomics_general_common.pops.tsv" \
      --genoFile "$RESULTS/genomics_general_common.haplo.geno.gz" \
      --genoFormat haplo --Threads 1 \
      --outFile "$RESULTS/genomics_general_${test}.csv" \
      > "$LOGS/genomics_general_${test}.stdout.log" \
      2> "$LOGS/genomics_general_${test}.stderr.log"
  done
  validate_extra+=(--official-b "$RESULTS/genomics_general_B.csv"
                   --official-c "$RESULTS/genomics_general_C.csv")
fi

"$ENV_DIR/bin/python" "$BASE/scripts/validate_ghost_results.py" \
  --results "$RESULTS" --window-snps "$WINDOW_SNPS" --step-snps "$STEP_SNPS" \
  "${validate_extra[@]}" \
  > "$LOGS/validate.stdout.log" 2> "$LOGS/validate.stderr.log"

if [[ "$REGION_ONLY" == NO ]]; then
  "$ENV_DIR/bin/python" "$BASE/scripts/plot_ghost_fd.py" \
    --windows "$WINDOW_FILE" --summary "$RESULTS/region_summary_ABC.tsv" \
    --region-label "$REGION_LABEL" --out-prefix "$RESULTS/${ANALYSIS_TAG}_ABC_D_fd" \
    > "$LOGS/plot.stdout.log" 2> "$LOGS/plot.stderr.log"
fi

printf '%s\n' "$([[ "$REGION_ONLY" == YES ]] && echo REGION_ONLY || echo WINDOWED)" \
  > "$RESULTS/analysis_mode.txt"

sha256sum "$CONFIG" "$BED" "$FAM" "$KG_VCF" "$TBN_VCF" "$DEN_VCF" "$NEA_VCF" \
  "$ANCESTRAL_FASTA" "$BASE/scripts/build_ghost_fd.py" "$BASE/scripts/make_gg_input.py" \
  "$BASE/scripts/validate_ghost_results.py" "$BASE/scripts/plot_ghost_fd.py" \
  "$BASE/scripts/run_ghost_fd.sh" > "$RESULTS/input_and_script_sha256.txt"

{
  "$ENV_DIR/bin/python" --version
  "$ENV_DIR/bin/bcftools" --version | head -n 2
  "$ENV_DIR/bin/python" -c 'import cyvcf2,matplotlib,numpy,pandas,scipy; print("cyvcf2="+cyvcf2.__version__); print("matplotlib="+matplotlib.__version__); print("numpy="+numpy.__version__); print("pandas="+pandas.__version__); print("scipy="+scipy.__version__)'
  printf 'genomics_general_commit='
  (cd "$EPAS_FD_BASE/software/genomics_general" && git rev-parse HEAD)
} > "$RESULTS/software_versions.txt"

echo "DONE: $REGION_LABEL A+B+C -> $RESULTS"
