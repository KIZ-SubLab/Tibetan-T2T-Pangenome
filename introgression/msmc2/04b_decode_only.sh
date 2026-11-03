#!/usr/bin/env bash
set -euo pipefail

: "${BASE:?Set BASE to the TMRCA analysis directory}"
ENV_PREFIX=${ENV_PREFIX:-/work/home/songjiaxu/anaconda3/envs/tmrca_fig5b}
export PATH="$ENV_PREFIX/bin:$PATH"

OUTDIR=${OUTDIR:-$BASE/msmc2_TBN2_hap2_vs_hap1}
THREADS=${SLURM_CPUS_PER_TASK:-4}
MULTI="$OUTDIR/input/TBN2.chr1.multihetsep.txt"
MODEL_PREFIX="$OUTDIR/model/TBN2.chr1.hap2_vs_hap1"
DECODE_ALL="$OUTDIR/decode/TBN2.chr1.hap2_vs_hap1.posterior.tsv"
DECODE_TARGET="$OUTDIR/decode/TBN2.target.hap2_vs_hap1.posterior.tsv"

MUTATION_RATE=$(awk '/mutationRate:/ {print $NF; exit}' "$MODEL_PREFIX.log")
LAST_LOOP=$(awk 'NF{line=$0} END{print line}' "$MODEL_PREFIX.loop.txt")
RECOMBINATION_RATE=$(printf '%s\n' "$LAST_LOOP" | cut -f1)
LAMBDA_VEC=$(printf '%s\n' "$LAST_LOOP" | cut -f4)
DECODER_LAMBDA_VEC=$(awk -F, -v m="$MUTATION_RATE" 'BEGIN{OFS=","} {for(i=1;i<=NF;i++){if(i>1)printf OFS; printf "%.17g", $i*m} print ""}' <<< "$LAMBDA_VEC")
NR_STATES=$(awk -F, '{print NF}' <<< "$DECODER_LAMBDA_VEC")

msmc2_decode -m "$MUTATION_RATE" -r "$RECOMBINATION_RATE" \
  -t "$NR_STATES" -l "$DECODER_LAMBDA_VEC" -I 0,1 --nrThreads "$THREADS" -s 1000 "$MULTI" \
  > "$DECODE_ALL" 2> "$OUTDIR/decode/decoder.stderr.log"
[[ -s "$DECODE_ALL" ]] || { echo "Decoder produced an empty output" >&2; exit 4; }

awk 'BEGIN{OFS="\t"} /^#/ {print; next} $1>=231395000 && $1<=232000000 {print}' \
  "$DECODE_ALL" > "$DECODE_TARGET"
[[ $(awk '!/^#/{n++}END{print n+0}' "$DECODE_TARGET") -gt 0 ]] || {
  echo "No decoded rows fall in the target region" >&2; exit 5;
}

{
  echo -e "parameter\tvalue\tsource"
  echo -e "sample\tTBN2\tconfig"
  echo -e "haplotype_pair\thap2(index1) vs hap1(index0); symmetric\tVCF GT allele order"
  echo -e "target_GRCh38\tchr1:231395000-232000000\tconfig"
  echo -e "MSMC2_scaled_mutation_rate\t${MUTATION_RATE}\tMSMC2 log"
  echo -e "MSMC2_fitted_recombination_rate\t${RECOMBINATION_RATE}\tlast EM iteration"
  echo -e "MSMC2_output_lambda_vector\t${LAMBDA_VEC}\tlast EM iteration; internal lambda divided by scaled mutation rate"
  echo -e "decoder_lambda_vector\t${DECODER_LAMBDA_VEC}\tMSMC2 output lambda multiplied by scaled mutation rate"
  echo -e "decoder_states\t${NR_STATES}\tlambda vector length"
  echo -e "decoder_stride_bp\t1000\tconfig"
} > "$OUTDIR/analysis_parameters.tsv"

awk 'BEGIN{n=0} !/^#/{n++}END{print "target_decoded_rows\t"n}' "$DECODE_TARGET" \
  > "$OUTDIR/qc/decoder_qc.tsv"
