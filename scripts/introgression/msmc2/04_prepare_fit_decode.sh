#!/usr/bin/env bash
set -euo pipefail

: "${BASE:?Set BASE to the TMRCA analysis directory}"
ENV_PREFIX=${ENV_PREFIX:-/work/home/songjiaxu/anaconda3/envs/tmrca_fig5b}
export PATH="$ENV_PREFIX/bin:$PATH"

PHASED_VCF=${PHASED_VCF:?}
SAMPLE=${SAMPLE:-TBN2}
OUTDIR=${OUTDIR:-$BASE/msmc2_TBN2_hap2_vs_hap1}
MSMC_TOOLS=${MSMC_TOOLS:-$BASE/software/msmc-tools}
FULL_MASK=${FULL_MASK:-$BASE/callable_mask/T2T012.chr1.GRCh38.GQ20.DP10-60.strict.callable.bed}
TARGET_START=${TARGET_START:-231395000}
TARGET_END=${TARGET_END:-232000000}
THREADS=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$OUTDIR/input" "$OUTDIR/model" "$OUTDIR/decode" "$OUTDIR/qc"

RAW_SAMPLE_VCF="$OUTDIR/input/${SAMPLE}.chr1.GRCh38.phased.biSNV.raw.vcf.gz"
SAMPLE_VCF="$OUTDIR/input/${SAMPLE}.chr1.GRCh38.phased.biSNV.vcf.gz"
DUPLICATE_BED="$OUTDIR/qc/duplicate_multiallelic_positions.bed"
MULTI="$OUTDIR/input/${SAMPLE}.chr1.multihetsep.txt"
MODEL_PREFIX="$OUTDIR/model/${SAMPLE}.chr1.hap2_vs_hap1"
DECODE_ALL="$OUTDIR/decode/${SAMPLE}.chr1.hap2_vs_hap1.posterior.tsv"
DECODE_TARGET="$OUTDIR/decode/${SAMPLE}.target.hap2_vs_hap1.posterior.tsv"

[[ -s "$FULL_MASK" ]] || { echo "Missing full callable mask: $FULL_MASK" >&2; exit 2; }

bcftools view -r chr1 -s "$SAMPLE" -c 1 -m2 -M2 -v snps \
  -e 'GT="mis"' "$PHASED_VCF" -Oz -o "$RAW_SAMPLE_VCF"
bcftools index -f -t "$RAW_SAMPLE_VCF"

# The cohort call set represents a few triallelic sites as two biallelic rows
# at the same coordinate. Exclude those coordinates rather than passing two
# incompatible observations to MSMC2.
bcftools query -f '%CHROM\t%POS\n' "$RAW_SAMPLE_VCF" \
  | awk 'prev==$2 {print $1, $2-1, $2} {prev=$2}' OFS='\t' \
  | sort -k1,1 -k2,2n -u > "$DUPLICATE_BED"
if [[ -s "$DUPLICATE_BED" ]]; then
  bcftools view -T ^"$DUPLICATE_BED" "$RAW_SAMPLE_VCF" -Oz -o "$SAMPLE_VCF"
else
  bcftools view "$RAW_SAMPLE_VCF" -Oz -o "$SAMPLE_VCF"
fi
bcftools index -f -t "$SAMPLE_VCF"

python "$MSMC_TOOLS/generate_multihetsep.py" \
  --mask "$FULL_MASK" --as_phased "$SAMPLE_VCF" > "$MULTI"

awk 'BEGIN{n=0;bad=0} /^#/ {next} {
       n++; if(length($4)!=2) bad++; if($2<=prev) unsorted++; prev=$2
     } END{print "records\t"n; print "bad_haplotype_columns\t"bad;
           print "nonincreasing_positions\t"unsorted}' "$MULTI" \
  > "$OUTDIR/qc/multihetsep_qc.tsv"

# Fit the pairwise model across all callable chr1 sequence. Haplotype indices
# 1 and 0 correspond to TBN2 hap2 and hap1; the pairwise model is symmetric.
msmc2 -I 0-1 -t "$THREADS" -o "$MODEL_PREFIX" "$MULTI"

MUTATION_RATE=$(awk '/mutationRate:/ {print $NF; exit}' "$MODEL_PREFIX.log")
LAST_LOOP=$(awk 'NF{line=$0} END{print line}' "$MODEL_PREFIX.loop.txt")
RECOMBINATION_RATE=$(printf '%s\n' "$LAST_LOOP" | cut -f1)
LAMBDA_VEC=$(printf '%s\n' "$LAST_LOOP" | cut -f4)
DECODER_LAMBDA_VEC=$(awk -F, -v m="$MUTATION_RATE" 'BEGIN{OFS=","} {for(i=1;i<=NF;i++){if(i>1)printf OFS; printf "%.17g", $i*m} print ""}' <<< "$LAMBDA_VEC")
NR_STATES=$(awk -F, '{print NF}' <<< "$DECODER_LAMBDA_VEC")

[[ -n "$MUTATION_RATE" && -n "$RECOMBINATION_RATE" && -n "$LAMBDA_VEC" ]] || {
  echo "Could not extract fitted decoder parameters" >&2; exit 3;
}

msmc2_decode -m "$MUTATION_RATE" -r "$RECOMBINATION_RATE" \
  -t "$NR_STATES" -l "$DECODER_LAMBDA_VEC" -I 0,1 --nrThreads "$THREADS" -s 1000 "$MULTI" \
  > "$DECODE_ALL" 2> "$OUTDIR/decode/decoder.stderr.log"
[[ -s "$DECODE_ALL" ]] || { echo "Decoder produced an empty output" >&2; exit 4; }

awk -v s="$TARGET_START" -v e="$TARGET_END" \
  'BEGIN{OFS="\t"} /^#/ {print; next} $1>=s && $1<=e {print}' \
  "$DECODE_ALL" > "$DECODE_TARGET"

{
  echo -e "parameter\tvalue\tsource"
  echo -e "sample\t${SAMPLE}\tconfig"
  echo -e "haplotype_pair\thap2(index1) vs hap1(index0); symmetric\tVCF GT allele order"
  echo -e "target_GRCh38\tchr1:${TARGET_START}-${TARGET_END}\tconfig"
  echo -e "MSMC2_scaled_mutation_rate\t${MUTATION_RATE}\tMSMC2 log"
  echo -e "MSMC2_fitted_recombination_rate\t${RECOMBINATION_RATE}\tlast EM iteration"
  echo -e "MSMC2_output_lambda_vector\t${LAMBDA_VEC}\tlast EM iteration; internal lambda divided by scaled mutation rate"
  echo -e "decoder_lambda_vector\t${DECODER_LAMBDA_VEC}\tMSMC2 output lambda multiplied by scaled mutation rate"
  echo -e "decoder_stride_bp\t1000\tconfig"
} > "$OUTDIR/analysis_parameters.tsv"

bcftools stats "$SAMPLE_VCF" > "$OUTDIR/qc/sample_vcf.stats.txt"
awk 'BEGIN{n=0;b=0}{n++;b+=$3-$2}END{print "intervals\tbases";print n"\t"b}' "$FULL_MASK" \
  > "$OUTDIR/qc/full_mask_coverage.tsv"
awk 'BEGIN{n=0} !/^#/{n++}END{print "target_decoded_rows\t"n}' "$DECODE_TARGET" \
  > "$OUTDIR/qc/decoder_qc.tsv"
