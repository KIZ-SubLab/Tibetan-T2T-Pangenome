#!/usr/bin/env bash
set -euo pipefail

: "${BASE:?Set BASE to the TMRCA analysis directory}"
ENV_PREFIX=${ENV_PREFIX:-/work/home/songjiaxu/anaconda3/envs/tmrca_fig5b}
export PATH="$ENV_PREFIX/bin:$PATH"

GVCF=${GVCF:-$BASE/gvcf/T2T012.chr1.chm13.g.vcf.gz}
CHAIN=${CHAIN:-$BASE/gvcf/chm13v2-grch38.chain}
STRICT_MASK=${STRICT_MASK:-$BASE/../TRACE_find.super.archic/resources/accessibility_mask/1000G_GRCh38_StrictMask_20160622/20160622.allChr.mask.bed}
TARGET_BED=${TARGET_BED:-$BASE/config/target.GRCh38.bed}
OUTDIR=${OUTDIR:-$BASE/callable_mask}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUTDIR"

CHM13_MASK="$OUTDIR/T2T012.chr1.chm13.GQ20.DP10-60.callable.bed"
GRCH38_RAW="$OUTDIR/T2T012.chr1.GRCh38.GQ20.DP10-60.callable.raw.bed"
GRCH38_MASK="$OUTDIR/T2T012.chr1.GRCh38.GQ20.DP10-60.callable.bed"
UNMAPPED="$OUTDIR/T2T012.chr1.chm13.GQ20.DP10-60.unmapped.bed"
FULL_STRICT="$OUTDIR/T2T012.chr1.GRCh38.GQ20.DP10-60.strict.callable.bed"
FINAL="$OUTDIR/T2T012.target.GRCh38.GQ20.DP10-60.strict.callable.bed"

if [[ ! -s "$CHM13_MASK" ]]; then
  python "$SCRIPT_DIR/02_make_callable_mask.py" \
    --gvcf "$GVCF" --sample T2T012 --contig chr1 \
    --min-gq 20 --min-depth 10 --max-depth 60 \
    --output "$CHM13_MASK" --stats "$OUTDIR/callable_build_stats.json"
fi

# BED liftover converts both ends of every reference-confidence block. This is
# intentionally separate from GATK LiftoverVcf, which drops gVCF INFO/END.
if [[ ! -s "$GRCH38_MASK" ]]; then
  liftOver -minMatch=0.95 "$CHM13_MASK" "$CHAIN" "$GRCH38_RAW" "$UNMAPPED"
  bedtools sort -i "$GRCH38_RAW" | bedtools merge -i - > "$GRCH38_MASK"
fi
bedtools intersect -a "$GRCH38_MASK" -b "$STRICT_MASK" \
  | bedtools sort -i - | bedtools merge -i - > "$FULL_STRICT"
bedtools intersect -a "$FULL_STRICT" -b "$TARGET_BED" \
  | bedtools sort -i - | bedtools merge -i - > "$FINAL"

awk 'BEGIN{OFS="\t"; print "mask","intervals","bases"}
     FNR==1{n=0;b=0} {n++; b+=($3-$2)}
     ENDFILE{print FILENAME,n,b}' \
  "$CHM13_MASK" "$GRCH38_MASK" "$UNMAPPED" "$FULL_STRICT" "$TARGET_BED" "$FINAL" \
  > "$OUTDIR/mask_coverage.tsv"

bgzip -f -k "$FINAL"
tabix -f -p bed "$FINAL.gz"
