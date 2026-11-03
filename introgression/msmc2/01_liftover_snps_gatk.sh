#!/usr/bin/env bash
set -euo pipefail

: "${BASE:?Set BASE to the TMRCA analysis directory}"
ENV_PREFIX=${ENV_PREFIX:-/work/home/songjiaxu/anaconda3/envs/tmrca_fig5b}
export PATH="$ENV_PREFIX/bin:$PATH"

GVCF=${GVCF:-$BASE/gvcf/T2T012.chr1.chm13.g.vcf.gz}
CHAIN=${CHAIN:-$BASE/gvcf/chm13v2-grch38.chain}
REFERENCE=${REFERENCE:-/work/share/acc2159rgq/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa}
OUTDIR=${OUTDIR:-$BASE/gvcf/liftover_snps}
mkdir -p "$OUTDIR"

SOURCE="$OUTDIR/T2T012.chr1.chm13.PASS.biSNV.vcf.gz"
TARGET="$OUTDIR/T2T012.chr1.GRCh38.PASS.biSNV.vcf.gz"
REJECT="$OUTDIR/T2T012.chr1.chm13.PASS.biSNV.rejected.vcf.gz"

# Split multiallelic records, discard the gVCF symbolic allele, and retain only
# genuine PASS biallelic SNP calls. This is an ordinary VCF suitable for GATK.
bcftools norm -m -any "$GVCF" -Ou \
  | bcftools view -v snps -f PASS -e 'ALT="<*>" || ALT="<NON_REF>"' -m2 -M2 -Oz -o "$SOURCE"
bcftools index -f -t "$SOURCE"

gatk --java-options "-Xmx8g" LiftoverVcf \
  -I "$SOURCE" \
  -O "$TARGET" \
  -CHAIN "$CHAIN" \
  -REJECT "$REJECT" \
  -R "$REFERENCE" \
  --CREATE_INDEX true \
  --LIFTOVER_MIN_MATCH 1.0 \
  >"$OUTDIR/gatk.stdout.log" 2>"$OUTDIR/gatk.stderr.log"

bcftools stats "$SOURCE" > "$OUTDIR/source.stats.txt"
bcftools stats "$TARGET" > "$OUTDIR/target.stats.txt"
bcftools stats "$REJECT" > "$OUTDIR/rejected.stats.txt"
bcftools norm -f "$REFERENCE" -c e "$TARGET" -Ou >/dev/null

{
  echo -e "file\trecords"
  echo -e "source_chm13\t$(bcftools view -H "$SOURCE" | wc -l)"
  echo -e "target_GRCh38\t$(bcftools view -H "$TARGET" | wc -l)"
  echo -e "rejected_chm13\t$(bcftools view -H "$REJECT" | wc -l)"
} > "$OUTDIR/liftover_counts.tsv"
