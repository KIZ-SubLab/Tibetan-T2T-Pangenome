#!/usr/bin/env bash
# vcfbub: keep only reference-reachable alleles (drop alleles > 100 kb, level > 0).
set -euo pipefail

SIF=/path/to/cactus.sif
OUT=vcfbub

mkdir -p "$OUT"

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    IN=mc/$CHR/graph/$CHR-mc.raw.vcf.gz

    singularity exec "$SIF" vcfbub --input "$IN" -l 0 -a 100000 \
        | bgzip -@ 8 > "$OUT/$CHR.vcfbub.vcf.gz"
    bcftools index -t "$OUT/$CHR.vcfbub.vcf.gz"
done
