#!/usr/bin/env bash
# vcfwave: split each chromosome VCF into 5 Mb windows, realign per window,
# then concatenate and sort per chromosome.
set -euo pipefail

WIN=5000000
THREADS=24
OUT=wave

for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    IN=vcfbub/$CHR.vcfbub.vcf.gz
    WIN_DIR=$OUT/windows/$CHR
    mkdir -p "$WIN_DIR"

    # 1. split into windows
    LEN=$(bcftools view -h "$IN" | awk -v c="$CHR" \
        'index($0, "##contig=<ID=" c ",") == 1 && match($0, /length=[0-9]+/) { print substr($0, RSTART + 7, RLENGTH - 7); exit }')

    for ((S = 1; S <= LEN; S += WIN)); do
        E=$((S + WIN - 1))
        ((E > LEN)) && E=$LEN
        W=$WIN_DIR/$(printf "%s_win_%04d" "$CHR" $(((S - 1) / WIN + 1)))
        bcftools view -r "$CHR:$S-$E" -Oz -o "$W.vcf.gz" "$IN"
        bcftools index -t "$W.vcf.gz"
    done

    # 2. vcfwave per window (serial: parallel instances inside one container deadlock)
    for W in "$WIN_DIR"/*.vcf.gz; do
        vcfwave -t "$THREADS" -I 1000 "$W" | bgzip -@ 8 > "${W%.vcf.gz}.wave.vcf.gz"
        bcftools sort -Oz -o "${W%.vcf.gz}.wave.sorted.vcf.gz" "${W%.vcf.gz}.wave.vcf.gz"
        mv "${W%.vcf.gz}.wave.sorted.vcf.gz" "${W%.vcf.gz}.wave.vcf.gz"
        bcftools index -t "${W%.vcf.gz}.wave.vcf.gz"
    done

    # 3. concat windows, then sort (vcfwave can emit a POS rollback at boundaries)
    mkdir -p "$OUT/final"
    ls -1 "$WIN_DIR"/*.wave.vcf.gz | sort -V | xargs bcftools concat -Oz -o "$OUT/final/$CHR.wave.vcf.gz"
    bcftools sort -Oz -o "$OUT/final/$CHR.sorted.vcf.gz" "$OUT/final/$CHR.wave.vcf.gz"
    mv "$OUT/final/$CHR.sorted.vcf.gz" "$OUT/final/$CHR.wave.vcf.gz"
    bcftools index -t "$OUT/final/$CHR.wave.vcf.gz"
done
