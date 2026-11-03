#!/usr/bin/env bash
# Assembly-based SV calling with minimap2 + SVIM-asm (v1.0.3)
set -euo pipefail

# ---- edit paths ----
REF=/path/to/reference.fa
PAT=/path/to/sample.paternal.fasta
MAT=/path/to/sample.maternal.fasta
SAMPLE=sample
THREADS=30
MASK=/path/to/cen_tel_mask.bed
# --------------------

# 1. Align both haplotypes to the reference
minimap2 -a -x asm5 --cs -r2k -t "${THREADS}" "${REF}" "${PAT}" > "${SAMPLE}.pat.sam"
minimap2 -a -x asm5 --cs -r2k -t "${THREADS}" "${REF}" "${MAT}" > "${SAMPLE}.mat.sam"
samtools sort -m4G -@25 -o "${SAMPLE}.pat.sorted.bam" "${SAMPLE}.pat.sam"
samtools sort -m4G -@25 -o "${SAMPLE}.mat.sorted.bam" "${SAMPLE}.mat.sam"
samtools index "${SAMPLE}.pat.sorted.bam"
samtools index "${SAMPLE}.mat.sorted.bam"

# 2. Diploid SV calling
svim-asm diploid "${SAMPLE}_workdir" \
    "${SAMPLE}.pat.sorted.bam" "${SAMPLE}.mat.sorted.bam" "${REF}"

# 3. Keep SVs >=50 bp passing FILTER
bcftools view -i 'SVLEN>=50 || SVLEN<=-50' "${SAMPLE}_workdir/variants.vcf" \
    | bcftools view -f PASS > "${SAMPLE}.svimasm.precise.vcf"

# 4. Exclude centromeres, telomeres and flanking regions
bcftools view -T "^${MASK}" "${SAMPLE}.svimasm.precise.vcf" > "${SAMPLE}.svimasm.clean.vcf"
