#!/usr/bin/env bash
# Read-based SV calling with minimap2 + Sniffles2 (v2.3.3)
set -euo pipefail

# ---- edit paths ----
REF=/path/to/reference.fa
READS=/path/to/sample.hifi_reads.fastq.gz
SAMPLE=sample
THREADS=32
MASK=/path/to/cen_tel_mask.bed
# --------------------

# 1. Alignment to the reference
minimap2 -ax map-hifi -t "${THREADS}" "${REF}" "${READS}" \
    | samtools sort -m 7g --threads "${THREADS}" -o "${SAMPLE}.hifi.bam"
samtools index "${SAMPLE}.hifi.bam"

# 2. SV calling
sniffles --threads "${THREADS}" \
    --input "${SAMPLE}.hifi.bam" \
    --reference "${REF}" \
    --vcf "${SAMPLE}.sniffles.vcf" \
    --snf "${SAMPLE}.sniffles.snf"

# 3. Keep SVs >=50 bp with precise breakpoints
bcftools view -i 'SVLEN>=50 || SVLEN<=-50' "${SAMPLE}.sniffles.vcf" \
    | awk '$1 ~ /#/ || ($8 !~ "IMPRECISE")' > "${SAMPLE}.sniffles.precise.vcf"

# 4. Exclude centromeres, telomeres and flanking regions
bcftools view -T "^${MASK}" "${SAMPLE}.sniffles.precise.vcf" > "${SAMPLE}.sniffles.clean.vcf"
