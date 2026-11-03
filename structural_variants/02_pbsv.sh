#!/usr/bin/env bash
# Read-based SV calling with pbmm2 + PBSV (v2.9.0)
set -euo pipefail

# ---- edit paths ----
REF=/path/to/reference.fa
READS=/path/to/sample.hifi_reads.fastq.gz
SAMPLE=sample
THREADS=32
MASK=/path/to/cen_tel_mask.bed
# --------------------

# 1. Alignment to the reference
pbmm2 align -j "${THREADS}" "${REF}" "${READS}" "${SAMPLE}.hifi.bam" \
    --sort --preset HIFI --sample "${SAMPLE}" \
    --rg "@RG\tID:${SAMPLE}\tSM:${SAMPLE}"

# 2. SV calling
pbsv discover "${SAMPLE}.hifi.bam" "${SAMPLE}.svsig.gz" --hifi
pbsv call -j "${THREADS}" "${REF}" "${SAMPLE}.svsig.gz" "${SAMPLE}.pbsv.vcf" --hifi

# 3. Keep SVs >=50 bp with precise breakpoints
bcftools view -i 'SVLEN>=50 || SVLEN<=-50' "${SAMPLE}.pbsv.vcf" \
    | awk '$1 ~ /#/ || ($8 !~ "IMPRECISE")' > "${SAMPLE}.pbsv.precise.vcf"

# 4. Exclude centromeres, telomeres and flanking regions
bcftools view -T "^${MASK}" "${SAMPLE}.pbsv.precise.vcf" > "${SAMPLE}.pbsv.clean.vcf"
