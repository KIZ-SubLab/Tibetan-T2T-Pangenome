#!/usr/bin/env bash
# Assembly-based SV calling with PAV (v1.1.2)
set -euo pipefail

# ---- edit paths ----
REF=/path/to/reference.fa
HAP1=/path/to/sample.hap1.fasta
HAP2=/path/to/sample.hap2.fasta
SAMPLE=sample
THREADS=32
MASK=/path/to/cen_tel_mask.bed
# --------------------

# 1. SV calling on both haplotype assemblies
pav --threads "${THREADS}" \
    --reference "${REF}" \
    --vcf "${SAMPLE}.pav.vcf" \
    "${HAP1}" "${HAP2}"

# 2. Keep SVs >=50 bp
bcftools view -i 'SVLEN>=50 || SVLEN<=-50' "${SAMPLE}.pav.vcf" > "${SAMPLE}.pav.precise.vcf"

# 3. Exclude centromeres, telomeres and flanking regions
bcftools view -T "^${MASK}" "${SAMPLE}.pav.precise.vcf" > "${SAMPLE}.pav.clean.vcf"
