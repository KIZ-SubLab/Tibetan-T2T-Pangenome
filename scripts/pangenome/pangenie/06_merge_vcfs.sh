#!/usr/bin/env bash
# Compress per-sample VCFs and merge them (in groups of 250).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

TYPING_DIR="${WORKDIR}/typing"
OUTDIR="${WORKDIR}/results/genotypes"
LIST="${OUTDIR}/vcf_list.txt"
MERGED="${OUTDIR}/all_samples.merged.vcf.gz"
mkdir -p "${OUTDIR}"

while IFS= read -r -d '' VCF; do
    bgzip -f "${VCF}"
    bcftools index -t "${VCF}.gz"
done < <(find "${TYPING_DIR}" -name '*_genotyping.vcf' -print0)

mapfile -t VCFS < <(find "${TYPING_DIR}" -name '*_genotyping.vcf.gz' | sort)
[[ ${#VCFS[@]} -gt 0 ]] || { echo "no per-sample VCFs found in ${TYPING_DIR}"; exit 1; }
printf '%s\n' "${VCFS[@]}" > "${LIST}"

GROUP_DIR="${OUTDIR}/groups_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${GROUP_DIR}"
split -l 250 -d -a 3 "${LIST}" "${GROUP_DIR}/group_"

for GROUP_LIST in "${GROUP_DIR}"/group_*; do
    mapfile -t GROUP < "${GROUP_LIST}"
    bcftools merge --threads "${THREADS}" -Oz \
        -o "${GROUP_LIST}.vcf.gz" "${GROUP[@]}"
    bcftools index -t "${GROUP_LIST}.vcf.gz"
done

mapfile -t GROUPS < <(find "${GROUP_DIR}" -name 'group_*.vcf.gz' | sort)
bcftools merge --threads "${THREADS}" -Oz -o "${MERGED}" "${GROUPS[@]}"
bcftools index -t "${MERGED}"

echo "output: ${MERGED}"
echo "samples: $(bcftools query -l "${MERGED}" | wc -l)"
echo "records: $(bcftools view -H "${MERGED}" | wc -l)"
