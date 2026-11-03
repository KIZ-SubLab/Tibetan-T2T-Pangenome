#!/usr/bin/env bash
# Convert bubble genotypes to biallelic variant genotypes (optional).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

MERGED=${1:-${WORKDIR}/results/genotypes/all_samples.merged.vcf.gz}
BIALLELIC="${WORKDIR}/results/input-vcf/callset-biallelic.vcf.gz"
CONVERTER="${PANGENIE_REPO}/pipelines/run-from-callset/scripts/convert-to-biallelic.py"
OUT="${WORKDIR}/results/genotypes/all_samples.biallelic.vcf.gz"

[[ -s "${MERGED}" ]] || { echo "missing ${MERGED}"; exit 1; }
[[ -s "${BIALLELIC}" ]] || { echo "missing ${BIALLELIC}"; exit 1; }
[[ -s "${CONVERTER}" ]] || { echo "missing converter: ${CONVERTER}"; exit 1; }

zcat "${MERGED}" | python3 "${CONVERTER}" "${BIALLELIC}" | bgzip -c > "${OUT}"
bcftools index -t "${OUT}"

echo "output: ${OUT}"
echo "records: $(bcftools view -H "${OUT}" | wc -l)"
