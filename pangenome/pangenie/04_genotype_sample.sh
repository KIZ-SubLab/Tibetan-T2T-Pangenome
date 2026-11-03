#!/usr/bin/env bash
# Genotype one sample with PanGenie.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

if [[ $# -ne 3 ]]; then
    echo "usage: bash 04_genotype_sample.sh SAMPLE R1.fastq.gz R2.fastq.gz"
    exit 1
fi

SAMPLE=$1
R1=$2
R2=$3
PREFIX="${WORKDIR}/results/index/pangenie"
OUTDIR="${WORKDIR}/typing/${SAMPLE}"
FASTA="${OUTDIR}/${SAMPLE}.fasta"

[[ -s "${PREFIX}_UniqueKmersMap.cereal" ]] || {
    echo "missing index; run 03_build_index.sh first"; exit 1;
}
mkdir -p "${OUTDIR}"

case "${R1}" in
    *.gz) CAT=(zcat) ;;
    *)    CAT=(cat) ;;
esac

"${CAT[@]}" "${R1}" "${R2}" | \
    awk 'NR%4==1{printf ">%s\n", substr($0,2)} NR%4==2{print}' > "${FASTA}"

pangenie_exec PanGenie \
    -f "${PREFIX}" \
    -i "${FASTA}" \
    -s "${SAMPLE}" \
    -o "${OUTDIR}/${SAMPLE}" \
    -j "${THREADS}" \
    -t "${THREADS}"

if [[ "${KEEP_FASTA:-0}" != "1" ]]; then
    rm -f "${FASTA}"
fi

echo "output: ${OUTDIR}/${SAMPLE}_genotyping.vcf"
