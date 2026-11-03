#!/usr/bin/env bash
# Build the PanGenie index.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

VCF="${WORKDIR}/results/pangenome/pangenome.vcf"
INDEX_DIR="${WORKDIR}/results/index"
PREFIX="${INDEX_DIR}/pangenie"

[[ -s "${VCF}" ]] || { echo "missing ${VCF}; run 02_prepare_pangenome.sh first"; exit 1; }
mkdir -p "${INDEX_DIR}"

pangenie_exec PanGenie-index \
    -v "${VCF}" \
    -r "${REFERENCE}" \
    -o "${PREFIX}" \
    -t "${THREADS}"

test -s "${PREFIX}_UniqueKmersMap.cereal"
test -s "${PREFIX}_path_segments.fasta"
echo "index: ${PREFIX}"
