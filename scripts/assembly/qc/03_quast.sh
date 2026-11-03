#!/usr/bin/env bash
# QC 03: Assembly evaluation with QUAST against CHM13v2 and GRCh38.
#
# Both haplotypes are evaluated together so that hap1 and hap2 are reported in
# consecutive columns of each report.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_QC}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_QC}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    hap1="${ASSEMBLY_DIR}/${child}/${family}.hap1.np2.fa"
    hap2="${ASSEMBLY_DIR}/${child}/${family}.hap2.np2.fa"

    quast.py -t "${THREADS}" "${hap1}" "${hap2}" \
        -r "${REF_CHM13}" -o "${family}.np2.chm13.quast" \
        --large --est-ref-size "${EST_REF_SIZE}" --no-icarus

    quast.py -t "${THREADS}" "${hap1}" "${hap2}" \
        -r "${REF_HG38}" -o "${family}.np2.hg38.quast" \
        --large --est-ref-size "${EST_REF_SIZE}" --no-icarus
done < "${SAMPLE_SHEET}"
