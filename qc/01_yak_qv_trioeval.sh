#!/usr/bin/env bash
# QC 01: Consensus quality value (QV) and trio inheritance of polished
# assemblies.
#
# QV is computed with yak qv using child k21 k-mers; the reported value is the
# third field of the "QV" line, as in the original analysis.
# SE/HE are computed with yak trioeval using parental k19 k-mers; the W and H
# rates correspond to SE and HE, respectively.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_QC}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_QC}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    k21="${NGS_QC_DIR}/${child}/${child}.k21.yak"
    pat_k19="${NGS_QC_DIR}/${paternal}/${paternal}.k19.yak"
    mat_k19="${NGS_QC_DIR}/${maternal}/${maternal}.k19.yak"

    for hap in hap1 hap2; do
        fa="${ASSEMBLY_DIR}/${child}/${family}.${hap}.np2.fa"

        yak qv -t "${THREADS}" -p "${k21}" "${fa}" \
            > "${family}.${hap}.np2.qv.k21.txt"
        yak trioeval -t "${THREADS}" -p "${pat_k19}" "${mat_k19}" "${fa}" \
            > "${family}.${hap}.np2.trioeval.k19.txt"
    done
done < "${SAMPLE_SHEET}"
