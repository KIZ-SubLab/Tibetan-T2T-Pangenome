#!/usr/bin/env bash
# QC 02: Detection of base- and structure-level assembly errors with
# inspector, using the filtered HiFi reads.
#
# Set RUN_INSPECTOR_CORRECT=true in config.sh to additionally run
# inspector-correct.py and re-evaluate the corrected assembly.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_QC}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_QC}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    hifi="${HIFI_QC_DIR}/${child}/${child}.hifi_reads.filt.fastq.gz"

    for hap in hap1 hap2; do
        fa="${ASSEMBLY_DIR}/${child}/${family}.${hap}.np2.fa"
        out="${family}.${hap}.np2.inspector"

        inspector.py -c "${fa}" -r "${hifi}" -o "${out}" \
            --datatype hifi -t "${THREADS}"

        if [[ "${RUN_INSPECTOR_CORRECT}" == "true" ]]; then
            inspector-correct.py -i "${out}/" --datatype pacbio-hifi \
                -o "${out}.corrected" -t "${THREADS}"
            inspector.py -c "${out}.corrected/contig_corrected.fa" \
                -r "${hifi}" -o "${out}.corrected.inspector" \
                --datatype hifi -t "${THREADS}"
        fi
    done
done < "${SAMPLE_SHEET}"
