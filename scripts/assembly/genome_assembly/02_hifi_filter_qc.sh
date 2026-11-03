#!/usr/bin/env bash
# Step 02: HiFi adapter/contaminant filtering and read statistics.
#
# HiFiAdapterFilt removes adapter and vector contamination. The filtered FASTQ
# is used for hifiasm assembly, assembly evaluation and nextPolish2 mapping.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_HIFI_QC}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_HIFI_QC}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    # HiFiAdapterFilt expects <prefix>.hifi_reads.bam in the working directory.
    ln -sf "${RAW_HIFI}/${child}/${child}.hifi_reads.bam" "${sample_dir}/"
    hifiadapterfilt.sh -l 44 -m 97 -t "${THREADS}" -p "${child}"

    seqkit stats -aT -j "${THREADS}" -t dna \
        "${child}.hifi_reads.filt.fastq.gz" > "${child}.hifi.seqkit.stats"
    samtools stats -@ "${THREADS}" \
        "${child}.hifi_reads.bam" > "${child}.hifi.samtools.stats"
done < "${SAMPLE_SHEET}"
