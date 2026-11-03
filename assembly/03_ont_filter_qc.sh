#!/usr/bin/env bash
# Step 03: ONT read filtering and statistics.
#
# Reads >= 100 kb are retained as ultra-long reads and supplied to hifiasm
# with --ul. Both the full pass reads and the ultra-long subset are summarised
# with NanoStat.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_ONT_QC}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_ONT_QC}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    zcat "${RAW_ONT}/${child}"/*pass.fastq.gz \
        | NanoFilt -l 100000 \
        | gzip > "${child}.pass_100k.fastq.gz"

    NanoStat --fastq "${RAW_ONT}/${child}"/*pass.fastq.gz \
        -t "${THREADS}" --tsv -o "${sample_dir}" -n "${child}.pass.stat"
    NanoStat --fastq "${child}.pass_100k.fastq.gz" \
        -t "${THREADS}" --tsv -o "${sample_dir}" -n "${child}.pass_100k.stat"
done < "${SAMPLE_SHEET}"
