#!/usr/bin/env bash
# Step 01: Illumina read trimming and k-mer counting.
#
# fastp removes adapters and low-quality bases. yak k-mer counts are built for
# downstream steps:
#   k19, parents -> yak trioeval (SE/HE)
#   k21, child   -> yak qv and nextPolish2
#   k31, parents -> hifiasm trio phasing
#   k31, child   -> nextPolish2
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_NGS_QC}"

count_kmers() {
    local sample="$1" k="$2" out_dir="$3"
    yak count -k"${k}" -b37 -t "${THREADS}" \
        -o "${out_dir}/${sample}.k${k}.yak" \
        <(zcat "${out_dir}/${sample}_1.clean.fq.gz") \
        <(zcat "${out_dir}/${sample}_2.clean.fq.gz")
}

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    for sample in "${paternal}" "${maternal}" "${child}"; do
        sample_dir="${OUT_NGS_QC}/${sample}"
        mkdir -p "${sample_dir}"

        fastp -i "${RAW_NGS}/${sample}/${sample}_1.fq.gz" \
              -o "${sample_dir}/${sample}_1.clean.fq.gz" \
              -I "${RAW_NGS}/${sample}/${sample}_2.fq.gz" \
              -O "${sample_dir}/${sample}_2.clean.fq.gz" \
              -j "${sample_dir}/${sample}.fastp.json" \
              -h "${sample_dir}/${sample}.fastp.html" \
              -u 30 -q 20 -w "${THREADS}"

        if [[ "${sample}" == "${child}" ]]; then
            count_kmers "${sample}" 21 "${sample_dir}"
        else
            count_kmers "${sample}" 19 "${sample_dir}"
        fi
        count_kmers "${sample}" 31 "${sample_dir}"
    done
done < "${SAMPLE_SHEET}"
