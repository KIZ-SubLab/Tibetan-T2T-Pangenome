#!/usr/bin/env bash
# Step 06: Long-read polishing of the phased assemblies with nextPolish2.
#
# HiFi reads are mapped to each haplotype with minimap2; nextPolish2 corrects
# base-level errors using child k-mer sets (k21 and k31). The polished contigs
# are the final assemblies used for downstream QC and analysis.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_POLISH}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_POLISH}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    hifi="${OUT_HIFI_QC}/${child}/${child}.hifi_reads.filt.fastq.gz"
    k21="${OUT_NGS_QC}/${child}/${child}.k21.yak"
    k31="${OUT_NGS_QC}/${child}/${child}.k31.yak"

    for hap in 1 2; do
        genome="${OUT_HIFIASM}/${child}/${family}_hifiasm_trio_ont.dip.hap${hap}_ctg.fa"
        bam="${family}.hap${hap}.lgs.bam"

        minimap2 -ax map-hifi -t "${THREADS}" "${genome}" "${hifi}" \
            | samtools sort -@ "${THREADS}" -m 1G -o "${bam}" -
        samtools index "${bam}"

        nextPolish2 -t "${NP2_THREADS}" "${bam}" "${genome}" "${k21}" "${k31}" \
            > "${family}.hap${hap}.np2.fa"
    done
done < "${SAMPLE_SHEET}"
