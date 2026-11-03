#!/usr/bin/env bash
# Step 04: Trio-binned, ultra-long-read-assisted assembly with hifiasm.
#
# -1 / -2 are parental k-mer sets used for haplotype phasing; --ul supplies
# ONT reads >= 100 kb. Each phased haplotype is converted from GFA to FASTA.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUT_HIFIASM}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_HIFIASM}/${child}"
    mkdir -p "${sample_dir}"
    cd "${sample_dir}"

    hifi="${OUT_HIFI_QC}/${child}/${child}.hifi_reads.filt.fastq.gz"
    ont="$(ls "${OUT_ONT_QC}/${child}/${child}.pass_100k.fastq.gz" | paste -sd, -)"
    pat_yak="${OUT_NGS_QC}/${paternal}/${paternal}.k31.yak"
    mat_yak="${OUT_NGS_QC}/${maternal}/${maternal}.k31.yak"

    hifiasm -1 "${pat_yak}" -2 "${mat_yak}" \
            -o "${family}_hifiasm_trio_ont" \
            -t "${THREADS}" \
            --ul "${ont}" \
            "${hifi}"

    for hap in 1 2; do
        gfa="${family}_hifiasm_trio_ont.dip.hap${hap}.p_ctg.gfa"
        fa="${family}_hifiasm_trio_ont.dip.hap${hap}_ctg.fa"
        awk '/^S/{print ">"$2"\n"$3}' "${gfa}" | fold -w 80 > "${fa}"
    done
done < "${SAMPLE_SHEET}"
