#!/usr/bin/env bash
# Step 05: Assembly statistics for each phased haplotype.
#
# seqkit reports contiguity statistics; yak trioeval reports parental k-mer
# inheritance (W rate ~ single-base/switch errors, H rate ~ heterozygosity).
# Dot plots against a reference are produced only when REFERENCE_FASTA and
# UNIMAP_DOTPLOT are set in config.sh.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_HIFIASM}/${child}"
    pat_yak="${OUT_NGS_QC}/${paternal}/${paternal}.k31.yak"
    mat_yak="${OUT_NGS_QC}/${maternal}/${maternal}.k31.yak"

    for hap in 1 2; do
        fa="${sample_dir}/${family}_hifiasm_trio_ont.dip.hap${hap}_ctg.fa"
        prefix="${fa%.fa}"

        seqkit stats -aT -j "${THREADS}" -t dna "${fa}" > "${prefix}.seqkit.stats"
        yak trioeval -t "${THREADS}" "${pat_yak}" "${mat_yak}" "${fa}" \
            > "${prefix}.trioeval"

        if [[ -n "${REFERENCE_FASTA}" && -n "${UNIMAP_DOTPLOT}" ]]; then
            sh "${UNIMAP_DOTPLOT}" "${REFERENCE_FASTA}" "${fa}" "${prefix}"
        fi
    done
done < "${SAMPLE_SHEET}"
