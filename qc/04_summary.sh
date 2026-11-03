#!/usr/bin/env bash
# QC 04: Collect QV, SE/HE, inspector and QUAST metrics into one table.
#
# Inputs are the outputs of QC steps 01-03. QUAST reports hap1 and hap2 in
# columns 2 and 3, respectively; "Genome fraction (%)" occupies columns 4-5.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

summary="${OUT_QC}/T2T_QC_summary.tsv"
printf 'family\thap\tQV_k21\tSE_k19\tHE_k19\tinspector_contigs\tinspector_total_length\tinspector_N50\tinspector_mapping_rate\tinspector_structural_errors\tinspector_QV\tN50_chm13\tNG50_chm13\tgenome_fraction_chm13\tN50_hg38\tNG50_hg38\tgenome_fraction_hg38\n' > "${summary}"

while read -r family paternal maternal child; do
    [[ -z "${family:-}" || "${family}" == \#* ]] && continue

    sample_dir="${OUT_QC}/${child}"

    for hap in hap1 hap2; do
        qv_file="${sample_dir}/${family}.${hap}.np2.qv.k21.txt"
        tri_file="${sample_dir}/${family}.${hap}.np2.trioeval.k19.txt"
        ins_file="${sample_dir}/${family}.${hap}.np2.inspector/summary_statistics"
        chm_report="${sample_dir}/${family}.np2.chm13.quast/report.txt"
        hg_report="${sample_dir}/${family}.np2.hg38.quast/report.txt"

        qv="$(awk '$1 == "QV" {v = $3} END {print v}' "${qv_file}")"
        se="$(awk '$1 == "W" {v = $4} END {print v}' "${tri_file}")"
        he="$(awk '$1 == "H" {v = $4} END {print v}' "${tri_file}")"

        ins_contigs="$(awk -F'\t' '$1 == "Number of contigs" {print $2}' "${ins_file}")"
        ins_len="$(awk -F'\t' '$1 == "Total length" {print $2}' "${ins_file}")"
        ins_n50="$(awk -F'\t' '$1 == "N50" {print $2}' "${ins_file}")"
        ins_map="$(awk -F'\t' '$1 == "Mapping rate /%" {print $2}' "${ins_file}")"
        ins_struct="$(awk -F'\t' '$1 == "Structural error" {print $2}' "${ins_file}")"
        ins_qv="$(awk -F'\t' '$1 == "QV" {print $2}' "${ins_file}")"

        if [[ "${hap}" == "hap1" ]]; then
            col=2
            gf_col=4
        else
            col=3
            gf_col=5
        fi

        chm_n50="$(awk -v c="${col}" '$1 == "N50" {print $c}' "${chm_report}")"
        chm_ng50="$(awk -v c="${col}" '$1 == "NG50" {print $c}' "${chm_report}")"
        chm_gf="$(awk -v c="${gf_col}" '$1 == "Genome" && $2 == "fraction" {print $c}' "${chm_report}")"

        hg_n50="$(awk -v c="${col}" '$1 == "N50" {print $c}' "${hg_report}")"
        hg_ng50="$(awk -v c="${col}" '$1 == "NG50" {print $c}' "${hg_report}")"
        hg_gf="$(awk -v c="${gf_col}" '$1 == "Genome" && $2 == "fraction" {print $c}' "${hg_report}")"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${family}" "${hap}" "${qv}" "${se}" "${he}" \
            "${ins_contigs}" "${ins_len}" "${ins_n50}" "${ins_map}" \
            "${ins_struct}" "${ins_qv}" \
            "${chm_n50}" "${chm_ng50}" "${chm_gf}" \
            "${hg_n50}" "${hg_ng50}" "${hg_gf}" >> "${summary}"
    done
done < "${SAMPLE_SHEET}"
