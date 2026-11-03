#!/usr/bin/env bash
# Shared configuration for assembly quality control (QC).
# Edit the paths below before running the numbered scripts.

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/project}"

# Final polished assemblies:
#   ${ASSEMBLY_DIR}/${child}/${family}.hap1.np2.fa
#   ${ASSEMBLY_DIR}/${child}/${family}.hap2.np2.fa
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${PROJECT_ROOT}/05_nextpolish2}"

# Sample sheet: family_id, paternal_id, maternal_id, child_id (tab or space
# separated). Lines starting with '#' are ignored.
SAMPLE_SHEET="${SAMPLE_SHEET:-${PROJECT_ROOT}/trio_info.tsv}"

# Output directory for QC results.
OUT_QC="${OUT_QC:-${PROJECT_ROOT}/QC_results}"

# Intermediate files produced by the assembly pipeline.
NGS_QC_DIR="${NGS_QC_DIR:-${PROJECT_ROOT}/01_NGS_QC}"
HIFI_QC_DIR="${HIFI_QC_DIR:-${PROJECT_ROOT}/02_hifi_filter_QC}"

# References used by QUAST.
REF_CHM13="${REF_CHM13:-/path/to/chm13v2.0.fa}"
REF_HG38="${REF_HG38:-/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa}"
EST_REF_SIZE="${EST_REF_SIZE:-3100000000}"

THREADS="${THREADS:-30}"

# Set to "true" to run inspector-correct.py after the initial inspector run.
RUN_INSPECTOR_CORRECT="${RUN_INSPECTOR_CORRECT:-false}"
