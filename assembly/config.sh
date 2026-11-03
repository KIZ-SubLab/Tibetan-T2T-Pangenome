#!/usr/bin/env bash
# Shared configuration for the genome assembly pipeline.
# Edit the paths below, then run the numbered scripts in order (01 -> 06).

# Project root. All outputs are written below this directory.
PROJECT_ROOT="${PROJECT_ROOT:-/path/to/project}"

# Raw data layout (one directory per sample):
#   ${RAW_NGS}/${sample}/${sample}_1.fq.gz
#   ${RAW_NGS}/${sample}/${sample}_2.fq.gz
#   ${RAW_HIFI}/${sample}/${sample}.hifi_reads.bam
#   ${RAW_ONT}/${sample}/${sample}.pass.fastq.gz
RAW_NGS="${RAW_NGS:-${PROJECT_ROOT}/data_tmp/NGS}"
RAW_HIFI="${RAW_HIFI:-${PROJECT_ROOT}/data_tmp/HIFI}"
RAW_ONT="${RAW_ONT:-${PROJECT_ROOT}/data_tmp/ONT}"

# Sample sheet: family_id, paternal_id, maternal_id, child_id (tab or space
# separated). Lines starting with '#' are ignored.
SAMPLE_SHEET="${SAMPLE_SHEET:-${PROJECT_ROOT}/trio_info.tsv}"

# Output directories.
OUT_NGS_QC="${PROJECT_ROOT}/01_NGS_QC"
OUT_HIFI_QC="${PROJECT_ROOT}/02_hifi_filter_QC"
OUT_ONT_QC="${PROJECT_ROOT}/03_ont_filter_QC"
OUT_HIFIASM="${PROJECT_ROOT}/04_hifiasm"
OUT_POLISH="${PROJECT_ROOT}/05_nextpolish2"

# Threads for each step. Adjust to the resources of the compute node.
THREADS="${THREADS:-32}"
NP2_THREADS="${NP2_THREADS:-5}"

# Optional: reference FASTA and a wrapper script producing a minimap-based
# dot plot (unimap). Leave empty to skip the dot-plot step.
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
UNIMAP_DOTPLOT="${UNIMAP_DOTPLOT:-}"
