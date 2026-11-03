#!/usr/bin/env bash
# =============================================================================
# config.example.sh — copy to config.sh and edit
# All paths are placeholders; no absolute paths are shipped with the scripts.
# =============================================================================

# ---- Part I : telomere detection and architecture -------------------------
export ASSEMBLY_DIR="/path/to/assemblies"          # <sample>.fasta
export VGP_PIPELINE_DIR="/path/to/vgp-assembly/pipeline"
export TELOMERE_WORK="/path/to/work/telomere"      # outputs of 01_

# ---- Part II : TVR / TAR1 -------------------------------------------------
export TELO_FASTA_DIR="${TELOMERE_WORK}/telo_fasta"
export TVR_WORK="/path/to/work/tvr"
export VPG_BED_DIR="/path/to/vgp/windows_bed"      # <sample>.windows.0.4.15kb.ends.bed
export RM_GFF_DIR="/path/to/repeatmasker/gff"      # <sample>.gff
export ASSEMBLY_FAI_DIR="/path/to/assembly/fai"    # <sample>.fai  (windowed GFF only)
export EXTRACT_FAI_DIR="/path/to/extract/fai"      # <sample>.fai  (windowed GFF only)

# ---- sample lists (one sample id per line) --------------------------------
export SAMPLE_LIST_ALL="/path/to/samples_all.txt"
export SAMPLE_LIST_A="/path/to/samples_group_a.txt"
export SAMPLE_LIST_B="/path/to/samples_group_b.txt"

# ---- detection parameters --------------------------------------------------
export VGP_THRESHOLD="0.4"        # repeat-density threshold
export VGP_ENDS="15000"           # terminal search window (bp)
export VGP_MERGE_DIST="100"       # bedtools merge -d
export EXPECTED_ENDS_PER_HAPLOID="46"

# ---- classification parameters --------------------------------------------
export CFB_MIN_LENGTH_BP="2000"   # B-type minimum length
export CFB_NONCANON_PCT="0.05"    # B-type non-canonical fraction
export CF_THRESHOLD="0.75"        # C vs F main-unit ratio cut-off

# ---- TVR parameters --------------------------------------------------------
export TVR_MIN_BLOCK="6"
export TVR_MAX_BLOCK="300"

# ---- TAR1 parameters -------------------------------------------------------
export TAR1_WINDOW="10000"        # towards chromosome interior only

echo "config loaded: VGP threshold=${VGP_THRESHOLD}, ends=${VGP_ENDS}, merge=${VGP_MERGE_DIST}"
