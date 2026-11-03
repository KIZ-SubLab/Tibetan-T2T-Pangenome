#!/usr/bin/env bash
# =============================================================================
# 01_telomere_detection.sh
# Telomere identification, recovery and sequence extraction.
#
# Pipeline:
#   assembly FASTA
#     -> VGP FindTelomereWindows (identity 99.9, repeat-density threshold,
#        terminal-terminal search window) -> raw telomeric windows
#     -> bedtools merge -d <MERGE_DIST>
#     -> keep windows within <ENDS> bp of a contig end
#     -> drop unplaced contigs (Unctg)
#     -> seqkit subseq -> per-chromosome-end telomeric FASTA
#
# The detection kernel (1-kbp windows, 200-bp step, canonical TTAGGG >= 10%)
# is implemented inside the VGP telomere.jar; this script wraps the official
# VGP driver so that detection parameters stay identical to the published run.
#
# Requires: VGP_PIPELINE_DIR (dir containing telomere/telomere_analysis.sh
#           and telomere/telomere.jar), java, bedtools, seqkit, samtools
# =============================================================================
set -euo pipefail

FASTA_DIR=""
OUT_DIR=""
SAMPLE_LIST=""
VGP_PIPELINE_DIR="${VGP_PIPELINE_DIR:-}"
THRESHOLD="0.4"
ENDS="15000"
MERGE_DIST="100"
EXCLUDE_CONTIG_REGEX="Unctg"
JOBS="1"

usage() {
    cat <<EOF
Usage: $(basename "$0") --fasta-dir DIR --sample-list FILE --output-dir DIR [options]

Required:
  --fasta-dir DIR        directory containing <sample>.fasta assemblies
  --sample-list FILE     one sample id per line
  --output-dir DIR       output directory

Options:
  --vgp-pipeline DIR     VGP assembly pipeline directory (or set VGP_PIPELINE_DIR)
  --threshold FLOAT      repeat-density threshold (default: 0.4)
  --ends INT             terminal search window in bp (default: 15000)
  --merge-dist INT       bedtools merge -d value (default: 100)
  --exclude REGEX        contig name regex to drop (default: Unctg)
  --jobs INT             parallel samples (default: 1)
  -h, --help

Outputs (under OUTPUT_DIR):
  windows/<sample>.windows.<threshold>.bed       merged telomeric windows
  telo_bed/<sample>.telo.bed                     chrom-end telomeric intervals
  telo_fasta/<sample>.telo.fasta                 extracted telomeric sequence
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta-dir)     FASTA_DIR="$2"; shift 2 ;;
        --sample-list)   SAMPLE_LIST="$2"; shift 2 ;;
        --output-dir)    OUT_DIR="$2"; shift 2 ;;
        --vgp-pipeline)  VGP_PIPELINE_DIR="$2"; shift 2 ;;
        --threshold)     THRESHOLD="$2"; shift 2 ;;
        --ends)          ENDS="$2"; shift 2 ;;
        --merge-dist)    MERGE_DIST="$2"; shift 2 ;;
        --exclude)       EXCLUDE_CONTIG_REGEX="$2"; shift 2 ;;
        --jobs)          JOBS="$2"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -z "$FASTA_DIR" || -z "$SAMPLE_LIST" || -z "$OUT_DIR" ]] && { usage; exit 1; }
[[ -z "$VGP_PIPELINE_DIR" ]] && { echo "ERROR: --vgp-pipeline (or VGP_PIPELINE_DIR) is required" >&2; exit 1; }
VGP_DRIVER="${VGP_PIPELINE_DIR}/telomere/telomere_analysis.sh"
[[ -f "$VGP_DRIVER" ]] || { echo "ERROR: VGP driver not found: $VGP_DRIVER" >&2; exit 1; }

mkdir -p "${OUT_DIR}"/{windows,telo_bed,telo_fasta,logs}

detect_one() {
    local sample="$1"
    local fasta="${FASTA_DIR}/${sample}.fasta"
    [[ -f "$fasta" ]] || { echo "[skip] missing ${fasta}" >&2; return 0; }

    local tmp
    tmp="$(mktemp -d)"
    cp "$fasta" "${tmp}/${sample}.fasta"

    # --- VGP detection (same parameters as the published run) ---
    ( cd "$tmp" && bash "$VGP_DRIVER" "$sample" "$THRESHOLD" "$ENDS" "${sample}.fasta" ) \
        > "${OUT_DIR}/logs/${sample}.detect.log" 2>&1

    local bed
    bed="$(find "$tmp" -name "*.windows.${THRESHOLD}.${ENDS}.ends.bed" | head -1)"
    if [[ -z "$bed" ]]; then
        echo "[warn] ${sample}: no telomeric window detected" >&2
        rm -rf "$tmp"; return 0
    fi

    # --- drop unplaced contigs ---
    grep -v "$EXCLUDE_CONTIG_REGEX" "$bed" > "${OUT_DIR}/telo_bed/${sample}.telo.bed"

    # --- extract telomeric sequence ---
    seqkit subseq --bed "${OUT_DIR}/telo_bed/${sample}.telo.bed" "$fasta" \
        > "${OUT_DIR}/telo_fasta/${sample}.telo.fasta" 2>/dev/null

    local n
    n=$(wc -l < "${OUT_DIR}/telo_bed/${sample}.telo.bed")
    echo "[done] ${sample}: ${n} chromosome-end telomeres"

    rm -rf "$tmp"
}

export -f detect_one
export FASTA_DIR OUT_DIR THRESHOLD ENDS VGP_DRIVER EXCLUDE_CONTIG_REGEX

if [[ "$JOBS" -gt 1 ]]; then
    grep -v '^\s*$' "$SAMPLE_LIST" | xargs -P "$JOBS" -I{} bash -c 'detect_one "$@"' _ {}
else
    while read -r sample; do
        [[ -z "$sample" ]] && continue
        detect_one "$sample"
    done < "$SAMPLE_LIST"
fi

echo "Telomere detection complete: ${OUT_DIR}"
