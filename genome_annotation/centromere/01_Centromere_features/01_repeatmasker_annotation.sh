#!/usr/bin/env bash
# ============================================================
# Repeat annotation and satellite-class extraction.
#
# Runs RepeatMasker on haplotype assemblies and derives the
# satellite classes used downstream (alpha, beta, gamma, HSat1-3).
#
# NOTE: the alpha-satellite monomer / HOR annotation itself is
#       provided by the APGp1 collaborator (CENANNO_DIR); this
#       script only produces the RepeatMasker layer.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

mkdir -p "$WORKDIR/repeatmasker"

# --- RepeatMasker (v4.1.6) ---
# -species human, soft-masked output retained for downstream parsing
for asm in "$ASSEMBLY_DIR"/*.fasta; do
    sample=$(basename "$asm" .fasta)
    "$REPEATMASKER" \
        -species human \
        -pa "$THREADS" \
        -xsmall \
        -dir "$WORKDIR/repeatmasker" \
        "$asm"
done

# --- Extract satellite classes from the .out table ---
# Classes retained: (alpha|beta|gamma)-satellite and HSat1/2/3.
# Class assignment follows the RepeatMasker 'class/family' field.
python3 "$(dirname "$0")/01b_parse_satellite_classes.py" \
    --outdir "$WORKDIR/repeatmasker" \
    --output "$WORKDIR/repeatmasker/satellite_classes.bed"
