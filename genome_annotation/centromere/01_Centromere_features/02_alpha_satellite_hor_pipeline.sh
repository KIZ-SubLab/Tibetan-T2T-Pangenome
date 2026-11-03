#!/usr/bin/env bash
# ============================================================
# Alpha-satellite monomer -> higher-order repeat (HOR) annotation.
#
# Steps
#   1. from the APGp1 cenanno table, drop non-alpha satellites
#      (BSat, GSat, HSat1A, HSat1B, HSat2, HSat3)
#   2. normalise the monomer id field  ("-" -> "/")
#   3. mon2stv.py            (APGp1 tool) -> STV / HOR structure
#   4. k-mer based HOR classification -> active HOR / inactive /
#      divergent monomers
#   5. emit merged HOR intervals (per haplotype, whole chromosome)
#
# External programs supplied by the APGp1 collaborators:
#   CENANNO_DIR/*.cenanno.bed   monomer annotation
#   MON2STV                     mon2stv.py
#   HOR_CLASSIFIER              k-mer HOR classifier
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

MON2STV=${MON2STV:-"$WORKDIR/bin/mon2stv.py"}
HOR_CLASSIFIER=${HOR_CLASSIFIER:-"$WORKDIR/bin/hor_kmer_classifier.py"}
EXCL='BSat|GSat|HSat1A|HSat1B|HSat2|HSat3'

OUT="$WORKDIR/hor"
mkdir -p "$OUT"/{input,out,out_bed6,HOR_merge}

for bed in "$CENANNO_DIR"/*.cenanno.bed; do
    sample=$(basename "$bed" .v0.9.full.cenanno.bed)

    # 1. keep alpha-satellite monomers only
    grep -Ev "$EXCL" "$bed" > "$OUT/input/${sample}.row.bed"

    # 2. normalise monomer id separators
    sed 's/-/\//g' "$OUT/input/${sample}.row.bed" > "$OUT/input/${sample}.input.bed"

    # 3. monomer -> STV
    python3 "$MON2STV" "$OUT/input/${sample}.input.bed" > "$OUT/out/${sample}.stv_raw.bed"

    # 4. k-mer based HOR classification
    python3 "$HOR_CLASSIFIER" "$OUT/out/${sample}.stv_raw.bed" \
        > "$OUT/out/${sample}.k_hor"

    # 5. merged HOR intervals (merge across the whole chromosome)
    bedtools merge -d 1000000000 -i "$OUT/out/${sample}.k_hor" \
        > "$OUT/HOR_merge/${sample}.HOR_merge.bed"
done

echo "HOR annotation written to $OUT"
