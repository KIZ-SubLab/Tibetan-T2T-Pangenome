#!/usr/bin/env bash
# ============================================================
# African background k-mer database (661 unrelated genomes).
#
# Per-sample meryl count from CRAM/FASTQ (see AFR_SAMPLE_LIST),
# then union-sum into one joint database, then greater-than 9
# (== total count >= 10) to suppress rare/erroneous k-mers.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

DB="$WORKDIR/kmerdb"
mkdir -p "$DB/afr"

# per-sample counting
while read -r sample r1 r2; do
    [ -z "${sample:-}" ] && continue
    "$MERYL" count k="$KMER_SIZE" threads="$THREADS" memory=4G \
        output "$DB/afr/${sample}.meryl" "$r1" ${r2:+"$r2"}
done < "$AFR_SAMPLE_LIST"

# joint database
shopt -s nullglob
afr=("$DB/afr/"*.meryl)
shopt -u nullglob
"$MERYL" union-sum output "$DB/AFR.merge.meryl" "${afr[@]}"

# greater-than 9 == total count >= 10
"$MERYL" greater-than "$AFR_MIN_COUNT" \
    "$DB/AFR.merge.meryl" \
    output "$DB/AFR.merge.filter.meryl"

"$MERYL" statistics "$DB/AFR.merge.filter.meryl"
