#!/usr/bin/env bash
# ============================================================
# Archaic-enriched canonical 15-mer databases.
#
# For each archaic source:
#   meryl count (k=15, canonical default)
#   meryl union-sum across libraries of the same source
#   meryl greater-than 2      -> retain total count >= 3
#   meryl difference  (archaic - AFR joint)  -> archaic-specific
#
# Finally the five archaic-specific databases are merged with
# meryl union-sum (counts are SUMMED across sources).
#
# Input reads are the public SRA/ENA runs listed in $ARCHAIC_MANIFEST
# (one line per run: source<TAB>run<TAB>r1[<TAB>r2]).
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

DB="$WORKDIR/kmerdb"
mkdir -p "$DB"/{raw,source,specific}

# ---- 1. per-run counting (canonical 15-mers) ----
while IFS=$'\t' read -r source run r1 r2; do
    [ -z "${source:-}" ] && continue
    "$MERYL" count k="$KMER_SIZE" threads="$THREADS" memory=32G \
        output "$DB/raw/${source}_${run}.meryl" \
        "$r1" ${r2:+"$r2"}
done < "$ARCHAIC_MANIFEST"

# ---- 2. union-sum per source, then low-frequency filter ----
for source in $ARCHAIC_DATASETS; do
    shopt -s nullglob
    libs=("$DB/raw/${source}_"*.meryl)
    shopt -u nullglob
    [ ${#libs[@]} -gt 0 ] || { echo "no libraries for $source" >&2; continue; }

    "$MERYL" union-sum output "$DB/source/${source}.meryl" "${libs[@]}"

    # greater-than 2  ==  total count >= 3
    "$MERYL" greater-than "$ARCHAIC_MIN_COUNT" \
        "$DB/source/${source}.meryl" \
        output "$DB/source/${source}.filtered.meryl"

    # ---- 3. subtract the African joint database ----
    "$MERYL" difference \
        "$DB/source/${source}.filtered.meryl" \
        "$AFR_DB" \
        output "$DB/specific/${source}_specific.meryl"
done

# ---- 4. merged archaic-specific database (union-sum) ----
shopt -s nullglob
spec=("$DB/specific/"*_specific.meryl)
shopt -u nullglob
"$MERYL" union-sum output "$DB/archaic_all_specific.meryl" "${spec[@]}"
"$MERYL" statistics "$DB/archaic_all_specific.meryl"
