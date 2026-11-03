#!/usr/bin/env bash
# ============================================================
# Centromeric sequence comparison (20-kb windows, 1-kb step).
#
# For every window of the query haplotype, the window sequence is
# mapped to the target haplotype with minimap2 asm5; only the
# highest-identity alignment is kept, then the identity
# distributions of the two groups are compared with a one-sided
# Mann-Whitney U test (implemented in 06b).
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

OUT="$WORKDIR/seqcomp"
mkdir -p "$OUT"/{windows,paf,identity}

for asm in "$ASSEMBLY_DIR"/*.fasta; do
    sample=$(basename "$asm" .fasta)
    # 1. 20-kb windows, 1-kb step
    bedtools makewindows -g "${asm}.fai" -w "$SEQ_WINDOW" -s "$SEQ_STEP" \
        > "$OUT/windows/${sample}.bed"
done

# 2. map query windows to the target haplotype
while read -r query query_fa target target_fa; do
    [ -z "${query:-}" ] && continue
    bedtools getfasta -fi "$query_fa" -bed "$OUT/windows/${query}.bed" \
        -fo "$OUT/windows/${query}.fa"
    "$MINIMAP2" -x asm5 --eqx --secondary=no -t "$THREADS" \
        "$target_fa" "$OUT/windows/${query}.fa" \
        > "$OUT/paf/${query}_vs_${target}.paf"

    # 3. keep highest-identity alignment per query window
    awk -v q="$SEQ_MAPQ_MIN" -v L="$SEQ_MIN_ALIGN_LEN" -v id="$SEQ_MIN_IDENTITY" \
        'BEGIN{OFS="\t"}
         $12 >= q && $11 >= L && ($10/$11*100) >= id {
             print $1, $6, $8, $9, $10, $11, ($10/$11*100), $12
         }' "$OUT/paf/${query}_vs_${target}.paf" \
        | sort -k1,1 -k7,7nr | awk '!seen[$1]++' \
        > "$OUT/identity/${query}_vs_${target}.tsv"
done < "$SEQ_COMPARISON_MANIFEST"   # query<TAB>query_fa<TAB>target<TAB>target_fa

python3 "$(dirname "$0")/06b_mannwhitney.py" \
    --indir "$OUT/identity" --out "$OUT/mannwhitney.tsv"
