#!/usr/bin/env bash
# ============================================================
# Full-flank orthologue -> 50-kb window phylogeny inputs.
#
# Four distinct steps (do NOT conflate them):
#   A. genome mapping      minimap2 -x asm5 (same-chromosome, PAF)
#   B. full-flank MSA      MAFFT over all 481 taxa
#   C. 50-kb projection    cut sequence from each taxon's own
#                          orthologue using the ruler-taxon row of the
#                          full-flank MSA as the coordinate ruler
#   D. per-window MSA      a NEW MAFFT run for every 50-kb window
#
# Taxon set: 480 human haplotypes + chimpanzee (gorilla excluded).
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

PHY="$WORKDIR/phylo"
mkdir -p "$PHY"/{paf,ortholog,fullflank_msa,windows}

# ---- A. orthologue mapping (one query per taxon) ----
while read -r taxon taxon_fa; do
    [ -z "${taxon:-}" ] && continue
    "$MINIMAP2" -x asm5 -t "$THREADS" "$REFERENCE_FASTA" "$taxon_fa" \
        > "$PHY/paf/${taxon}.paf"

    # same-chromosome + highest-score hit; PAF strand is kept so the
    # extracted sequence can be oriented to the forward strand
    awk 'BEGIN{OFS="\t"} {
        split($6, t, "#"); if (t[length(t)] != CHROM) next
        print $6, $8, $9, $5, $10, $11, $12
    }' CHROM="$TARGET_CHROM" "$PHY/paf/${taxon}.paf" \
        | sort -k5,5nr | head -1 > "$PHY/ortholog/${taxon}.hit"
done < "$TAXON_MANIFEST"   # taxon<TAB>assembly_fasta

# ---- B. full-flank MSA (481 taxa) ----
"$MAFFT" --thread "$THREADS" --retree 2 --maxiterate 0 --inputorder \
    "$PHY/ortholog/fullflank_unaligned.fa" \
    > "$PHY/fullflank_msa/fullflank.mafft.fa"

# ---- C + D. 50-kb windows ----
# WINDOW_MANIFEST: window_id<TAB>side<TAB>start_in_fullflank_msa<TAB>end
while IFS=$'\t' read -r wid side start end; do
    [ -z "${wid:-}" ] && continue
    outdir="$PHY/windows/${side}/${wid}"
    mkdir -p "$outdir"

    # C. project the column block back onto each taxon's own sequence
    python3 "$(dirname "$0")/07b_project_windows.py" \
        --msa "$PHY/fullflank_msa/fullflank.mafft.fa" \
        --ruler-taxon "$RULER_TAXON" \
        --start "$start" --end "$end" \
        --output "$outdir/${wid}.unaligned.fa"

    # D. per-window realignment
    "$MAFFT" --thread "$MAFFT_THREADS" --retree 2 --maxiterate 0 --inputorder \
        "$outdir/${wid}.unaligned.fa" > "$outdir/${wid}.mafft.fa"
done < "$WINDOW_MANIFEST"
