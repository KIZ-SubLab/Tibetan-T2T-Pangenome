#!/usr/bin/env bash
# ============================================================
# ARCkmerFinder scan of the 16 Tibetan haplotypes.
#
# ARCkmerFinder writes 2-kbp non-overlapping windows across the
# assembly; for each window the value in column 4 is the number
# of k-mer occurrences in that window that are present in the
# archaic-specific database (presence lookup, not DB-count weighted).
#
# Output of the tool:
#   <assembly>_in_2000_windows.bed
#   windows_across_genome_with_zero_and_nonzero_matching_kmers.bed
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

SCAN="$WORKDIR/arcscan"
mkdir -p "$SCAN"

# ARCkmerFinder config template is supplied with the tool; only the
# two fields that matter are overridden here.
for asm in "$ASSEMBLY_DIR"/*.fasta; do
    sample=$(basename "$asm" .fasta)
    run="$SCAN/$sample"
    mkdir -p "$run"
    cp "$ARCKMERFINDER_SNAKEFILE" "$run/Snakefile"
    cat > "$run/config.yaml" << CONF
kmer_counting:
  assembly: $asm
  existing_meryl_database_path: $WORKDIR/kmerdb/archaic_all_specific.meryl
  window_size: $ARC_WINDOW
  remove_intermediate_files: True
CONF
    ( cd "$run" && snakemake --cores "$THREADS" --latency-wait 60 )
done

# ---- significant windows (uniform threshold) ----
mkdir -p "$SCAN/sig_windows"
for bed in "$SCAN"/*/windows_across_genome_with_zero_and_nonzero_matching_kmers.bed; do
    sample=$(basename "$(dirname "$bed")")
    # keep contig as-is (TBN native coordinates)
    awk -F'\t' -v t="$ARC_SIGNAL_MIN" 'BEGIN{OFS="\t"} $4>=t' "$bed" \
        | sort -k1,1 -k2,2n > "$SCAN/sig_windows/${sample}.sig.bed"
done

# ---- intersect with haplotype-specific HOR intervals (>=1 bp) ----
mkdir -p "$SCAN/cen_sig_windows"
while read -r sample hor_bed; do
    [ -z "${sample:-}" ] && continue
    bedtools intersect \
        -a "$SCAN/sig_windows/${sample}.sig.bed" \
        -b "$hor_bed" \
        -wa -u \
        > "$SCAN/cen_sig_windows/${sample}.cen_sig.bed"
done < "$HOR_MANIFEST"   # sample<TAB>path/to/<sample>.HOR_merge.bed

echo "significant centromeric windows: $SCAN/cen_sig_windows"
