#!/usr/bin/env bash
# ============================================================
# hg38 -> Tibetan haplotype projection of SNV-based archaic segments.
#
# Pipeline:
#   bedtools getfasta  (hg38 sequence of each segment)
#   minimap2 -x asm5   (segment -> haplotype assembly)
#   awk filter         (same chromosome, MAPQ, query coverage)
#   best hit per query
#   bedtools merge
#
# One or more BEDs per haplotype (Neanderthal / Denisovan / unresolved)
# are projected independently and then unioned by 04_introgression_fraction.py.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/../config.env"

OUT="$WORKDIR/snv_projected"
TMP="$WORKDIR/tmp_snv"
mkdir -p "$OUT" "$TMP"

while read -r hap hap_fasta ancestral_bed; do
    [ -z "${hap:-}" ] && continue

    # 1. hg38 sequence for each segment
    tail -n +2 "$ancestral_bed" | cut -f1-3 | sort -k1,1 -k2,2n > "$TMP/${hap}.snv.bed"
    bedtools getfasta -fi "$HG38_FASTA" -bed "$TMP/${hap}.snv.bed" \
        -fo "$TMP/${hap}.snv.fa"

    # 2. map asm5
    "$MINIMAP2" -x asm5 -t "$THREADS" "$hap_fasta" "$TMP/${hap}.snv.fa" \
        > "$TMP/${hap}.paf"

    # 3. same-chromosome + MAPQ + query coverage; 4. best hit per query
    #    PAF fields: 1 qname 3 qstart 4 qend 6 tname 8 tstart 9 tend
    #                10 matches 11 blocklen 12 mapq
    awk -v mq="$SNV_MAPQ_MIN" -v cov="$SNV_COV_MIN" 'BEGIN{OFS="\t"} {
        if ($12 < mq) next
        if ($11 == 0 || $10 / $11 < cov) next
        split($1, a, ":")
        hg = a[1]; t = $6; sub(/^.*#/, "", t)
        if (hg != t) next
        print $6, $8, $9, $1, int(100*$10/$11), int(100*$10/($4-$3)), $12, $10
    }' "$TMP/${hap}.paf" | sort -k1,1 -k2,2n > "$TMP/${hap}.raw.bed"

    sort -k4,4 -k7,7nr -k6,6nr -k5,5nr -k8,8nr "$TMP/${hap}.raw.bed" \
        | awk '!seen[$4]++' | sort -k1,1 -k2,2n > "$TMP/${hap}.best.bed"

    # 5. merge overlapping / book-ended segments
    cut -f1-3 "$TMP/${hap}.best.bed" | bedtools merge -i - > "$OUT/${hap}.bed"
    echo "$hap: $(wc -l < "$OUT/${hap}.bed") merged segments"
done < "$SNV_SEGMENT_MANIFEST"   # hap<TAB>hap_fasta<TAB>ancestral_bed(hg38)
