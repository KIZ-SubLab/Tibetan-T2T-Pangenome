#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_telomere_summary.py
Telomere recovery statistics: count, length, coverage, total length per haploid.

Input : per-sample telomeric BED files produced by 01_telomere_detection.sh
Output: per-sample and pooled summary TSVs (no plotting)

Definitions
-----------
telomere        : one telomeric interval (one chromosome end) reported in the BED
telomere length : interval end - start (bp)
coverage        : detected telomeres / (n_haplotypes * EXPECTED_ENDS)
EXPECTED_ENDS   : theoretical chromosome ends per haploid
                  default 46 = 22 autosomes x2 + (X + Y) or (X x2)
                  i.e. every chromosome contributes two ends

Example
-------
python 02_telomere_summary.py \
    --bed-dir telo_bed --sample-list samples.txt \
    --expected-ends 46 --output telomere_summary.tsv
"""
import argparse
import os
import sys
from collections import defaultdict


def read_bed(path):
    """Yield (chrom, start, end) from a BED file."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                c = line.split()
            if len(c) < 3:
                continue
            try:
                yield c[0], int(c[1]), int(c[2])
            except ValueError:
                continue


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed-dir", required=True,
                    help="directory with <sample>.telo.bed files")
    ap.add_argument("--sample-list", required=True,
                    help="one sample id per line; matches <sample>.telo.bed")
    ap.add_argument("--expected-ends", type=int, default=46,
                    help="theoretical chromosome ends per haploid (default: 46)")
    ap.add_argument("--bed-suffix", default=".telo.bed",
                    help="BED file suffix (default: .telo.bed)")
    ap.add_argument("--per-sample-out", default=None,
                    help="per-sample TSV path (default: <bed-dir>/../telomere_per_sample.tsv)")
    ap.add_argument("--output", default=None,
                    help="pooled summary TSV path (default: <bed-dir>/../telomere_summary.tsv)")
    args = ap.parse_args()

    samples = [s.strip() for s in open(args.sample_list) if s.strip()]
    if not samples:
        sys.exit("ERROR: empty sample list")

    base = os.path.dirname(os.path.abspath(args.bed_dir))
    per_sample_out = args.per_sample_out or os.path.join(base, "telomere_per_sample.tsv")
    pooled_out = args.output or os.path.join(base, "telomere_summary.tsv")

    rows = []
    total_n = 0
    total_bp = 0
    per_chrom = defaultdict(int)

    for sample in samples:
        bed = os.path.join(args.bed_dir, f"{sample}{args.bed_suffix}")
        if not os.path.exists(bed):
            print(f"[warn] missing BED for {sample}: {bed}", file=sys.stderr)
            rows.append((sample, 0, 0, 0.0, 0.0))
            continue
        n = 0
        bp = 0
        for chrom, start, end in read_bed(bed):
            n += 1
            bp += max(0, end - start)
            per_chrom[chrom] += 1
        frac = n / args.expected_ends if args.expected_ends else 0.0
        rows.append((sample, n, bp, frac, bp / 1000.0))
        total_n += n
        total_bp += bp

    with open(per_sample_out, "w") as out:
        out.write("sample\tn_telomeres\ttotal_telomere_bp\tfraction_of_expected_ends\ttotal_telomere_kbp\n")
        for s, n, bp, frac, kbp in rows:
            out.write(f"{s}\t{n}\t{bp}\t{frac:.4f}\t{kbp:.2f}\n")

    n_hap = len(samples)
    with open(pooled_out, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"n_haplotypes\t{n_hap}\n")
        out.write(f"expected_ends_per_haploid\t{args.expected_ends}\n")
        out.write(f"n_telomeres_recovered\t{total_n}\n")
        out.write(f"expected_telomeres\t{n_hap * args.expected_ends}\n")
        out.write(f"recovery_percent\t{total_n / (n_hap * args.expected_ends) * 100:.2f}\n"
                  if n_hap * args.expected_ends else "recovery_percent\t0\n")
        out.write(f"mean_telomeres_per_haploid\t{total_n / n_hap:.2f}\n")
        out.write(f"mean_total_telomere_kbp_per_haploid\t{total_bp / n_hap / 1000:.2f}\n")
        out.write(f"mean_telomere_length_bp\t{total_bp / total_n:.1f}\n" if total_n else "")
        for chrom in sorted(per_chrom):
            out.write(f"n_telomeres_{chrom}\t{per_chrom[chrom]}\n")

    print(f"[done] per-sample -> {per_sample_out}")
    print(f"[done] pooled     -> {pooled_out}")
    if n_hap * args.expected_ends:
        print(f"  recovered {total_n} telomeres = "
              f"{total_n / (n_hap * args.expected_ends) * 100:.2f}% of "
              f"{n_hap * args.expected_ends} expected ends")
    print(f"  mean total telomere length {total_bp / n_hap / 1000:.1f} kbp per haploid")


if __name__ == "__main__":
    main()
