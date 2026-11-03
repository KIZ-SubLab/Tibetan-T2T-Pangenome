#!/usr/bin/env python3
"""
Satellite composition of centromeric / pericentromeric intervals.

For every haplotype:
  * intersect the satellite-class BED with the centromere intervals
  * merge overlapping intervals per class (no double counting)
  * report bp per class, total annotated satellite bp, and the
    proportion of the centromere interval covered

The denominator is the centromere interval length actually used by
the project (CENTROMERE_BED), NOT the whole assembly.
"""
import argparse
import subprocess
import tempfile
import os
from collections import defaultdict


def merge_by_class(sat_bed, cen_bed, chrom, start, end):
    """Return {class: merged_bp} inside [start,end)."""
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
        tmp.write(f"{chrom}\t{start}\t{end}\n")
        region = tmp.name
    try:
        out = subprocess.run(
            ["bedtools", "intersect", "-a", sat_bed, "-b", region, "-wa", "-u"],
            capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(region)

    per_class = defaultdict(list)
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 4:
            continue
        per_class[p[3]].append((int(p[1]), int(p[2])))

    merged_bp = {}
    for cls, iv in per_class.items():
        iv.sort()
        merged = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        merged_bp[cls] = sum(e - s for s, e in merged)
    return merged_bp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--satellite-bed", required=True)
    ap.add_argument("--centromere-bed", required=True,
                    help="BED4: contig start end sample")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    with open(a.centromere_bed) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            contig, start, end, sample = p[0], int(p[1]), int(p[2]), p[3]
            bp = merge_by_class(a.satellite_bed, a.centromere_bed,
                                contig, start, end)
            total = sum(bp.values())
            span = end - start
            for cls, v in sorted(bp.items()):
                rows.append((sample, contig.split("#")[-1], start, end,
                             cls, v, span, v / span if span else 0.0))
            rows.append((sample, contig.split("#")[-1], start, end,
                         "TOTAL", total, span, total / span if span else 0.0))

    with open(a.out, "w") as out:
        out.write("sample\tchrom\tstart\tend\tsatellite_class\tbp\t"
                  "cen_span_bp\tfraction_of_cen_span\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
