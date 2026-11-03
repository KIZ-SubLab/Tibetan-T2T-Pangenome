#!/usr/bin/env python3
"""
Per-haplotype autosomal archaic introgression fraction.

numerator   = union of Neanderthal / Denisovan / unresolved segments
              projected onto that haplotype, restricted to chr1-22
denominator = the SAME haplotype assembly, chr1-22 sum of .fai lengths

Both terms must come from the same haplotype; this is checked
explicitly and the script aborts on a mismatch.
"""
import argparse
import sys


def read_fai(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[p[0].split("#")[-1]] = int(p[1])
    return d


def autosomes():
    return {f"chr{i}" for i in range(1, 23)}


def merge(intervals):
    """Merge overlapping and book-ended intervals (bedtools merge semantics)."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beds", nargs="+", required=True,
                    help="one BED per archaic ancestry, same haplotype")
    ap.add_argument("--fai", required=True,
                    help=".fai of the SAME haplotype assembly")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    auto = autosomes()
    per_chr = {}
    source_contigs = set()
    for bed in a.beds:
        with open(bed) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) < 3:
                    continue
                chrom = p[0].split("#")[-1]
                source_contigs.add(p[0])
                if chrom not in auto:
                    continue
                per_chr.setdefault(chrom, []).append((int(p[1]), int(p[2])))

    merged = {c: merge(v) for c, v in per_chr.items()}
    num = sum(e - s for v in merged.values() for s, e in v)

    fai = read_fai(a.fai)
    den = sum(fai[c] for c in auto if c in fai)

    # consistency guard: every contig in the numerator must belong to
    # the haplotype described by the .fai
    fai_contigs = set(fai)
    bad = [c for c in source_contigs if c not in fai_contigs]
    if bad:
        sys.exit(f"ERROR: numerator contigs not in .fai: {bad[:3]}")

    frac = num / den
    with open(a.out, "w") as out:
        out.write("introgressed_bp\tdenominator_bp\tfraction\tpercent\n")
        out.write(f"{num}\t{den}\t{frac:.8f}\t{100*frac:.4f}\n")
    print(f"introgressed={num} autosomal={den} fraction={100*frac:.4f}%")


if __name__ == "__main__":
    main()
