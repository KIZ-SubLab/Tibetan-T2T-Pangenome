#!/usr/bin/env python3
"""
Per-haplotype, per-chromosome HOR array length and among-chromosome tests.

HOR length extraction (unchanged from the original pipeline)
------------------------------------------------------------
length(chr) = sum of (end - start) over the merged HOR intervals of that
              chromosome (merged with `bedtools merge -d 1000000000`),
              computed per haplotype.

Statistical analysis (submitted Methods)
----------------------------------------
HOR lengths are grouped BY CHROMOSOME (all haplotypes pooled within a
chromosome), then compared among chromosomes with:

  1. Kruskal-Wallis rank test   (scipy.stats.kruskal)
  2. Dunn's post-hoc pairwise test with a tie-corrected rank variance

Multiple-comparison correction
------------------------------
Benjamini-Hochberg is applied to the Dunn p-values in this implementation
and is reported EXPLICITLY as such. It could NOT be confirmed from the
historical production scripts which correction (if any) was used, so this
script does not claim that BH was the original choice. Raw p-values are
also reported so either can be quoted.

PROVENANCE / LIMITATION
-----------------------
This is an independent implementation of the statistics described in the
submitted Methods. It was NOT the program that produced the original
results: the historical production scripts (4.HOR_length/2.py and
4.1_all_APG_HOR_length/2.py) computed a per-chromosome Welch two-sample
t-test (scipy.stats.ttest_ind) and emitted 4.HOR_length/p_value.out.
The Kruskal-Wallis / Dunn statistics below have NOT been validated to
reproduce the published p-values.
"""
import argparse
import sys
from collections import defaultdict, Counter

from scipy.stats import kruskal, rankdata, norm


# ----------------------------------------------------------------------
# HOR length extraction (original logic)
# ----------------------------------------------------------------------
def hor_lengths(merge_bed):
    """Return {chrom: total HOR bp} for one merged HOR BED."""
    d = {}
    with open(merge_bed) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 3:
                continue
            chrom = p[0].split("#")[-1]
            d[chrom] = d.get(chrom, 0) + (int(p[2]) - int(p[1]))
    return d


# ----------------------------------------------------------------------
# Dunn's post-hoc test
# ----------------------------------------------------------------------
def dunn(groups, p_adjust="bh"):
    """Dunn's test for all group pairs.

    groups : {label: [values]}
    Returns rows: (g1, g2, n1, n2, mean_rank1, mean_rank2, z, p_raw, p_adj)

    Tie-corrected rank variance (Dunn 1964):
        A = N(N+1)/12 - sum(t^3 - t) / (12(N-1))
        z = (mean_rank_i - mean_rank_j) / sqrt(A * (1/n_i + 1/n_j))
    """
    from itertools import combinations

    pooled, labels = [], []
    for g, vals in groups.items():
        pooled.extend(vals)
        labels.extend([g] * len(vals))

    n_total = len(pooled)
    if n_total < 3 or len(groups) < 2:
        return []

    ranks = rankdata(pooled)
    tie_term = sum(c ** 3 - c for c in Counter(pooled).values())
    var = n_total * (n_total + 1) / 12.0 - tie_term / (12.0 * (n_total - 1))

    mean_rank, sizes = {}, {}
    for g, vals in groups.items():
        sizes[g] = len(vals)
    sums = defaultdict(float)
    for lab, r in zip(labels, ranks):
        sums[lab] += r
    for g in groups:
        mean_rank[g] = sums[g] / sizes[g]

    rows = []
    for a, b in combinations(sorted(groups), 2):
        se = (var * (1.0 / sizes[a] + 1.0 / sizes[b])) ** 0.5
        z = (mean_rank[a] - mean_rank[b]) / se if se > 0 else 0.0
        p = 2.0 * (1.0 - norm.cdf(abs(z)))
        rows.append([a, b, sizes[a], sizes[b], mean_rank[a], mean_rank[b], z, p])

    if p_adjust == "bh" and rows:
        order = sorted(range(len(rows)), key=lambda i: rows[i][7])
        m = len(rows)
        prev = 1.0
        for rank_pos in range(m - 1, -1, -1):
            i = order[rank_pos]
            prev = min(prev, rows[i][7] * m / (rank_pos + 1))
            rows[i].append(prev)
    else:
        for r in rows:
            r.append(r[7])

    return [tuple(r) for r in rows]


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="TSV: sample<TAB>population<TAB>merged_hor_bed")
    ap.add_argument("--out", required=True,
                    help="output TSV prefix; writes <out>.kw.tsv and <out>.dunn.tsv")
    ap.add_argument("--p-adjust", default="bh", choices=["bh", "none"],
                    help="multiple-comparison correction for Dunn (default bh)")
    a = ap.parse_args()

    # {chrom: [length across haplotypes]}
    per_chr = defaultdict(list)
    pops = defaultdict(list)
    with open(a.manifest) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            sample, pop, bed = p[0], p[1], p[2]
            for chrom, bp in hor_lengths(bed).items():
                per_chr[chrom].append(float(bp))
                pops[chrom].append(pop)

    # ---- Kruskal-Wallis among chromosomes ----
    with open(f"{a.out}.kw.tsv", "w") as out:
        out.write("n_chromosomes\tn_haplotype_observations\tH_statistic\tp_value\n")
        usable = {c: v for c, v in per_chr.items() if len(v) >= 2}
        if len(usable) >= 2:
            h, p = kruskal(*usable.values())
            total = sum(len(v) for v in usable.values())
            out.write(f"{len(usable)}\t{total}\t{h:.6e}\t{p:.6e}\n")
            print(f"Kruskal-Wallis: H={h:.6e} p={p:.6e} "
                  f"(k={len(usable)}, N={total})", file=sys.stderr)
        else:
            out.write("0\t0\tNA\tNA\n")

    # ---- Dunn's post-hoc ----
    rows = dunn(usable, p_adjust=a.p_adjust)
    with open(f"{a.out}.dunn.tsv", "w") as out:
        out.write("chrom1\tchrom2\tn1\tn2\tmean_rank1\tmean_rank2\t"
                  "z\tp_raw\tp_adjusted\n")
        for r in rows:
            out.write("\t".join(
                f"{x:.6e}" if isinstance(x, float) else str(x) for x in r) + "\n")

    print(f"wrote {a.out}.kw.tsv and {a.out}.dunn.tsv ({len(rows)} pairs)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
