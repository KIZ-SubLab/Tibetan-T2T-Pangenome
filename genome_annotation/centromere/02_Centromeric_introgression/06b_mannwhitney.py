#!/usr/bin/env python3
"""One-sided Mann-Whitney U test on per-window sequence identity."""
import argparse
import glob
import os

from scipy.stats import mannwhitneyu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--group-a", default="candidate")
    ap.add_argument("--group-b", default="control")
    a = ap.parse_args()

    groups = {a.group_a: [], a.group_b: []}
    for f in sorted(glob.glob(os.path.join(a.indir, "*.tsv"))):
        name = os.path.basename(f).replace(".tsv", "")
        grp = a.group_a if name.startswith(a.group_a) else a.group_b
        with open(f) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 7:
                    groups[grp].append(float(p[6]))

    with open(a.out, "w") as out:
        out.write("group_a\tgroup_b\tn_a\tn_b\tU\tp_one_sided\n")
        if not groups[a.group_a] or not groups[a.group_b]:
            out.write(f"{a.group_a}\t{a.group_b}\t0\t0\tNA\tNA\n")
            return
        u, p = mannwhitneyu(groups[a.group_a], groups[a.group_b],
                            alternative="greater")
        out.write(f"{a.group_a}\t{a.group_b}\t{len(groups[a.group_a])}\t"
                  f"{len(groups[a.group_b])}\t{u:.1f}\t{p:.3e}\n")


if __name__ == "__main__":
    main()
