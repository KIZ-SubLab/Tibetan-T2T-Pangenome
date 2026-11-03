#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_TVR_population_comparison.py
Chromosome-arm TVR repertoire similarity and between-group composition test.

Two analyses, both driven by CLI arguments (no group-specific code):

1. Within-group chromosome-arm repertoire similarity
       For each arm, a haplotype-pooled set of unique TVR hexamers is built.
       J(A, B) = |A n B| / |A u B| for every pair of arms (48 x 48 when all
       arms are present; p-p, q-q and p-q comparisons are all included).

2. Between-group composition test
       For every chromosome arm x TVR hexamer:
           a = number of group-A haplotypes carrying that hexamer (>=1 copy)
           b = number of group-B haplotypes carrying that hexamer
           two-sided Fisher exact test on
               [[a, A_total - a], [b, B_total - b]]
           where A_total / B_total are the summed haplotype-presence counts of
           that arm. The statistical unit is therefore a haplotype's
           presence/absence of a hexamer, NOT the number of copies.
       Benjamini-Hochberg FDR is applied across all tests.

Only class == 'tvr' 6-mers are used; canonical / rotated / ambiguous slots are
excluded upstream by 05_TVR_calling.py.

Example
-------
python 07_TVR_population_comparison.py \
    --group-a-dir hex_A --group-a-samples A.txt --group-a-label A \
    --group-b-dir hex_B --group-b-samples B.txt --group-b-label B \
    --output-dir comparison/
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import fisher_exact


def parse_arm(chrom_field, p_arm_limit=20000):
    import re
    if "#" in chrom_field:
        chrom_field = chrom_field.split("#")[-1]
    m = re.match(r"(chr(?:\d+|X|Y))[_:]", chrom_field, re.I)
    if not m:
        return None
    chrom = m.group(1)
    mi = re.search(r"_(\d+)-(\d+)", chrom_field)
    if not mi:
        return None
    return f"{chrom}_{'p' if int(mi.group(1)) < p_arm_limit else 'q'}"


def arm_sort_key(arm):
    import re
    m = re.match(r"chr(\d+)_([pq])", arm, re.I)
    if m:
        return (int(m.group(1)), 0 if m.group(2).lower() == "p" else 1)
    m = re.match(r"chr([XY])_([pq])", arm, re.I)
    if m:
        return (23 if m.group(1).upper() == "X" else 24,
                0 if m.group(2).lower() == "p" else 1)
    return (99, 0)


def load_presence(hex_dir, samples, p_arm_limit=20000):
    """-> {arm: {hexamer: set(haplotype)}}"""
    data = defaultdict(lambda: defaultdict(set))
    used = []
    for sample in samples:
        path = os.path.join(hex_dir, f"{sample}_hexamers.tsv")
        if not os.path.exists(path):
            print(f"[warn] missing {path}", file=sys.stderr)
            continue
        used.append(sample)
        with open(path) as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 7 or c[6] != "tvr":
                    continue
                arm = parse_arm(c[2], p_arm_limit)
                if arm:
                    data[arm][c[5].upper()].add(sample)
    return data, used


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum(1.0, ranked * n / np.arange(1, n + 1))
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group-a-dir", required=True)
    ap.add_argument("--group-a-samples", required=True)
    ap.add_argument("--group-a-label", default="A")
    ap.add_argument("--group-b-dir", required=True)
    ap.add_argument("--group-b-samples", required=True)
    ap.add_argument("--group-b-label", default="B")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--fdr-threshold", type=float, default=0.05)
    ap.add_argument("--min-presence", type=int, default=5,
                    help="require a+b >= this to test a hexamer (default: 5)")
    ap.add_argument("--p-arm-limit", type=int, default=20000)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    a_samples = [s.strip() for s in open(args.group_a_samples) if s.strip()]
    b_samples = [s.strip() for s in open(args.group_b_samples) if s.strip()]

    A, used_a = load_presence(args.group_a_dir, a_samples, args.p_arm_limit)
    B, used_b = load_presence(args.group_b_dir, b_samples, args.p_arm_limit)
    print(f"{args.group_a_label}: {len(used_a)} haplotypes, {len(A)} arms")
    print(f"{args.group_b_label}: {len(used_b)} haplotypes, {len(B)} arms")

    # ---------------- 1. within-group arm-arm Jaccard ----------------
    for label, data, used in ((args.group_a_label, A, used_a),
                              (args.group_b_label, B, used_b)):
        arms = sorted(data.keys(), key=arm_sort_key)
        out = os.path.join(args.output_dir, f"arm_similarity_{label}.tsv")
        with open(out, "w") as fh:
            fh.write("arm_a\tarm_b\tjaccard\tn_shared\tn_a_only\tn_b_only\n")
            for i, a1 in enumerate(arms):
                s1 = set(data[a1])
                for a2 in arms[i:]:
                    s2 = set(data[a2])
                    inter, union = s1 & s2, s1 | s2
                    jac = len(inter) / len(union) if union else 0.0
                    fh.write(f"{a1}\t{a2}\t{jac:.4f}\t{len(inter)}\t"
                             f"{len(s1 - s2)}\t{len(s2 - s1)}\n")
        print(f"[done] arm-arm Jaccard ({label}, {len(arms)} arms) -> {out}")

    # ---------------- 2. between-group composition test ----------------
    common_arms = sorted(set(A) & set(B), key=arm_sort_key)
    tests = []
    for arm in common_arms:
        a_hex, b_hex = A[arm], B[arm]
        a_total = sum(len(v) for v in a_hex.values())   # haplotype-presence total
        b_total = sum(len(v) for v in b_hex.values())
        for h in set(a_hex) | set(b_hex):
            a_cnt = len(a_hex.get(h, ()))
            b_cnt = len(b_hex.get(h, ()))
            if a_cnt + b_cnt < args.min_presence:
                continue
            try:
                odds, p = fisher_exact([[a_cnt, a_total - a_cnt],
                                        [b_cnt, b_total - b_cnt]])
            except Exception:
                continue
            tests.append({
                "arm": arm, "hexamer": h,
                "a_count": a_cnt, "b_count": b_cnt,
                "a_other": a_total - a_cnt, "b_other": b_total - b_cnt,
                "odds": odds, "p": p,
                "a_freq": a_cnt / a_total if a_total else 0.0,
                "b_freq": b_cnt / b_total if b_total else 0.0,
            })

    fdr = bh_fdr([t["p"] for t in tests]) if tests else np.array([])
    out = os.path.join(args.output_dir, "composition_fisher.tsv")
    with open(out, "w") as fh:
        fh.write(f"arm\thexamer\t{args.group_a_label}_count\t{args.group_b_label}_count\t"
                 f"{args.group_a_label}_other\t{args.group_b_label}_other\t"
                 f"{args.group_a_label}_freq\t{args.group_b_label}_freq\t"
                 "odds_ratio\tp_value\tFDR\tsignificant\n")
        for t, f in zip(tests, fdr):
            sig = "yes" if f < args.fdr_threshold else "no"
            fh.write(f"{t['arm']}\t{t['hexamer']}\t{t['a_count']}\t{t['b_count']}\t"
                     f"{t['a_other']}\t{t['b_other']}\t{t['a_freq']:.4f}\t"
                     f"{t['b_freq']:.4f}\t{t['odds']:.4g}\t{t['p']:.6e}\t"
                     f"{f:.6e}\t{sig}\n")

    n_sig = int(sum(1 for f in fdr if f < args.fdr_threshold))
    print(f"[done] composition test: {len(tests)} Fisher tests, "
          f"{n_sig} significant at FDR < {args.fdr_threshold} -> {out}")
    if tests:
        print(f"  min p = {min(t['p'] for t in tests):.3e}, min FDR = {min(fdr):.3e}")
        sig_arms = sorted({t["arm"] for t, f in zip(tests, fdr)
                           if f < args.fdr_threshold}, key=arm_sort_key)
        if sig_arms:
            print(f"  arms with significant hexamers: {sig_arms}")


if __name__ == "__main__":
    main()
