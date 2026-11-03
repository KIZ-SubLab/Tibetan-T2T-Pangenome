#!/usr/bin/env python3
"""Validate DISC1 A+B+C statistics and official B/C window reproduction."""

import argparse
import math
import os

import numpy as np
import pandas as pd


COMPS = {
    "A_relative": ("den", "nea", "tbn", "yri"),
    "B_Denisovan": ("yri", "tbn", "den", "anc"),
    "C_Neanderthal": ("yri", "tbn", "nea", "anc"),
}


def check(value, message):
    if not value:
        raise AssertionError(message)


def equal_or_both_nan(left: float, right: float, atol: float = 1e-12) -> bool:
    return (
        math.isclose(left, right, abs_tol=atol)
        if math.isfinite(left) and math.isfinite(right)
        else math.isnan(left) and math.isnan(right)
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", required=True)
    p.add_argument("--window-snps", type=int, required=True)
    p.add_argument("--step-snps", type=int, required=True)
    p.add_argument("--official-b")
    p.add_argument("--official-c")
    p.add_argument("--region-only", action="store_true")
    args = p.parse_args()
    sites = pd.read_csv(os.path.join(args.results, "sites_common_ancestral_polarized.tsv.gz"), sep="\t")
    summary = pd.read_csv(os.path.join(args.results, "region_summary_ABC.tsv"), sep="\t")
    windows = None if args.region_only else pd.read_csv(os.path.join(
        args.results, f"windows_{args.window_snps}SNP_step{args.step_snps}.tsv"), sep="\t")
    check(sites.pos.is_monotonic_increasing, "Sites are not sorted")
    check(not sites.pos.duplicated().any(), "Duplicate site positions")
    for group in ["tbn", "yri", "den", "nea"]:
        check(sites[f"p_{group}"].between(0, 1).all(), f"ALT frequency outside [0,1]: {group}")
        check(sites[f"derived_{group}"].between(0, 1).all(),
              f"Derived frequency outside [0,1]: {group}")
    for code, groups in COMPS.items():
        p1, p2, p3, p4 = [(0.0 if g == "anc" else sites[f"derived_{g}"]) for g in groups]
        algebra = (p2-p1)*(p3-p4)
        check(np.allclose(algebra, sites[f"{code}_num"], atol=1e-12),
              f"Numerator algebra failed: {code}")
        numerator = float(sites[f"{code}_num"].sum())
        d_denominator = float(sites[f"{code}_d_den"].sum())
        fd_denominator = float(sites[f"{code}_fd_den"].sum())
        d = numerator / d_denominator if d_denominator != 0 else math.nan
        fd = numerator / fd_denominator if fd_denominator != 0 else math.nan
        row = summary.loc[summary.comparison.eq(code)].iloc[0]
        check(equal_or_both_nan(d, row.D_estimate), f"Region D failed: {code}")
        check(equal_or_both_nan(fd, row.fd_estimate), f"Region fd failed: {code}")

    finite_region_p = 0
    for _, row in summary.iterrows():
        for stat in ["D", "fd"]:
            pvalue = row[f"{stat}_p_two_sided"]
            if math.isfinite(pvalue):
                finite_region_p += 1
                check(0 <= pvalue <= 1,
                      f"Region P outside [0,1]: {row.comparison} {stat}")
                check(row[f"{stat}_n_blocks"] >= 3,
                      f"Fewer than three region blocks: {row.comparison} {stat}")
            else:
                # An undefined ratio, fewer than three valid delete blocks, or
                # a zero jackknife SE makes a Z/P value mathematically undefined.
                check(not math.isfinite(row[f"{stat}_z"]),
                      f"Finite Z paired with missing P: {row.comparison} {stat}")
    if args.region_only:
        check(finite_region_p >= 4,
              "Region-only case has fewer than four estimable A/B/C D/fd P values")

    if not args.region_only:
        for code in COMPS:
            for suffix in ["D_p_two_sided", "D_q_BH", "fd_p_two_sided", "fd_q_BH_exploratory"]:
                col = f"{code}_{suffix}"
                check(windows[col].dropna().between(0, 1).all(), f"P/q outside [0,1]: {col}")
        check(windows.n_snps.eq(args.window_snps).all(), "Incorrect SNP count in a window")
        check(windows.first_site_index.diff().dropna().eq(args.step_snps).all(), "Incorrect window step")

        check(args.official_b and args.official_c,
              "Official B/C files are required outside region-only mode")
        for code, official_path in [("B_Denisovan", args.official_b),
                                    ("C_Neanderthal", args.official_c)]:
            official = pd.read_csv(official_path)
            check(len(official) == len(windows), f"Official window count differs: {code}")
            check(np.array_equal(official.start, windows.start), f"Official starts differ: {code}")
            check(np.array_equal(official.end, windows.end), f"Official ends differ: {code}")
            # genomics_general can emit NA when a window has a zero denominator.
            check(np.allclose(official.D, windows[f"{code}_D_estimate"],
                              atol=5.1e-5, equal_nan=True),
                  f"Official D differs: {code}")
            check(np.allclose(official.fd, windows[f"{code}_fd_estimate"],
                              atol=5.1e-5, equal_nan=True),
                  f"Official fd differs: {code}")

    with open(os.path.join(args.results, "validation.txt"), "w") as out:
        out.write("PASS: common sites sorted and unique\n")
        out.write("PASS: ALT and derived frequencies in [0,1]\n")
        out.write("PASS: A+B+C ABBA-BABA numerator algebra\n")
        out.write("PASS: A+B+C whole-region D/fd exact reproduction\n")
        out.write("PASS: region-level jackknife P values valid whenever mathematically estimable\n")
        if args.region_only:
            out.write("PASS: region-only mode intentionally omitted SNP windows\n")
            out.write("PASS: region-only result contains all three A/B/C comparisons\n")
            out.write(f"mode=REGION_ONLY sites={len(sites)} windows=0\n")
        else:
            out.write("PASS: window SNP count, step and P/q ranges\n")
            out.write("PASS: B/C D and fd reproduce official genomics_general windows\n")
            out.write(f"mode=WINDOWED sites={len(sites)} windows={len(windows)}\n")


if __name__ == "__main__":
    main()
