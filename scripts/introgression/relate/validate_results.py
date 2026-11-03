#!/usr/bin/env python3
"""Structural and numerical validation for one RELATE/CLUES run."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True)
    p.add_argument("--windows", required=True, nargs="+")
    p.add_argument("--expected-df", required=True, type=int)
    p.add_argument("--expected-haplotypes", required=True, type=int)
    args = p.parse_args()

    root = Path(args.output_dir)
    results = root / "results"
    required = [
        results / "focal_summary.tsv",
        results / "tree" / "focal_tree.png",
        results / "tree" / "focal_tree.svg",
        results / "tree" / "focal_tree.pdf",
        results / "branch_samples" / "focal.newick",
        results / "branch_samples" / "focal.mut",
        results / "focal_times.txt",
        results / "derived_haplotypes.txt",
        results / "local_genealogy" / "region.anc",
        results / "local_genealogy" / "region.mut",
    ]
    rows = []

    def record(check, status, detail):
        rows.append((check, status, str(detail)))

    for path in required:
        ok = path.is_file() and path.stat().st_size > 0
        record(f"file:{path.relative_to(root)}", "PASS" if ok else "FAIL", path.stat().st_size if path.is_file() else "missing")

    derived = results / "derived_haplotypes.txt"
    if derived.is_file():
        vals = np.loadtxt(derived, ndmin=1)
        ok = len(vals) == args.expected_haplotypes and set(np.unique(vals)).issubset({0, 1})
        record("derived_haplotype_vector", "PASS" if ok else "FAIL", f"n={len(vals)} sum={vals.sum():.0f}")

    times = results / "focal_times.txt"
    if times.is_file():
        nlines = sum(1 for line in times.open(encoding="utf-8") if line.strip())
        record("clues_times_even_lines", "PASS" if nlines > 0 and nlines % 2 == 0 else "FAIL", nlines)

    for label in args.windows:
        prefix = results / f"clues2_{label}" / "focal"
        inference = Path(str(prefix) + "_inference.txt")
        ci = Path(str(prefix) + "_CI.txt")
        freqs = Path(str(prefix) + "_freqs.txt")
        post = Path(str(prefix) + "_post.txt")
        figure_prefix = results / "figures" / f"focal_CLUES2_{label}"
        for path in [inference, ci, freqs, post]:
            ok = path.is_file() and path.stat().st_size > 0
            record(f"file:{path.relative_to(root)}", "PASS" if ok else "FAIL", path.stat().st_size if path.is_file() else "missing")
        for ext in ["png", "svg", "pdf", "tiff", "source_data.tsv", "display_smooth.tsv"]:
            path = Path(str(figure_prefix) + "." + ext)
            ok = path.is_file() and path.stat().st_size > 0
            record(f"file:{path.relative_to(root)}", "PASS" if ok else "FAIL", path.stat().st_size if path.is_file() else "missing")
        if freqs.is_file() and post.is_file():
            grid = np.loadtxt(freqs, delimiter=",", ndmin=1)
            matrix = np.loadtxt(post, delimiter=",", ndmin=2)
            if matrix.shape[0] != len(grid) and matrix.shape[1] == len(grid):
                matrix = matrix.T
            shape_ok = matrix.shape[0] == len(grid) == args.expected_df and matrix.shape[1] > 1
            record(f"{label}:posterior_dimensions", "PASS" if shape_ok else "FAIL", f"grid={len(grid)} post={matrix.shape}")
            finite_ok = np.all(np.isfinite(grid)) and np.all(np.isfinite(matrix))
            record(f"{label}:posterior_finite", "PASS" if finite_ok else "FAIL", finite_ok)
        if inference.is_file():
            table = pd.read_csv(inference, sep="\t")
            needed = {"logLR", "-log10(p-value)", "SelectionMLE1"}
            ok = len(table) == 1 and needed.issubset(table.columns) and np.all(np.isfinite(table[list(needed)].to_numpy(float)))
            record(f"{label}:inference_fields", "PASS" if ok else "FAIL", ",".join(table.columns))

    report = results / "validation_report.tsv"
    report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["check", "status", "detail"]).to_csv(report, sep="\t", index=False)
    if any(status == "FAIL" for _, status, _ in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

