#!/usr/bin/env python3
"""Validate tables and exported artifacts for the PBS/fD locus figure."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--archaic-code", choices=["DEN", "NEA", "GHOST"], required=True)
    parser.add_argument("--archaic-label", required=True)
    parser.add_argument("--pbs-ymax", type=float)
    args = parser.parse_args()
    chrom = args.chrom if args.chrom.startswith("chr") else f"chr{args.chrom}"
    root = args.run_root
    pbs = pd.read_csv(root / "data/PBS_sites.region.tsv", sep="\t")
    fd_name = "GHOST_ABC_fD_windows.region.tsv" if args.archaic_code == "GHOST" else f"{args.archaic_code}_fD_windows.region.tsv"
    fd = pd.read_csv(root / "data" / fd_name, sep="\t")
    summary = json.loads((root / "metadata/analysis_summary.json").read_text())
    spec = json.loads((root / "metadata/figure_spec.json").read_text())
    checks = []

    def check(condition: bool, label: str):
        if not condition:
            raise AssertionError(label)
        checks.append(f"PASS: {label}")

    check(len(pbs) == summary["pbs_sites"], "PBS row count agrees with metadata")
    check(not pbs.duplicated(["CHR", "POS"]).any(), "PBS coordinates are unique")
    check((pbs.CHR == chrom).all() and pbs.POS.between(args.start, args.end).all(), "PBS coordinates are inside locus")
    check(np.isfinite(pbs.PBS).all() and (pbs.PBS > 0).all(), "PBS values are finite and positive")
    check(np.isclose(pbs.PBS.max(), summary["pbs_max"]), "PBS color normalization maximum is traceable")
    check(len(fd) == summary["fd_windows_overlapping_locus"], "fD row count agrees with metadata")
    check((fd.n_snps == 100).all(), "all fD windows contain 100 SNPs")
    if args.archaic_code == "GHOST":
        check(np.isfinite(fd.A_relative_fd_estimate).all(),
              "all displayed test A relative fD values are finite")
        check("only test A" in spec["background_bars"] and "A_relative_fd_estimate" in spec["background_bars"],
              "ghost background is explicitly restricted to test A")
        check("not a directly observed ghost fD" in spec["background_bars"],
              "ghost mode explicitly avoids claiming a directly observed ghost fD")
    else:
        check(np.isfinite(fd.fD).all(), "all displayed fD values are finite")
        check(spec["background_bars"].startswith(f"fD(EAS,TBN;{args.archaic_label},YRI)"), "displayed fD definition is explicit")
    check(spec["pbs_y_axis"]["maximum"] >= float(pbs.PBS.max()), "PBS y-axis does not clip source points")
    if args.pbs_ymax is not None:
        check(np.isclose(spec["pbs_y_axis"]["maximum"], args.pbs_ymax), "user-defined PBS y-axis maximum is preserved")

    prefix = root / "figures" / f"{chrom}_{args.start}_{args.end}_PBS_{args.archaic_code}_fD"
    for suffix in [".svg", ".pdf", ".png", ".tiff"]:
        path = prefix.with_suffix(suffix)
        check(path.exists() and path.stat().st_size > 1000, f"{suffix[1:].upper()} exists and is non-empty")
    svg_text = prefix.with_suffix(".svg").read_text(errors="replace")
    check("<text" in svg_text, "SVG contains editable text nodes")
    if args.archaic_code == "GHOST":
        check("ghost-candidate test A" in svg_text and "Test A" in svg_text,
              "SVG labels ghost-candidate mode and test A")
        check("test B" not in svg_text and "test C" not in svg_text,
              "SVG does not display tests B or C")
    else:
        check(args.archaic_label in svg_text, f"SVG labels {args.archaic_label} fD")
    with Image.open(prefix.with_suffix(".png")) as image:
        check(image.width >= 2000 and image.height >= 1600, "PNG is at least 2000 x 1600 pixels")
    with Image.open(prefix.with_suffix(".tiff")) as image:
        check(image.info.get("dpi", (0, 0))[0] >= 590, "TIFF is exported at approximately 600 dpi")

    report = "\n".join(checks) + f"\nPBS_sites={len(pbs)}\nfD_windows={len(fd)}\n"
    (root / "qa/final_artifact_checks.txt").write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
