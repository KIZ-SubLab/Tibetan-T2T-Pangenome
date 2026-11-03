#!/usr/bin/env python3
"""Prepare auditable PBS, Denisovan fD, and GENCODE tables for one locus."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def gtf_attrs(text: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+) "([^"]*)"', text)}


def load_gencode(path: Path, chrom: str, start: int, end: int):
    genes: dict[str, dict] = {}
    transcripts: dict[str, dict] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != chrom:
                continue
            feature = fields[2]
            feature_start, feature_end = int(fields[3]), int(fields[4])
            if feature_start > end or feature_end < start:
                continue
            attrs = gtf_attrs(fields[8])
            gene_id = attrs.get("gene_id", "").split(".")[0]
            if not gene_id:
                continue
            if feature == "gene":
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_name": attrs.get("gene_name", gene_id),
                    "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "unknown")),
                    "start": feature_start,
                    "end": feature_end,
                    "strand": fields[6],
                }
            elif feature in {"transcript", "exon"}:
                transcript_id = attrs.get("transcript_id", "").split(".")[0]
                if not transcript_id:
                    continue
                if feature == "transcript":
                    transcripts[transcript_id] = {
                        "transcript_id": transcript_id,
                        "gene_id": gene_id,
                        "gene_name": attrs.get("gene_name", gene_id),
                        "start": feature_start,
                        "end": feature_end,
                        "strand": fields[6],
                        "tags": sorted(set(re.findall(r'tag "([^"]+)"', fields[8]))),
                    }
                else:
                    exons[transcript_id].append((max(start, feature_start), min(end, feature_end)))

    by_gene: dict[str, list[dict]] = defaultdict(list)
    for transcript_id, meta in transcripts.items():
        if transcript_id in exons:
            by_gene[meta["gene_id"]].append(meta)

    def rank(meta: dict):
        tags = set(meta["tags"])
        priority = (
            0 if "MANE_Select" in tags else
            1 if "Ensembl_canonical" in tags else
            2 if any(tag.startswith("appris_principal") for tag in tags) else 3
        )
        return priority, -(meta["end"] - meta["start"]), meta["transcript_id"]

    gene_rows, exon_rows = [], []
    for gene_id, gene in genes.items():
        candidates = by_gene.get(gene_id, [])
        if not candidates:
            continue
        best = sorted(candidates, key=rank)[0]
        gene_rows.append({**gene, "transcript_id": best["transcript_id"], "selection_rank": rank(best)[0]})
        for exon_start, exon_end in sorted(exons[best["transcript_id"]]):
            exon_rows.append({
                "transcript_id": best["transcript_id"],
                "gene_id": gene_id,
                "start": exon_start,
                "end": exon_end,
            })

    gene_columns = [
        "gene_id", "gene_name", "gene_type", "start", "end", "strand",
        "transcript_id", "selection_rank",
    ]
    exon_columns = ["transcript_id", "gene_id", "start", "end"]
    gene_df = pd.DataFrame(gene_rows, columns=gene_columns)
    exon_df = pd.DataFrame(exon_rows, columns=exon_columns)
    return assign_lanes(gene_df, start, end), exon_df


def assign_lanes(frame: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    result = frame.copy()
    if result.empty:
        result["lane"] = pd.Series(dtype=int)
        return result
    span = end - start
    lane_ends: list[int] = []
    assignments: list[tuple[int, int]] = []
    for row in result.sort_values(["start", "end", "gene_name"]).itertuples():
        padded_start = max(start, row.start) - int(span * 0.02)
        padded_end = min(end, row.end) + max(int(span * 0.02), len(str(row.gene_name)) * int(span * 0.006))
        lane = next((idx for idx, last in enumerate(lane_ends) if padded_start > last), None)
        if lane is None:
            lane = len(lane_ends)
            lane_ends.append(padded_end)
        else:
            lane_ends[lane] = padded_end
        assignments.append((row.Index, lane))
    result["lane"] = 0
    for index, lane in assignments:
        result.loc[index, "lane"] = lane
    return result


def read_pbs(path: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    required = ["CHR", "POS", "LOG_R", "LOG_MU", "PBS"]
    selected = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=500_000):
        if list(chunk.columns) != required:
            raise RuntimeError(f"Unexpected PBS columns: {list(chunk.columns)}; expected {required}")
        hit = chunk[(chunk.CHR == chrom) & chunk.POS.between(start, end)]
        if not hit.empty:
            selected.append(hit)
    if not selected:
        raise RuntimeError(f"No PBS sites in {chrom}:{start}-{end}")
    frame = pd.concat(selected, ignore_index=True).sort_values("POS").reset_index(drop=True)
    if frame.duplicated(["CHR", "POS"]).any():
        raise RuntimeError("Duplicate PBS CHR/POS keys in requested interval")
    for column in ["POS", "LOG_R", "LOG_MU", "PBS"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if (~np.isfinite(frame.PBS)).any() or (frame.PBS <= 0).any():
        raise RuntimeError("PBS must be finite and strictly positive")
    frame["regional_PBS_percentile"] = frame.PBS.rank(method="average", pct=True)
    return frame


def read_fd_windows(path: Path, chrom: str, start: int, end: int, archaic_code: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", na_values=["NA"])
    if archaic_code == "GHOST":
        value_columns = [
            "A_relative_fd_estimate", "B_Denisovan_fd_estimate",
            "C_Neanderthal_fd_estimate",
        ]
        required = {"chrom", "start", "end", "midpoint", "n_snps", *value_columns}
    else:
        value_columns = ["fD", "fD_ancestral"]
        required = {"chrom", "start", "end", "midpoint", "n_snps", *value_columns}
    missing = required.difference(frame.columns)
    if missing:
        raise RuntimeError(f"Missing fD columns: {sorted(missing)}")
    frame = frame[(frame.chrom == chrom) & (frame.end >= start) & (frame.start <= end)].copy()
    if frame.empty:
        raise RuntimeError(f"No fD windows overlap {chrom}:{start}-{end}")
    frame = frame.sort_values("midpoint").reset_index(drop=True)
    if frame.midpoint.duplicated().any():
        raise RuntimeError("Duplicate fD window midpoints")
    for column in value_columns:
        if frame[column].isna().any() or (~np.isfinite(frame[column])).any():
            raise RuntimeError(f"Selected fD windows contain non-finite {column}")
    if (frame.n_snps != 100).any():
        raise RuntimeError("Expected exactly 100 SNPs in every selected fD window")
    mid = frame.midpoint.to_numpy(dtype=float)
    gaps = np.diff(mid)
    widths = np.empty(len(frame), dtype=float)
    if len(frame) == 1:
        widths[:] = min(end - start + 1, frame.loc[0, "end"] - frame.loc[0, "start"] + 1)
    else:
        widths[0] = gaps[0]
        widths[-1] = gaps[-1]
        if len(frame) > 2:
            widths[1:-1] = np.minimum(gaps[:-1], gaps[1:])
    frame["display_width"] = np.maximum(1.0, widths * 0.84)
    frame["display_left"] = frame.midpoint - frame.display_width / 2
    frame["display_right"] = frame.midpoint + frame.display_width / 2
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--pbs-file", type=Path, required=True)
    parser.add_argument("--fd-windows", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--archaic-label", required=True)
    parser.add_argument("--archaic-code", choices=["DEN", "NEA", "GHOST"], required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    args = parser.parse_args()
    chrom = args.chrom if args.chrom.startswith("chr") else f"chr{args.chrom}"
    if args.start >= args.end:
        raise ValueError("--start must be smaller than --end")
    data = args.run_root / "data"
    metadata = args.run_root / "metadata"
    data.mkdir(parents=True, exist_ok=True)
    metadata.mkdir(parents=True, exist_ok=True)

    pbs = read_pbs(args.pbs_file, chrom, args.start, args.end)
    fd = read_fd_windows(args.fd_windows, chrom, args.start, args.end, args.archaic_code)
    genes, exons = load_gencode(args.gtf, chrom, args.start, args.end)
    pbs.to_csv(data / "PBS_sites.region.tsv", sep="\t", index=False, float_format="%.10g")
    fd_name = "GHOST_ABC_fD_windows.region.tsv" if args.archaic_code == "GHOST" else f"{args.archaic_code}_fD_windows.region.tsv"
    fd.to_csv(data / fd_name, sep="\t", index=False, float_format="%.10g")
    genes.to_csv(data / "gencode_v50_basic.representative_genes.tsv", sep="\t", index=False)
    exons.to_csv(data / "gencode_v50_basic.representative_exons.tsv", sep="\t", index=False)

    summary = {
        "region": f"{chrom}:{args.start}-{args.end}",
        "coordinate_system": "GRCh38; locus bounds are 1-based inclusive",
        "pbs_source": str(args.pbs_file),
        "pbs_sites": int(len(pbs)),
        "pbs_min": float(pbs.PBS.min()),
        "pbs_max": float(pbs.PBS.max()),
        "pbs_quantiles": {str(q): float(pbs.PBS.quantile(q)) for q in [0.5, 0.9, 0.95, 0.99]},
        "fd_source": str(args.fd_windows),
        "archaic_label": args.archaic_label,
        "fd_definition": (
            "Ghost candidate test A is displayed: D/fD(Denisovan,Neanderthal,TBN,YRI), "
            "a relative Denisovan-vs-Neanderthal affinity statistic. Tests B and C are retained "
            "in the source results; no directly sampled ghost fD exists."
            if args.archaic_code == "GHOST" else
            f"fD(EAS,TBN;{args.archaic_label},YRI), generated by the integrated fd workflow"
        ),
        "fd_windows_overlapping_locus": int(len(fd)),
        "fd_window_scheme": "100 SNPs per window; 50-SNP step",
        "fd_first_window_start": int(fd.start.min()),
        "fd_last_window_end": int(fd.end.max()),
        "fd_min": float(fd.A_relative_fd_estimate.min()) if args.archaic_code == "GHOST" else float(fd.fD.min()),
        "fd_max": float(fd.A_relative_fd_estimate.max()) if args.archaic_code == "GHOST" else float(fd.fD.max()),
        "gene_models": int(len(genes)),
        "gencode": "GENCODE Human release 50 basic annotation, GRCh38.p14",
    }
    (metadata / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
