#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Join SV326 direct cCRE overlaps to ENCODE4/S SCREEN tissue or biosample z-scores.

Input:
  - sv326_encode4_ccre_hits.tsv from annotate_sv326_with_encode4_ccre.py
  - ENCODE4 z-score matrices, ENCFF*.txt.gz with experiment and biosample headers

The join is accession-based through EH38E/EH38D cCRE IDs, so it is compatible
with cCRE coordinates lifted to T2T/CHM13 as long as the original ENCODE cCRE
accessions are preserved.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CCRE_HITS = SCRIPT_DIR / "SV326_encode4_ccre_out" / "sv326_encode4_ccre_hits.tsv"
DEFAULT_ZSCORE_DIR = Path("data/zscores")
DEFAULT_OUTDIR = SCRIPT_DIR / "SV326_encode4_ccre_tissue_out"

ACCESSION_RE = re.compile(r"^EH(?:38|37)[DE]\d+$")

TISSUE_KEYWORDS = [
    ("vascular_endothelial", ["huvec", "endothelial", "artery", "aorta", "vascular", "vein", "vena", "umbilical vein"]),
    ("lung", ["lung", "pulmonary", "bronch", "alveol"]),
    ("heart", ["heart", "cardiac", "cardiomy", "myocardium", "ventricle", "atrium"]),
    ("blood_immune", ["blood", "eryth", "monocyte", "macrophage", "dendritic", "neutrophil", "immune", "t cell", "b cell", "lymph", "nk cell", "jurkat", "k562"]),
    ("placenta", ["placenta", "trophoblast", "chorion"]),
    ("brain_neural", ["brain", "neuron", "astrocyte", "cortex", "hippocampus", "cingulate", "neural"]),
    ("skeletal_muscle", ["skeletal muscle", "gastrocnemius", "myoblast", "myotube", "muscle"]),
    ("kidney", ["kidney", "renal"]),
    ("liver", ["liver", "hepatic", "hepatocyte", "hepg2"]),
]


def clean_label(value: object) -> str:
    text = str(value).strip().replace('"', "").replace("'", "")
    return text if text else "NA"


def tissue_group(label: object) -> str:
    text = clean_label(label).lower().replace("_", " ")
    for group, keywords in TISSUE_KEYWORDS:
        if any(keyword in text for keyword in keywords):
            return group
    return "other"


def is_accession(value: object) -> bool:
    return bool(ACCESSION_RE.match(str(value).strip()))


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_float(text: str) -> float:
    try:
        value = float(text)
    except ValueError:
        return float("nan")
    return value if math.isfinite(value) else float("nan")


def load_direct_hits(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"cCRE hits file not found: {path}")
    hits = pd.read_csv(path, sep="\t", dtype=str)
    required = {"SV_ID", "relation", "ccre_id", "d_accession", "e_accession", "ccre_class"}
    missing = sorted(required - set(hits.columns))
    if missing:
        raise ValueError(f"cCRE hits file missing columns: {missing}")
    direct = hits.loc[hits["relation"].astype(str).str.lower() == "direct"].copy()
    if direct.empty:
        raise ValueError("No direct cCRE hits found. Run the SV326 cCRE overlap step first.")
    for col in ["PBS", "P_BH_PBSpos", "AF_TBN", "AF_HAN", "AF_EUR", "AF_TBN_minus_AF_HAN", "distance_bp", "overlap_bp"]:
        if col in direct.columns:
            direct[col] = pd.to_numeric(direct[col], errors="coerce")
    direct["direct_hit_key"] = direct["SV_ID"].astype(str) + "|" + direct["ccre_id"].astype(str)
    for col in ["ccre_id", "d_accession", "e_accession"]:
        direct[col] = direct[col].fillna("").astype(str).str.strip()
    return direct


def target_accessions_from_hits(direct_hits: pd.DataFrame) -> Set[str]:
    targets: Set[str] = set()
    for col in ["ccre_id", "e_accession", "d_accession"]:
        targets.update(x for x in direct_hits[col].dropna().astype(str).str.strip() if is_accession(x))
    return targets


def discover_zscore_files(zscore_dir: Path, explicit_files: Sequence[Path]) -> List[Path]:
    if explicit_files:
        return [p if p.is_absolute() else zscore_dir / p for p in explicit_files]
    files = sorted(zscore_dir.glob("ENCFF*.txt.gz"))
    files += sorted(zscore_dir.glob("*.zscore*.txt.gz"))
    files += sorted(zscore_dir.glob("*Z*score*.txt.gz"))
    return list(dict.fromkeys(files))


def find_target_accession(parts: Sequence[str], targets: Set[str], max_scan_columns: int) -> Tuple[Optional[str], Optional[int]]:
    for i, value in enumerate(parts[:max_scan_columns]):
        token = value.strip()
        if token in targets:
            return token, i
    return None, None


def infer_assay_from_filename(path: Path) -> str:
    known = {"ENCFF833CSA": "DNase", "ENCFF560UZM": "H3K27ac", "ENCFF508PCP": "H3K4me3"}
    for accession, assay in known.items():
        if accession in path.name:
            return assay
    name = path.name.lower()
    for assay in ["dnase", "h3k4me3", "h3k27ac", "ctcf"]:
        if assay in name:
            return assay.upper() if assay == "ctcf" else assay
    return path.stem.replace(".txt", "")


def scan_zscore_matrix(
    path: Path,
    targets: Set[str],
    z_threshold: float,
    max_scan_columns: int,
    keep_all_scores: bool,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    stats = {
        "zscore_file": str(path),
        "assay": infer_assay_from_filename(path),
        "header_rows": 0,
        "header_columns": 0,
        "data_rows_scanned": 0,
        "target_rows_found": 0,
        "records_written": 0,
    }

    with open_text(path) as handle:
        experiment_header: Optional[List[str]] = None
        biosample_header: Optional[List[str]] = None
        pending_data_parts: Optional[List[str]] = None

        for line in handle:
            if line.strip() and not line.startswith("#"):
                experiment_header = line.rstrip("\n").split("\t")
                break
        if experiment_header is None:
            return rows, stats

        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] == experiment_header[0] and not is_accession(parts[0]):
                biosample_header = parts
                stats["header_rows"] = 2
            else:
                biosample_header = experiment_header
                pending_data_parts = parts
                stats["header_rows"] = 1
            break

        if biosample_header is None:
            return rows, stats
        stats["header_columns"] = len(experiment_header)
        if len(biosample_header) != len(experiment_header):
            raise ValueError(f"Inconsistent matrix header widths: {path}")

        def process(parts: List[str]) -> None:
            if len(parts) != len(experiment_header):
                raise ValueError(f"Matrix row width differs from header: {path}")
            stats["data_rows_scanned"] += 1
            accession, accession_col = find_target_accession(parts, targets, max_scan_columns)
            if accession is None:
                return
            stats["target_rows_found"] += 1
            for i, raw in enumerate(parts):
                if i == accession_col:
                    continue
                z = parse_float(raw)
                if not np.isfinite(z):
                    continue
                active = z >= z_threshold
                if not active and not keep_all_scores:
                    continue
                experiment = clean_label(experiment_header[i] if i < len(experiment_header) else f"column_{i + 1}")
                biosample = clean_label(biosample_header[i] if i < len(biosample_header) else experiment)
                if is_accession(biosample) or biosample.lower() in {"start", "end", "score"}:
                    continue
                rows.append(
                    {
                        "ccre_accession": accession,
                        "assay": stats["assay"],
                        "zscore_file": path.name,
                        "experiment_accession": experiment,
                        "biosample": biosample,
                        "tissue_group": tissue_group(biosample),
                        "zscore": z,
                        "active_by_z_threshold": active,
                    }
                )

        if pending_data_parts is not None:
            process(pending_data_parts)
        for line in handle:
            if line.strip() and not line.startswith("#"):
                process(line.rstrip("\n").split("\t"))

    stats["records_written"] = len(rows)
    return rows, stats


def load_zscores(args: argparse.Namespace, targets: Set[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    files = discover_zscore_files(args.zscore_dir, args.zscore_files)
    if not files:
        raise FileNotFoundError(f"No z-score matrix files found in {args.zscore_dir}")
    records: List[Dict[str, object]] = []
    diagnostics: List[Dict[str, object]] = []
    for path in files:
        if not path.exists():
            raise FileNotFoundError(f"z-score matrix not found: {path}")
        print(f"[INFO] scan z-score matrix: {path}", flush=True)
        recs, stats = scan_zscore_matrix(path, targets, args.z_threshold, args.max_scan_columns, args.keep_all_scores)
        records.extend(recs)
        diagnostics.append(stats)
    return pd.DataFrame(records), pd.DataFrame(diagnostics)


def join_activity_to_hits(direct_hits: pd.DataFrame, activity: pd.DataFrame) -> pd.DataFrame:
    if activity.empty:
        return pd.DataFrame()
    base_cols = [
        c
        for c in [
            "direct_hit_key",
            "SV_ID",
            "SVTYPE",
            "sv_size_from_id",
            "ccre_id",
            "d_accession",
            "e_accession",
            "ccre_class",
            "ccre_group",
            "PBS",
            "P_BH_PBSpos",
            "AF_TBN",
            "AF_HAN",
            "AF_EUR",
            "AF_TBN_minus_AF_HAN",
            "in_360_overlap",
        ]
        if c in direct_hits.columns
    ]
    maps = []
    for col in ["ccre_id", "e_accession", "d_accession"]:
        tmp = direct_hits[base_cols].copy()
        tmp["ccre_accession"] = direct_hits[col].astype(str).values
        tmp = tmp.loc[tmp["ccre_accession"].map(is_accession)]
        maps.append(tmp)
    mapping = pd.concat(maps, ignore_index=True).drop_duplicates()
    return mapping.merge(activity, on="ccre_accession", how="inner").drop_duplicates()


def collapse_values(values: Iterable[object], max_items: int = 30) -> str:
    seen: List[str] = []
    for value in values:
        text = str(value).strip()
        if not text or text == "nan":
            continue
        if text not in seen:
            seen.append(text)
        if len(seen) >= max_items:
            break
    return ";".join(seen)


def numeric_max(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce")
    return float(values.max()) if values.notna().any() else float("nan")


def make_outputs(direct_hits: pd.DataFrame, joined: pd.DataFrame, outdir: Path) -> None:
    active = joined.loc[joined["active_by_z_threshold"]] if not joined.empty else pd.DataFrame()

    if not active.empty:
        tissue_summary = (
            active.groupby(["tissue_group", "biosample"], dropna=False)
            .agg(
                selected_sv_count=("SV_ID", "nunique"),
                direct_ccre_count=("ccre_id", "nunique"),
                hit_count=("direct_hit_key", "nunique"),
                max_zscore=("zscore", "max"),
                median_zscore=("zscore", "median"),
            )
            .reset_index()
            .sort_values(["selected_sv_count", "direct_ccre_count", "max_zscore"], ascending=[False, False, False])
        )
        tissue_group_summary = (
            active.groupby("tissue_group", dropna=False)
            .agg(
                selected_sv_count=("SV_ID", "nunique"),
                direct_ccre_count=("ccre_id", "nunique"),
                active_biosample_count=("biosample", "nunique"),
                hit_count=("direct_hit_key", "nunique"),
                max_zscore=("zscore", "max"),
                median_zscore=("zscore", "median"),
            )
            .reset_index()
            .sort_values(["selected_sv_count", "direct_ccre_count", "max_zscore"], ascending=[False, False, False])
        )
        sv_tissue_long = (
            active.groupby(["SV_ID", "tissue_group"], dropna=False)
            .agg(
                active_direct_ccre_count=("ccre_id", "nunique"),
                active_biosample_count=("biosample", "nunique"),
                max_zscore=("zscore", "max"),
                biosamples=("biosample", lambda x: collapse_values(x, 25)),
                ccre_classes=("ccre_class", lambda x: collapse_values(x, 10)),
            )
            .reset_index()
            .sort_values(["SV_ID", "active_direct_ccre_count", "max_zscore"], ascending=[True, False, False])
        )
    else:
        tissue_summary = pd.DataFrame()
        tissue_group_summary = pd.DataFrame()
        sv_tissue_long = pd.DataFrame()

    direct_rows = []
    for sv_id, sub in direct_hits.groupby("SV_ID", sort=False):
        sub_active = active.loc[active["SV_ID"] == sv_id] if not active.empty else pd.DataFrame()
        direct_rows.append(
            {
                "SV_ID": sv_id,
                "SVTYPE": collapse_values(sub.get("SVTYPE", []), 3),
                "PBS": numeric_max(sub.get("PBS", pd.Series(dtype=float))),
                "P_BH_PBSpos": numeric_max(sub.get("P_BH_PBSpos", pd.Series(dtype=float))),
                "AF_TBN": numeric_max(sub.get("AF_TBN", pd.Series(dtype=float))),
                "AF_HAN": numeric_max(sub.get("AF_HAN", pd.Series(dtype=float))),
                "AF_EUR": numeric_max(sub.get("AF_EUR", pd.Series(dtype=float))),
                "AF_TBN_minus_AF_HAN": numeric_max(sub.get("AF_TBN_minus_AF_HAN", pd.Series(dtype=float))),
                "direct_cCRE_count": sub["ccre_id"].nunique(),
                "direct_cCRE_classes": collapse_values(sub["ccre_class"], 15),
                "active_tissue_group_count": sub_active["tissue_group"].nunique() if not sub_active.empty else 0,
                "active_biosample_count": sub_active["biosample"].nunique() if not sub_active.empty else 0,
                "active_tissue_groups": collapse_values(sub_active["tissue_group"] if not sub_active.empty else [], 25),
                "max_active_zscore": numeric_max(sub_active["zscore"]) if not sub_active.empty else np.nan,
                "top_active_biosamples": collapse_values(
                    sub_active.sort_values("zscore", ascending=False)["biosample"] if not sub_active.empty else [], 20
                ),
            }
        )
    sv_summary = pd.DataFrame(direct_rows).sort_values(
        ["PBS", "direct_cCRE_count", "active_tissue_group_count", "max_active_zscore"],
        ascending=[False, False, False, False],
    )

    joined.to_csv(outdir / "sv326_direct_ccre_tissue_zscore_long.tsv", sep="\t", index=False)
    tissue_summary.to_csv(outdir / "sv326_direct_ccre_biosample_summary.tsv", sep="\t", index=False)
    tissue_group_summary.to_csv(outdir / "sv326_direct_ccre_tissue_group_summary.tsv", sep="\t", index=False)
    sv_tissue_long.to_csv(outdir / "sv326_direct_sv_tissue_group_long.tsv", sep="\t", index=False)
    sv_summary.to_csv(outdir / "sv326_direct_sv_tissue_summary.tsv", sep="\t", index=False)

    if not sv_tissue_long.empty:
        matrix = sv_tissue_long.pivot_table(
            index="SV_ID",
            columns="tissue_group",
            values="max_zscore",
            aggfunc="max",
        )
        matrix.to_csv(outdir / "sv326_direct_sv_by_tissue_max_zscore_matrix.tsv", sep="\t")
    else:
        pd.DataFrame().to_csv(outdir / "sv326_direct_sv_by_tissue_max_zscore_matrix.tsv", sep="\t", index=False)

    report_lines = [
        "# SV326 Direct cCRE Tissue Z-Score Report",
        "",
        f"- Direct SV-cCRE hit rows: {len(direct_hits)}",
        f"- SVs with direct cCRE overlap: {direct_hits['SV_ID'].nunique()}",
        f"- Active joined z-score rows: {len(active)}",
        f"- SVs with at least one active tissue/biosample cCRE: {active['SV_ID'].nunique() if not active.empty else 0}",
        f"- Tissue groups detected: {tissue_group_summary['tissue_group'].nunique() if not tissue_group_summary.empty else 0}",
        "",
        "## Output Files",
        "- `sv326_direct_ccre_tissue_zscore_long.tsv`: long table of direct cCRE z-score activity.",
        "- `sv326_direct_ccre_biosample_summary.tsv`: biosample-level summary.",
        "- `sv326_direct_ccre_tissue_group_summary.tsv`: tissue-group-level summary.",
        "- `sv326_direct_sv_tissue_group_long.tsv`: SV by tissue-group summary.",
        "- `sv326_direct_sv_by_tissue_max_zscore_matrix.tsv`: SV by tissue-group max z-score matrix.",
        "- `sv326_direct_sv_tissue_summary.tsv`: one row per SV with direct cCRE and tissue summary.",
        "",
        "## Notes",
        "- Activity is defined as z-score >= threshold. Default threshold is 1.64.",
        "- Tissue groups are keyword summaries from ENCODE biosample labels; keep the raw biosample table for manual review.",
        "- This is ENCODE/S SCREEN reference activity, not HUVEC-specific experimental activity.",
    ]
    if not tissue_group_summary.empty:
        report_lines.extend(["", "## Top Tissue Groups"])
        for _, row in tissue_group_summary.head(12).iterrows():
            report_lines.append(
                f"- {row['tissue_group']}: {int(row['selected_sv_count'])} SVs, "
                f"{int(row['direct_ccre_count'])} cCREs, {int(row['active_biosample_count'])} biosamples"
            )
    (outdir / "sv326_direct_ccre_tissue_zscore_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ccre-hits", type=Path, default=DEFAULT_CCRE_HITS)
    parser.add_argument("--zscore-dir", type=Path, default=DEFAULT_ZSCORE_DIR)
    parser.add_argument("--zscore-files", nargs="*", type=Path, default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--z-threshold", type=float, default=1.64)
    parser.add_argument("--max-scan-columns", type=int, default=12)
    parser.add_argument("--keep-all-scores", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    direct_hits = load_direct_hits(args.ccre_hits)
    direct_hits.to_csv(args.outdir / "sv326_direct_ccre_hits_used.tsv", sep="\t", index=False)
    targets = target_accessions_from_hits(direct_hits)
    print(f"[INFO] direct hit rows={len(direct_hits)}; target accessions={len(targets)}", flush=True)
    activity, diagnostics = load_zscores(args, targets)
    diagnostics.to_csv(args.outdir / "sv326_direct_ccre_zscore_scan_diagnostics.tsv", sep="\t", index=False)
    joined = join_activity_to_hits(direct_hits, activity)
    make_outputs(direct_hits, joined, args.outdir)
    print(f"Wrote {args.outdir}")


if __name__ == "__main__":
    main()
