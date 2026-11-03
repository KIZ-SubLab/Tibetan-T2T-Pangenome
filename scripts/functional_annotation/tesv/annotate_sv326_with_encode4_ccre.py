#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate the SV326 list with human ENCODE4/S SCREEN cCRE regulatory elements.

The SV list can use either the new SV326 column name `svID_5caller` or the old
`SV_ID` name. cCRE coordinates must be in the same genome coordinate system as
the SV IDs. For the T2T-CHM13 analysis, use the ENCODE4 GRCh38 cCRE BED lifted
to CHM13v2.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SV_CSV = SCRIPT_DIR / "SV326.csv"
DEFAULT_CCRE_BED = Path("data/GRCh38-cCREs.liftedToCHM13v2.bed")
DEFAULT_OUTDIR = SCRIPT_DIR / "SV326_encode4_ccre_out"
SV_RE = re.compile(r"^(chr[^-]+)-(\d+)-(\d+)-([A-Za-z]+)(?:-(\d+))?$")


def pick_first_existing(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    available = set(columns)
    for col in candidates:
        if col in available:
            return col
    return None


def parse_sv_id(sv_id: object) -> Dict[str, object]:
    text = str(sv_id).strip()
    m = SV_RE.match(text)
    if not m:
        return {
            "chrom": np.nan,
            "start_1based": np.nan,
            "end_1based": np.nan,
            "sv_type_from_id": np.nan,
            "sv_size_from_id": np.nan,
            "parse_ok": False,
        }
    chrom, start, end, sv_type, size = m.groups()
    start_i = int(start)
    end_i = int(end)
    if start_i < 1 or end_i < start_i:
        raise ValueError(f"Invalid 1-based inclusive SV coordinates: {text}")
    size_i = int(size) if size is not None else max(1, abs(end_i - start_i) + 1)
    return {
        "chrom": chrom,
        "start_1based": start_i,
        "end_1based": end_i,
        "sv_type_from_id": sv_type.upper(),
        "sv_size_from_id": size_i,
        "parse_ok": True,
    }


def sv_interval(row: pd.Series, flank_bp: int = 0, ins_overlap_mode: str = "span") -> Tuple[str, int, int, str]:
    chrom = str(row["chrom"])
    start = int(row["start_1based"])
    end = int(row["end_1based"])
    sv_type = str(row.get("SVTYPE", row.get("sv_type_from_id", ""))).upper()

    if sv_type == "INS":
        if ins_overlap_mode == "span" and end > start:
            start0 = max(0, start - 1 - flank_bp)
            end0 = max(start0 + 1, end + flank_bp)
            return chrom, start0, end0, "INS_span"
        start0 = max(0, start - 1 - flank_bp)
        end0 = max(start0 + 1, start + flank_bp)
        return chrom, start0, end0, "INS_point"

    start0 = max(0, start - 1 - flank_bp)
    end0 = max(start0 + 1, end + flank_bp)
    return chrom, start0, end0, "SV_span"


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if start_a < end_b and start_b < end_a:
        return 0
    if end_a <= start_b:
        return int(start_b - end_a)
    return int(start_a - end_b)


def overlap_bp(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b))


def normalize_ccre_class(text: object) -> str:
    value = str(text).strip()
    if value in {"PLS", "pELS", "dELS", "CA-CTCF", "CA-H3K4me3", "CA-TF", "CA"}:
        return value
    return value if value and value.lower() != "nan" else "NA"


def ccre_group(ccre_class: object) -> str:
    value = normalize_ccre_class(ccre_class)
    if value == "PLS":
        return "promoter_like"
    if value in {"pELS", "dELS"}:
        return "enhancer_like"
    if "CTCF" in value:
        return "ctcf_like"
    if value.startswith("CA"):
        return "chromatin_accessible"
    return "other"


def read_sv_table(path: Path) -> tuple[pd.DataFrame, str]:
    if not path.exists():
        raise FileNotFoundError(f"SV table not found: {path}")
    sv = pd.read_csv(path, encoding="utf-8", encoding_errors="replace")
    id_col = pick_first_existing(sv.columns, ["SV_ID", "svID_5caller", "sv_id", "svID", "ID"])
    if id_col is None:
        raise ValueError(f"Cannot find an SV ID column. Available columns: {list(sv.columns)}")

    sv = sv.copy()
    sv["SV_ID"] = sv[id_col].astype(str).str.strip()
    if sv["SV_ID"].duplicated().any():
        raise ValueError("Duplicate SV IDs; resolve before annotation")
    parsed = pd.DataFrame([parse_sv_id(x) for x in sv["SV_ID"]])
    sv = pd.concat([sv, parsed], axis=1)
    if int(sv["parse_ok"].sum()) != len(sv):
        bad = sv.loc[~sv["parse_ok"], "SV_ID"].head(10).astype(str).tolist()
        raise ValueError(f"Unparseable SV_ID examples: {bad}")

    numeric_candidates = [
        "PBS",
        "P_BH_PBSpos",
        "FST",
        "AF_TBN",
        "AF_TIB",
        "AF_HAN",
        "AF_EUR",
        "AF_DIFF",
        "A_POS_COUNT",
        "W_POS_COUNT",
    ]
    for col in numeric_candidates:
        if col in sv.columns:
            sv[col] = pd.to_numeric(sv[col], errors="coerce")

    if "AF_TBN" in sv.columns and "AF_HAN" in sv.columns and "AF_TBN_minus_AF_HAN" not in sv.columns:
        sv["AF_TBN_minus_AF_HAN"] = sv["AF_TBN"] - sv["AF_HAN"]
    if "AF_TIB" in sv.columns and "AF_HAN" in sv.columns and "AF_TIB_minus_AF_HAN" not in sv.columns:
        sv["AF_TIB_minus_AF_HAN"] = sv["AF_TIB"] - sv["AF_HAN"]
    return sv, id_col


def read_ccre_bed(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"cCRE BED not found: {path}")
    names = ["ccre_chrom", "ccre_start", "ccre_end", "d_accession", "e_accession", "ccre_class"]
    ccre = pd.read_csv(path, sep="\t", header=None, names=names, usecols=list(range(6)), comment="#")
    ccre["ccre_chrom"] = ccre["ccre_chrom"].astype(str)
    for col in ("ccre_start", "ccre_end"):
        ccre[col] = pd.to_numeric(ccre[col], errors="raise")
        if not np.isfinite(ccre[col]).all() or (ccre[col] != np.floor(ccre[col])).any():
            raise ValueError("cCRE BED coordinates must be finite integers")
    if (ccre.ccre_start < 0).any() or (ccre.ccre_end <= ccre.ccre_start).any():
        raise ValueError("Invalid cCRE BED interval")
    ccre["ccre_start"] = ccre["ccre_start"].astype(int)
    ccre["ccre_end"] = ccre["ccre_end"].astype(int)
    ccre["ccre_class"] = ccre["ccre_class"].map(normalize_ccre_class)
    ccre["ccre_group"] = ccre["ccre_class"].map(ccre_group)
    ccre["ccre_id"] = ccre["e_accession"].fillna(ccre["d_accession"]).astype(str)
    return ccre


def overlap_rows(ccre_chr: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    return ccre_chr[(ccre_chr["ccre_start"] < end) & (ccre_chr["ccre_end"] > start)].copy()


def class_counts(classes: pd.Series) -> str:
    if classes.empty:
        return ""
    counts = classes.astype(str).value_counts()
    return ";".join(f"{k}:{v}" for k, v in counts.items())


def annotate_sv(
    row: pd.Series,
    ccre_by_chrom: Dict[str, pd.DataFrame],
    nearby_bp: int,
    ins_overlap_mode: str,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    chrom, start, end, rule = sv_interval(row, 0, ins_overlap_mode)
    # A gap exactly equal to the window is included; direct overlap stays >=1 bp.
    _chrom, near_start, near_end, _rule = sv_interval(row, nearby_bp + 1, ins_overlap_mode)
    ccre_chr = ccre_by_chrom.get(chrom, pd.DataFrame())
    direct = overlap_rows(ccre_chr, start, end) if not ccre_chr.empty else pd.DataFrame()
    nearby_all = overlap_rows(ccre_chr, near_start, near_end) if not ccre_chr.empty else pd.DataFrame()
    nearby = nearby_all.loc[~nearby_all.index.isin(direct.index)].copy() if not nearby_all.empty and nearby_bp > 0 else pd.DataFrame()

    input_meta_cols = [
        c
        for c in [
            "PBS",
            "P_BH_PBSpos",
            "FST",
            "AF_TBN",
            "AF_TIB",
            "AF_HAN",
            "AF_EUR",
            "AF_TBN_minus_AF_HAN",
            "AF_TIB_minus_AF_HAN",
            "AF_DIFF",
            "in_360_overlap",
        ]
        if c in row.index
    ]

    hits: List[Dict[str, object]] = []
    for relation, hit_df in [("direct", direct), ("nearby", nearby)]:
        for _, hit in hit_df.iterrows():
            distance = interval_distance(start, end, int(hit["ccre_start"]), int(hit["ccre_end"]))
            ov = overlap_bp(start, end, int(hit["ccre_start"]), int(hit["ccre_end"]))
            record = {
                "SV_ID": row["SV_ID"],
                "relation": relation,
                "distance_bp": distance,
                "overlap_bp": ov,
                "SVTYPE": row.get("SVTYPE", row.get("sv_type_from_id", "")),
                "sv_size_from_id": row.get("sv_size_from_id", np.nan),
                "sv_chrom": chrom,
                "sv_start_0based": start,
                "sv_end_0based": end,
                "sv_overlap_rule": rule,
                "ccre_chrom": hit["ccre_chrom"],
                "ccre_start": int(hit["ccre_start"]),
                "ccre_end": int(hit["ccre_end"]),
                "ccre_id": hit["ccre_id"],
                "d_accession": hit["d_accession"],
                "e_accession": hit["e_accession"],
                "ccre_class": hit["ccre_class"],
                "ccre_group": hit["ccre_group"],
            }
            for col in input_meta_cols:
                record[col] = row.get(col, np.nan)
            hits.append(record)

    direct_classes = direct["ccre_class"] if not direct.empty else pd.Series(dtype=str)
    nearby_classes = nearby["ccre_class"] if not nearby.empty else pd.Series(dtype=str)
    direct_groups = set(direct["ccre_group"].astype(str)) if not direct.empty else set()
    nearby_groups = set(nearby["ccre_group"].astype(str)) if not nearby.empty else set()

    summary = {
        "SV_ID": row["SV_ID"],
        "SVTYPE": row.get("SVTYPE", row.get("sv_type_from_id", "")),
        "sv_size_from_id": row.get("sv_size_from_id", np.nan),
        "chrom": chrom,
        "start_1based": row["start_1based"],
        "end_1based": row["end_1based"],
        "sv_overlap_rule": rule,
        "direct_encode4_ccre": "YES" if len(direct) > 0 else "NO",
        "nearby_encode4_ccre_within_window": "YES" if len(nearby) > 0 else "NO",
        "direct_cCRE_count": len(direct),
        "nearby_cCRE_count": len(nearby),
        "direct_cCRE_classes": class_counts(direct_classes),
        "nearby_cCRE_classes": class_counts(nearby_classes),
        "direct_promoter_like": "YES" if "promoter_like" in direct_groups else "NO",
        "direct_enhancer_like": "YES" if "enhancer_like" in direct_groups else "NO",
        "direct_ctcf_like": "YES" if "ctcf_like" in direct_groups else "NO",
        "nearby_promoter_like": "YES" if "promoter_like" in nearby_groups else "NO",
        "nearby_enhancer_like": "YES" if "enhancer_like" in nearby_groups else "NO",
        "nearby_ctcf_like": "YES" if "ctcf_like" in nearby_groups else "NO",
        "nearest_cCRE_distance_bp": np.nan,
        "nearest_cCRE_id": "",
        "nearest_cCRE_class": "",
        "nearest_cCRE_relation": "",
    }
    for col in input_meta_cols:
        summary[col] = row.get(col, np.nan)

    if hits:
        nearest = sorted(hits, key=lambda x: (x["distance_bp"], 0 if x["relation"] == "direct" else 1))[0]
        summary["nearest_cCRE_distance_bp"] = nearest["distance_bp"]
        summary["nearest_cCRE_id"] = nearest["ccre_id"]
        summary["nearest_cCRE_class"] = nearest["ccre_class"]
        summary["nearest_cCRE_relation"] = nearest["relation"]
    return hits, summary


def counter_from_class_strings(values: pd.Series) -> Dict[str, int]:
    counter: Dict[str, int] = {}
    for text in values.dropna().astype(str):
        for token in text.split(";"):
            if not token:
                continue
            cls, count = token.rsplit(":", 1)
            counter[cls] = counter.get(cls, 0) + int(count)
    return counter


def write_report(
    outdir: Path,
    args: argparse.Namespace,
    sv: pd.DataFrame,
    id_col: str,
    ccre: pd.DataFrame,
    summary: pd.DataFrame,
    hits: pd.DataFrame,
) -> None:
    direct_n = int((summary["direct_cCRE_count"] > 0).sum())
    nearby_only_n = int(((summary["direct_cCRE_count"] == 0) & (summary["nearby_cCRE_count"] > 0)).sum())
    no_nearby_n = int(((summary["direct_cCRE_count"] == 0) & (summary["nearby_cCRE_count"] == 0)).sum())
    direct_classes = counter_from_class_strings(summary["direct_cCRE_classes"])
    nearby_classes = counter_from_class_strings(summary["nearby_cCRE_classes"])

    lines = [
        "# SV326 ENCODE4 cCRE Annotation Report",
        "",
        "## Inputs",
        f"- SV table: `{args.sv_csv}`",
        f"- SV ID column: `{id_col}`",
        f"- cCRE BED: `{args.ccre_bed}`",
        f"- cCRE coordinate system: `{args.ccre_coordinate_system}`",
        f"- nearby window: {args.nearby_bp} bp",
        f"- INS overlap mode: `{args.ins_overlap_mode}`",
        "",
        "## Main Counts",
        f"- SVs analyzed: {len(sv)}",
        f"- cCRE records loaded: {len(ccre)}",
        f"- direct SV-cCRE overlaps: {direct_n} SVs",
        f"- nearby-only cCREs within {args.nearby_bp} bp: {nearby_only_n} SVs",
        f"- no cCRE within {args.nearby_bp} bp: {no_nearby_n} SVs",
        f"- total direct/nearby hit rows: {len(hits)}",
        "",
        "## Direct cCRE Classes",
    ]
    if direct_classes:
        for cls, count in sorted(direct_classes.items(), key=lambda x: (-x[1], x[0])):
            lines.append(f"- {cls}: {count}")
    else:
        lines.append("- none")

    lines.append("")
    lines.append("## Nearby cCRE Classes")
    if nearby_classes:
        for cls, count in sorted(nearby_classes.items(), key=lambda x: (-x[1], x[0])):
            lines.append(f"- {cls}: {count}")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Output Files",
            "- `sv326_encode4_ccre_hits.tsv`: one row per SV-cCRE direct or nearby hit.",
            "- `sv326_encode4_ccre_summary.tsv`: one row per SV with YES/NO cCRE annotation.",
            "- `sv326_encode4_ccre_report.md`: this summary report.",
            "",
            "## Interpretation Notes",
            "- `direct_encode4_ccre=YES` means the SV interval directly overlaps a known ENCODE4/S SCREEN cCRE in the selected coordinate system.",
            "- `nearby_encode4_ccre_within_window=YES` means no direct overlap is required; the cCRE is within the configured flanking window.",
            "- For T2T SVs, use a cCRE BED lifted from GRCh38 to CHM13v2. Mixing GRCh38 cCRE coordinates with CHM13/T2T SV coordinates will give incorrect overlaps.",
            "- The 6-column ENCODE cCRE registry BED is a catalog of known regulatory elements. It is not HUVEC-specific activity evidence by itself.",
        ]
    )
    (outdir / "sv326_encode4_ccre_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-csv", type=Path, default=DEFAULT_SV_CSV)
    parser.add_argument("--ccre-bed", type=Path, default=DEFAULT_CCRE_BED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--nearby-bp", type=int, default=5000)
    parser.add_argument("--ins-overlap-mode", choices=["span", "point"], default="span")
    parser.add_argument("--ccre-coordinate-system", default="CHM13v2_lifted_from_GRCh38")
    args = parser.parse_args()
    if args.nearby_bp < 0:
        parser.error("--nearby-bp must be nonnegative")

    args.outdir.mkdir(parents=True, exist_ok=True)
    sv, id_col = read_sv_table(args.sv_csv)
    ccre = read_ccre_bed(args.ccre_bed)
    ccre_by_chrom = {chrom: sub.copy() for chrom, sub in ccre.groupby("ccre_chrom", sort=False)}

    all_hits: List[Dict[str, object]] = []
    summaries: List[Dict[str, object]] = []
    for _, row in sv.iterrows():
        hits, summary = annotate_sv(row, ccre_by_chrom, args.nearby_bp, args.ins_overlap_mode)
        all_hits.extend(hits)
        summaries.append(summary)

    hit_df = pd.DataFrame(all_hits)
    summary_df = pd.DataFrame(summaries)

    hit_path = args.outdir / "sv326_encode4_ccre_hits.tsv"
    summary_path = args.outdir / "sv326_encode4_ccre_summary.tsv"
    hit_df.to_csv(hit_path, sep="\t", index=False)
    summary_df.to_csv(summary_path, sep="\t", index=False)
    write_report(args.outdir, args, sv, id_col, ccre, summary_df, hit_df)

    print(f"Wrote {hit_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {args.outdir / 'sv326_encode4_ccre_report.md'}")


if __name__ == "__main__":
    main()
