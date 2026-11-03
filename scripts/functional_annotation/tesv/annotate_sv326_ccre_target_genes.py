#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate SV326 direct cCRE overlaps with ENCODE4/S SCREEN target-gene links.

The join is accession-based. This is important because the cCRE coordinates may
have been lifted from GRCh38 to CHM13/T2T, while SCREEN gene-link files remain
keyed by stable EH38E cCRE accessions.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CCRE_HITS = SCRIPT_DIR / "sv326_encode4_ccre_hits.tsv"
DEFAULT_GENE_LINK_DIR = Path("data/Human_Gene_Links")
DEFAULT_OUTDIR = SCRIPT_DIR / "SV326_encode4_ccre_target_gene_out"

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


def clean_text(value: object) -> str:
    text = str(value).strip().replace('"', "").replace("'", "")
    return text if text and text.lower() != "nan" else "NA"


def tissue_group(label: object) -> str:
    text = clean_text(label).lower().replace("_", " ")
    for group, keywords in TISSUE_KEYWORDS:
        if any(keyword in text for keyword in keywords):
            return group
    return "other"


def is_accession(value: object) -> bool:
    return bool(ACCESSION_RE.match(str(value).strip()))


def parse_float(value: object) -> float:
    try:
        parsed = float(str(value).strip())
    except Exception:
        return float("nan")
    return parsed if math.isfinite(parsed) else float("nan")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def collapse_values(values: Iterable[object], max_items: int = 80) -> str:
    seen: List[str] = []
    for value in values:
        text = clean_text(value)
        if text == "NA" or text in seen:
            continue
        seen.append(text)
        if len(seen) >= max_items:
            break
    return ";".join(seen) if seen else "NA"


def load_direct_hits(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"cCRE hits file not found: {path}")
    hits = pd.read_csv(path, sep="\t", dtype=str)
    required = {"SV_ID", "relation", "ccre_id", "d_accession", "e_accession", "ccre_class"}
    missing = sorted(required - set(hits.columns))
    if missing:
        raise ValueError(f"cCRE hits file missing required columns: {missing}")
    direct = hits.loc[hits["relation"].astype(str).str.lower() == "direct"].copy()
    if direct.empty:
        raise ValueError("No direct cCRE hits found. Use sv326_encode4_ccre_hits.tsv from the overlap step.")
    for col in [
        "PBS",
        "P_BH_PBSpos",
        "AF_TBN",
        "AF_HAN",
        "AF_EUR",
        "AF_TBN_minus_AF_HAN",
        "distance_bp",
        "overlap_bp",
        "ccre_start",
        "ccre_end",
    ]:
        if col in direct.columns:
            direct[col] = pd.to_numeric(direct[col], errors="coerce")
        else:
            direct[col] = np.nan
    for col in ["ccre_id", "d_accession", "e_accession"]:
        direct[col] = direct[col].fillna("").astype(str).str.strip()
    direct["direct_hit_key"] = direct["SV_ID"].astype(str) + "|" + direct["ccre_id"].astype(str)
    return direct


def target_accessions_from_hits(direct_hits: pd.DataFrame) -> Set[str]:
    targets: Set[str] = set()
    for col in ["ccre_id", "e_accession", "d_accession"]:
        targets.update(x for x in direct_hits[col].dropna().astype(str).str.strip() if is_accession(x))
    return targets


def link_type_from_path(path: Path) -> str:
    name = path.name
    if "CRISPR" in name:
        return "CRISPR"
    if "3D" in name or "Chromatin" in name:
        return "3D-Chromatin"
    if "eQTL" in name or "eQTLs" in name:
        return "eQTL"
    return name


def split_eqtl_context(context: str) -> tuple[str, str]:
    context = clean_text(context)
    if context.startswith("GTEx:"):
        parts = context.split(maxsplit=1)
        dataset = parts[0]
        tissue = parts[1] if len(parts) > 1 else "NA"
        return dataset, tissue
    return context, context


def parse_eqtl_fields(parts: List[str]) -> tuple[str, str, float, float]:
    """Parse SCREEN eQTL links, supporting both 8-column and 9-column exports."""
    if len(parts) >= 9:
        dataset = clean_text(parts[5])
        tissue = clean_text(parts[6])
        effect = parse_float(parts[7])
        pvalue = parse_float(parts[8])
        return dataset, tissue, effect, pvalue
    dataset, tissue = split_eqtl_context(parts[5] if len(parts) > 5 else "")
    effect = parse_float(parts[6] if len(parts) > 6 else "")
    pvalue = parse_float(parts[7] if len(parts) > 7 else "")
    return dataset, tissue, effect, pvalue


def read_gene_link_file(path: Path, targets: Set[str]) -> tuple[pd.DataFrame, Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    link_type = link_type_from_path(path)
    scanned = 0
    target_rows = 0

    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            scanned += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            accession = parts[0].strip()
            if accession not in targets:
                continue
            target_rows += 1

            gene_id = clean_text(parts[1] if len(parts) > 1 else "")
            gene_symbol = clean_text(parts[2] if len(parts) > 2 else gene_id)
            gene_biotype = clean_text(parts[3] if len(parts) > 3 else "")

            evidence_id = "NA"
            method = "NA"
            dataset = "NA"
            biosample_or_tissue = "NA"
            effect_or_score = float("nan")
            p_or_q = float("nan")

            if link_type == "3D-Chromatin":
                method = clean_text(parts[4] if len(parts) > 4 else "")
                dataset = clean_text(parts[5] if len(parts) > 5 else "")
                biosample_or_tissue = clean_text(parts[6] if len(parts) > 6 else "")
                effect_or_score = parse_float(parts[7] if len(parts) > 7 else "")
                p_or_q = parse_float(parts[8] if len(parts) > 8 else "")
            elif link_type == "CRISPR":
                evidence_id = clean_text(parts[4] if len(parts) > 4 else "")
                method = clean_text(parts[5] if len(parts) > 5 else "")
                dataset = clean_text(parts[6] if len(parts) > 6 else "")
                biosample_or_tissue = clean_text(parts[7] if len(parts) > 7 else "")
                effect_or_score = parse_float(parts[8] if len(parts) > 8 else "")
                p_or_q = parse_float(parts[9] if len(parts) > 9 else "")
            elif link_type == "eQTL":
                evidence_id = clean_text(parts[4] if len(parts) > 4 else "")
                dataset, biosample_or_tissue, effect_or_score, p_or_q = parse_eqtl_fields(parts)
                method = "GTEx eQTL"
            else:
                evidence_id = clean_text(parts[4] if len(parts) > 4 else "")
                method = clean_text(parts[5] if len(parts) > 5 else "")
                dataset = clean_text(parts[6] if len(parts) > 6 else "")
                biosample_or_tissue = clean_text(parts[7] if len(parts) > 7 else "")
                effect_or_score = parse_float(parts[8] if len(parts) > 8 else "")
                p_or_q = parse_float(parts[9] if len(parts) > 9 else "")

            rows.append(
                {
                    "ccre_accession": accession,
                    "link_type": link_type,
                    "source_file": path.name,
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "gene_biotype": gene_biotype,
                    "evidence_id": evidence_id,
                    "method": method,
                    "dataset": dataset,
                    "biosample_or_tissue": biosample_or_tissue,
                    "tissue_group": tissue_group(biosample_or_tissue),
                    "effect_or_score": effect_or_score,
                    "p_or_q": p_or_q,
                }
            )

    stats = {
        "gene_link_file": str(path),
        "link_type": link_type,
        "data_rows_scanned": scanned,
        "target_rows_found": target_rows,
        "records_written": len(rows),
    }
    return pd.DataFrame(rows), stats


def discover_gene_link_files(gene_link_dir: Path, explicit_files: Sequence[Path]) -> List[Path]:
    if explicit_files:
        return [p if p.is_absolute() else gene_link_dir / p for p in explicit_files]
    files = [
        gene_link_dir / "V4-hg38.Gene-Links.3D-Chromatin.txt",
        gene_link_dir / "V4-hg38.Gene-Links.CRISPR.txt",
        gene_link_dir / "V4-hg38.Gene-Links.eQTLs.txt",
    ]
    return files


def load_gene_links(gene_link_dir: Path, gene_link_files: Sequence[Path], targets: Set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables: List[pd.DataFrame] = []
    diagnostics: List[Dict[str, object]] = []
    for path in discover_gene_link_files(gene_link_dir, gene_link_files):
        if not path.exists():
            raise FileNotFoundError(f"Requested gene-link file missing: {path}")
        print(f"[INFO] scan gene-link file: {path}", flush=True)
        table, stats = read_gene_link_file(path, targets)
        stats["status"] = "ok"
        diagnostics.append(stats)
        if not table.empty:
            tables.append(table)
    if not tables:
        return pd.DataFrame(), pd.DataFrame(diagnostics)
    links = pd.concat(tables, ignore_index=True).drop_duplicates()
    return links, pd.DataFrame(diagnostics)


def join_links_to_hits(direct_hits: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    if links.empty:
        return pd.DataFrame()
    base_cols = [
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
    base_cols = [c for c in base_cols if c in direct_hits.columns]
    maps = []
    for col in ["ccre_id", "e_accession", "d_accession"]:
        tmp = direct_hits[base_cols].copy()
        tmp["ccre_accession"] = direct_hits[col].astype(str).values
        tmp = tmp.loc[tmp["ccre_accession"].map(is_accession)]
        maps.append(tmp)
    mapping = pd.concat(maps, ignore_index=True).drop_duplicates()
    joined = mapping.merge(links, on="ccre_accession", how="inner")
    return joined.drop_duplicates()


def max_abs(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    if not numeric.notna().any():
        return float("nan")
    return float(np.nanmax(np.abs(numeric.to_numpy(dtype=float))))


def make_outputs(direct_hits: pd.DataFrame, link_hits: pd.DataFrame, outdir: Path) -> None:
    direct_hits.to_csv(outdir / "sv326_direct_ccre_hits_used.tsv", sep="\t", index=False)
    link_hits.to_csv(outdir / "sv326_direct_ccre_gene_link_long.tsv", sep="\t", index=False)

    if link_hits.empty:
        for filename in [
            "sv326_direct_ccre_target_gene_summary.tsv",
            "sv326_direct_sv_target_gene_summary.tsv",
            "sv326_target_gene_summary.tsv",
            "sv326_target_gene_link_type_summary.tsv",
        ]:
            pd.DataFrame(columns=["SV_ID", "gene_symbol", "link_type"]).to_csv(outdir / filename, sep="\t", index=False)
        return

    ccre_gene = (
        link_hits.groupby(
            [
                "SV_ID",
                "ccre_id",
                "ccre_class",
                "ccre_group",
                "gene_symbol",
                "gene_id",
                "gene_biotype",
            ],
            dropna=False,
        )
        .agg(
            link_types=("link_type", lambda x: collapse_values(x, 10)),
            link_evidence_count=("link_type", "size"),
            assay_or_method_count=("method", "nunique"),
            methods=("method", lambda x: collapse_values(x, 30)),
            datasets=("dataset", lambda x: collapse_values(x, 30)),
            biosamples_or_tissues=("biosample_or_tissue", lambda x: collapse_values(x, 60)),
            tissue_groups=("tissue_group", lambda x: collapse_values(x, 20)),
            max_abs_effect_or_score=("effect_or_score", max_abs),
            min_p_or_q=("p_or_q", "min"),
            PBS=("PBS", "max"),
            P_BH_PBSpos=("P_BH_PBSpos", "min"),
            AF_TBN=("AF_TBN", "max"),
            AF_HAN=("AF_HAN", "max"),
            AF_EUR=("AF_EUR", "max"),
            AF_TBN_minus_AF_HAN=("AF_TBN_minus_AF_HAN", "max"),
        )
        .reset_index()
        .sort_values(
            ["SV_ID", "link_evidence_count", "min_p_or_q", "max_abs_effect_or_score"],
            ascending=[True, False, True, False],
        )
    )
    ccre_gene.to_csv(outdir / "sv326_direct_ccre_target_gene_summary.tsv", sep="\t", index=False)

    sv_summary = (
        link_hits.groupby("SV_ID", dropna=False)
        .agg(
            SVTYPE=("SVTYPE", lambda x: collapse_values(x, 3)),
            PBS=("PBS", "max"),
            P_BH_PBSpos=("P_BH_PBSpos", "min"),
            AF_TBN=("AF_TBN", "max"),
            AF_HAN=("AF_HAN", "max"),
            AF_EUR=("AF_EUR", "max"),
            AF_TBN_minus_AF_HAN=("AF_TBN_minus_AF_HAN", "max"),
            direct_cCRE_count=("ccre_id", "nunique"),
            direct_cCRE_classes=("ccre_class", lambda x: collapse_values(x, 20)),
            linked_gene_count=("gene_symbol", "nunique"),
            linked_genes=("gene_symbol", lambda x: collapse_values(x, 120)),
            protein_coding_gene_count=("gene_biotype", lambda x: int((x.astype(str) == "protein_coding").sum())),
            link_types=("link_type", lambda x: collapse_values(x, 10)),
            gene_link_record_count=("link_type", "size"),
            methods=("method", lambda x: collapse_values(x, 40)),
            tissue_groups=("tissue_group", lambda x: collapse_values(x, 20)),
            biosamples_or_tissues=("biosample_or_tissue", lambda x: collapse_values(x, 80)),
            min_p_or_q=("p_or_q", "min"),
            max_abs_effect_or_score=("effect_or_score", max_abs),
        )
        .reset_index()
        .sort_values(
            ["PBS", "AF_TBN_minus_AF_HAN", "linked_gene_count", "gene_link_record_count"],
            ascending=[False, False, False, False],
        )
    )
    sv_summary.to_csv(outdir / "sv326_direct_sv_target_gene_summary.tsv", sep="\t", index=False)

    gene_summary = (
        link_hits.groupby(["gene_symbol", "gene_id", "gene_biotype"], dropna=False)
        .agg(
            selected_sv_count=("SV_ID", "nunique"),
            direct_cCRE_count=("ccre_id", "nunique"),
            gene_link_record_count=("link_type", "size"),
            link_types=("link_type", lambda x: collapse_values(x, 10)),
            SV_IDs=("SV_ID", lambda x: collapse_values(x, 80)),
            cCRE_classes=("ccre_class", lambda x: collapse_values(x, 20)),
            methods=("method", lambda x: collapse_values(x, 40)),
            tissue_groups=("tissue_group", lambda x: collapse_values(x, 20)),
            biosamples_or_tissues=("biosample_or_tissue", lambda x: collapse_values(x, 80)),
            min_p_or_q=("p_or_q", "min"),
            max_abs_effect_or_score=("effect_or_score", max_abs),
            max_PBS=("PBS", "max"),
            max_AF_TBN_minus_AF_HAN=("AF_TBN_minus_AF_HAN", "max"),
        )
        .reset_index()
        .sort_values(
            ["selected_sv_count", "direct_cCRE_count", "gene_link_record_count", "max_PBS"],
            ascending=[False, False, False, False],
        )
    )
    gene_summary.to_csv(outdir / "sv326_target_gene_summary.tsv", sep="\t", index=False)

    link_type_summary = (
        link_hits.groupby("link_type", dropna=False)
        .agg(
            selected_sv_count=("SV_ID", "nunique"),
            direct_cCRE_count=("ccre_id", "nunique"),
            linked_gene_count=("gene_symbol", "nunique"),
            gene_link_record_count=("link_type", "size"),
            min_p_or_q=("p_or_q", "min"),
            max_abs_effect_or_score=("effect_or_score", max_abs),
        )
        .reset_index()
        .sort_values(["selected_sv_count", "linked_gene_count"], ascending=False)
    )
    link_type_summary.to_csv(outdir / "sv326_target_gene_link_type_summary.tsv", sep="\t", index=False)


def write_report(outdir: Path, direct_hits: pd.DataFrame, links: pd.DataFrame, link_hits: pd.DataFrame, diagnostics: pd.DataFrame) -> None:
    sv_summary_path = outdir / "sv326_direct_sv_target_gene_summary.tsv"
    gene_summary_path = outdir / "sv326_target_gene_summary.tsv"
    link_type_path = outdir / "sv326_target_gene_link_type_summary.tsv"

    linked_sv_count = int(link_hits["SV_ID"].nunique()) if not link_hits.empty else 0
    linked_gene_count = int(link_hits["gene_symbol"].nunique()) if not link_hits.empty else 0
    linked_ccre_count = int(link_hits["ccre_id"].nunique()) if not link_hits.empty else 0

    top_sv_text = "No linked SVs found."
    if sv_summary_path.exists() and sv_summary_path.stat().st_size > 0:
        sv_summary = pd.read_csv(sv_summary_path, sep="\t")
        if not sv_summary.empty:
            lines = []
            for _, row in sv_summary.head(12).iterrows():
                lines.append(
                    f"- {row['SV_ID']}: PBS={row.get('PBS', np.nan):.3g}, "
                    f"AF_TBN-HAN={row.get('AF_TBN_minus_AF_HAN', np.nan):.3g}, "
                    f"genes={row.get('linked_gene_count', 'NA')}, "
                    f"link_types={row.get('link_types', 'NA')}, "
                    f"top_genes={row.get('linked_genes', 'NA')}"
                )
            top_sv_text = "\n".join(lines)

    top_gene_text = "No linked genes found."
    if gene_summary_path.exists() and gene_summary_path.stat().st_size > 0:
        gene_summary = pd.read_csv(gene_summary_path, sep="\t")
        if not gene_summary.empty:
            top_gene_text = "\n".join(
                f"- {row['gene_symbol']}: {row['selected_sv_count']} SVs, "
                f"{row['direct_cCRE_count']} cCREs, link_types={row['link_types']}"
                for _, row in gene_summary.head(15).iterrows()
            )

    link_type_text = "No link-type summary available."
    if link_type_path.exists() and link_type_path.stat().st_size > 0:
        link_type_summary = pd.read_csv(link_type_path, sep="\t")
        if not link_type_summary.empty:
            link_type_text = "\n".join(
                f"- {row['link_type']}: {row['selected_sv_count']} SVs, "
                f"{row['direct_cCRE_count']} cCREs, {row['linked_gene_count']} genes, "
                f"{row['gene_link_record_count']} records"
                for _, row in link_type_summary.iterrows()
            )

    report = f"""# SV326 cCRE Target-Gene Link Report

## Run Summary

- Direct SV-cCRE hit rows used: {len(direct_hits)}
- SVs with direct cCRE overlap: {direct_hits['SV_ID'].nunique()}
- cCRE accessions queried in gene-link files: {len(target_accessions_from_hits(direct_hits))}
- Gene-link records loaded before SV join: {len(links)}
- Joined SV-cCRE-gene evidence rows: {len(link_hits)}
- SVs with at least one target-gene link: {linked_sv_count}
- Direct cCREs with at least one target-gene link: {linked_ccre_count}
- Unique linked genes: {linked_gene_count}

## Link-Type Summary

{link_type_text}

## Top SVs by Selection Signal

{top_sv_text}

## Top Recurrently Linked Genes

{top_gene_text}

## Output Files

- `sv326_direct_ccre_gene_link_long.tsv`: detailed SV-cCRE-gene evidence rows.
- `sv326_direct_ccre_target_gene_summary.tsv`: one row per SV-cCRE-gene.
- `sv326_direct_sv_target_gene_summary.tsv`: one row per SV with linked target genes.
- `sv326_target_gene_summary.tsv`: one row per linked gene.
- `sv326_target_gene_link_type_summary.tsv`: summary by 3D-Chromatin, CRISPR, and eQTL evidence.
- `sv326_gene_link_scan_diagnostics.tsv`: input scan diagnostics.

## Notes

- Gene links are from ENCODE4/S SCREEN files in GRCh38 accession space.
- The join is by stable EH38E/EH38D cCRE accessions, not by coordinates.
- These links suggest candidate regulatory target genes; they should be integrated with HUVEC RNA, ATAC, and Hi-C evidence before final biological interpretation.
"""
    (outdir / "sv326_ccre_target_gene_report.md").write_text(report, encoding="utf-8")
    diagnostics.to_csv(outdir / "sv326_gene_link_scan_diagnostics.tsv", sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ccre-hits", type=Path, default=DEFAULT_CCRE_HITS)
    parser.add_argument("--gene-link-dir", type=Path, default=DEFAULT_GENE_LINK_DIR)
    parser.add_argument("--gene-link-files", type=Path, nargs="*", default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    direct_hits = load_direct_hits(args.ccre_hits)
    targets = target_accessions_from_hits(direct_hits)
    if not targets:
        raise ValueError("No EH38E/EH38D accessions found in direct cCRE hits.")
    print(f"[INFO] direct hit rows={len(direct_hits)}; target accessions={len(targets)}", flush=True)
    links, diagnostics = load_gene_links(args.gene_link_dir, args.gene_link_files, targets)
    print(f"[INFO] gene-link records loaded={len(links)}", flush=True)
    link_hits = join_links_to_hits(direct_hits, links)
    print(f"[INFO] joined SV-cCRE-gene rows={len(link_hits)}", flush=True)
    make_outputs(direct_hits, link_hits, args.outdir)
    write_report(args.outdir, direct_hits, links, link_hits, diagnostics)
    print(f"Wrote {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
