#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare GWAS phenotype-gene associations.

Reads the GWAS Catalog "All associations" TSV and produces a long-format
phenotype-gene table for enrichment analysis:

  - retain associations with INTERGENIC = 0 and a non-empty MAPPED_GENE
  - GWAS significance: PVALUE_MLOG >= 7.30102999566 (P <= 5e-8); if PVALUE_MLOG
    is missing, fall back to P-VALUE <= 5e-8
  - expand MAPPED_TRAIT_URI into one row per ontology term
  - split MAPPED_GENE on commas
  - match gene symbols exactly to the GENCODE v38 protein-coding universe
  - deduplicate on (TRAIT_URI, TRAIT_LABEL, GENE)
  - require at least --n-min distinct protein-coding genes per phenotype

If preprocessed phenotype-gene pairs are available, supply them directly to
04_enrichment.py / 05_union_modelB.py instead of running this script.

Usage:
  python3 01_prepare_gwas.py --gwas GWAS.tsv --uri-labels labels.tsv \
      --background gencode_v38_protein_coding_genes.txt --out pairs.tsv
"""

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Build phenotype-gene pairs from a GWAS Catalog associations TSV")
    ap.add_argument("--gwas", required=True,
                    help="GWAS Catalog All associations TSV")
    ap.add_argument("--uri-labels", required=True,
                    help="TSV with columns TRAIT_URI and TRAIT_LABEL")
    ap.add_argument("--background", required=True,
                    help="GENCODE v38 protein-coding gene symbols (one per line)")
    ap.add_argument("--out", required=True,
                    help="Output phenotype-gene pairs TSV")
    ap.add_argument("--n-min", type=int, default=20,
                    help="Minimum distinct genes per phenotype (default: 20)")
    ap.add_argument("--mlog-cutoff", type=float, default=7.30102999566,
                    help="Minimum PVALUE_MLOG (default: -log10(5e-8))")
    ap.add_argument("--p-cutoff", type=float, default=5e-8,
                    help="Fallback P-value cutoff when PVALUE_MLOG is missing")
    args = ap.parse_args()

    pc_genes = set(pd.read_csv(args.background, sep="\t", header=None,
                               dtype=str)[0])
    labels = pd.read_csv(args.uri_labels, sep="\t", dtype=str,
                         keep_default_na=False)
    uri_to_label = dict(zip(labels["TRAIT_URI"], labels["TRAIT_LABEL"]))

    print("Preparing phenotype-gene associations")
    df = pd.read_csv(args.gwas, sep="\t", dtype=str, keep_default_na=False)
    n_records = len(df)

    # Intragenic associations with a mapped gene
    df = df[(df["INTERGENIC"] == "0") & (df["MAPPED_GENE"].str.strip() != "")]

    # GWAS significance threshold
    mlog = pd.to_numeric(df["PVALUE_MLOG"], errors="coerce")
    pval = pd.to_numeric(df["P-VALUE"], errors="coerce")
    df = df[mlog.ge(args.mlog_cutoff) | (mlog.isna() & pval.le(args.p_cutoff))]

    # One row per mapped gene
    df = df.copy()
    df["_GENES"] = df["MAPPED_GENE"].map(
        lambda v: [x.strip() for x in str(v).split(",") if x.strip()])
    df = df.explode("_GENES").dropna(subset=["_GENES"])
    df = df.rename(columns={"_GENES": "GENE"})

    # One row per ontology trait term
    df = df.copy()
    df["_URIS"] = df["MAPPED_TRAIT_URI"].map(
        lambda v: [x.strip() for x in str(v).split(",") if x.strip()])
    df = df.explode("_URIS").dropna(subset=["_URIS"])
    df = df.rename(columns={"_URIS": "TRAIT_URI"})
    df["TRAIT_LABEL"] = df["TRAIT_URI"].map(uri_to_label).fillna("NA")

    # Restrict to the protein-coding universe
    df = df[df["GENE"].isin(pc_genes)]

    # Deduplicate phenotype-gene pairs
    pairs = df[["TRAIT_URI", "TRAIT_LABEL", "GENE"]].drop_duplicates()

    # Minimum number of genes per phenotype
    n_genes = pairs.groupby("TRAIT_URI")["GENE"].nunique()
    keep = set(n_genes[n_genes >= args.n_min].index)
    pairs = pairs[pairs["TRAIT_URI"].isin(keep)].sort_values(
        ["TRAIT_URI", "GENE"]).reset_index(drop=True)

    pairs.to_csv(args.out, sep="\t", index=False)
    print(f"GWAS records: {n_records:,} -> {len(pairs):,} phenotype-gene pairs, "
          f"{pairs['TRAIT_URI'].nunique():,} phenotypes "
          f"(>= {args.n_min} genes)")
    print(f"Output: {args.out}")


if __name__ == "__main__":
    main()
