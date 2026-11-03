#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run gene-set enrichment by background permutation.

The target gene set is compared with either GWAS phenotype gene sets (--mode gwas)
or a hypoxia gene set (--mode hypoxia). Random gene sets of the same size as the
target are drawn without replacement from the supplied background.

Statistics:
  expected overlap  = permutation mean
  fold enrichment   = observed / expected
  empirical P       = (1 + number of null overlaps >= observed) / (N_PERM + 1)
  FDR               = Benjamini-Hochberg over all phenotypes tested (GWAS mode)

In GWAS mode with --restrict-to-background yes, phenotype genes are first
restricted to the background and phenotypes are then re-filtered at --n-min
genes; this is required when the background is a matched subset of the genome
(e.g. SV +/- 5 kb genes). Hypoxia mode performs a single test and applies no
multiple-testing correction.

Usage:
  python3 04_enrichment.py --mode gwas --target target.txt \
      --background background.txt --pairs pairs.tsv \
      --n-perm 100000 --seed 42 --outdir out/gwas
  python3 04_enrichment.py --mode hypoxia --target target.txt \
      --background background.txt --hypoxia hypoxia.txt \
      --n-perm 100000 --seed 42 --outdir out/hypoxia
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg adjusted P values."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p, kind="mergesort")
    adjusted = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def read_list(path):
    return sorted(set(pd.read_csv(path, sep="\t", header=None, dtype=str)[0]))


def main():
    ap = argparse.ArgumentParser(
        description="Background permutation enrichment for a target gene set")
    ap.add_argument("--mode", choices=["gwas", "hypoxia"], required=True,
                    help="Enrichment mode")
    ap.add_argument("--target", required=True, help="Target gene list")
    ap.add_argument("--background", required=True,
                    help="Permutation background gene list")
    ap.add_argument("--pairs", help="Phenotype-gene pairs TSV (GWAS mode)")
    ap.add_argument("--hypoxia", help="Hypoxia gene list (hypoxia mode)")
    ap.add_argument("--n-perm", type=int, default=100000,
                    help="Number of permutations (default: 100000)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--restrict-to-background", choices=["yes", "no"],
                    default="no",
                    help="GWAS mode: restrict phenotype genes to the background "
                         "and re-apply the minimum gene count")
    ap.add_argument("--n-min", type=int, default=20,
                    help="Minimum genes per phenotype when restricting")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    target = read_list(args.target)
    background = read_list(args.background)
    background_set = set(background)
    assert set(target) <= background_set, \
        "Target genes must be a subset of the background"
    n_target, n_background = len(target), len(background)
    index = {gene: i for i, gene in enumerate(background)}
    rng = np.random.default_rng(args.seed)
    print(f"Mode: {args.mode} | target: {n_target} | background: "
          f"{n_background:,} | permutations: {args.n_perm:,} | seed: {args.seed}")

    if args.mode == "gwas":
        assert args.pairs, "--pairs is required in GWAS mode"
        pairs = pd.read_csv(args.pairs, sep="\t", dtype=str,
                            keep_default_na=False)
        if args.restrict_to_background == "yes":
            pairs = pairs[pairs["GENE"].isin(background_set)]
            gene_counts = pairs.groupby("TRAIT_URI")["GENE"].nunique()
            keep = set(gene_counts[gene_counts >= args.n_min].index)
            pairs = pairs[pairs["TRAIT_URI"].isin(keep)].copy()
            print(f"Phenotypes after background restriction "
                  f"(>= {args.n_min} genes): {pairs['TRAIT_URI'].nunique():,}")

        phenotypes = sorted(pairs["TRAIT_URI"].unique())
        phenotype_index = {p: i for i, p in enumerate(phenotypes)}
        rows = pairs["TRAIT_URI"].map(phenotype_index).values
        cols = pairs["GENE"].map(lambda g: index.get(g, -1)).values
        valid = cols >= 0
        matrix = csr_matrix(
            (np.ones(valid.sum(), dtype=np.int32), (rows[valid], cols[valid])),
            shape=(len(phenotypes), n_background))
        n_gwas = np.asarray(matrix.sum(axis=1)).ravel()
        labels = pairs.drop_duplicates("TRAIT_URI").set_index(
            "TRAIT_URI")["TRAIT_LABEL"]
        target_vector = np.zeros(n_background, dtype=np.int32)
        target_vector[[index[g] for g in target]] = 1
        observed = np.asarray(matrix @ target_vector).ravel()

        null = np.zeros((len(phenotypes), args.n_perm), dtype=np.uint16)
        for j in range(args.n_perm):
            sampled = rng.choice(n_background, size=n_target, replace=False)
            vector = np.zeros(n_background, dtype=np.int32)
            vector[sampled] = 1
            null[:, j] = np.asarray(matrix @ vector).ravel()

        expected = null.mean(axis=1)
        n_ge = (null >= observed[:, None]).sum(axis=1)
        empirical_p = (1 + n_ge) / (args.n_perm + 1)
        fdr = benjamini_hochberg(empirical_p)
        with np.errstate(divide="ignore", invalid="ignore"):
            fold = np.where(expected > 0, observed / expected, np.nan)
        output = pd.DataFrame({
            "TRAIT_URI": phenotypes,
            "TRAIT_LABEL": [labels.get(p, "") for p in phenotypes],
            "N_GWAS_GENES": n_gwas.astype(int),
            "N_TARGET": n_target,
            "OBSERVED_OVERLAP": observed.astype(int),
            "EXPECTED_OVERLAP": np.round(expected, 4),
            "FOLD_ENRICHMENT": np.round(fold, 4),
            "EMPIRICAL_P": empirical_p,
            "FDR": np.round(fdr, 6),
            "SIGNIFICANT": ((fdr < 0.05) & (fold > 1)).astype(int),
        }).sort_values(["FDR", "FOLD_ENRICHMENT"], ascending=[True, False])
        output.to_csv(os.path.join(args.outdir, "GWAS_enrichment_full.tsv"),
                      sep="\t", index=False, na_rep="NA")
        output[output["SIGNIFICANT"] == 1].to_csv(
            os.path.join(args.outdir, "GWAS_significant_phenotypes.tsv"),
            sep="\t", index=False, na_rep="NA")
        print(f"Phenotypes tested: {len(output):,} | significant: "
              f"{int(output['SIGNIFICANT'].sum())}")

    else:
        assert args.hypoxia, "--hypoxia is required in hypoxia mode"
        hypoxia = sorted(set(pd.read_csv(args.hypoxia, sep="\t", header=None,
                                         dtype=str)[0]) & background_set)
        hypoxia_vector = np.zeros(n_background, dtype=np.int32)
        hypoxia_vector[[index[g] for g in hypoxia]] = 1
        target_vector = np.zeros(n_background, dtype=np.int32)
        target_vector[[index[g] for g in target]] = 1
        observed = int((target_vector * hypoxia_vector).sum())

        counts, n_ge = {}, 0
        for _ in range(args.n_perm):
            sampled = rng.choice(n_background, size=n_target, replace=False)
            overlap = int(hypoxia_vector[sampled].sum())
            counts[overlap] = counts.get(overlap, 0) + 1
            n_ge += int(overlap >= observed)
        expected = sum(a * b for a, b in counts.items()) / args.n_perm
        empirical_p = (1 + n_ge) / (args.n_perm + 1)

        pd.DataFrame([{
            "MODE": "hypoxia",
            "N_TARGET": n_target,
            "N_BACKGROUND": n_background,
            "N_HYPOXIA_IN_BACKGROUND": len(hypoxia),
            "OBSERVED_OVERLAP": observed,
            "EXPECTED_OVERLAP": round(expected, 4),
            "FOLD_ENRICHMENT": round(observed / expected, 4) if expected > 0 else np.nan,
            "EMPIRICAL_P": empirical_p,
            "N_PERM": args.n_perm,
            "SEED": args.seed,
            "OVERLAP_GENES": ",".join(sorted(set(target) & set(hypoxia))),
        }]).to_csv(os.path.join(args.outdir, "Hypoxia_enrichment_statistics.tsv"),
                   sep="\t", index=False)
        pd.DataFrame({"OVERLAP_COUNT": sorted(counts),
                      "N_PERMUTATIONS": [counts[x] for x in sorted(counts)]}
                     ).to_csv(os.path.join(
                         args.outdir, "Hypoxia_permutation_null_distribution.tsv"),
                     sep="\t", index=False)
        print(f"Observed overlap: {observed} | expected: {expected:.4f} | "
              f"fold: {observed / expected:.4f} | P: {empirical_p:.5e}")


if __name__ == "__main__":
    main()
