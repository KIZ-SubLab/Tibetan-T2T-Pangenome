#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Union gene-set enrichment with stratified Model B permutations.

Each permutation:
  1. samples --n-sv genes (default 167) without replacement from the SV +/- 5 kb
     matched background;
  2. excludes those genes from the genome-wide background and samples --n-wg
     additional genes (default 130);
  3. combines both parts into a random union of the same size as the observed set;
  4. computes overlaps with all GWAS phenotypes and with the hypoxia gene set
     from the same random gene set.

Statistics:
  expected overlap = permutation mean
  fold enrichment  = observed / expected
  empirical P      = (1 + number of null overlaps >= observed) / (N_PERM + 1)
  GWAS FDR         = Benjamini-Hochberg over all phenotypes tested
  Hypoxia          = single prespecified test, no multiple-testing correction

Usage:
  python3 05_union_modelB.py \
      --union union297.txt --background-b SV5kb_background.txt \
      --genome gencode_v38_protein_coding_genes.txt --pairs pairs.tsv \
      --hypoxia hypoxia_harmonized.txt \
      --n-perm 100000 --seed 42 --outdir out/union_modelB
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
        description="Union gene-set enrichment with stratified Model B permutations")
    ap.add_argument("--union", required=True, help="Observed union gene list")
    ap.add_argument("--background-b", required=True,
                    help="SV +/- 5 kb matched background gene list")
    ap.add_argument("--genome", required=True,
                    help="Genome-wide protein-coding background gene list")
    ap.add_argument("--pairs", required=True,
                    help="Phenotype-gene pairs TSV")
    ap.add_argument("--hypoxia", required=True, help="Hypoxia gene list")
    ap.add_argument("--n-perm", type=int, default=100000,
                    help="Number of permutations (default: 100000)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--n-sv", type=int, default=167,
                    help="Genes sampled from the SV background per permutation")
    ap.add_argument("--n-wg", type=int, default=130,
                    help="Genes sampled from the remaining genome-wide background "
                         "per permutation")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    union = read_list(args.union)
    sv_background = read_list(args.background_b)
    genome = read_list(args.genome)
    hypoxia = sorted(set(read_list(args.hypoxia)) & set(genome))
    sv_set, genome_set = set(sv_background), set(genome)

    assert sv_set <= genome_set, "SV background must be a subset of the genome background"
    assert set(union) <= genome_set, "Union genes must belong to the genome background"
    assert args.n_sv + args.n_wg == len(union), (
        f"--n-sv + --n-wg ({args.n_sv} + {args.n_wg}) must equal the union size "
        f"({len(union)})")
    print(f"Union: {len(union)} | SV background: {len(sv_background):,} | "
          f"genome background: {len(genome):,} | Model B: n_sv={args.n_sv}, "
          f"n_wg={args.n_wg} | permutations: {args.n_perm:,} | seed: {args.seed}")

    index = {gene: i for i, gene in enumerate(genome)}
    genome_index = np.arange(len(genome), dtype=np.int64)
    sv_index = np.array([index[g] for g in sv_background], dtype=np.int64)

    # GWAS phenotype matrix (phenotypes x genome background)
    pairs = pd.read_csv(args.pairs, sep="\t", dtype=str, keep_default_na=False)
    phenotypes = sorted(pairs["TRAIT_URI"].unique())
    phenotype_index = {p: i for i, p in enumerate(phenotypes)}
    rows = pairs["TRAIT_URI"].map(phenotype_index).values
    cols = pairs["GENE"].map(lambda g: index.get(g, -1)).values
    valid = cols >= 0
    matrix = csr_matrix(
        (np.ones(valid.sum(), dtype=np.int32), (rows[valid], cols[valid])),
        shape=(len(phenotypes), len(genome)))
    n_gwas = np.asarray(matrix.sum(axis=1)).ravel()
    labels = pairs.drop_duplicates("TRAIT_URI").set_index(
        "TRAIT_URI")["TRAIT_LABEL"]

    union_vector = np.zeros(len(genome), dtype=np.int32)
    union_vector[[index[g] for g in union]] = 1
    observed = np.asarray(matrix @ union_vector).ravel()

    hypoxia_vector = np.zeros(len(genome), dtype=np.int32)
    hypoxia_vector[[index[g] for g in hypoxia]] = 1
    observed_hypoxia = int((union_vector * hypoxia_vector).sum())
    overlap_genes = sorted(set(union) & set(hypoxia))
    print(f"Phenotypes: {len(phenotypes):,} | hypoxia genes in background: "
          f"{len(hypoxia)} | observed hypoxia overlap: {observed_hypoxia}")

    print("Running stratified permutations")
    rng = np.random.default_rng(args.seed)
    null = np.zeros((len(phenotypes), args.n_perm), dtype=np.uint16)
    hypoxia_counts, hypoxia_n_ge = {}, 0
    for j in range(args.n_perm):
        # 1) sample from the SV +/- 5 kb background
        sampled_sv = sv_index[rng.choice(len(sv_index), size=args.n_sv,
                                         replace=False)]
        # 2) exclude them from the genome background and sample the remainder
        mask = np.ones(len(genome_index), dtype=bool)
        mask[sampled_sv] = False
        sampled_wg = rng.choice(genome_index[mask], size=args.n_wg,
                                replace=False)
        # 3) combine into unique random union genes
        random_union = np.concatenate([sampled_sv, sampled_wg])
        assert len(np.unique(random_union)) == len(union)
        # 4) same random set for GWAS and hypoxia overlaps
        vector = np.zeros(len(genome), dtype=np.int32)
        vector[random_union] = 1
        null[:, j] = np.asarray(matrix @ vector).ravel()
        overlap = int(hypoxia_vector[random_union].sum())
        hypoxia_counts[overlap] = hypoxia_counts.get(overlap, 0) + 1
        hypoxia_n_ge += int(overlap >= observed_hypoxia)
        if (j + 1) % 20000 == 0:
            print(f"  {j + 1:,}/{args.n_perm:,} permutations completed")

    # GWAS summary
    expected = null.mean(axis=1)
    n_ge = (null >= observed[:, None]).sum(axis=1)
    empirical_p = (1 + n_ge) / (args.n_perm + 1)
    fdr = benjamini_hochberg(empirical_p)
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = np.where(expected > 0, observed / expected, np.nan)
    gwas = pd.DataFrame({
        "TRAIT_URI": phenotypes,
        "TRAIT_LABEL": [labels.get(p, "") for p in phenotypes],
        "N_GWAS_GENES": n_gwas.astype(int),
        "UNION_SIZE": len(union),
        "OBSERVED_OVERLAP": observed.astype(int),
        "EXPECTED_OVERLAP": np.round(expected, 4),
        "FOLD_ENRICHMENT": np.round(fold, 4),
        "EMPIRICAL_P": empirical_p,
        "FDR": np.round(fdr, 6),
        "SIGNIFICANT": ((fdr < 0.05) & (fold > 1)).astype(int),
    }).sort_values(["FDR", "FOLD_ENRICHMENT"], ascending=[True, False])
    gwas.to_csv(os.path.join(args.outdir, "GWAS_enrichment_full.tsv"),
                sep="\t", index=False, na_rep="NA")
    gwas[gwas["SIGNIFICANT"] == 1].to_csv(
        os.path.join(args.outdir, "GWAS_significant_phenotypes.tsv"),
        sep="\t", index=False, na_rep="NA")

    # Hypoxia summary (single prespecified test)
    expected_hypoxia = sum(a * b for a, b in hypoxia_counts.items()) / args.n_perm
    hypoxia_p = (1 + hypoxia_n_ge) / (args.n_perm + 1)
    pd.DataFrame([{
        "MODEL": f"Model_B_WG{args.n_wg}_SV{args.n_sv}",
        "SAMPLING": f"SV{args.n_sv}+WG{args.n_wg}",
        "UNION_SIZE": len(union),
        "HYPOXIA_IN_GENOME": len(hypoxia),
        "N_PERM": args.n_perm,
        "SEED": args.seed,
        "OBSERVED_OVERLAP": observed_hypoxia,
        "EXPECTED_OVERLAP": round(expected_hypoxia, 4),
        "FOLD_ENRICHMENT": round(observed_hypoxia / expected_hypoxia, 4)
        if expected_hypoxia > 0 else np.nan,
        "EMPIRICAL_P": hypoxia_p,
        "OVERLAP_GENES": ",".join(overlap_genes),
    }]).to_csv(os.path.join(args.outdir, "Hypoxia_enrichment_statistics.tsv"),
               sep="\t", index=False)
    pd.DataFrame({"OVERLAP_COUNT": sorted(hypoxia_counts),
                  "N_PERMUTATIONS": [hypoxia_counts[x]
                                     for x in sorted(hypoxia_counts)]}
                 ).to_csv(os.path.join(
                     args.outdir, "Hypoxia_permutation_null_distribution.tsv"),
                 sep="\t", index=False)
    with open(os.path.join(args.outdir, "Hypoxia_overlap_genes.txt"), "w") as handle:
        for gene in overlap_genes:
            handle.write(gene + "\n")

    print(f"GWAS significant phenotypes: {int(gwas['SIGNIFICANT'].sum())} | "
          f"hypoxia observed: {observed_hypoxia} | expected: "
          f"{expected_hypoxia:.4f} | fold: "
          f"{observed_hypoxia / expected_hypoxia:.4f} | P: {hypoxia_p:.5e}")
    print("Analysis completed")


if __name__ == "__main__":
    main()
