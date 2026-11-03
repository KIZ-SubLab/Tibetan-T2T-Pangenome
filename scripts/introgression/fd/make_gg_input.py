#!/usr/bin/env python3
"""Create deterministic haploid-count genomics_general input for B/C checks."""

import argparse
import gzip
import os

import numpy as np
import pandas as pd


def alleles(ref, alt, freq, n_called, n_max):
    n_called = int(n_called)
    dosage_float = float(freq) * n_called
    dosage = int(round(dosage_float))
    if not np.isclose(dosage_float, dosage, atol=1e-6):
        raise ValueError(f"Non-integer reconstructed dosage: {dosage_float}")
    return [alt] * dosage + [ref] * (n_called-dosage) + ["N"] * (n_max-n_called)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sites", required=True)
    p.add_argument("--geno", required=True)
    p.add_argument("--pops", required=True)
    args = p.parse_args()
    d = pd.read_csv(args.sites, sep="\t")
    maxima = {g: int(d[f"n_{g.lower()}_chr"].max()) for g in ["YRI", "TBN", "DEN", "NEA"]}
    maxima["ANC"] = 1
    names = {g: [f"{g}_h{i+1:03d}" for i in range(n)] for g, n in maxima.items()}
    order = ["YRI", "TBN", "DEN", "NEA", "ANC"]
    ordered_names = sum((names[g] for g in order), [])
    os.makedirs(os.path.dirname(args.geno), exist_ok=True)
    with open(args.pops, "w") as out:
        for group in order:
            for sample in names[group]:
                out.write(f"{sample}\t{group}\n")
    with gzip.open(args.geno, "wt") as out:
        out.write("#CHROM\tPOS\t" + "\t".join(ordered_names) + "\n")
        for r in d.itertuples(index=False):
            row = []
            for group in ["yri", "tbn", "den", "nea"]:
                row += alleles(r.ref, r.alt, getattr(r, f"p_{group}"),
                               getattr(r, f"n_{group}_chr"), maxima[group.upper()])
            row += [r.ancestral]
            out.write(f"{r.chrom}\t{r.pos}\t" + "\t".join(row) + "\n")


if __name__ == "__main__":
    main()
