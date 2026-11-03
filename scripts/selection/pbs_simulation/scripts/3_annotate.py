#!/usr/bin/env python
"""Annotate each observed site with its own POINT recombination and mutation rate.

R_point  : cM/Mb of the gmap segment containing POS
MU_point : MU_PER_BP_PER_GEN of the 25-kb mutation-map bin containing POS

LOG_R  = log10(R_point + 1e-4)
LOG_MU = log10(MU_point + 1e-12)

Outputs
-------
<out>         LOG_R <tab> LOG_MU <tab> PBS     (one row per input site)
<out>.index   CHROM <tab> POS                  (same order)

Example
-------
python 03_annotate_point_rates.py --sites chr22.sites.pbs.tsv.gz --chrom chr22 \
    --recomb-map maps/chr22_noMask.scaled.gmap \
    --mutation-map maps/mutation_rate_25kb.tsv.gz \
    --out queries/chr22.queries.tsv
"""
import argparse
import gzip

import numpy as np

ER, EM = 1e-4, 1e-12


def load_gmap(path):
    pos, rate = [], []
    with open(path) as f:
        f.readline()
        for line in f:
            a = line.split()
            pos.append(int(a[0]))
            rate.append(float(a[1]))
    return np.array(pos, dtype=np.int64), np.array(rate)


def load_mutation_map(path, chrom):
    start, rate = [], []
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            a = line.rstrip("\n").split("\t")
            if a[0] != chrom:
                continue
            start.append(int(a[1]))
            rate.append(float(a[9]))
    return np.array(start, dtype=np.int64), np.array(rate)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True, help="POS<TAB>PBS (gzip or plain)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--recomb-map", required=True)
    ap.add_argument("--mutation-map", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    gp, gr = load_gmap(args.recomb_map)
    ms, mu = load_mutation_map(args.mutation_map, args.chrom)

    opener = gzip.open if args.sites.endswith(".gz") else open
    pos, pbs = [], []
    with opener(args.sites, "rt") as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            pos.append(int(a[0]))
            pbs.append(float(a[1]))
    pos = np.array(pos, dtype=np.int64)
    pbs = np.array(pbs)

    i = np.clip(np.searchsorted(gp, pos, side="right") - 1, 0, len(gp) - 1)
    j = np.clip(np.searchsorted(ms, pos, side="right") - 1, 0, len(ms) - 1)
    log_r = np.log10(gr[i] + ER)
    log_mu = np.log10(mu[j] + EM)

    with open(args.out, "w") as fo, open(args.out + ".index", "w") as fx:
        for k in range(len(pos)):
            fo.write(f"{log_r[k]:.15g}\t{log_mu[k]:.15g}\t{pbs[k]:.15g}\n")
            fx.write(f"{args.chrom}\t{pos[k]}\n")
    print(f"{args.chrom}: {len(pos)} sites -> {args.out}")


if __name__ == "__main__":
    main()
