#!/usr/bin/env python
"""Per-site Weir-Cockerham FST and PBS for one chromosome.

scikit-allel's weir_cockerham_fst returns the (a, b, c) variance components as a
tuple of arrays; the per-variant FST is a / (a + b + c).  Applying np.ravel() to
the tuple and treating the result as per-site FST silently mis-assigns values
and must not be done.

Output (gzip, no header): POS <tab> PBS

Example
-------
python 02_site_pbs.py --vcf input.vcf.gz --chrom chr22 \
    --tib TIB.samples --han HAN.samples --eur EUR.samples \
    --out chr22.sites.pbs.tsv.gz
"""
import argparse
import gzip
import subprocess
import tempfile
import os
import shutil

import numpy as np
import allel


def clean_halfcalls(vcf_in, vcf_out):
    """Set any genotype containing a missing allele to ./. (keeps FORMAT fields)."""
    with gzip.open(vcf_in, "rt") as f, gzip.open(vcf_out, "wt") as fo:
        for line in f:
            if line.startswith("#"):
                fo.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            out = []
            for g in cols[9:]:
                gt = g.split(":")[0].replace("|", "/").split("/")
                rest = ":".join(g.split(":")[1:])
                out.append(("./." + (":" + rest if rest else ""))
                           if any(x == "." for x in gt) else g)
            fo.write("\t".join(cols[:9] + out) + "\n")


def read_samples(vcf_path):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError("no #CHROM header")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--tib", required=True, help="one sample ID per line")
    ap.add_argument("--han", required=True)
    ap.add_argument("--eur", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-clean", action="store_true",
                    help="skip half-call cleaning (e.g. simulation output)")
    args = ap.parse_args()

    tmpdir = tempfile.mkdtemp(prefix="site_pbs_", dir=".")
    try:
        clean = os.path.join(tmpdir, "clean.vcf.gz")
        r = subprocess.run(["bcftools", "view", "-r", args.chrom, "-Oz",
                            "-o", clean, args.vcf], capture_output=True, text=True)
        if r.returncode != 0:  # not bgzipped/indexed
            subprocess.run(["bcftools", "view", "-Oz", "-o", clean, args.vcf], check=True)
        if not args.no_clean:
            clean_halfcalls(clean, clean + ".clean.gz")
            os.replace(clean + ".clean.gz", clean)

        samples = read_samples(clean)
        idx = {s: i for i, s in enumerate(samples)}
        pops = []
        for path in (args.tib, args.han, args.eur):
            with open(path) as f:
                pops.append([idx[s.strip()] for s in f if s.strip() in idx])

        cs = allel.read_vcf(clean, fields=["variants/POS", "calldata/GT"])
        pos_arr = cs["variants/POS"]
        gt = allel.GenotypeArray(cs["calldata/GT"])

        # Pairwise FST. NOTE: the function returns (a, b, c) components.
        fst = {}
        for name, sub in (("TIB_HAN", [pops[0], pops[1]]),
                          ("TIB_EUR", [pops[0], pops[2]]),
                          ("HAN_EUR", [pops[1], pops[2]])):
            a, b, c = allel.weir_cockerham_fst(gt, subpops=sub)
            den = a + b + c
            with np.errstate(invalid="ignore", divide="ignore"):
                f = np.where(den > 0, a / np.where(den == 0, 1.0, den), 0.0)
            f = f[:, 0]                      # biallelic: allele columns agree
            fst[name] = np.where(np.isfinite(f) & (f > 0), f, 0.0)

        t_th = -np.log1p(-np.clip(fst["TIB_HAN"], 0, 0.999999))
        t_te = -np.log1p(-np.clip(fst["TIB_EUR"], 0, 0.999999))
        t_he = -np.log1p(-np.clip(fst["HAN_EUR"], 0, 0.999999))
        pbs = (t_th + t_te - t_he) / 2.0

        n = 0
        with gzip.open(args.out, "wt") as fo:
            for k, p in enumerate(pos_arr):
                if not np.isfinite(pbs[k]):
                    continue
                fo.write(f"{p}\t{pbs[k]:.10g}\n")
                n += 1
        print(f"{args.chrom}: {n} sites -> {args.out}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
