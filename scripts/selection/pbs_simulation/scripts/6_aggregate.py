#!/usr/bin/env python
"""Aggregate shard query results into empirical P values and BH-adjusted values.

Inputs
------
--queries   LOG_R <tab> LOG_MU <tab> PBS        (same order as queried)
--index     CHROM <tab> POS                     (same order)
--shards    glob or list of per-shard outputs   (TOL_USED TOL N_MATCH N_GE)

Output
------
CHROM POS PBS LOG_R LOG_MU N_MATCH_TOTAL N_GE_TOTAL P_SIM_PLUS1 P_BH_PBSpos MATCH_VALID

Statistical definition
----------------------
N_MATCH_TOTAL = sum over shards
N_GE_TOTAL    = sum over shards
P_SIM_PLUS1   = (N_GE + 1) / (N_MATCH + 1)
MATCH_VALID   = N_MATCH_TOTAL >= --min-match
BH            = Benjamini-Hochberg, applied only to PBS > 0 and valid sites
                and restricted to the test family given by --bh-mode
                  separate : each row of --family is corrected on its own
                  pooled   : all rows of --family share one correction

Example
-------
python 06_aggregate_bh.py --queries queries.tsv --index queries.tsv.index \
    --shards 'results/shard_*.out' --out SNP.sites.p_results.BH.tsv
"""
import argparse
import glob

import numpy as np


def bh_adjust(p, mask):
    """Benjamini-Hochberg q-values; NaN outside the tested set."""
    q = np.full(len(p), np.nan)
    idx = np.where(mask)[0]
    n = idx.size
    if n:
        order = idx[np.argsort(p[idx], kind="mergesort")]
        r = p[order] * n / np.arange(1, n + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        q[order] = np.minimum(r, 1.0)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--shards", required=True,
                    help="glob pattern or comma-separated list of shard outputs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-match", type=int, default=10000)
    ap.add_argument("--family", default=None,
                    help="optional file with one dataset label per query; "
                         "default: all queries form a single family")
    ap.add_argument("--bh-mode", choices=["separate", "pooled"], default="separate")
    args = ap.parse_args()

    q = np.loadtxt(args.queries)
    log_r, log_mu, pbs = q[:, 0], q[:, 1], q[:, 2]
    nq = len(q)

    nm = np.zeros(nq, dtype=np.int64)
    ng = np.zeros(nq, dtype=np.int64)
    paths = (sorted(glob.glob(args.shards)) if any(ch in args.shards for ch in "*?[")
             else args.shards.split(","))
    for path in paths:
        d = np.loadtxt(path)
        assert len(d) == nq, f"{path}: {len(d)} rows != {nq} queries"
        nm += d[:, 2].astype(np.int64)
        ng += d[:, 3].astype(np.int64)

    valid = nm >= args.min_match
    p_sim = np.where(valid, (ng + 1.0) / (nm + 1.0), np.nan)

    if args.family:
        labels = np.array([l.strip() for l in open(args.family)])
        assert len(labels) == nq
        groups = [(lab, labels == lab) for lab in dict.fromkeys(labels)]
    else:
        groups = [("all", np.ones(nq, dtype=bool))]

    p_bh = np.full(nq, np.nan)
    if args.bh_mode == "pooled":
        mask = valid & (pbs > 0)
        p_bh = bh_adjust(p_sim, mask)
    else:
        for _, m in groups:
            sel = np.where(m)[0]
            p_bh[sel] = bh_adjust(p_sim[sel], valid[sel] & (pbs[sel] > 0))

    idx = [l.rstrip("\n").split("\t") for l in open(args.index)]
    with open(args.out, "w") as f:
        f.write("CHROM\tPOS\tPBS\tLOG_R\tLOG_MU\tN_MATCH_TOTAL\tN_GE_TOTAL\t"
                "P_SIM_PLUS1\tP_BH_PBSpos\tMATCH_VALID\n")
        for i in range(nq):
            f.write(f"{idx[i][0]}\t{idx[i][1]}\t{pbs[i]:.8g}\t{log_r[i]:.8g}\t"
                    f"{log_mu[i]:.8g}\t{nm[i]}\t{ng[i]}\t"
                    + ("%.6g" % p_sim[i] if valid[i] else "NA") + "\t"
                    + ("%.6g" % p_bh[i] if np.isfinite(p_bh[i]) else "NA") + "\t"
                    + ("1" if valid[i] else "0") + "\n")

    n_valid = int(valid.sum())
    n_tested = int((valid & (pbs > 0)).sum())
    print(f"queries={nq:,} valid={n_valid:,} tested(PBS>0)={n_tested:,} "
          f"P_BH<0.05={int(np.nansum(p_bh < 0.05)):,} -> {args.out}")


if __name__ == "__main__":
    main()
