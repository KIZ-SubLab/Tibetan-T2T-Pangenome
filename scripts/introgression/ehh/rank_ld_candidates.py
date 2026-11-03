#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
from collections import defaultdict


def read_meta(path):
    with open(path, newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


ap = argparse.ArgumentParser()
ap.add_argument("--ld", required=True)
ap.add_argument("--ids", required=True)
ap.add_argument("--snp-meta", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

ids = [x.strip() for x in open(a.ids) if x.strip() and not x.startswith("#")]
idset = set(ids)
snp = read_meta(a.snp_meta)
candidates = defaultdict(dict)

with gzip.open(a.ld, "rt") as fh:
    header = fh.readline().split()
    col = {x: i for i, x in enumerate(header)}
    for line in fh:
        f = line.split()
        va, vb = f[col["SNP_A"]], f[col["SNP_B"]]
        if va in idset and vb not in idset:
            sid, cand = va, vb
            svpos, snppos = int(f[col["BP_A"]]), int(f[col["BP_B"]])
        elif vb in idset and va not in idset:
            sid, cand = vb, va
            svpos, snppos = int(f[col["BP_B"]]), int(f[col["BP_A"]])
        else:
            continue
        if cand not in snp:
            continue
        try:
            value = float(f[col["R2"]])
        except ValueError:
            continue
        if math.isfinite(value):
            candidates[sid][cand] = (value, abs(snppos - svpos), snppos)

fields = ["sv_id", "rank", "lead_snp_id", "chrom_chm13", "pos_chm13", "ref_chm13", "alt_chm13", "ac", "an", "af", "r2", "distance_bp"]
with open(a.out, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for sid in ids:
        ranked = sorted(candidates[sid].items(), key=lambda x: (-x[1][0], x[1][1], x[1][2], x[0]))
        for rank, (cand, (value, dist, _)) in enumerate(ranked, 1):
            m = snp[cand]
            w.writerow({"sv_id": sid, "rank": rank, "lead_snp_id": cand, "chrom_chm13": m["chrom"], "pos_chm13": m["pos"], "ref_chm13": m["ref"], "alt_chm13": m["alt"], "ac": m["ac"], "an": m["an"], "af": m["af"], "r2": f"{value:.12g}", "distance_bp": dist})

