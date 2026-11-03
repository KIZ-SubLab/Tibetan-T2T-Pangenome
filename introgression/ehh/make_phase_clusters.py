#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
from collections import defaultdict


ap = argparse.ArgumentParser()
ap.add_argument("--proxy", required=True)
ap.add_argument("--bcftools", required=True)
ap.add_argument("--example-vcf", required=True)
ap.add_argument("--hg38-fai", required=True)
ap.add_argument("--padding", type=int, default=5_000_000)
ap.add_argument("--out", required=True)
a = ap.parse_args()

header = subprocess.check_output([a.bcftools, "view", "-h", a.example_vcf], text=True)
chm_len = {}
for line in header.splitlines():
    m = re.match(r"##contig=<ID=([^,>]+),length=([0-9]+)", line)
    if m:
        chm_len[m.group(1)] = int(m.group(2))
hg_len = {}
with open(a.hg38_fai) as fh:
    for line in fh:
        f = line.split("\t")
        hg_len[f[0]] = int(f[1])

by_chrom = defaultdict(list)
with open(a.proxy, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["status"] == "ok":
            by_chrom[row["lead_snp_chrom_grch38"]].append(row)

clusters = []
for chrom in sorted(by_chrom, key=lambda x: int(x.removeprefix("chr"))):
    ordered = sorted(by_chrom[chrom], key=lambda r: int(r["lead_snp_pos_grch38"]))
    current = []
    current_end = -1
    for row in ordered:
        pos = int(row["lead_snp_pos_grch38"])
        start, end = max(1, pos - a.padding), min(hg_len[chrom], pos + a.padding)
        if current and start > current_end + 1:
            clusters.append((chrom, current))
            current = []
            current_end = -1
        current.append(row)
        current_end = max(current_end, end)
    if current:
        clusters.append((chrom, current))

fields = ["index", "cluster_id", "chrom", "chm13_start", "chm13_end",
          "grch38_start", "grch38_end", "proxy_ids", "sv_ids"]
with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    per_chr = defaultdict(int)
    for index, (chrom, rows) in enumerate(clusters, 1):
        per_chr[chrom] += 1
        cluster_id = f"{chrom}.cluster{per_chr[chrom]:02d}"
        cp = [int(r["lead_snp_pos_chm13"]) for r in rows]
        hp = [int(r["lead_snp_pos_grch38"]) for r in rows]
        cchroms = {r["lead_snp_chrom_chm13"] for r in rows}
        if cchroms != {chrom}:
            raise SystemExit(f"{cluster_id}: unexpected cross-chromosome lift: {cchroms}")
        w.writerow({
            "index": index, "cluster_id": cluster_id, "chrom": chrom,
            "chm13_start": max(1, min(cp) - a.padding),
            "chm13_end": min(chm_len[chrom], max(cp) + a.padding),
            "grch38_start": max(1, min(hp) - a.padding),
            "grch38_end": min(hg_len[chrom], max(hp) + a.padding),
            "proxy_ids": ",".join(r["lead_snp_id"] for r in rows),
            "sv_ids": ",".join(r["sv_id"] for r in rows),
        })
print(f"phase_clusters={len(clusters)}")

