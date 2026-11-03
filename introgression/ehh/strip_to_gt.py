#!/usr/bin/env python3
import argparse
import sys


ap = argparse.ArgumentParser()
ap.add_argument("--sample-map", required=True)
a = ap.parse_args()
mapping = {}
with open(a.sample_map) as fh:
    for line in fh:
        old, new = line.rstrip().split("\t")[:2]
        mapping[old] = new

for raw in sys.stdin:
    if raw.startswith("##"):
        if raw.startswith("##FORMAT=") or raw.startswith("##INFO="):
            continue
        sys.stdout.write(raw)
    elif raw.startswith("#CHROM"):
        sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        sys.stdout.write('##source=DeepVariant-GLnexus GT-only input for SHAPEIT5\n')
        f = raw.rstrip().split("\t")
        samples = [mapping.get(x, x) for x in f[9:]]
        sys.stdout.write("\t".join(f[:9] + samples) + "\n")
    else:
        f = raw.rstrip().split("\t")
        fmt = f[8].split(":")
        gi = fmt.index("GT")
        gts = []
        for sample in f[9:]:
            x = sample.split(":")
            gt = x[gi] if gi < len(x) else "./."
            if gt in {"0/.", "./0", "1/.", "./1", "0|.", ".|0", "1|.", ".|1"}:
                gt = "./."
            gts.append(gt)
        sys.stdout.write("\t".join([f[0], f[1], f[2], f[3], f[4], f[5], f[6], ".", "GT"] + gts) + "\n")

