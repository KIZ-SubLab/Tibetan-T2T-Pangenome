#!/usr/bin/env python3
import csv
import gzip
import sys


path = sys.argv[1]
regions = []
for item in sys.argv[2:]:
    chrom, bounds = item.split(":", 1)
    start, end = map(int, bounds.split("-", 1))
    regions.append((chrom, start, end))

stats = {region: {"count": 0, "max": float("-inf")} for region in regions}
chrom_stats = {chrom: {"count": 0, "min": None, "max": None} for chrom, _, _ in regions}
opener = gzip.open if path.endswith(".gz") else open
with opener(path, "rt", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    print("header=" + repr(reader.fieldnames))
    for row in reader:
        chrom = row["CHR"]
        if chrom not in chrom_stats:
            continue
        pos = int(row["POS"])
        cs = chrom_stats[chrom]
        cs["count"] += 1
        cs["min"] = pos if cs["min"] is None else min(cs["min"], pos)
        cs["max"] = pos if cs["max"] is None else max(cs["max"], pos)
        for region in regions:
            rchrom, start, end = region
            if chrom == rchrom and start <= pos <= end:
                stats[region]["count"] += 1
                stats[region]["max"] = max(stats[region]["max"], float(row["PBS"]))

for region in regions:
    chrom, start, end = region
    value = stats[region]
    print(f"{chrom}:{start}-{end}\tcount={value['count']}\tPBS_max={value['max']}")
for chrom, value in chrom_stats.items():
    print(f"{chrom}\ttotal={value['count']}\trange={value['min']}-{value['max']}")
