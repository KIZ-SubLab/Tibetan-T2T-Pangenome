#!/usr/bin/env python3
import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True)
ap.add_argument("--required-an", type=int, required=True)
ap.add_argument("--output", required=True)
a = ap.parse_args()

with open(a.input, newline="") as src, open(a.output, "w", newline="") as dst:
    reader = csv.DictReader(src, delimiter="\t")
    writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        if int(row["an"]) == a.required_an:
            writer.writerow(row)
