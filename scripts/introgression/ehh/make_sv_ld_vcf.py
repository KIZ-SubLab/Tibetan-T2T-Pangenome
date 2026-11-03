#!/usr/bin/env python3
"""Replace long sequence-resolved SV alleles with A/C for PLINK LD only."""
import sys

seen_format = False
for raw in sys.stdin:
    if raw.startswith("##"):
        if raw.startswith("##FORMAT="):
            continue
        if raw.startswith("##INFO=") or raw.startswith("##FILTER="):
            continue
        sys.stdout.write(raw)
        continue
    if raw.startswith("#CHROM"):
        if not seen_format:
            sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            sys.stdout.write('##source=SV alleles recoded A/C for PLINK LD; genotypes unchanged\n')
            seen_format = True
        cols = raw.rstrip("\n").split("\t")
        sys.stdout.write("\t".join(cols[:9] + cols[9:]) + "\n")
        continue
    cols = raw.rstrip("\n").split("\t")
    fmt = cols[8].split(":")
    try:
        gt_idx = fmt.index("GT")
    except ValueError as exc:
        raise SystemExit(f"GT absent at {cols[0]}:{cols[1]}") from exc
    gts = []
    for sample in cols[9:]:
        fields = sample.split(":")
        gts.append(fields[gt_idx] if gt_idx < len(fields) else ".")
    out = [cols[0], cols[1], cols[2], "A", "C", ".", "PASS", ".", "GT"] + gts
    sys.stdout.write("\t".join(out) + "\n")

