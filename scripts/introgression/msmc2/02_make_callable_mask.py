#!/usr/bin/env python3
"""Create a merged callable BED from one-sample gVCF reference blocks/calls."""

import argparse
import json
from collections import Counter

import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gvcf", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--stats", required=True)
    p.add_argument("--contig", default="chr1")
    p.add_argument("--min-gq", type=float, default=20)
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--max-depth", type=int, default=60)
    return p.parse_args()


def main():
    args = parse_args()
    stats = Counter()
    total_bases = Counter()
    current = None

    def write_or_merge(out, chrom, start, end):
        nonlocal current
        if current is None:
            current = [chrom, start, end]
        elif chrom == current[0] and start <= current[2]:
            current[2] = max(current[2], end)
        else:
            out.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")
            current = [chrom, start, end]

    with pysam.VariantFile(args.gvcf) as vcf, open(args.output, "w") as out:
        if args.sample not in vcf.header.samples:
            raise ValueError(f"sample {args.sample!r} is absent from the gVCF")

        for rec in vcf.fetch(args.contig):
            stats["records_seen"] += 1
            length = max(1, rec.stop - rec.start)
            total_bases["record_span_seen"] += length
            call = rec.samples[args.sample]
            gt = call.get("GT")
            if not gt or any(a is None for a in gt):
                stats["fail_missing_gt"] += 1
                total_bases["fail_missing_gt"] += length
                continue

            gq = call.get("GQ")
            if gq is None or gq < args.min_gq:
                stats["fail_gq"] += 1
                total_bases["fail_gq"] += length
                continue

            depth = call.get("MIN_DP")
            if depth is None:
                depth = call.get("DP")
            if depth is None:
                stats["fail_missing_depth"] += 1
                total_bases["fail_missing_depth"] += length
                continue
            if depth < args.min_depth:
                stats["fail_low_depth"] += 1
                total_bases["fail_low_depth"] += length
                continue
            if depth > args.max_depth:
                stats["fail_high_depth"] += 1
                total_bases["fail_high_depth"] += length
                continue

            # gVCF reference blocks use an empty FILTER value; true variant calls
            # must be PASS. RefCall records are deliberately excluded.
            filters = set(rec.filter.keys())
            has_real_alt = any(a not in (None, ".", "<*>", "<NON_REF>") for a in (rec.alts or ()))
            if has_real_alt and "PASS" not in filters:
                stats["fail_variant_filter"] += 1
                total_bases["fail_variant_filter"] += length
                continue

            stats["records_callable"] += 1
            total_bases["callable_record_span"] += length
            write_or_merge(out, rec.contig, rec.start, rec.stop)

        if current is not None:
            out.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")

    payload = {
        "input": args.gvcf,
        "sample": args.sample,
        "contig": args.contig,
        "thresholds": {
            "min_gq": args.min_gq,
            "min_depth": args.min_depth,
            "max_depth": args.max_depth,
            "variant_filter": "PASS; reference blocks permit empty FILTER",
        },
        "record_counts": dict(stats),
        "base_span_sums_before_union": dict(total_bases),
        "coordinate_system": "BED: zero-based, half-open",
    }
    with open(args.stats, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main()
