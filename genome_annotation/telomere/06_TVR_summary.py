#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_TVR_summary.py
Per-chromosome-arm TVR statistics for an arbitrary sample set (no plotting).

Inputs
------
--hexamer-dir : directory of <sample>_hexamers.tsv  (05_TVR_calling.py output)
--blocks-dir  : directory of <sample>_tvr_blocks.tsv
--sample-list : one sample id per line (sample IDs are never hard-coded)

Metrics (one row per chromosome arm)
------------------------------------
tvr_load
    sum(TVR block length) / analysed telomeric interval length
    - numerator  : sum of block lengths, only blocks with 6 <= length <= 300 bp
    - denominator: length of the telomeric interval the block was called in,
                   taken from the BED-style interval in the `chrom` column
                   (end - start + 1)
    - canonical / rotated / ambiguous slots are NOT part of any block, so they
      never contribute to the numerator.

variant length class (per 6-mer)
    hamming = number of mismatching bases between the 6-mer and the strand's
    canonical repeat (TTAGGG for G-rich, CCCTAA for C-rich); the strand is
    inferred per sample x arm from canonical 6-mer counts.
        1bp (hamming 1), 2bp (2), 3-4bp (3-4), 5-6bp (5-6)

substitution vs indel-associated (per 6-mer, only class == tvr)
    indel_associated : the enclosing TVR block has n_hexamers > 0 and
                       length != n_hexamers*6 (length anomaly), OR the slot is
                       not contiguous with the previous TVR slot on the same
                       arm (start != previous end)
    substitution     : every other TVR 6-mer
    No sequence alignment / edit-distance operation is used.

Example
-------
python 06_TVR_summary.py \
    --hexamer-dir hex --blocks-dir blocks --sample-list group_a.txt \
    --label A --output tvr_summary_A.tsv
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict

CANONICAL_G = "TTAGGG"
CANONICAL_C = "CCCTAA"
HAM_BINS = ("1bp", "2bp", "3-4bp", "5-6bp")
INTERVAL_RE = re.compile(r"[_:](\d+)-(\d+)")


def parse_arm(chrom_field, p_arm_limit=20000):
    """chr5_182094801-182102400 -> chr5_q ; chr2_1-9600 -> chr2_p"""
    if "#" in chrom_field:
        chrom_field = chrom_field.split("#")[-1]
    m = re.match(r"(chr(?:\d+|X|Y))[_:]", chrom_field, re.I)
    if not m:
        m = re.match(r"(chr(?:\d+|X|Y))$", chrom_field, re.I)
        return m.group(1).lower() if m else None
    chrom = m.group(1)
    mi = INTERVAL_RE.search(chrom_field)
    if not mi:
        return None
    start = int(mi.group(1))
    return f"{chrom}_{'p' if start < p_arm_limit else 'q'}"


def interval_length(chrom_field):
    mi = INTERVAL_RE.search(chrom_field)
    if not mi:
        return None
    return int(mi.group(2)) - int(mi.group(1)) + 1


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def hamming_bin(hd):
    if hd <= 0:
        return None
    if hd == 1:
        return "1bp"
    if hd == 2:
        return "2bp"
    if hd <= 4:
        return "3-4bp"
    return "5-6bp"


def infer_strand(hexamer_path):
    """Per arm: G-rich if TTAGGG canonical slots outnumber CCCTAA slots."""
    counts = defaultdict(Counter)
    with open(hexamer_path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 7 or c[6] != "canonical":
                continue
            arm = parse_arm(c[2])
            h = c[5].upper()
            if arm and h in (CANONICAL_G, CANONICAL_C):
                counts[arm][h] += 1
    return {a: ("G-rich" if c.get(CANONICAL_G, 0) >= c.get(CANONICAL_C, 0) else "C-rich")
            for a, c in counts.items()}


def load_blocks(path):
    """-> {(arm): [(start, end, length, n_hex)]}"""
    blocks = defaultdict(list)
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 7:
                continue
            arm = parse_arm(c[2])
            if not arm:
                continue
            blocks[arm].append((int(c[3]), int(c[4]), int(c[5]), int(c[6])))
    return blocks


def summarise_sample(hexamer_path, blocks_path, p_arm_limit):
    """Return {arm: dict} for one sample."""
    strands = infer_strand(hexamer_path)
    blocks = load_blocks(blocks_path)

    acc = defaultdict(lambda: {
        "interval_bp": 0, "block_bp": 0, "n_tvr": 0,
        "hb": Counter(), "mt": Counter(), "blocks": 0,
    })
    seen_interval = {}
    prev_end = {}
    indel_keys = set()
    for arm, bl in blocks.items():
        for s, e, ln, nh in bl:
            if nh > 0 and ln != nh * 6:
                indel_keys.add((arm, s, e))

    with open(hexamer_path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 7:
                continue
            chrom, start_s, end_s, hexamer = c[2], c[3], c[4], c[5].upper()
            cls = c[6]
            arm = parse_arm(chrom, p_arm_limit)
            if not arm:
                continue
            if arm not in seen_interval:
                L = interval_length(chrom)
                if L:
                    seen_interval[arm] = L
                    acc[arm]["interval_bp"] = L
            if cls != "tvr":
                continue

            start, end = int(start_s), int(end_s)
            strand = strands.get(arm, "G-rich")
            canonical = CANONICAL_G if strand == "G-rich" else CANONICAL_C
            hd = hamming(hexamer, canonical)
            hb = hamming_bin(hd)

            is_indel = (arm, start, end) in indel_keys
            pe = prev_end.get(arm)
            if not is_indel and pe is not None and start != pe:
                is_indel = True
            prev_end[arm] = end

            a = acc[arm]
            a["n_tvr"] += 1
            if hb:
                a["hb"][hb] += 1
            a["mt"]["indel_associated" if is_indel else "substitution"] += 1

    for arm, bl in blocks.items():
        kept = [ln for (s, e, ln, nh) in bl if 6 <= ln <= 300]
        acc[arm]["block_bp"] += sum(kept)
        acc[arm]["blocks"] += len(kept)

    return acc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hexamer-dir", required=True)
    ap.add_argument("--blocks-dir", required=True)
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--label", default="population", help="group label for output rows")
    ap.add_argument("--p-arm-limit", type=int, default=20000)
    args = ap.parse_args()

    samples = [s.strip() for s in open(args.sample_list) if s.strip()]
    per_arm_loads = defaultdict(list)
    per_arm_totals = defaultdict(lambda: Counter())
    overall_loads = []

    for sample in samples:
        hp = os.path.join(args.hexamer_dir, f"{sample}_hexamers.tsv")
        bp = os.path.join(args.blocks_dir, f"{sample}_tvr_blocks.tsv")
        if not (os.path.exists(hp) and os.path.exists(bp)):
            print(f"[warn] missing input for {sample}", file=sys.stderr)
            continue
        acc = summarise_sample(hp, bp, args.p_arm_limit)
        for arm, a in acc.items():
            if a["interval_bp"] <= 0:
                continue
            load = a["block_bp"] / a["interval_bp"]
            per_arm_loads[arm].append(load)
            overall_loads.append(load)
            c = per_arm_totals[arm]
            c["n_hap"] += 1
            c["n_tvr"] += a["n_tvr"]
            c.update(a["hb"])
            c.update(a["mt"])

    def pct(num, den):
        return num / den * 100 if den else 0.0

    with open(args.output, "w") as out:
        out.write("population\tarm\tn_hap\ttvr_load_mean\ttvr_load_sd\t"
                  "pct_1bp\tpct_2bp\tpct_3_4bp\tpct_5_6bp\t"
                  "pct_substitution\tpct_indel_associated\tn_tvr\n")
        for arm in sorted(per_arm_loads):
            loads = per_arm_loads[arm]
            mean = sum(loads) / len(loads)
            sd = (sum((x - mean) ** 2 for x in loads) / len(loads)) ** 0.5
            c = per_arm_totals[arm]
            n = c["n_tvr"]
            out.write(
                f"{args.label}\t{arm}\t{c['n_hap']}\t{mean:.6f}\t{sd:.6f}\t"
                f"{pct(c['1bp'], n):.2f}\t{pct(c['2bp'], n):.2f}\t"
                f"{pct(c['3-4bp'], n):.2f}\t{pct(c['5-6bp'], n):.2f}\t"
                f"{pct(c['substitution'], n):.2f}\t{pct(c['indel_associated'], n):.2f}\t{n}\n")

    if overall_loads:
        gm = sum(overall_loads) / len(overall_loads)
        print(f"[done] {args.label}: {len(per_arm_loads)} arms, "
              f"{len(overall_loads)} hap-arm observations")
        print(f"  pooled overall mean TVR load = {gm:.4f} ({gm * 100:.2f}%)")
    print(f"  -> {args.output}")


if __name__ == "__main__":
    main()
