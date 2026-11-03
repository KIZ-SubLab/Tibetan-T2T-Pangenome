#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_Btype_variants.py
Non-canonical motif spectrum of B-type telomeres.

Reads the per-telomere motif count table produced by 03_CFB_classification.py
(columns: sample, telomere, length_bp, type, n_TTAGGG, n_CTAGGG, ...).

Outputs
-------
1. Motif frequency table (non-canonical motifs only), per telomere type.
2. Substitution spectrum: each non-canonical 6-mer is aligned to its
   strand-matched canonical repeat (G-rich -> TTAGGG, C-rich -> CCCTAA) and
   every mismatching position is recorded as ref>alt at position 1..6.
   No gaps / no alignment algorithm is used (position-wise comparison of
   equal-length 6-mers).

Example
-------
python 04_Btype_variants.py \
    --input telomere_types.tsv --type B \
    --motif-output Btype_motif_frequency.tsv \
    --substitution-output Btype_substitutions.tsv
"""
import argparse
import sys
from collections import Counter, defaultdict

CANONICAL_G = "TTAGGG"
CANONICAL_C = "CCCTAA"
CANONICAL = {CANONICAL_G, CANONICAL_C}


def is_crich(motif):
    """C-rich motifs are the reverse-complementary family (start with CC)."""
    return motif.startswith("CC")


def strand_canonical(motif):
    return CANONICAL_C if is_crich(motif) else CANONICAL_G


def read_counts(path):
    """Read the motif count table written by 03_CFB_classification.py."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        motif_cols = [h[2:] for h in header if h.startswith("n_")]
        rows = []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(header):
                continue
            rows.append({
                "sample": c[idx["sample"]],
                "telomere": c[idx["telomere"]],
                "length_bp": int(c[idx["length_bp"]]),
                "type": c[idx["type"]],
                "counts": {m: int(c[idx["n_" + m]]) for m in motif_cols},
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="telomere type table from 03_CFB_classification.py")
    ap.add_argument("--type", default="B", help="telomere type to summarise (default: B)")
    ap.add_argument("--motif-output", required=True, help="motif frequency TSV")
    ap.add_argument("--substitution-output", required=True,
                    help="substitution spectrum TSV")
    args = ap.parse_args()

    rows = [r for r in read_counts(args.input) if r["type"] == args.type]
    if not rows:
        sys.exit(f"ERROR: no telomeres of type {args.type} in {args.input}")

    # --- motif frequency ---
    motif_total = Counter()
    per_telomere_hits = Counter()
    for r in rows:
        for m, n in r["counts"].items():
            if m in CANONICAL or n == 0:
                continue
            motif_total[m] += n
            per_telomere_hits[m] += 1

    grand_total = sum(motif_total.values())
    with open(args.motif_output, "w") as out:
        out.write("motif\tcount\tfraction_of_noncanonical_motifs\t"
                  "n_telomeres_with_motif\tstrand\n")
        for m, n in motif_total.most_common():
            out.write(f"{m}\t{n}\t{n / grand_total:.4f}\t{per_telomere_hits[m]}\t"
                      f"{'C-rich' if is_crich(m) else 'G-rich'}\n")

    # --- substitution spectrum ---
    # weight each motif by its total occurrence (motif-level spectrum)
    pos_change = defaultdict(int)       # (pos, ref>alt) -> weighted count
    pos_total = Counter()               # pos -> weighted count
    strand_change = defaultdict(int)    # (strand, pos, ref>alt) -> count

    for r in rows:
        for m, n in r["counts"].items():
            if m in CANONICAL or n == 0:
                continue
            ref = strand_canonical(m)
            strand = "C-rich" if is_crich(m) else "G-rich"
            for i in range(6):
                if m[i] != ref[i]:
                    key = (i + 1, f"{ref[i]}>{m[i]}")
                    pos_change[key] += n
                    pos_total[i + 1] += n
                    strand_change[(strand, i + 1, f"{ref[i]}>{m[i]}")] += n

    with open(args.substitution_output, "w") as out:
        out.write("position\tchange\tcount\tfraction_at_position\n")
        for (pos, change), n in sorted(pos_change.items(),
                                       key=lambda x: (-x[1], x[0])):
            out.write(f"{pos}\t{change}\t{n}\t{n / pos_total[pos]:.4f}\n")

    # --- console summary ---
    print(f"[done] {len(rows)} {args.type}-type telomeres")
    if motif_total:
        top_m, top_n = motif_total.most_common(1)[0]
        print(f"  most frequent non-canonical motif: {top_m} "
              f"({top_n} copies, {top_n / grand_total * 100:.2f}%)")
    print("  motif frequency  -> " + args.motif_output)
    print("  substitution     -> " + args.substitution_output)
    top_subs = sorted(pos_change.items(), key=lambda x: -x[1])[:5]
    if top_subs:
        print("  top substitutions: " +
              ", ".join(f"pos{p} {c} ({n})" for (p, c), n in top_subs))


if __name__ == "__main__":
    main()
