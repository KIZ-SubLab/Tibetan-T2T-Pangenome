#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_CFB_classification.py
Telomere architecture classification: Continuous (C), Fragmented (F), Blocky (B).

Rules (current final version)
-----------------------------
Motif counting
    Each telomeric sequence is scanned for the 10 telomeric 6-mers
    (TTAGGG, CTAGGG, TCAGGG, GTAGGG, TGAGGG and their C-rich complements)
    using non-overlapping string counts, as in the original implementation.
    total_motifs = sum of all 10 counts;  length_bp = total_motifs * 6.

B-type (blocky)
    total_motifs >= floor(min_length_bp / 6)          (default 2000 bp -> 333)
    AND at least one non-canonical motif (any motif other than TTAGGG/CCCTAA)
        has count > noncanon_pct * total_motifs       (default 5%)
    B-type telomeres are reported with their dominant non-canonical motif.

AC-type (candidate C or F)
    all remaining telomeres.

C vs F split within AC-type
    main_unit_ratio = max(count(TTAGGG), count(CCCTAA)) * 6 / length_bp
    C-type (continuous) if main_unit_ratio >  cf_threshold   (default 0.75)
    F-type (fragmented) otherwise.

    The threshold 0.75 is the value fixed after KDE (0.40), k-means elbow and
    Gaussian mixture model exploration; only the final fixed cut-off is applied
    here, no clustering is re-run.

Example
-------
python 03_CFB_classification.py \
    --fasta-dir telo_fasta --sample-list samples.txt \
    --output telomere_types.tsv
"""
import argparse
import os
import sys
from collections import Counter

CANONICAL = ("TTAGGG", "CCCTAA")
MOTIFS = ("TTAGGG", "CTAGGG", "TCAGGG", "GTAGGG", "TGAGGG",
          "CCCTAA", "CCCTGA", "CCCTAG", "CCCTCA", "CCCTAC")


def iter_fasta(path):
    """Yield (record_id, sequence) from a FASTA file."""
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def count_motifs(seq, motifs=MOTIFS):
    """Non-overlapping occurrence counts for each motif (original behaviour)."""
    return {m: seq.count(m) for m in motifs}


def classify(counts, length_bp, min_length_bp, noncanon_pct, cf_threshold):
    """Return (type, major_noncanonical_motif, main_unit_ratio, detail)."""
    total = sum(counts.values())
    min_motifs = int(min_length_bp / 6)
    main_unit = max(counts.get("TTAGGG", 0), counts.get("CCCTAA", 0))
    ratio = (main_unit * 6 / length_bp) if length_bp else 0.0

    # --- B-type ---
    if total >= min_motifs:
        threshold = total * noncanon_pct
        over = [(m, counts[m]) for m in MOTIFS
                if m not in CANONICAL and counts[m] > threshold]
        if over:
            # dominant non-canonical motif = first motif exceeding the cut
            major = max(over, key=lambda x: x[1])
            return "B", major[0], ratio, {
                "total_motifs": total,
                "noncanon_motifs_over_threshold": ";".join(f"{m}:{c}" for m, c in over),
            }

    # --- AC-type: split C vs F ---
    ctype = "C" if ratio > cf_threshold else "F"
    return ctype, "", ratio, {"total_motifs": total,
                              "noncanon_motifs_over_threshold": ""}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta-dir", required=True,
                    help="directory with <sample>.telo.fasta (multi-sequence)")
    ap.add_argument("--sample-list", required=True, help="one sample id per line")
    ap.add_argument("--output", required=True, help="output TSV")
    ap.add_argument("--counts-output", default=None,
                    help="optional per-telomere motif count table")
    ap.add_argument("--min-length-bp", type=float, default=2000.0,
                    help="minimum length for B-type (default: 2000)")
    ap.add_argument("--noncanon-pct", type=float, default=0.05,
                    help="non-canonical motif fraction for B-type (default: 0.05)")
    ap.add_argument("--cf-threshold", type=float, default=0.75,
                    help="main-unit ratio cut-off for C vs F (default: 0.75)")
    ap.add_argument("--fasta-suffix", default=".telo.fasta")
    args = ap.parse_args()

    samples = [s.strip() for s in open(args.sample_list) if s.strip()]
    rows = []
    type_counter = Counter()
    major_counter = Counter()

    for sample in samples:
        fasta = os.path.join(args.fasta_dir, f"{sample}{args.fasta_suffix}")
        if not os.path.exists(fasta):
            print(f"[warn] missing FASTA for {sample}", file=sys.stderr)
            continue
        for rec, seq in iter_fasta(fasta):
            counts = count_motifs(seq)
            length_bp = len(seq)
            ttype, major, ratio, detail = classify(
                counts, length_bp, args.min_length_bp,
                args.noncanon_pct, args.cf_threshold)
            type_counter[ttype] += 1
            if major:
                major_counter[major] += 1
            rows.append({
                "sample": sample, "telomere": rec, "length_bp": length_bp,
                "total_motifs": detail["total_motifs"], "main_unit_ratio": ratio,
                "type": ttype, "major_noncanonical": major,
                "noncanon_over_threshold": detail["noncanon_motifs_over_threshold"],
                "counts": counts,
            })

    with open(args.output, "w") as out:
        out.write("sample\ttelomere\tlength_bp\ttotal_motifs\tmain_unit_ratio\t"
                  "type\tmajor_noncanonical\t"
                  + "\t".join(f"n_{m}" for m in MOTIFS) + "\n")
        for r in rows:
            out.write(
                f"{r['sample']}\t{r['telomere']}\t{r['length_bp']}\t{r['total_motifs']}\t"
                f"{r['main_unit_ratio']:.4f}\t{r['type']}\t{r['major_noncanonical']}\t"
                + "\t".join(str(r["counts"][m]) for m in MOTIFS) + "\n")

    if args.counts_output:
        with open(args.counts_output, "w") as out:
            out.write("sample\ttelomere\tlength_bp\ttype\t"
                      + "\t".join(MOTIFS) + "\n")
            for r in rows:
                out.write(f"{r['sample']}\t{r['telomere']}\t{r['length_bp']}\t{r['type']}\t"
                          + "\t".join(str(r["counts"][m]) for m in MOTIFS) + "\n")

    total = len(rows)
    print(f"[done] {total} telomeres -> {args.output}")
    for t in ("C", "F", "B"):
        n = type_counter.get(t, 0)
        print(f"  {t}-type: {n} ({n / total * 100:.2f}%)" if total else f"  {t}-type: 0")
    if major_counter:
        top = major_counter.most_common(5)
        print("  top non-canonical motifs in B-type: "
              + ", ".join(f"{m}({c})" for m, c in top))


if __name__ == "__main__":
    main()
