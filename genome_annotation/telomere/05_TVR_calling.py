#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_TVR_calling.py
Frame-aware TVR (telomere variant repeat) caller with orientation normalisation.

Orientation
-----------
Telomeric FASTA is extracted in assembly orientation, so p-arm sequence runs
distal -> proximal while q-arm runs proximal -> distal. To describe every arm
in a single coordinate system (0 bp = subtelomere-telomere boundary,
increasing = proximal -> distal) the caller:

    p-arm : reverse-complement the sequence before calling
    q-arm : keep the sequence as extracted

The arm is inferred from the FASTA record header (interval start < P_ARM_LIMIT
=> p-arm); use --orientation-mode none to disable normalisation.

Calling
-------
    for each of 2 strands (G-rich: TTAGGG, C-rich: CCCTAA)
        for each of 6 frame offsets
            scan consecutive 6-bp slots and label each slot:
                canonical  : equals the strand canonical 6-mer
                rotated    : one of the 5 phase-shifted rotations of it
                ambiguous  : contains N
                tvr        : any other 6-mer (true telomere variant repeat)
    choose the (strand, offset) maximising
                score = 10*n_canonical - 2*n_rotated - 3*n_ambiguous - 0.2*n_tvr
            (ties broken by more canonical, then fewer ambiguous, fewer rotated)

TVR blocks
----------
Consecutive TVR slots (contiguous in coordinate) are merged into one block.
A block is retained when min_tvr <= length <= max_tvr (default 6-300 bp).

Outputs
-------
--hexamer-tsv : sample, hap, chrom, start, end, hexamer, class
--blocks-tsv  : sample, hap, chrom, start, end, length, n_hexamers,
                dominant_hexamer, sequence

Example
-------
python 05_TVR_calling.py \
    --input sample.telo.fasta \
    --hexamer-tsv sample_hexamers.tsv \
    --blocks-tsv sample_tvr_blocks.tsv
"""
import argparse
import re
import sys
from collections import Counter, defaultdict


# -----------------------------------------------------------------------------
# sequence utilities
# -----------------------------------------------------------------------------
def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def canonical_for_strand(strand):
    return "TTAGGG" if strand == "G-rich" else "CCCTAA"


def rotations_of(hexamer):
    return {hexamer[i:] + hexamer[:i] for i in range(6)}


def revcomp_hexamer(h):
    return revcomp(h)


# -----------------------------------------------------------------------------
# header parsing / orientation
# -----------------------------------------------------------------------------
INTERVAL_RE = re.compile(r"(chr(?:[0-9]+|X|Y))[_:](\d+)-(\d+)", re.I)


def parse_header(header, p_arm_limit):
    """Return (sample, hap, chrom_field, start, end, is_p_arm)."""
    parts = header.split("#")
    if len(parts) >= 3:
        sample, hap = parts[0], parts[1]
        chrom_field = "#".join(parts[2:])
    elif len(parts) == 2:
        sample, hap, chrom_field = parts[0], "", parts[1]
    else:
        sample, hap, chrom_field = header, "", header

    m = INTERVAL_RE.search(chrom_field)
    if not m:
        return sample, hap, chrom_field, None, None, None
    start, end = int(m.group(2)), int(m.group(3))
    return sample, hap, chrom_field, start, end, start < p_arm_limit


# -----------------------------------------------------------------------------
# calling
# -----------------------------------------------------------------------------
def scan_frame(seq, strand, offset):
    canonical = canonical_for_strand(strand)
    rotated = rotations_of(canonical) - {canonical}
    records = []
    for i in range(offset, len(seq) - 5, 6):
        h = seq[i:i + 6]
        if len(h) < 6:
            continue
        if "N" in h:
            cls = "ambiguous"
        elif h == canonical:
            cls = "canonical"
        elif h in rotated:
            cls = "rotated"
        else:
            cls = "tvr"
        records.append({"start": i, "end": i + 6, "hexamer": h, "class": cls})
    return records


def score_frame(records):
    c = Counter(r["class"] for r in records)
    n_can = c.get("canonical", 0)
    n_rot = c.get("rotated", 0)
    n_amb = c.get("ambiguous", 0)
    n_tvr = c.get("tvr", 0)
    return 10 * n_can - 2 * n_rot - 3 * n_amb - 0.2 * n_tvr, c


def choose_best(seq):
    """Return (strand, offset, records, counts)."""
    best_key = None
    best = None
    for strand in ("G-rich", "C-rich"):
        for offset in range(6):
            records = scan_frame(seq, strand, offset)
            score, counts = score_frame(records)
            key = (score, counts.get("canonical", 0),
                   -counts.get("ambiguous", 0), -counts.get("rotated", 0))
            if best_key is None or key > best_key:
                best_key = key
                best = (strand, offset, records, counts)
    return best


def call_blocks(records):
    """Merge consecutive tvr slots into blocks (non-TVR classes break blocks)."""
    blocks = []
    current = []

    def flush():
        nonlocal current
        if not current:
            return
        start = current[0]["start"]
        end = current[-1]["end"]
        seq = "".join(r["hexamer"] for r in current)
        dominant = Counter(r["hexamer"] for r in current).most_common(1)[0][0]
        blocks.append({
            "start": start, "end": end, "length": end - start,
            "n_hexamers": len(current), "sequence": seq, "dominant_hexamer": dominant,
        })
        current = []

    for rec in records:
        if rec["class"] == "tvr":
            if current and rec["start"] != current[-1]["end"]:
                flush()
            current.append(rec)
        else:
            flush()
    flush()
    return blocks


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def iter_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="telomeric FASTA")
    ap.add_argument("--hexamer-tsv", required=True, help="output per-6-mer TSV")
    ap.add_argument("--blocks-tsv", required=True, help="output TVR block TSV")
    ap.add_argument("--orientation-mode", choices=["auto", "none"], default="auto",
                    help="auto = reverse-complement p-arm (default); none = as extracted")
    ap.add_argument("--p-arm-limit", type=int, default=20000,
                    help="interval start below which a record is treated as p-arm "
                         "(default: 20000)")
    ap.add_argument("--min-tvr", type=int, default=6, help="min TVR block length (default: 6)")
    ap.add_argument("--max-tvr", type=int, default=300, help="max TVR block length (default: 300)")
    ap.add_argument("--include-rotated-as-tvr", action="store_true",
                    help="treat rotated canonical slots as TVR (default: off)")
    args = ap.parse_args()

    n_rec = n_tvr_slots = n_blocks = 0
    frame_stats = defaultdict(int)

    with open(args.hexamer_tsv, "w") as fh_out, open(args.blocks_tsv, "w") as fb_out:
        fh_out.write("sample\thap\tchrom\tstart\tend\thexamer\tclass\n")
        fb_out.write("sample\thap\tchrom\tstart\tend\tlength\tn_hexamers\t"
                     "dominant_hexamer\tsequence\n")

        for header, seq in iter_fasta(args.input):
            sample, hap, chrom, start, end, is_p = parse_header(header, args.p_arm_limit)
            if seq is None or not seq:
                continue
            seq = seq.upper()

            normalised = False
            if args.orientation_mode == "auto" and is_p:
                seq = revcomp(seq)
                normalised = True

            strand, offset, records, counts = choose_best(seq)

            if args.include_rotated_as_tvr:
                for r in records:
                    if r["class"] == "rotated":
                        r["class"] = "tvr"
                records = [r for r in records]

            for r in records:
                h = r["hexamer"]
                if normalised:          # sequence reversed -> report reverse complement
                    h = revcomp_hexamer(h)
                fh_out.write(f"{sample}\t{hap}\t{chrom}\t{r['start']}\t{r['end']}\t"
                             f"{h}\t{r['class']}\n")

            for b in call_blocks(records):
                if not (args.min_tvr <= b["length"] <= args.max_tvr):
                    continue
                sequence = b["sequence"]
                if normalised:
                    sequence = revcomp(sequence)
                fb_out.write(f"{sample}\t{hap}\t{chrom}\t{b['start']}\t{b['end']}\t"
                             f"{b['length']}\t{b['n_hexamers']}\t{b['dominant_hexamer']}\t"
                             f"{sequence}\n")
                n_blocks += 1

            n_rec += 1
            n_tvr_slots += counts.get("tvr", 0)
            frame_stats[(strand, offset)] += 1

    print(f"[done] records={n_rec}, tvr_slots={n_tvr_slots}, blocks={n_blocks}")
    print(f"  hexamers -> {args.hexamer_tsv}")
    print(f"  blocks   -> {args.blocks_tsv}")
    top = sorted(frame_stats.items(), key=lambda x: -x[1])[:3]
    print("  most frequent (strand, offset): " +
          ", ".join(f"{s}/off{o}({n})" for (s, o), n in top))


if __name__ == "__main__":
    main()
