#!/usr/bin/env python3
"""
Project a column block of the full-flank MSA back onto each taxon's
own unaligned orthologue sequence.

The ruler taxon has no gaps in the full-flank alignment (its non-gap
count equals its orthologue length), so the k-th non-gap base of the
ruler identifies which base of every other taxon belongs in the window.
"""
import argparse


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msa", required=True)
    ap.add_argument("--ruler-taxon", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    msa = read_fasta(a.msa)
    ruler = msa[a.ruler_taxon]

    # k-th non-gap base of the ruler -> alignment column
    cols = [i for i, c in enumerate(ruler) if c not in "-Nn"]
    block = cols[a.start:a.end]

    with open(a.output, "w") as out:
        for taxon, seq in msa.items():
            out.write(f">{taxon}\n")
            s = "".join(seq[i] for i in block).replace("-", "")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")


if __name__ == "__main__":
    main()
