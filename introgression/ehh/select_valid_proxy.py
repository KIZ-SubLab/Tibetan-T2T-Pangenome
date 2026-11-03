#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def comp(x):
    return x.translate(COMP).upper()


class Fasta:
    def __init__(self, path):
        self.path = path
        self.fp = open(path, "rb")
        self.index = {}
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, line_bases, line_bytes, *_ = line.rstrip().split("\t")
                self.index[name] = tuple(map(int, (length, offset, line_bases, line_bytes)))

    def base(self, name, pos1):
        length, offset, line_bases, line_bytes = self.index[name]
        if pos1 < 1 or pos1 > length:
            return ""
        p0 = pos1 - 1
        byte = offset + (p0 // line_bases) * line_bytes + (p0 % line_bases)
        self.fp.seek(byte)
        return self.fp.read(1).decode().upper()


ap = argparse.ArgumentParser()
ap.add_argument("--candidates", required=True)
ap.add_argument("--ids", required=True)
ap.add_argument("--sv-meta", required=True)
ap.add_argument("--mapped", required=True)
ap.add_argument("--unmapped", required=True)
ap.add_argument("--hg38-ref", required=True)
ap.add_argument("--ancestral-dir", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

ids = [x.strip() for x in open(a.ids) if x.strip() and not x.startswith("#")]
with open(a.sv_meta, newline="") as fh:
    svmeta = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
by_sv = defaultdict(list)
with open(a.candidates, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        by_sv[row["sv_id"]].append(row)

mapped = {}
with open(a.mapped) as fh:
    for line in fh:
        if line.strip() and not line.startswith("#"):
            chrom, start, end, vid, score, strand = line.rstrip().split("\t")[:6]
            mapped[vid] = (chrom, int(start) + 1, strand)

hg = Fasta(a.hg38_ref)
anc = {}


def anc_base(chrom, pos):
    n = chrom.removeprefix("chr")
    if n not in anc:
        path = str(Path(a.ancestral_dir) / f"homo_sapiens_ancestor_{n}.fa")
        anc[n] = Fasta(path)
    seq = next(iter(anc[n].index))
    return anc[n].base(seq, pos)


fields = [
    "status", "sv_id", "sv_chrom_chm13", "sv_pos_chm13", "sv_ref", "sv_alt", "selected_rank",
    "lead_snp_id", "lead_snp_chrom_chm13", "lead_snp_pos_chm13", "lead_snp_ref_chm13", "lead_snp_alt_chm13",
    "lead_snp_ac", "lead_snp_an", "lead_snp_af", "r2", "distance_bp", "lead_snp_chrom_grch38",
    "lead_snp_pos_grch38", "chain_strand", "lead_snp_ref_grch38", "lead_snp_alt_grch38",
    "genotype_allele_index_swap", "ancestral_allele_grch38", "derived_allele_grch38",
    "ancestral_allele_chm13", "derived_allele_chm13", "selection_rule"
]
with open(a.out, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for sid in ids:
        sm = svmeta[sid]
        answer = {"status": "no_mappable_polarizable_candidate", "sv_id": sid, "sv_chrom_chm13": sm["chrom"], "sv_pos_chm13": sm["pos"], "sv_ref": sm["ref"], "sv_alt": sm["alt"], "selection_rule": "AN=70 required; max_r2_then_nearest_then_lowest_position_then_id; first candidate mappable to GRCh38 and polarizable; no minimum r2"}
        for row in by_sv[sid]:
            vid = row["lead_snp_id"]
            if vid not in mapped:
                continue
            hchrom, hpos, strand = mapped[vid]
            r0, a0 = row["ref_chm13"].upper(), row["alt_chm13"].upper()
            ro = r0 if strand == "+" else comp(r0)
            ao = a0 if strand == "+" else comp(a0)
            href = hg.base(hchrom, hpos)
            if href == ro:
                halt, swap = ao, "0"
            elif href == ao:
                halt, swap = ro, "1"
            else:
                continue
            aa = anc_base(hchrom, hpos)
            if aa == href:
                da = halt
            elif aa == halt:
                da = href
            else:
                continue
            answer.update({
                "status": "ok", "selected_rank": row["rank"], "lead_snp_id": vid,
                "lead_snp_chrom_chm13": row["chrom_chm13"], "lead_snp_pos_chm13": row["pos_chm13"],
                "lead_snp_ref_chm13": r0, "lead_snp_alt_chm13": a0, "lead_snp_ac": row["ac"],
                "lead_snp_an": row["an"], "lead_snp_af": row["af"], "r2": row["r2"],
                "distance_bp": row["distance_bp"], "lead_snp_chrom_grch38": hchrom,
                "lead_snp_pos_grch38": hpos, "chain_strand": strand, "lead_snp_ref_grch38": href,
                "lead_snp_alt_grch38": halt, "genotype_allele_index_swap": swap,
                "ancestral_allele_grch38": aa, "derived_allele_grch38": da,
                "ancestral_allele_chm13": aa if strand == "+" else comp(aa),
                "derived_allele_chm13": da if strand == "+" else comp(da),
            })
            break
        w.writerow(answer)
