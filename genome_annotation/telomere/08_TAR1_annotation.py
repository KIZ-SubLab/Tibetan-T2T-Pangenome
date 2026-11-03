#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_TAR1_annotation.py
TAR1 annotation relative to the VGP telomere boundary (single code path for
any sample group; CZA and CHA differ only through CLI arguments).

Boundary
--------
Derived from the VGP telomeric windows (identity 99.9, threshold 0.4, 15 kb):

    p-arm : proximal end of the p-arm telomeric array  = max(end) of windows
    q-arm : proximal start of the q-arm telomeric array = min(start) of windows

Search window (towards chromosome interior only)
------------------------------------------------
    p-arm : [boundary, boundary + WINDOW]
    q-arm : [boundary - WINDOW, boundary]
    default WINDOW = 10,000 bp

Coordinate systems
------------------
RepeatMasker may have been run either on the whole assembly (GFF coordinates
are genomic) or on extracted subtelomeric windows (GFF coordinates are window
coordinates). This is handled by one generic mapping step:

    --gff-coordinate-mode genomic
        use GFF coordinates directly
    --gff-coordinate-mode windowed
        window coordinate x -> genomic coordinate = offset + x
        where offset is computed from the assembly .fai (chromosome length)
        and the extracted-window .fai (window length):
            chrN_start : offset = 0                    (window = [1, L])
            chrN_end   : offset = chr_len - L          (window = [chr_len-L+1, chr_len])

Effective haplotypes / arm classification
-----------------------------------------
n_effective(arm) = haplotypes that have a valid VGP boundary AND a RepeatMasker
                   record for that arm. Missing haplotypes are excluded from the
                   denominator and are NEVER scored as TAR1-negative.

    highest-positive threshold = ceil(n_effective / 2)
    arm is TAR1+ when n_positive >= ceil(n_effective / 2)   (>= 50%, inclusive)

Example
-------
python 08_TAR1_annotation.py \
    --vpg-bed-dir vpg_beds --gff-dir rm_gff --sample-list samples.txt \
    --output-dir tar1_out
# windowed GFF (RepeatMasker run on 150-kb extracts):
python 08_TAR1_annotation.py \
    --vpg-bed-dir vpg_beds --gff-dir rm_gff --sample-list samples.txt \
    --gff-coordinate-mode windowed \
    --assembly-fai-dir asm_fai --extract-fai-dir extract_fai \
    --output-dir tar1_out
"""
import argparse
import os
import re
import sys
from collections import defaultdict

CHROM_RE = re.compile(r"chr(?:\d+|X|Y)$", re.I)
SEQID_GENOMIC = re.compile(r"(?:.*#)?(chr(?:\d+|X|Y))$", re.I)
SEQID_WINDOW = re.compile(r"(chr(?:\d+|X|Y))_(start|end)$", re.I)


# -----------------------------------------------------------------------------
# VGP boundary
# -----------------------------------------------------------------------------
def load_boundaries(bed_path):
    """-> {arm: boundary_genomic_coord}; unplaced contigs are ignored."""
    if not os.path.exists(bed_path):
        return None
    p_end = {}
    q_start = {}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                c = line.split()
            if len(c) < 3:
                continue
            chrom = c[0].split("#")[-1]
            if not CHROM_RE.match(chrom):
                continue
            try:
                start, end = int(c[1]), int(c[2])
            except ValueError:
                continue
            if start < 20000:                       # p-arm telomeric array
                p_end[chrom] = max(p_end.get(chrom, 0), end)
            else:                                   # q-arm telomeric array
                q_start[chrom] = min(q_start.get(chrom, 10 ** 12), start)
    arms = {f"{c}_p": v for c, v in p_end.items()}
    arms.update({f"{c}_q": v for c, v in q_start.items()})
    return arms


# -----------------------------------------------------------------------------
# coordinate mapping for windowed GFF
# -----------------------------------------------------------------------------
def load_fai_lengths(path, key_regex=None):
    """-> {key: length}; key is the first column, optionally regex-extracted."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 2:
                continue
            key = c[0]
            if key_regex:
                m = key_regex.search(key)
                if not m:
                    continue
                key = m.group(1)
            try:
                out[key] = int(c[1])
            except ValueError:
                continue
    return out


def build_window_offset(assembly_fai, extract_fai, sample):
    """Return {seqid: genomic_offset} for chrN_start / chrN_end windows."""
    a = load_fai_lengths(os.path.join(assembly_fai, f"{sample}.fai"),
                         key_regex=CHROM_RE)
    e = load_fai_lengths(os.path.join(extract_fai, f"{sample}.fai"))
    offset = {}
    for seqid, L in e.items():
        m = SEQID_WINDOW.match(seqid)
        if not m:
            continue
        chrom, side = m.group(1), m.group(2).lower()
        if chrom not in a:
            continue
        offset[seqid] = 0 if side == "start" else a[chrom] - L
    return offset


# -----------------------------------------------------------------------------
# TAR1 features
# -----------------------------------------------------------------------------
def load_tar1(gff_path, mode, offset_map=None, p_arm_limit=20000):
    """-> {arm: [(genomic_start, genomic_end), ...]} or None if no GFF."""
    if not os.path.exists(gff_path):
        return None
    feats = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "TAR1" not in line:
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            try:
                gstart, gend = int(c[3]), int(c[4])
            except ValueError:
                continue
            seqid = c[0]
            if mode == "genomic":
                m = SEQID_GENOMIC.search(seqid)
                if not m:
                    continue
                chrom = m.group(1)
                feats[f"{chrom}_{'p' if gstart < p_arm_limit else 'q'}"].append(
                    (gstart, gend))
            else:
                m = SEQID_WINDOW.match(seqid)
                if not m or offset_map is None or seqid not in offset_map:
                    continue
                chrom, side = m.group(1), m.group(2).lower()
                off = offset_map[seqid]
                feats[f"{chrom}_{'p' if side == 'start' else 'q'}"].append(
                    (gstart + off, gend + off))
    return feats


def overlap(feats, window_start, window_end):
    return sum(1 for s, e in feats if s <= window_end and e >= window_start)


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def arm_sort_key(arm):
    m = re.match(r"chr(\d+)_([pq])", arm, re.I)
    if m:
        return (int(m.group(1)), 0 if m.group(2).lower() == "p" else 1)
    m = re.match(r"chr([XY])_([pq])", arm, re.I)
    if m:
        return (23 if m.group(1).upper() == "X" else 24,
                0 if m.group(2).lower() == "p" else 1)
    return (99, 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vpg-bed-dir", required=True,
                    help="dir with <sample> VGP 0.4 windows bed (one file or subdir per sample)")
    ap.add_argument("--vpg-bed-pattern", default="{sample}.windows.0.4.15kb.ends.bed",
                    help="filename pattern inside --vpg-bed-dir (default: %(default)s)")
    ap.add_argument("--gff-dir", required=True, help="dir with <sample>.gff RepeatMasker output")
    ap.add_argument("--gff-suffix", default=".gff")
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--window", type=int, default=10000,
                    help="search window from boundary towards interior (default: 10000)")
    ap.add_argument("--gff-coordinate-mode", choices=["genomic", "windowed"],
                    default="genomic")
    ap.add_argument("--assembly-fai-dir", default=None,
                    help="dir with <sample>.fai (required for windowed mode)")
    ap.add_argument("--extract-fai-dir", default=None,
                    help="dir with <sample>.fai of extracted windows (required for windowed)")
    ap.add_argument("--group-label", default="group")
    args = ap.parse_args()

    if args.gff_coordinate_mode == "windowed" and not (
            args.assembly_fai_dir and args.extract_fai_dir):
        sys.exit("ERROR: windowed mode requires --assembly-fai-dir and --extract-fai-dir")

    os.makedirs(args.output_dir, exist_ok=True)
    samples = [s.strip() for s in open(args.sample_list) if s.strip()]

    per_hap = []
    arm_pos = defaultdict(int)
    arm_eff = defaultdict(int)

    for sample in samples:
        bed = os.path.join(args.vpg_bed_dir, args.vpg_bed_pattern.format(sample=sample))
        if not os.path.exists(bed):
            print(f"[warn] no VGP bed for {sample}", file=sys.stderr)
            continue
        boundaries = load_boundaries(bed)
        if not boundaries:
            continue

        offset_map = None
        if args.gff_coordinate_mode == "windowed":
            offset_map = build_window_offset(args.assembly_fai_dir,
                                             args.extract_fai_dir, sample)
            if not offset_map:
                print(f"[warn] no window mapping for {sample}", file=sys.stderr)
                continue

        feats = load_tar1(os.path.join(args.gff_dir, f"{sample}{args.gff_suffix}"),
                          args.gff_coordinate_mode, offset_map)
        if feats is None:
            print(f"[warn] no GFF for {sample}", file=sys.stderr)
            continue

        for arm in sorted(boundaries, key=arm_sort_key):
            b = boundaries[arm]
            if arm.endswith("_p"):
                ws, we = b, b + args.window
            else:
                ws, we = b - args.window, b
            n_feat = overlap(feats.get(arm, []), ws, we)
            present = n_feat > 0
            per_hap.append((sample, arm, b, ws, we, present, n_feat))
            arm_eff[arm] += 1
            if present:
                arm_pos[arm] += 1

    out_hap = os.path.join(args.output_dir, f"tar1_per_hap_arm_{args.group_label}.tsv")
    with open(out_hap, "w") as fh:
        fh.write("sample\tarm\tVGP_boundary\twindow_start\twindow_end\t"
                 "TAR1_present\tn_TAR1_features\n")
        for r in per_hap:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t"
                     f"{'yes' if r[5] else 'no'}\t{r[6]}\n")

    out_arm = os.path.join(args.output_dir, f"tar1_arm_summary_{args.group_label}.tsv")
    with open(out_arm, "w") as fh:
        fh.write("group\tarm\tpositive\teffective\tfraction\tthreshold\tclass\n")
        for arm in sorted(arm_eff, key=arm_sort_key):
            n_pos, n_eff = arm_pos[arm], arm_eff[arm]
            thr = (n_eff + 1) // 2                   # ceil(n_eff / 2) -> >= 50%
            cls = "TAR1+" if n_pos >= thr else "TAR1-"
            frac = n_pos / n_eff if n_eff else 0.0
            fh.write(f"{args.group_label}\t{arm}\t{n_pos}\t{n_eff}\t{frac:.4f}\t"
                     f"{thr}\t{cls}\n")

    n_pos_arm = sum(1 for a in arm_eff if arm_pos[a] >= (arm_eff[a] + 1) // 2)
    print(f"[done] {args.group_label}: {len(per_hap)} hap x arm records, "
          f"{len(arm_eff)} arms")
    print(f"  TAR1+ arms: {n_pos_arm} / {len(arm_eff)}")
    print(f"  per-hap -> {out_hap}")
    print(f"  arms    -> {out_arm}")


if __name__ == "__main__":
    main()
