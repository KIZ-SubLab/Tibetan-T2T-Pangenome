#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_gt(gt):
    return (int(gt[0]), int(gt[2])) if len(gt) == 3 and gt[1] == "|" and gt[0] in "01" and gt[2] in "01" else None


def load_vcf(bcf, vcf):
    samples = subprocess.check_output([bcf, "query", "-l", str(vcf)], text=True).splitlines()
    text = subprocess.check_output([bcf, "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n", str(vcf)], text=True)
    records, invalid = [], 0
    for line in text.splitlines():
        f = line.split("\t")
        parsed = [parse_gt(x) for x in f[5:]]
        if len(f) != 5 + len(samples) or len(f[3]) != 1 or len(f[4]) != 1 or any(x is None for x in parsed):
            invalid += 1
            continue
        records.append({"chrom": f[0], "pos": int(f[1]), "id": f[2], "ref": f[3].upper(),
                        "alt": f[4].upper(), "hap": np.array([x for pair in parsed for x in pair], dtype=np.int8)})
    return samples, records, invalid


def extend_classes(classes, alleles):
    _, inv = np.unique(np.column_stack((classes, alleles)), axis=0, return_inverse=True)
    return inv


def class_ehh(classes):
    n = len(classes)
    if n < 2:
        return np.nan
    counts = np.bincount(classes)
    return float(np.sum(counts * (counts - 1)) / (n * (n - 1)))


def directional_ehh(matrix, core, carriers):
    out = np.full(matrix.shape[1], np.nan)
    out[core] = 1.0
    classes = np.zeros(len(carriers), dtype=np.int32)
    for j in range(core - 1, -1, -1):
        classes = extend_classes(classes, matrix[carriers, j]); out[j] = class_ehh(classes)
    classes = np.zeros(len(carriers), dtype=np.int32)
    for j in range(core + 1, matrix.shape[1]):
        classes = extend_classes(classes, matrix[carriers, j]); out[j] = class_ehh(classes)
    return out


def unique_positions(records, core_id):
    by_pos = defaultdict(list)
    for rec in records:
        by_pos[rec["pos"]].append(rec)
    kept, dropped = [], 0
    for pos in sorted(by_pos):
        group = by_pos[pos]
        hit = [r for r in group if r["id"] == core_id]
        kept.append(hit[0] if hit else group[0]); dropped += len(group) - 1
    return kept, dropped


def sval(value):
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    return f"{value:.12g}" if isinstance(value, (float, np.floating)) else value


ap = argparse.ArgumentParser()
ap.add_argument("--manifest", required=True)
ap.add_argument("--proxy", required=True)
ap.add_argument("--clusters", required=True)
ap.add_argument("--phased-dir", required=True)
ap.add_argument("--bcftools", required=True)
ap.add_argument("--out-root", required=True)
ap.add_argument("--window", type=int, default=300000)
a = ap.parse_args()

manifest = rows(a.manifest)
proxy = {r["sv_id"]: r for r in rows(a.proxy) if r["status"] == "ok"}
cluster_rows = rows(a.clusters)
proxy_to_cluster = {}
for c in cluster_rows:
    for pid in c["proxy_ids"].split(","):
        proxy_to_cluster[pid] = c

cache = {}
summaries = []
for case in manifest:
    group, case_id, sid = case["group"], case["case_id"], case["original_sv_id"]
    row = proxy.get(sid)
    if row is None:
        summaries.append({"group": group, "case_id": case_id, "original_sv_id": sid,
                          "status": "NO_VALID_PROXY"})
        continue
    core_id = row["lead_snp_id"]
    cluster = proxy_to_cluster[core_id]
    cid = cluster["cluster_id"]
    if cid not in cache:
        vcf = Path(a.phased_dir) / f"{cid}.CZA35.GRCh38.phased.vcf.gz"
        cache[cid] = load_vcf(a.bcftools, vcf)
    samples, all_records, invalid = cache[cid]
    chrom, core_pos = row["lead_snp_chrom_grch38"], int(row["lead_snp_pos_grch38"])
    selected = [r for r in all_records if r["chrom"] == chrom and abs(r["pos"] - core_pos) <= a.window]
    selected, duplicate_drops = unique_positions(selected, core_id)
    selected.sort(key=lambda r: r["pos"])
    core_hits = [i for i, r in enumerate(selected) if r["id"] == core_id]
    if len(core_hits) != 1:
        raise RuntimeError(f"{group}/{case_id}: expected one core {core_id}; found {len(core_hits)}")
    core = core_hits[0]; core_rec = selected[core]
    aa, da = row["ancestral_allele_grch38"], row["derived_allele_grch38"]
    if aa == core_rec["ref"]:
        aa_index = 0
    elif aa == core_rec["alt"]:
        aa_index = 1
    else:
        raise RuntimeError(f"{group}/{case_id}: AA {aa} does not match {core_rec['ref']}/{core_rec['alt']}")
    da_index = 1 - aa_index
    matrix = np.stack([r["hap"] for r in selected], axis=1)
    aa_haps = np.flatnonzero(matrix[:, core] == aa_index)
    da_haps = np.flatnonzero(matrix[:, core] == da_index)
    pos = np.array([r["pos"] for r in selected]); rel = pos - core_pos
    case_dir = Path(a.out_root) / group / case_id
    case_dir.mkdir(parents=True, exist_ok=True)

    two_group = len(aa_haps) >= 2 and len(da_haps) >= 2
    if two_group:
        status = "OK_TWO_GROUP_AA_DA_EHH"
        ehh_aa = directional_ehh(matrix, core, aa_haps)
        ehh_da = directional_ehh(matrix, core, da_haps)
        ehh_major = np.full(len(selected), np.nan)
        major_allele = ""
    else:
        status = "INSUFFICIENT_FOR_TWO_GROUP_EHH_SINGLE_GROUP_MAJOR_ALLELE"
        major_index = 0 if np.sum(matrix[:, core] == 0) >= np.sum(matrix[:, core] == 1) else 1
        major_haps = np.flatnonzero(matrix[:, core] == major_index)
        major_allele = core_rec["ref"] if major_index == 0 else core_rec["alt"]
        ehh_major = directional_ehh(matrix, core, major_haps)
        ehh_aa = np.full(len(selected), np.nan); ehh_da = np.full(len(selected), np.nan)

    source = case_dir / f"{case_id}.GRCh38.EHH_source.tsv"
    with open(source, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chrom", "position", "distance_from_core_bp", "marker_id", "ref", "alt", "EHH_AA", "EHH_DA", "EHH_major"])
        for j, rec in enumerate(selected):
            w.writerow([chrom, rec["pos"], rel[j], rec["id"], rec["ref"], rec["alt"],
                        sval(ehh_aa[j]), sval(ehh_da[j]), sval(ehh_major[j])])

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = rel / 1000
    if two_group:
        ax.plot(x, ehh_aa, color="#D79232", lw=1.35, label=f"AA ancestral ({aa}); n={len(aa_haps)}")
        ax.plot(x, ehh_da, color="#327DA8", lw=1.35, label=f"DA derived ({da}); n={len(da_haps)}")
        plot_label = "AA_DA"
    else:
        ax.plot(x, ehh_major, color="#327DA8", lw=1.35, label=f"major allele ({major_allele}); n={len(major_haps)}")
        ax.text(0.01, 0.02, "Two-group EHH unavailable: one core-allele group has <2 haplotypes",
                transform=ax.transAxes, fontsize=7, color="#8A1C1C")
        plot_label = "single_group_major_allele"
    ax.axvline(0, color="#333333", lw=0.8, ls="--")
    ax.set_xlim(-a.window / 1000, a.window / 1000); ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Distance from proxy SNP (kb)"); ax.set_ylabel("EHH")
    ax.set_title(f"{group}: {case_id}\nGRCh38 proxy {core_id}; r²={float(row['r2']):.3f}")
    ax.legend(frameon=False, fontsize=8); ax.grid(axis="y", color="#DDDDDD", lw=0.5)
    fig.tight_layout()
    stem = case_dir / f"{case_id}.GRCh38.{plot_label}_EHH"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=300 if ext == "png" else None)
    plt.close(fig)

    metadata = {
        "group": group, "case_id": case_id, "original_sv_id": sid, "hg38_sv_id": case["hg38_sv_id"],
        "status": status, "proxy_id": core_id, "proxy_r2": row["r2"], "proxy_ld_rank": row["selected_rank"],
        "proxy_chrom_grch38": chrom, "proxy_pos_grch38": core_pos, "ancestral_allele": aa,
        "derived_allele": da, "AA_haplotype_count": len(aa_haps), "DA_haplotype_count": len(da_haps),
        "major_allele_if_single_group": major_allele, "sample_count": len(samples),
        "haplotype_count": 2 * len(samples), "markers_used": len(selected), "window_bp": a.window,
        "leftmost_distance_bp": int(rel[0]), "rightmost_distance_bp": int(rel[-1]),
        "duplicate_position_markers_dropped": duplicate_drops, "invalid_phased_records_excluded": invalid,
        "phase_cluster": cid,
    }
    with open(case_dir / "case_metadata.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(metadata), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerow(metadata)
    summaries.append(metadata)

fields = []
for row in summaries:
    for key in row:
        if key not in fields:
            fields.append(key)
for group in ("DEN", "NEA", "GHOST"):
    subset = [r for r in summaries if r["group"] == group]
    with open(Path(a.out_root) / group / "EHH_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(subset)
with open(Path(a.out_root) / "results" / "all_groups_EHH_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader(); w.writerows(summaries)
if len(summaries) != len(manifest):
    raise SystemExit(f"summary rows {len(summaries)} != manifest rows {len(manifest)}")
print(f"EHH_COMPLETE cases={len(summaries)}")

