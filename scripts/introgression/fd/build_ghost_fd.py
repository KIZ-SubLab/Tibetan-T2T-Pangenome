#!/usr/bin/env python3
"""DISC1 ghost-introgression A+B+C D/fd analysis on a common SNP set.

A: (P1,P2,P3,P4) = (Denisovan, Neanderthal, TBN, YRI)
B: (P1,P2,P3,P4) = (YRI, TBN, Denisovan, ancestral)
C: (P1,P2,P3,P4) = (YRI, TBN, Neanderthal, ancestral)

TBN variant-only SNPs define the candidate universe.  All three comparisons use
the same sites, requiring exact/called data in both archaics, YRI, TBN and a
resolved GRCh38 ancestral allele.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
from collections import Counter

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from scipy.stats import norm


COMPARISONS = {
    "A_relative": {
        "label": "D/fd(Denisovan,Neanderthal,TBN,YRI)",
        "groups": ("den", "nea", "tbn", "yri"),
        "positive": "positive: TBN is relatively closer to Neanderthal than Denisovan",
        "negative": "negative: TBN is relatively closer to Denisovan than Neanderthal",
    },
    "B_Denisovan": {
        "label": "D/fd(YRI,TBN,Denisovan,ancestral)",
        "groups": ("yri", "tbn", "den", "anc"),
        "positive": "positive: excess TBN-Denisovan derived-allele sharing",
        "negative": "negative: no targeted TBN-Denisovan excess; direction favors YRI-Denisovan",
    },
    "C_Neanderthal": {
        "label": "D/fd(YRI,TBN,Neanderthal,ancestral)",
        "groups": ("yri", "tbn", "nea", "anc"),
        "positive": "positive: excess TBN-Neanderthal derived-allele sharing",
        "negative": "negative: no targeted TBN-Neanderthal excess; direction favors YRI-Neanderthal",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--kg-vcf", required=True)
    p.add_argument("--fam", required=True)
    p.add_argument("--tbn-vcf", required=True)
    p.add_argument("--den-vcf", required=True)
    p.add_argument("--nea-vcf", required=True)
    p.add_argument("--bed", required=True)
    p.add_argument("--ancestral-fasta", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--window-snps", type=int, default=100)
    p.add_argument("--step-snps", type=int, default=50)
    p.add_argument("--window-jackknife-block-snps", type=int, default=10)
    p.add_argument("--region-jackknife-bp", type=int, default=50000)
    p.add_argument(
        "--region-only", action="store_true",
        help=("compute whole-region A/B/C D and fD with physical-block "
              "jackknife, but do not require or emit SNP windows"),
    )
    return p.parse_args()


def read_bed(path: str) -> tuple[str, int, int]:
    with open(path) as fh:
        fields = fh.readline().split()
    if len(fields) < 3:
        raise ValueError("BED requires at least three columns")
    chrom, start0, end0 = fields[0], int(fields[1]), int(fields[2])
    if end0 <= start0:
        raise ValueError("BED end must exceed start")
    return chrom, start0, end0


def read_yri(fam: str) -> list[str]:
    samples = []
    with open(fam) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.split()
            if len(f) >= 2 and f[0] == "YRI":
                samples.append(f[1])
    if not samples:
        raise RuntimeError("No YRI samples found in FAM")
    return samples


def alt_list(record) -> list[str]:
    return [str(x).upper() for x in (record.ALT or []) if x not in (None, ".")]


def target_af(record, target_alt: str,
              sample_indices: list[int] | None = None) -> tuple[float, int]:
    alts = alt_list(record)
    if target_alt not in alts:
        return math.nan, 0
    target_index = alts.index(target_alt) + 1
    genotypes = record.genotypes
    if sample_indices is None:
        sample_indices = range(len(genotypes))
    dosage = 0
    called = 0
    for i in sample_indices:
        for allele in genotypes[i][:-1]:
            if allele is None or allele < 0:
                continue
            called += 1
            dosage += int(allele == target_index)
    return ((dosage / called, called) if called else (math.nan, 0))


def invariant_ref_af(record) -> tuple[float, int]:
    called = 0
    for gt in record.genotypes:
        for allele in gt[:-1]:
            if allele is None or allele < 0:
                continue
            called += 1
            if allele != 0:
                return math.nan, 0
    return ((0.0, called) if called else (math.nan, 0))


def scalar_format(record, field: str) -> float:
    try:
        value = record.format(field)
        if value is None or len(value) == 0:
            return math.nan
        x = np.asarray(value).reshape(-1)[0]
        return float(x) if np.isfinite(x) else math.nan
    except Exception:
        return math.nan


def load_ancestral_region(path: str, start1: int, end1: int) -> tuple[str, str]:
    opener = gzip.open if path.endswith(".gz") else open
    header = None
    chunks: list[str] = []
    cursor1 = 1
    with opener(path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError("Expected a single-sequence ancestral FASTA")
                header = line[1:].split()[0]
                continue
            line = line.upper()
            line_start = cursor1
            line_end = cursor1 + len(line) - 1
            ov_start = max(start1, line_start)
            ov_end = min(end1, line_end)
            if ov_start <= ov_end:
                chunks.append(line[ov_start-line_start:ov_end-line_start+1])
            cursor1 = line_end + 1
            if cursor1 > end1:
                break
    seq = "".join(chunks)
    expected = end1 - start1 + 1
    if len(seq) != expected:
        raise ValueError(f"Ancestral interval length {len(seq)} != {expected}")
    return header or "unknown", seq


def collect_archaic(vcf_path: str, region: str, candidates: set[int]) -> tuple[dict, list[str]]:
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)
    records: dict[int, list] = {}
    for rec in vcf(region):
        if rec.POS in candidates:
            records.setdefault(rec.POS, []).append(rec)
    vcf.close()
    return records, samples


def match_archaic(records: list | None, ref: str, alt: str) -> tuple[float, int, str, float, float]:
    if not records:
        return math.nan, 0, "no_record", math.nan, math.nan
    exact = [r for r in records if r.REF.upper() == ref and alt_list(r) == [alt]]
    invariant = [r for r in records if r.REF.upper() == ref and not alt_list(r)]
    if len(exact) == 1:
        rec = exact[0]
        af, nchr = target_af(rec, alt)
        status = "exact_alt_match"
    elif len(exact) == 0 and len(invariant) == 1:
        rec = invariant[0]
        af, nchr = invariant_ref_af(rec)
        status = "called_invariant_ref"
    else:
        return math.nan, 0, "allele_mismatch_or_duplicate", math.nan, math.nan
    return af, nchr, status, scalar_format(rec, "DP"), scalar_format(rec, "GQ")


def site_terms(p1: float, p2: float, p3: float, p4: float) -> tuple[float, float, float, float]:
    abba = ((1-p1)*p2*p3*(1-p4)) + (p1*(1-p2)*(1-p3)*p4)
    baba = (p1*(1-p2)*p3*(1-p4)) + ((1-p1)*p2*(1-p3)*p4)
    numerator = abba - baba
    pdyn = max(p2, p3)
    fd_den = (pdyn - p1) * (pdyn - p4)
    return abba, baba, numerator, fd_den


def ratio(frame: pd.DataFrame, numerator: str, denominator: str) -> float:
    den = float(frame[denominator].sum())
    return float(frame[numerator].sum()) / den if den != 0 else math.nan


def inference(theta: float, estimates: list[float]) -> dict[str, float]:
    b = len(estimates)
    if b < 3 or not np.isfinite(theta):
        return {"estimate": theta, "se": math.nan, "z": math.nan,
                "p_two_sided": math.nan, "minus_log10_p": math.nan,
                "n_blocks": b, "ci95_low": math.nan, "ci95_high": math.nan}
    arr = np.asarray(estimates, dtype=float)
    mean_leave = arr.mean()
    se = math.sqrt((b - 1) / b * np.sum((arr - mean_leave) ** 2))
    z = theta / se if se > 0 else math.nan
    if np.isfinite(z):
        logp = math.log(2.0) + float(norm.logsf(abs(z)))
        pval = math.exp(max(logp, math.log(np.finfo(float).tiny)))
        minus_log10 = -logp / math.log(10.0)
    else:
        pval = minus_log10 = math.nan
    return {"estimate": theta, "se": se, "z": z,
            "p_two_sided": pval, "minus_log10_p": minus_log10,
            "n_blocks": b, "ci95_low": theta - 1.96*se,
            "ci95_high": theta + 1.96*se}


def physical_jackknife(frame: pd.DataFrame, start1: int, block_bp: int,
                       numerator: str, denominator: str) -> dict[str, float]:
    theta = ratio(frame, numerator, denominator)
    blocks = ((frame["pos"] - start1) // block_bp).astype(int).to_numpy()
    estimates = []
    for block in np.unique(blocks):
        keep = frame.iloc[blocks != block]
        value = ratio(keep, numerator, denominator)
        if np.isfinite(value):
            estimates.append(value)
    return inference(theta, estimates)


def site_block_jackknife(frame: pd.DataFrame, block_snps: int,
                         numerator: str, denominator: str) -> dict[str, float]:
    theta = ratio(frame, numerator, denominator)
    blocks = np.arange(len(frame), dtype=int) // block_snps
    estimates = []
    for block in np.unique(blocks):
        keep = frame.iloc[blocks != block]
        value = ratio(keep, numerator, denominator)
        if np.isfinite(value):
            estimates.append(value)
    return inference(theta, estimates)


def bh(pvalues: pd.Series) -> np.ndarray:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    q = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if not len(valid):
        return q
    order = valid[np.argsort(p[valid])]
    adjusted = p[order] * len(valid) / np.arange(1, len(valid) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q[order] = np.minimum(adjusted, 1.0)
    return q


def add_terms(frame: pd.DataFrame, code: str, groups: tuple[str, str, str, str]) -> None:
    def get(r, group: str) -> float:
        return 0.0 if group == "anc" else float(r[f"derived_{group}"])
    terms = frame.apply(lambda r: site_terms(*(get(r, g) for g in groups)), axis=1)
    frame[[f"{code}_abba", f"{code}_baba", f"{code}_num", f"{code}_fd_den"]] = pd.DataFrame(
        terms.tolist(), index=frame.index)
    frame[f"{code}_d_den"] = frame[f"{code}_abba"] + frame[f"{code}_baba"]


def main() -> None:
    args = parse_args()
    if min(args.window_snps, args.step_snps, args.window_jackknife_block_snps,
           args.region_jackknife_bp) <= 0:
        raise ValueError("Window, step and jackknife sizes must be positive")
    if math.ceil(args.window_snps / args.window_jackknife_block_snps) < 3:
        raise ValueError("Window jackknife needs at least three blocks")
    os.makedirs(args.outdir, exist_ok=True)
    chrom, start0, end0 = read_bed(args.bed)
    start1, end1 = start0 + 1, end0
    region = f"{chrom}:{start1}-{end1}"
    yri_samples = read_yri(args.fam)
    qc = Counter()

    tbn_vcf = VCF(args.tbn_vcf)
    tbn_samples = list(tbn_vcf.samples)
    candidates: dict[int, dict] = {}
    duplicate_tbn = set()
    for rec in tbn_vcf(region):
        qc["tbn_records_in_bed"] += 1
        alts = alt_list(rec)
        if len(rec.REF) != 1 or len(alts) != 1 or len(alts[0]) != 1:
            qc["drop_tbn_not_biallelic_snp"] += 1
            continue
        af, nchr = target_af(rec, alts[0])
        if not np.isfinite(af) or not nchr:
            qc["drop_tbn_no_called_genotype"] += 1
            continue
        if rec.POS in candidates:
            duplicate_tbn.add(rec.POS)
            candidates.pop(rec.POS, None)
            continue
        candidates[rec.POS] = {"chrom": chrom, "pos": rec.POS,
                               "ref": rec.REF.upper(), "alt": alts[0],
                               "p_tbn": af, "n_tbn_chr": nchr}
    tbn_vcf.close()
    qc["drop_tbn_duplicate_position"] = len(duplicate_tbn)
    qc["tbn_candidate_sites"] = len(candidates)

    kg = VCF(args.kg_vcf, samples=yri_samples)
    kg_sample_to_i = {s: i for i, s in enumerate(kg.samples)}
    missing_yri = [s for s in yri_samples if s not in kg_sample_to_i]
    if missing_yri:
        raise RuntimeError(f"YRI absent from 1KGP VCF: {missing_yri[:10]}")
    yri_i = [kg_sample_to_i[s] for s in yri_samples]
    kg_by_pos = {}
    kg_dups = set()
    for rec in kg(region):
        if rec.POS not in candidates:
            continue
        if rec.POS in kg_by_pos:
            kg_dups.add(rec.POS)
        else:
            kg_by_pos[rec.POS] = rec

    den_by_pos, den_samples = collect_archaic(args.den_vcf, region, set(candidates))
    nea_by_pos, nea_samples = collect_archaic(args.nea_vcf, region, set(candidates))
    rows = []
    for pos in sorted(candidates):
        c = candidates[pos]
        ref, alt = c["ref"], c["alt"]
        if pos in kg_dups:
            qc["drop_kg_duplicate_position"] += 1
            continue
        kg_rec = kg_by_pos.get(pos)
        if kg_rec is None:
            p_yri, n_yri = 0.0, 2 * len(yri_samples)
            kg_status = "absent_assumed_all_ref"
        else:
            if kg_rec.REF.upper() != ref or alt_list(kg_rec) != [alt]:
                qc["drop_kg_allele_mismatch"] += 1
                continue
            p_yri, n_yri = target_af(kg_rec, alt, yri_i)
            if not np.isfinite(p_yri):
                qc["drop_kg_missing_yri"] += 1
                continue
            kg_status = "exact_alt_match"
        qc[f"kg_{kg_status}"] += 1

        den = match_archaic(den_by_pos.get(pos), ref, alt)
        nea = match_archaic(nea_by_pos.get(pos), ref, alt)
        if not np.isfinite(den[0]) or not den[1]:
            qc[f"drop_den_{den[2]}"] += 1
            continue
        if not np.isfinite(nea[0]) or not nea[1]:
            qc[f"drop_nea_{nea[2]}"] += 1
            continue
        qc[f"den_{den[2]}"] += 1
        qc[f"nea_{nea[2]}"] += 1
        rows.append({**c, "p_yri": p_yri, "n_yri_chr": n_yri,
                     "p_den": den[0], "n_den_chr": den[1], "den_status": den[2],
                     "den_dp": den[3], "den_gq": den[4],
                     "p_nea": nea[0], "n_nea_chr": nea[1], "nea_status": nea[2],
                     "nea_dp": nea[3], "nea_gq": nea[4], "kg_status": kg_status})
    kg.close()
    raw = pd.DataFrame(rows).sort_values("pos")
    if raw.empty:
        raise RuntimeError("No common four-population sites remain")
    qc["common_called_before_ancestral"] = len(raw)

    anc_header, anc_seq = load_ancestral_region(args.ancestral_fasta, start1, end1)
    raw["ancestral"] = [anc_seq[int(p)-start1] for p in raw["pos"]]
    raw["ancestor_status"] = np.where(raw.ancestral.eq(raw.ref), "REF",
                              np.where(raw.ancestral.eq(raw.alt), "ALT", "UNRESOLVED"))
    qc["ancestor_ref"] = int(raw.ancestor_status.eq("REF").sum())
    qc["ancestor_alt"] = int(raw.ancestor_status.eq("ALT").sum())
    qc["ancestor_unresolved"] = int(raw.ancestor_status.eq("UNRESOLVED").sum())
    raw.to_csv(os.path.join(args.outdir, "sites_all_common_TBN_candidates.tsv.gz"),
               sep="\t", index=False, compression="gzip", na_rep="NA")

    pol = raw[raw.ancestor_status.isin(["REF", "ALT"])].copy()
    flip = pol.ancestor_status.eq("ALT")
    for group in ["tbn", "yri", "den", "nea"]:
        pol[f"derived_{group}"] = np.where(flip, 1-pol[f"p_{group}"], pol[f"p_{group}"])

    # Require biallelic information for every A/B/C comparison so all tests use
    # exactly the same SNPs and official B/C validation has identical windows.
    a = pol[["derived_den", "derived_nea", "derived_tbn", "derived_yri"]]
    valid_a = (a.max(axis=1) > 0) & (a.min(axis=1) < 1)
    valid_b = pol[["derived_yri", "derived_tbn", "derived_den"]].max(axis=1) > 0
    valid_c = pol[["derived_yri", "derived_tbn", "derived_nea"]].max(axis=1) > 0
    valid = valid_a & valid_b & valid_c
    qc["drop_not_biallelic_in_all_tests"] = int((~valid).sum())
    pol = pol.loc[valid].copy().reset_index(drop=True)
    qc["final_common_polarized_sites"] = len(pol)
    if len(pol) < args.window_snps and not args.region_only:
        raise RuntimeError(f"Only {len(pol)} sites; fewer than one complete window")

    for code, spec in COMPARISONS.items():
        add_terms(pol, code, spec["groups"])
    pol.to_csv(os.path.join(args.outdir, "sites_common_ancestral_polarized.tsv.gz"),
               sep="\t", index=False, compression="gzip", na_rep="NA")

    summary_rows = []
    for code, spec in COMPARISONS.items():
        d = physical_jackknife(pol, start1, args.region_jackknife_bp,
                               f"{code}_num", f"{code}_d_den")
        fd = physical_jackknife(pol, start1, args.region_jackknife_bp,
                                f"{code}_num", f"{code}_fd_den")
        interpretation = spec["positive"] if d["estimate"] > 0 else spec["negative"]
        summary_rows.append({"comparison": code, "quartet": spec["label"],
                             "n_sites": len(pol), "interpretation": interpretation,
                             "fd_direction_valid": bool(d["estimate"] >= 0),
                             **{f"D_{k}": v for k, v in d.items()},
                             **{f"fd_{k}": v for k, v in fd.items()}})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(args.outdir, "region_summary_ABC.tsv"),
                   sep="\t", index=False, na_rep="NA")

    if args.region_only:
        qc["complete_windows"] = 0
        qc["analysis_mode_region_only"] = 1
    else:
        windows = []
        for i in range(0, len(pol)-args.window_snps+1, args.step_snps):
            w = pol.iloc[i:i+args.window_snps]
            row = {"chrom": chrom, "window_index": len(windows)+1,
                   "first_site_index": i+1, "last_site_index": i+len(w),
                   "start": int(w.pos.iloc[0]), "end": int(w.pos.iloc[-1]),
                   "midpoint": (int(w.pos.iloc[0])+int(w.pos.iloc[-1]))/2,
                   "n_snps": len(w),
                   "window_jackknife_block_snps": args.window_jackknife_block_snps}
            for code in COMPARISONS:
                d = site_block_jackknife(w, args.window_jackknife_block_snps,
                                         f"{code}_num", f"{code}_d_den")
                fd = site_block_jackknife(w, args.window_jackknife_block_snps,
                                          f"{code}_num", f"{code}_fd_den")
                row.update({f"{code}_sum_ABBA": float(w[f"{code}_abba"].sum()),
                            f"{code}_sum_BABA": float(w[f"{code}_baba"].sum()),
                            **{f"{code}_D_{k}": v for k, v in d.items()},
                            **{f"{code}_fd_{k}": v for k, v in fd.items()},
                            f"{code}_fd_direction_valid": bool(d["estimate"] >= 0)})
            windows.append(row)
        win = pd.DataFrame(windows)
        for code in COMPARISONS:
            win[f"{code}_D_q_BH"] = bh(win[f"{code}_D_p_two_sided"])
            win[f"{code}_fd_q_BH_exploratory"] = bh(win[f"{code}_fd_p_two_sided"])
        window_name = f"windows_{args.window_snps}SNP_step{args.step_snps}.tsv"
        win.to_csv(os.path.join(args.outdir, window_name), sep="\t", index=False,
                   na_rep="NA")
        qc["complete_windows"] = len(win)
        qc["analysis_mode_region_only"] = 0

    sensitivity = []
    filters = {
        "all_common_sites": np.ones(len(pol), dtype=bool),
        "both_archaics_DP_ge_5": (pol.den_dp.ge(5) & pol.nea_dp.ge(5)).to_numpy(),
        "both_archaics_GQ_ge_30": (pol.den_gq.ge(30) & pol.nea_gq.ge(30)).to_numpy(),
        "both_archaics_DP_ge_5_GQ_ge_30": (
            pol.den_dp.ge(5) & pol.nea_dp.ge(5) & pol.den_gq.ge(30) & pol.nea_gq.ge(30)).to_numpy(),
        "KG_exact_alt_only": pol.kg_status.eq("exact_alt_match").to_numpy(),
    }
    for filter_name, mask in filters.items():
        sub = pol.loc[mask]
        for code, spec in COMPARISONS.items():
            d = physical_jackknife(sub, start1, args.region_jackknife_bp,
                                   f"{code}_num", f"{code}_d_den")
            fd = physical_jackknife(sub, start1, args.region_jackknife_bp,
                                    f"{code}_num", f"{code}_fd_den")
            sensitivity.append({"filter": filter_name, "comparison": code,
                                "quartet": spec["label"], "n_sites": len(sub),
                                **{f"D_{k}": v for k, v in d.items()},
                                **{f"fd_{k}": v for k, v in fd.items()}})
    pd.DataFrame(sensitivity).to_csv(os.path.join(args.outdir, "quality_sensitivity_ABC.tsv"),
                                     sep="\t", index=False, na_rep="NA")

    with open(os.path.join(args.outdir, "site_QC.tsv"), "w") as out:
        out.write("metric\tvalue\n")
        for key in sorted(qc):
            out.write(f"{key}\t{qc[key]}\n")
        out.write(f"region\t{region}\nancestral_fasta_header\t{anc_header}\n")
        out.write(f"tbn_samples\t{len(tbn_samples)}\nyri_samples\t{len(yri_samples)}\n")
        out.write(f"den_samples\t{','.join(den_samples)}\nnea_samples\t{','.join(nea_samples)}\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
