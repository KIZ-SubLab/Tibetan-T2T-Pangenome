#!/usr/bin/env python3
"""Matched, without-replacement SV-set enrichment against independent regions."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

SV_RE = re.compile(r"^(chr[^-]+)-(\d+)-(\d+)-(DEL|INS)-(\d+)$")
LENGTH_BINS = [0, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000, np.inf]
MAF_BINS = [0, .01, .05, .1, .2, .3, .4, .5]


def save(frame, path):
    frame.to_csv(path, sep="\t", index=False)


def parse_sv(ids, mode):
    rows = []
    for value in ids:
        m = SV_RE.fullmatch(str(value).strip())
        if not m:
            raise ValueError(f"Unrecognized SV ID: {value}")
        chrom, start, end, typ, length = m.groups()
        start, end, length = int(start), int(end), int(length)
        if start < 1 or end < start or length < 1:
            raise ValueError(f"Invalid coordinates/length: {value}")
        # Keep the earlier population pipeline's 1-based inclusive ID convention.
        start0 = start - 1
        end0 = start if typ == "INS" and mode == "point" else max(start, end)
        rows.append((str(value).strip(), chrom, start0, end0, typ, length))
    return pd.DataFrame(rows, columns=["SV_ID", "chrom", "start", "end", "SVTYPE", "length"])


def studies(args, out):
    selected = pd.read_csv(args.selected)
    col = next(c for c in ("svID_5caller", "SV_ID", "KEY") if c in selected)
    selected = selected.loc[pd.to_numeric(selected.PBS, errors="coerce").notna()].copy()
    selected_ids = selected[col].astype(str).str.strip()
    if len(selected) != args.expected_selected or selected_ids.duplicated().any():
        raise ValueError("Selected count/uniqueness failed; VNTRs without PBS are excluded")
    universe = pd.read_csv(args.universe)
    ucol = next(c for c in ("KEY", "SV_ID", "svID_5caller") if c in universe)
    sv = parse_sv(universe[ucol], args.ins_mode)
    if sv.SV_ID.duplicated().any():
        raise ValueError("Duplicate universe SV IDs; resolve before sampling")
    missing = selected_ids.loc[~selected_ids.isin(sv.SV_ID)]
    save(pd.DataFrame({"missing_selected_SV_ID": missing}), out / "missing_selected.tsv")
    if len(missing):
        raise ValueError("Selected IDs absent from universe. No nearest-length/fuzzy matching allowed.")
    if args.maf_column:
        maf = pd.to_numeric(universe[args.maf_column], errors="coerce")
        definition = f"user-supplied universe column {args.maf_column}"
    else:
        af = universe[["AF_TIB", "AF_HAN"]].apply(pd.to_numeric, errors="coerce")
        if af.isna().any().any() or ((af < 0) | (af > 1)).any().any():
            raise ValueError("Invalid AF; missing frequencies cannot be assigned to a bin")
        avg = af.AF_TIB * args.tib_af_weight + af.AF_HAN * (1 - args.tib_af_weight)
        maf = np.minimum(avg, 1 - avg)
        definition = f"AF-based MAF proxy: min(w*AF_TIB+(1-w)*AF_HAN,1-AF); w={args.tib_af_weight}"
    if maf.isna().any() or ((maf < 0) | (maf > .5)).any():
        raise ValueError("MAF must be finite and in [0,0.5] for every universe SV")
    sv["MAF"] = np.asarray(maf)
    sv["is_TESV"] = sv.SV_ID.isin(selected_ids)
    sv["interval_length"] = sv.end - sv.start
    sv["length_bin"] = pd.cut(sv.length, LENGTH_BINS, labels=False, include_lowest=True)
    sv["MAF_bin"] = pd.cut(sv.MAF, MAF_BINS, labels=False, include_lowest=True)
    keys = ["SVTYPE", "length_bin", "MAF_bin"]
    if args.match_chrom:
        keys.append("chrom")
    if args.match_span:
        sv["span_bin"] = pd.cut(sv.interval_length, [0, 1, 50, 200, 1000, 5000, np.inf], labels=False)
        keys.append("span_bin")
    sv["stratum"] = sv[keys].astype(str).agg("|".join, axis=1)
    save(pd.DataFrame({"selected_SV_ID": selected_ids, "universe_SV_ID": selected_ids,
                       "match_method": "exact_id"}), out / "selected_to_universe.tsv")
    save(sv, out / "universe_covariates.tsv")
    return sv, definition


def overlap(svs, bed):
    """At least one shared bp in 0-based half-open intervals, counted once per SV."""
    hits = np.zeros(len(svs), dtype=bool)
    for chrom, group in svs.groupby("chrom", sort=False):
        reg = bed.loc[bed.chrom.eq(chrom)].sort_values("start")
        if reg.empty:
            continue
        starts = reg.start.to_numpy(dtype=np.int64)
        max_ends = np.maximum.accumulate(reg.end.to_numpy(dtype=np.int64))
        idx = np.searchsorted(starts, group.end.to_numpy(), side="left") - 1
        valid = idx >= 0
        found = np.zeros(len(group), dtype=bool)
        found[valid] = max_ends[idx[valid]] > group.start.to_numpy()[valid]
        hits[group.index.to_numpy()] = found
    return hits


def annotations(sv, manifest, out):
    meta = pd.read_csv(manifest, sep="\t").fillna("")
    required = {"endpoint", "bed", "status", "family", "scope", "assembly", "definition", "padding_bp"}
    if not required.issubset(meta):
        raise ValueError(f"Endpoint manifest needs {required}")
    if meta.endpoint.duplicated().any() or not meta.family.isin(["primary", "secondary"]).all():
        raise ValueError("Invalid endpoint names/families")
    columns = []
    for row in meta.itertuples(index=False):
        if row.status != "READY":
            continue
        if row.scope != "genomewide_independent_of_TESVs" or row.assembly != "CHM13v2.0":
            raise ValueError(f"Unusable region scope/assembly for {row.endpoint}")
        path = Path(row.bed)
        if not path.is_absolute():
            path = Path(manifest).parent / path
        if path.stat().st_size:
            bed = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                              names=["chrom", "start", "end"], dtype={"chrom": str})
            for c in ("start", "end"):
                nums = pd.to_numeric(bed[c], errors="raise")
                if nums.isna().any() or (nums != np.floor(nums)).any():
                    raise ValueError(f"Noninteger BED coordinates: {path}")
                bed[c] = nums.astype(np.int64)
            if (bed.start < 0).any() or (bed.end <= bed.start).any():
                raise ValueError(f"Invalid BED interval: {path}")
            pad = int(row.padding_bp)
            if pad < 0:
                raise ValueError("Negative padding")
            # Positive distance tests include exact endpoint gaps of padding_bp.
            # Direct overlap (pad=0) still requires at least one shared bp.
            expansion = pad + 1 if pad else 0
            bed.start = (bed.start - expansion).clip(lower=0)
            bed.end += expansion
        else:
            bed = pd.DataFrame(columns=["chrom", "start", "end"])
        sv[row.endpoint] = overlap(sv, bed)
        columns.append(row.endpoint)
    save(meta, out / "endpoint_status.tsv")
    save(sv, out / "all_SV_feature_annotation.tsv.gz")
    if not columns:
        raise ValueError("No READY genomewide endpoints")
    return columns, meta


def matching(sv, out):
    groups, records = [], []
    for key, selected in sv.loc[sv.is_TESV].groupby("stratum", sort=True):
        controls = sv.index[~sv.is_TESV & sv.stratum.eq(key)].to_numpy()
        n = len(selected)
        records.append((key, n, len(controls), len(controls) / n, len(controls) >= n))
        groups.append((controls, n))
    qc = pd.DataFrame(records, columns=["stratum", "selected_n", "background_n", "controls_per_TESV", "feasible"])
    save(qc, out / "matching_strata.tsv")
    if not qc.feasible.all():
        raise ValueError("Insufficient controls in exact strata; see matching_strata.tsv. No silent relaxation/replacement.")
    return groups


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return result


def simulate(sv, columns, groups, n_perm, seed):
    rng = np.random.default_rng(seed)
    matrix = sv[columns].to_numpy(dtype=np.int32)
    cov = np.column_stack([np.log2(sv.length), sv.MAF, np.log2(sv.interval_length)])
    null = np.zeros((n_perm, len(columns)), dtype=np.int32)
    balance = np.zeros((n_perm, 3))
    first = None
    for b in range(n_perm):
        idx = np.concatenate([rng.choice(pool, n, replace=False) for pool, n in groups])
        if b == 0:
            first = idx
            assert len(np.unique(idx)) == len(idx) and not sv.loc[idx, "is_TESV"].any()
        null[b] = matrix[idx].sum(axis=0)
        balance[b] = cov[idx].mean(axis=0)
    return null, balance, first




def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selected", type=Path, required=True)
    p.add_argument("--universe", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--expected-selected", type=int, default=326)
    p.add_argument("--n-perm", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260928)
    p.add_argument("--ins-mode", choices=["span", "point"], default="span")
    p.add_argument("--maf-column", default="")
    p.add_argument("--tib-af-weight", type=float, default=.5)
    p.add_argument("--match-chrom", action="store_true")
    p.add_argument("--match-span", action="store_true")
    p.add_argument("--require-complete", action="store_true")
    args = p.parse_args()
    if args.n_perm < 1 or not 0 <= args.tib_af_weight <= 1:
        p.error("Invalid permutations or AF weight")
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    sv, maf_definition = studies(args, out)
    groups = matching(sv, out)
    columns, meta = annotations(sv, args.manifest, out)
    null, balance, first = simulate(sv, columns, groups, args.n_perm, args.seed)
    null = pd.DataFrame(null, columns=columns)
    save(null, out / "permutation_counts.tsv.gz")
    save(sv.loc[first], out / "first_matched_background_set.tsv")
    result = meta.copy()
    n = int(sv.is_TESV.sum())
    for c in columns:
        obs = int(sv.loc[sv.is_TESV, c].sum())
        mean = float(null[c].mean())
        stats = dict(selected_n=n, background_n=int((~sv.is_TESV).sum()), observed_SVs=obs,
                     observed_fraction=obs / n, null_mean=mean, null_fraction=mean / n,
                     all_background_overlap_n=int(sv.loc[~sv.is_TESV, c].sum()),
                     all_background_overlap_fraction=float(sv.loc[~sv.is_TESV, c].mean()),
                     null_2_5_percentile=float(null[c].quantile(.025)),
                     null_97_5_percentile=float(null[c].quantile(.975)),
                     fold_enrichment=obs / mean if mean else np.nan,
                     P_enrichment=(1 + int((null[c] >= obs).sum())) / (args.n_perm + 1),
                     P_depletion=(1 + int((null[c] <= obs).sum())) / (args.n_perm + 1),
                     n_perm=args.n_perm)
        for key, value in stats.items():
            result.loc[result.endpoint.eq(c), key] = value
    # An anchor's ascertainment and an effect conditional on assay evaluability
    # are distinct questions; export both labels, without adding post-hoc tests.
    # Missing planned tests count as P=1 in family correction, but stay NA in output.
    for _, ix in result.groupby("family").groups.items():
        result.loc[ix, "FDR_enrichment"] = bh(result.loc[ix, "P_enrichment"].fillna(1))
    result.loc[~result.status.eq("READY"), "FDR_enrichment"] = np.nan
    save(result, out / "enrichment_summary.tsv")
    cov = np.column_stack([np.log2(sv.length), sv.MAF, np.log2(sv.interval_length)])
    qc = []
    for j, name in enumerate(["log2_SV_length", "MAF", "log2_interval_length"]):
        s = cov[sv.is_TESV, j]
        b = cov[~sv.is_TESV, j]
        sd = np.sqrt((s.var() + b.var()) / 2)
        smd = (s.mean() - balance[:, j].mean()) / sd if sd else 0
        qc.append((name, s.mean(), b.mean(), balance[:, j].mean(),
                   *np.quantile(balance[:, j], [.025, .975]), smd, abs(smd) > .1))
    save(pd.DataFrame(qc, columns=["covariate", "TESV_mean", "all_background_mean", "matched_mean",
                                  "null_mean_low", "null_mean_high", "standardized_mean_difference",
                                  "balance_warning_abs_SMD_gt_0_1"]), out / "matching_balance.tsv")
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    config.update(MAF_definition=maf_definition, universe_n=len(sv), selected_n=n,
                  background_n=int((~sv.is_TESV).sum()), length_bins=[str(x) for x in LENGTH_BINS], maf_bins=MAF_BINS)
    (out / "configuration.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    state = "COMPLETE" if meta.status.eq("READY").all() else "PARTIAL"
    report = ["# TESV matched-background enrichment", f"\nStatus: {state}",
              f"\nUniverse: {len(sv)}; TESVs: {n}; non-TESV background: {len(sv)-n}.",
              f"\nMAF: {maf_definition}.",
              f"\n{args.n_perm} independent matched sets, without replacement within each set; seed={args.seed}.",
              "\nSV type, prespecified SV-length bin and MAF bin are exactly matched. Continuous values are not exactly matched; inspect matching_balance.tsv.",
              "\nOne-sided enrichment P=(1+null counts >= observed)/(B+1). BH is calculated within each manifest family; missing tests contribute P=1 to correction.",
              "\nEach SV contributes at most one hit per endpoint. Empty available annotations are zero hits; unavailable endpoints are NA.",
              "\nThis tests association of selection status with annotated regions, not causal SV effects. Hi-C effects are descriptive with 2+2 donors, not FDR-significant differential loops.",
              "\nCheck INS span versus point and chromosome/span matching sensitivity; unmatched genomic confounders (TSS, repeats, coverage) remain."]
    (out / "report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"[INFO] {state}: {out / 'enrichment_summary.tsv'}", flush=True)
    if args.require_complete and state != "COMPLETE":
        raise SystemExit("PARTIAL: missing endpoints; available results were saved")


if __name__ == "__main__":
    main()
