#!/usr/bin/env python3
"""Build genomewide functional regions WITHOUT conditioning on selected SVs."""
import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = [f"{person}-{time}" for person in ["A-1", "A-3", "W-1", "W-3"] for time in [1, 4]]
ENDPOINTS = [
    ("ATAC_open_direct", "primary", "atac_open", 0),
    ("ATAC_population_DOR_direct", "primary", "atac_dor", 0),
    ("HiC_population_effect_direct", "primary", "hic_effect", 0),
    ("cCRE_hypoxia_gene_direct", "primary", "ccre_hypoxia", 0),
    ("ATAC_open_within_1kb", "secondary", "atac_open", 1000),
    ("ATAC_population_DOR_within_1kb", "secondary", "atac_dor", 1000),
    ("HiC_any_anchor_direct", "secondary", "hic_any", 0),
    ("HiC_evaluable_anchor_direct", "secondary", "hic_valid", 0),
    ("cCRE_any_direct", "secondary", "ccre_all", 0),
]


def require(path):
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        raise FileNotFoundError(f"Missing input: {path}")
    return Path(path)


def table(path):
    return pd.read_csv(require(path), sep="\t")


def reject_sv_conditioning(frame):
    forbidden = {"SV_ID", "TESV", "selected_SV_ID", "svID_5caller", "PBS"}
    if forbidden.intersection(frame.columns):
        raise ValueError("SV-filtered table rejected: supply full genomewide input")


def bed_write(frame, path):
    frame = frame[["chrom", "start", "end"]].drop_duplicates().sort_values(["chrom", "start", "end"])
    frame.to_csv(path, sep="\t", index=False, header=False)


def atac_regions(args, out):
    require(args.atac_rpkm)
    parts = []
    for chunk in pd.read_csv(args.atac_rpkm, sep="\t", chunksize=30000):
        samples = [c for c in chunk if re.fullmatch(r"[AW]-\d+-[1-5]", c)]
        if len(samples) != 50:
            raise ValueError(f"Expected 50 ATAC samples, found {len(samples)}")
        vals = chunk[samples].apply(pd.to_numeric, errors="raise")
        if not np.isfinite(vals.to_numpy()).all() or (vals < 0).any().any():
            raise ValueError("Missing/invalid RPKM values")
        chunk = chunk.loc[vals.max(axis=1) > args.peak_cutoff].copy()
        chrom = next(c for c in ["chr", "chrom"] if c in chunk)
        fid = next(c for c in ["peak_id", "feature_id"] if c in chunk)
        parts.append(chunk[[chrom, "start", "end", fid]].rename(columns={chrom: "chrom", fid: "feature_id"}))
    peaks = pd.concat(parts, ignore_index=True)
    if peaks.feature_id.duplicated().any():
        raise ValueError("Nonunique ATAC peak IDs")
    bed_write(peaks, out / "atac_open.bed")
    return peaks


def dor_regions(args, peaks, out):
    de = table(args.atac_population)
    reject_sv_conditioning(de)
    if de.feature_id.duplicated().any():
        raise ValueError("Duplicate differential feature IDs")
    mask = (de["adj.P.Val"] < args.fdr) & (de.logFC.abs() >= args.min_logfc)
    merged = peaks.merge(de[["feature_id", "logFC", "adj.P.Val"]], on="feature_id", how="inner", validate="one_to_one")
    if len(merged) == 0:
        raise ValueError("No common RPKM and differential peak IDs")
    # Only peaks in the genomewide tested feature universe can be DORs.
    dor = peaks.loc[peaks.feature_id.isin(de.loc[mask, "feature_id"])]
    bed_write(dor, out / "atac_dor.bed")
    pd.DataFrame([{"open_peaks": len(peaks), "tested_open_peaks": len(merged), "DOR_peaks": len(dor)}]).to_csv(out / "atac_qc.tsv", sep="\t", index=False)


def load_loops(path):
    loops = table(path)
    reject_sv_conditioning(loops)
    needed = {"interaction_id", "chrom", "anchor1_start", "anchor1_end", "anchor2_start", "anchor2_end"}
    if not needed.issubset(loops) or loops.interaction_id.duplicated().any():
        raise ValueError("Need unique genomewide interaction IDs and both anchor coordinates")
    for i in (1, 2):
        if ((loops[f"anchor{i}_end"] - loops[f"anchor{i}_start"]) != 20000).any() or (loops[f"anchor{i}_start"] % 20000 != 0).any():
            raise ValueError("Expected 20-kb bin-aligned Hi-C anchors")
    return loops


def anchors(loops):
    return pd.concat([loops[["chrom", f"anchor{i}_start", f"anchor{i}_end"]].rename(
        columns={f"anchor{i}_start": "start", f"anchor{i}_end": "end"}) for i in (1, 2)], ignore_index=True)


def quantify_hic(loops, metadata, out, block_rows):
    """Read sparse row blocks, preserving the existing balanced single-pixel O/E."""
    import cooler
    meta = table(metadata)
    if meta["sample"].duplicated().any() or set(meta["sample"]) != set(SAMPLES):
        raise ValueError(f"Hi-C metadata must contain exactly {SAMPLES}")
    result = loops.copy()
    for row in meta.itertuples(index=False):
        print(f"[INFO] genomewide Hi-C O/E: {row.sample}", flush=True)
        path = str(row.mcool)
        require(path.split("::")[0])
        clr = cooler.Cooler(path if "::" in path else path + "::/resolutions/20000")
        if clr.binsize != 20000:
            raise ValueError("Incorrect Hi-C resolution")
        exp = table(row.expected_cis)
        if exp.duplicated(["region1", "dist"]).any():
            raise ValueError("Expected-cis must have unique region1/dist keys")
        expected = {(str(c), int(d)): float(v) for c, d, v in
                    exp[["region1", "dist", "balanced.avg.smoothed"]].itertuples(index=False, name=None)
                    if pd.notna(d) and np.isfinite(v) and v > 0}
        values = np.full(len(loops), np.nan)
        for chrom, group in loops.groupby("chrom", sort=False):
            if chrom not in clr.chromnames:
                continue
            offset, end = clr.extent(chrom)
            ii = group.anchor1_start.to_numpy(dtype=np.int64) // 20000
            jj = group.anchor2_start.to_numpy(dtype=np.int64) // 20000
            dist = np.abs(jj - ii)
            weights = clr.bins().fetch(chrom)["weight"].to_numpy()
            if (np.maximum(ii, jj) >= len(weights)).any() or (np.minimum(ii, jj) < 0).any():
                raise ValueError(f"Anchor outside chromosome: {chrom}")
            expected_values = np.array([expected.get((chrom, int(d)), np.nan) for d in dist])
            observed = np.full(len(group), np.nan)
            for block in np.unique(ii // block_rows):
                which = np.flatnonzero(ii // block_rows == block)
                lo = int(block * block_rows)
                hi = min(lo + block_rows, end - offset)
                sparse = clr.matrix(balance=True, sparse=True)[offset + lo:offset + hi, offset:end].tocsr()
                observed[which] = np.asarray(sparse[ii[which] - lo, jj[which]]).ravel()
            valid = ((dist >= 2) & np.isfinite(observed) & (observed > 0) &
                     np.isfinite(expected_values) & (expected_values > 0) &
                     np.isfinite(weights[ii]) & np.isfinite(weights[jj]))
            idx = group.index.to_numpy()[valid]
            values[idx] = np.log2(observed[valid] / expected_values[valid])
        result[row.sample] = values
    result.to_csv(out / "genomewide_loop_log2OE.tsv.gz", sep="\t", index=False)
    return result


def hic_effect_regions(args, loops, out):
    if args.hic_oe:
        oe = table(args.hic_oe)
        reject_sv_conditioning(oe)
        if oe.interaction_id.duplicated().any() or set(oe.interaction_id) != set(loops.interaction_id):
            raise ValueError("O/E table must contain ALL candidate loops exactly once, including missing signals")
        oe = loops.merge(oe[["interaction_id"] + SAMPLES], on="interaction_id", validate="one_to_one")
    else:
        require(args.hic_metadata)
        oe = quantify_hic(loops, args.hic_metadata, out, args.hic_block_rows)
    vals = oe[SAMPLES].apply(pd.to_numeric, errors="raise")
    complete = np.isfinite(vals).all(axis=1)
    for time in (1, 4):
        oe[f"delta_{time}"] = (vals[[f"A-1-{time}", f"A-3-{time}"]].mean(axis=1, skipna=False) -
                                vals[[f"W-1-{time}", f"W-3-{time}"]].mean(axis=1, skipna=False))
    effect = complete & ((oe.delta_1.abs() >= args.hic_delta) | (oe.delta_4.abs() >= args.hic_delta))
    bed_write(anchors(oe.loc[complete]), out / "hic_valid.bed")
    bed_write(anchors(oe.loc[effect]), out / "hic_effect.bed")
    oe["complete_8"] = complete
    oe["population_effect"] = effect
    oe.to_csv(out / "genomewide_loop_population_effect.tsv.gz", sep="\t", index=False)
    pd.DataFrame([{"all_loops": len(oe), "complete_loops": int(complete.sum()), "effect_loops": int(effect.sum()),
                   "definition": "complete 8 samples AND abs(Tibetan-Han mean log2OE)>=cutoff at 0h OR 3d"}]).to_csv(out / "hic_qc.tsv", sep="\t", index=False)


def norm_gene(value):
    value = str(value).strip()
    return value.split(".")[0] if value.startswith("ENSG") else value


def responsive_accessions(args, out):
    rna = table(args.rna_time)
    reject_sv_conditioning(rna)
    cols = ["time6h", "time1d", "time3d", "time5d"]
    responsive = (rna["adj.P.Val"] < args.fdr) & (rna[cols].abs().max(axis=1) >= args.min_logfc)
    keys = {}
    ambiguous = set()
    for i, row in rna.iterrows():
        for c in ["feature_id", "gene_id"]:
            if c not in rna or pd.isna(row[c]):
                continue
            key = norm_gene(row[c])
            if key in keys and keys[key] != i:
                ambiguous.add(key)
            keys[key] = i
    for key in ambiguous:
        keys.pop(key, None)
    if not args.gene_link:
        raise FileNotFoundError("Provide ALL three raw ENCODE4 V4-hg38.Gene-Links files with --gene-link")
    for path in args.gene_link:
        require(path)
    accessions, qc = set(), []
    for path in args.gene_link:
        opener = gzip.open if str(path).endswith(".gz") else open
        scanned = matched = hits = conflicts = 0
        print(f"[INFO] scanning global gene links: {path}", flush=True)
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if not re.fullmatch(r"EH38[DE]\d+", fields[0]):
                    if fields[0].lower() in {"ccre", "ccre_id", "ccre_accession", "accession"}:
                        continue
                    raise ValueError(f"Not raw ENCODE gene-link format: {path}: {line[:100]}")
                if len(fields) < 4:
                    raise ValueError(f"Malformed gene link in {path}")
                scanned += 1
                targets = {keys[k] for k in [norm_gene(fields[1]), norm_gene(fields[2])] if k in keys}
                if len(targets) > 1:
                    conflicts += 1
                    continue
                if targets:
                    matched += 1
                    if responsive.iloc[next(iter(targets))]:
                        accessions.add(fields[0])
                        hits += 1
        qc.append(dict(file=str(path), scanned=scanned, matched_RNA_rows=matched,
                       responsive_link_rows=hits, conflicting_gene_ID_symbol=conflicts))
    pd.DataFrame(qc).to_csv(out / "gene_link_scan_qc.tsv", sep="\t", index=False)
    if not sum(r["matched_RNA_rows"] for r in qc):
        raise ValueError("No ENCODE genes mapped to RNA IDs/symbols; explicit gene mapping needed")
    rna.loc[responsive].to_csv(out / "hypoxia_time_responsive_genes.tsv", sep="\t", index=False)
    return accessions


def ccre_regions(args, out, accessions=None):
    require(args.ccre_bed)
    name = "ccre_all" if accessions is None else "ccre_hypoxia"
    with (out / f"{name}.bed").open("w", encoding="utf-8") as handle:
        for chunk in pd.read_csv(args.ccre_bed, sep="\t", header=None, chunksize=150000):
            if chunk.shape[1] < 6:
                raise ValueError("Expected lifted ENCODE registry BED with rDHS/cCRE IDs in columns 4/5")
            if accessions is not None:
                chunk = chunk.loc[chunk[3].isin(accessions) | chunk[4].isin(accessions)]
            chunk[[0, 1, 2]].to_csv(handle, sep="\t", index=False, header=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["atac-rpkm", "atac-population", "hic-loops", "hic-metadata", "hic-oe", "rna-time", "ccre-bed"]:
        p.add_argument("--" + name, type=Path)
    p.add_argument("--gene-link", type=Path, action="append", default=[])
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--peak-cutoff", type=float, default=2)
    p.add_argument("--fdr", type=float, default=.05)
    p.add_argument("--min-logfc", type=float, default=.5)
    p.add_argument("--hic-delta", type=float, default=1)
    p.add_argument("--hic-block-rows", type=int, default=128)
    args = p.parse_args()
    if not 0 < args.fdr <= 1 or args.hic_block_rows < 1 or args.min_logfc < 0 or args.hic_delta < 0:
        p.error("Invalid cutoff or block size")
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    inventory = []
    for key, value in vars(args).items():
        paths = value if isinstance(value, list) else [value]
        for path in paths:
            if isinstance(path, Path) and key != "outdir":
                ok = path.is_file()
                inventory.append(dict(input=key, path=str(path), exists=ok,
                                      size_bytes=path.stat().st_size if ok else None,
                                      mtime_ns=path.stat().st_mtime_ns if ok else None))
    pd.DataFrame(inventory).to_csv(out / "input_inventory.tsv", sep="\t", index=False)
    status, definitions = {}, {}
    def attempt(names, fn, definition):
        try:
            val = fn()
            for name in names:
                status[name] = "READY"
            return val
        except (FileNotFoundError, ImportError) as e:
            for name in names:
                status[name] = "MISSING: " + str(e)
            print(f"[WARN] {names}: {e}", flush=True)
        finally:
            for name in names:
                definitions[name] = definition
    peaks = attempt(["atac_open"], lambda: atac_regions(args, out), f"max RPKM across 50 samples > {args.peak_cutoff}")
    if peaks is not None:
        attempt(["atac_dor"], lambda: dor_regions(args, peaks, out), f"open peak AND population FDR < {args.fdr} AND abs(logFC) >= {args.min_logfc}")
    loops = attempt(["hic_any"], lambda: load_loops(args.hic_loops), "all genomewide candidate loop anchors; no TESV filtering")
    if loops is not None:
        bed_write(anchors(loops), out / "hic_any.bed")
        attempt(["hic_valid", "hic_effect"], lambda: hic_effect_regions(args, loops, out),
                f"8 valid balanced log2OE samples; effect subset: abs(Tibetan-Han mean) >= {args.hic_delta} at 0h OR 3d; descriptive")
    attempt(["ccre_all"], lambda: ccre_regions(args, out), "all hg38 ENCODE4 cCREs lifted to CHM13v2.0")
    def gene_regions():
        accessions = responsive_accessions(args, out)
        ccre_regions(args, out, accessions)
    attempt(["ccre_hypoxia"], gene_regions,
            f"global ENCODE cCRE gene links; RNA time omnibus FDR < {args.fdr} AND max abs(time logFC) >= {args.min_logfc}; not population DEG")
    rows = []
    for endpoint, family, stem, pad in ENDPOINTS:
        s = status.get(stem, "MISSING: upstream input unavailable")
        rows.append(dict(endpoint=endpoint, family=family, bed=f"{stem}.bed" if s == "READY" else "",
                         status="READY" if s == "READY" else "MISSING", detail=s, padding_bp=pad,
                         scope="genomewide_independent_of_TESVs", assembly="CHM13v2.0",
                         definition=definitions.get(stem, "upstream unavailable")))
    pd.DataFrame(rows).to_csv(out / "region_manifest.tsv", sep="\t", index=False)
    (out / "region_configuration.json").write_text(json.dumps(vars(args), default=str, indent=2), encoding="utf-8")
    print(f"[INFO] Region manifest: {out / 'region_manifest.tsv'}", flush=True)


if __name__ == "__main__":
    main()
