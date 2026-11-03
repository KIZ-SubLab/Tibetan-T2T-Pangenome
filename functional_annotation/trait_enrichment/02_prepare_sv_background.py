#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the SV-associated protein-coding gene background.

Gene-level intervals from the CHM13 v2.0 RefSeq Liftoff v5.2 annotation are
intersected with the merged SV +/- 5 kb intervals (at least 1 bp overlap).
Gene symbols are then harmonized to the GENCODE v38 protein-coding universe by
exact match, by Ensembl-ID-verified symbol updates, and optionally by
unambiguous mapping through a single stable Ensembl gene ID.

Coordinates: input intervals are 1-based inclusive; BED output is 0-based
half-open (START - 1, END). Insertion loci are represented as a single
breakpoint position in the input interval file.

Usage:
  python3 02_prepare_sv_background.py \
      --sv-regions sv_5kb_regions.tsv --chm13-genes chm13_gene_level.tsv \
      --gencode-pc gencode_v38_protein_coding_genes.txt \
      --gencode-id-map gencode_v38_gene_id_map.tsv \
      --outdir results/sv [--sv-bed sv5kb_regions.bed] [--allow-network no]

Outputs (in --outdir):
  SV5kb_regions.bed
  SV5kb_background_GENCODEv38_proteinCoding.txt
  SV5kb_background_symbol_harmonization.tsv
"""

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

# Symbol updates confirmed through stable Ensembl gene identifiers
KNOWN_MAPPINGS = {"ZFTA": ("C11orf95", "ENSG00000188070"),
                  "SMIM45": ("LINC00634", "ENSG00000205704"),
                  "PALS2": ("MPP6", "ENSG00000105926")}


def ensembl_xref(symbol, timeout=10):
    """Return Ensembl gene IDs for a symbol via the Ensembl REST API."""
    url = ("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/"
           f"{urllib.parse.quote(symbol)}?content-type=application/json")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            docs = json.load(response)
        return sorted({d["id"] for d in docs if d.get("type") == "gene"
                       and str(d.get("id", "")).startswith("ENSG")})
    except Exception:
        return []


def main():
    ap = argparse.ArgumentParser(
        description="Build the SV +/- 5 kb matched protein-coding gene background")
    ap.add_argument("--sv-regions", required=True,
                    help="Merged SV +/- 5 kb intervals (1-based inclusive TSV)")
    ap.add_argument("--chm13-genes", required=True,
                    help="CHM13 gene-level annotation TSV")
    ap.add_argument("--gencode-pc", required=True,
                    help="GENCODE v38 protein-coding gene symbols")
    ap.add_argument("--gencode-id-map", required=True,
                    help="TSV without header: gene_id, gene_name, gene_type")
    ap.add_argument("--outdir", required=True,
                    help="Output directory")
    ap.add_argument("--sv-bed", default=None,
                    help="Precomputed SV +/- 5 kb intervals in BED format")
    ap.add_argument("--allow-network", choices=["yes", "no"], default="no",
                    help="Allow Ensembl REST lookups for ambiguous symbols")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pc_genes = set(pd.read_csv(args.gencode_pc, sep="\t", header=None,
                               dtype=str)[0])
    id_map = pd.read_csv(args.gencode_id_map, sep="\t", header=None, dtype=str,
                         names=["gene_id", "gene_name", "gene_type"])
    id_to_symbol = dict(zip(id_map["gene_id"], id_map["gene_name"]))
    id_to_type = dict(zip(id_map["gene_id"], id_map["gene_type"]))

    print("Building SV-associated gene background")

    # SV intervals -> BED
    sv_bed = os.path.join(args.outdir, "SV5kb_regions.bed")
    if args.sv_bed and os.path.isfile(args.sv_bed):
        sv_bed_path = args.sv_bed
        print(f"Using precomputed SV intervals: {sv_bed_path}")
    else:
        sv = pd.read_csv(args.sv_regions, sep="\t", dtype=str)
        sv["START"] = pd.to_numeric(sv["START_1BASED"]) - 1  # 1-based -> 0-based
        sv["END"] = pd.to_numeric(sv["END_1BASED"])
        sv[["CHROM", "START", "END", "SV_ID", "SVTYPE", "SV_POS", "SV_END"]] \
            .to_csv(sv_bed, sep="\t", index=False, header=False)
        sv_bed_path = sv_bed
        print(f"SV intervals: {len(sv):,} -> {sv_bed}")

    # CHM13 gene-level annotation -> BED
    annotation = pd.read_csv(args.chm13_genes, sep="\t", dtype=str,
                             keep_default_na=False)
    gene_bed = os.path.join(args.outdir, "chm13_gene_level.bed")
    pd.DataFrame({
        "CHROM": annotation["CHROM"],
        "START": pd.to_numeric(annotation["START_1BASED"]) - 1,
        "END": pd.to_numeric(annotation["END_1BASED"]),
        "GENE_ID": annotation["GENE_ID"],
        "GENE_NAME": annotation["GENE_NAME"],
        "GENE_BIOTYPE": annotation["GENE_BIOTYPE"]}
    ).to_csv(gene_bed, sep="\t", index=False, header=False)
    print(f"CHM13 gene features: {len(annotation):,}")

    # Interval intersection (>= 1 bp overlap)
    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-a", sv_bed_path, "-b", gene_bed,
             "-wa", "-wb"], capture_output=True, text=True, check=True)
    except FileNotFoundError:
        sys.exit("Error: bedtools not found (bedtools intersect is required)")

    records = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        records.append({"SV_ID": fields[3], "SVTYPE": fields[4],
                        "GENE_NAME": fields[10], "GENE_BIOTYPE": fields[12]})
    mapping = pd.DataFrame(records)
    n_sv = int(pd.read_csv(sv_bed_path, sep="\t", header=None,
                           dtype=str)[3].nunique())
    chm13_genes = sorted(set(x for x in mapping["GENE_NAME"] if x and x != "NA"))
    print(f"SV-gene pairs: {len(mapping):,} | SVs with a gene: "
          f"{mapping['SV_ID'].nunique():,} / {n_sv:,} | "
          f"unique CHM13 genes: {len(chm13_genes):,}")

    # Harmonize to GENCODE v38 protein-coding symbols
    biotype = mapping.drop_duplicates("GENE_NAME").set_index(
        "GENE_NAME")["GENE_BIOTYPE"].to_dict()
    exact = sorted(set(chm13_genes) & pc_genes)
    unmatched = [g for g in chm13_genes
                 if g not in pc_genes and biotype.get(g) == "protein_coding"]
    harmonization, harmonized = [], []
    for symbol in unmatched:
        if symbol in KNOWN_MAPPINGS:
            updated, ensembl_id = KNOWN_MAPPINGS[symbol]
            harmonized.append(updated)
            harmonization.append({
                "INPUT_SYMBOL": symbol, "GENCODE_V38_SYMBOL": updated,
                "ENSEMBL_ID": ensembl_id, "STATUS": "mapped_symbol_update",
                "EVIDENCE": f"symbol update {symbol} -> {updated}"})
            continue
        if args.allow_network == "yes":
            ensembl_ids = ensembl_xref(symbol)
            time.sleep(0.1)
            if len(ensembl_ids) == 1 \
                    and id_to_type.get(ensembl_ids[0]) == "protein_coding" \
                    and id_to_symbol.get(ensembl_ids[0]) in pc_genes:
                updated = id_to_symbol[ensembl_ids[0]]
                harmonized.append(updated)
                harmonization.append({
                    "INPUT_SYMBOL": symbol, "GENCODE_V38_SYMBOL": updated,
                    "ENSEMBL_ID": ensembl_ids[0], "STATUS": "mapped_ensembl_id",
                    "EVIDENCE": f"unique Ensembl gene ID {ensembl_ids[0]}"})
                continue
        harmonization.append({
            "INPUT_SYMBOL": symbol, "GENCODE_V38_SYMBOL": "", "ENSEMBL_ID": "",
            "STATUS": "excluded",
            "EVIDENCE": "no unique mapping to GENCODE v38 protein-coding"})

    background = sorted(set(exact) | set(harmonized))
    assert set(background) <= pc_genes

    pd.DataFrame(harmonization).to_csv(os.path.join(
        args.outdir, "SV5kb_background_symbol_harmonization.tsv"),
        sep="\t", index=False)
    output_path = os.path.join(
        args.outdir, "SV5kb_background_GENCODEv38_proteinCoding.txt")
    with open(output_path, "w") as handle:
        for gene in background:
            handle.write(gene + "\n")
    print(f"Exact matches: {len(exact):,} | harmonized: {len(harmonized):,} | "
          f"background genes: {len(background):,}")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    main()
