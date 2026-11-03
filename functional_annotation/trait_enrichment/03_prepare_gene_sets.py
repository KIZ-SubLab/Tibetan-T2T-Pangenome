#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare candidate gene sets for enrichment analysis.

All gene sets are harmonized to the GENCODE v38 protein-coding universe by:
  1. exact symbol match
  2. confirmed symbol updates verified through stable Ensembl gene IDs
  3. optional unambiguous mapping through a single Ensembl gene ID
     (--allow-network yes)
No fuzzy matching is applied.

Outputs (in --outdir):
  selection_148.txt
  tesv_full_proteinCoding.txt              # TE-SV genes in the protein-coding universe
  tesv_in_SV5kb_background.txt             # TE-SV genes within the SV +/- 5 kb background
  tesv_excluded_outside_background.tsv     # TE-SV genes outside that background
  union297.txt                             # selection union TE-SV (within background)
  union297_composition.tsv
  hypoxia_harmonized.txt                   # hypoxia genes in the protein-coding universe
  hypoxia_in_SV5kb_background.txt          # hypoxia genes within the SV background
  symbol_harmonization.tsv

Usage:
  python3 03_prepare_gene_sets.py \
      --selection-input selection_symbols.txt --tesv-input tesv_symbols.txt \
      --hypoxia-input hypoxia_symbols.txt \
      --gencode-pc gencode_v38_protein_coding_genes.txt \
      --gencode-id-map gencode_v38_gene_id_map.tsv \
      --sv-background SV5kb_background_GENCODEv38_proteinCoding.txt \
      --outdir results/gene_sets [--allow-network no]
"""

import argparse
import json
import os
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


def harmonize(symbols, pc_genes, id_to_symbol, id_to_type, allow_network, label):
    """Map raw symbols to unique GENCODE v38 protein-coding symbols."""
    exact = sorted(set(symbols) & pc_genes)
    remaining = [s for s in sorted(set(symbols)) if s not in pc_genes]
    kept, rows = set(exact), []
    for symbol in remaining:
        if symbol in KNOWN_MAPPINGS:
            updated = KNOWN_MAPPINGS[symbol][0]
            kept.add(updated)
            rows.append({"SET": label, "INPUT_SYMBOL": symbol,
                         "GENCODE_V38_SYMBOL": updated,
                         "STATUS": "mapped_symbol_update"})
            continue
        mapped = None
        if allow_network == "yes":
            ensembl_ids = ensembl_xref(symbol)
            time.sleep(0.1)
            if len(ensembl_ids) == 1 \
                    and id_to_type.get(ensembl_ids[0]) == "protein_coding" \
                    and id_to_symbol.get(ensembl_ids[0]) in pc_genes:
                mapped = id_to_symbol[ensembl_ids[0]]
        if mapped:
            kept.add(mapped)
            rows.append({"SET": label, "INPUT_SYMBOL": symbol,
                         "GENCODE_V38_SYMBOL": mapped,
                         "STATUS": "mapped_ensembl_id"})
        else:
            rows.append({"SET": label, "INPUT_SYMBOL": symbol,
                         "GENCODE_V38_SYMBOL": "", "STATUS": "excluded"})
    result = sorted(kept)
    assert set(result) <= pc_genes
    print(f"{label}: {len(set(symbols))} input symbols -> {len(result)} "
          f"protein-coding genes")
    return result, pd.DataFrame(rows)


def write_list(path, genes):
    with open(path, "w") as handle:
        for gene in genes:
            handle.write(gene + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Prepare selection, TE-SV, union and hypoxia gene sets")
    ap.add_argument("--selection-input", required=True,
                    help="Selection candidate gene symbols (one per line)")
    ap.add_argument("--tesv-input", required=True,
                    help="TE-SV candidate gene symbols (one per line)")
    ap.add_argument("--hypoxia-input", required=True,
                    help="Hypoxia-related gene symbols (one per line)")
    ap.add_argument("--gencode-pc", required=True,
                    help="GENCODE v38 protein-coding gene symbols")
    ap.add_argument("--gencode-id-map", required=True,
                    help="TSV without header: gene_id, gene_name, gene_type")
    ap.add_argument("--sv-background", required=True,
                    help="SV +/- 5 kb matched background gene list")
    ap.add_argument("--outdir", required=True,
                    help="Output directory")
    ap.add_argument("--allow-network", choices=["yes", "no"], default="no",
                    help="Allow Ensembl REST lookups for ambiguous symbols")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pc_genes = set(pd.read_csv(args.gencode_pc, sep="\t", header=None,
                               dtype=str)[0])
    background = set(pd.read_csv(args.sv_background, sep="\t", header=None,
                                 dtype=str)[0])
    id_map = pd.read_csv(args.gencode_id_map, sep="\t", header=None, dtype=str,
                         names=["gene_id", "gene_name", "gene_type"])
    id_to_symbol = dict(zip(id_map["gene_id"], id_map["gene_name"]))
    id_to_type = dict(zip(id_map["gene_id"], id_map["gene_type"]))

    def read_symbols(path):
        return pd.read_csv(path, sep="\t", header=None, dtype=str,
                           keep_default_na=False)[0].str.strip().tolist()

    print("Preparing candidate gene sets")

    selection, h1 = harmonize(read_symbols(args.selection_input), pc_genes,
                              id_to_symbol, id_to_type, args.allow_network,
                              "SELECTION")
    tesv_full, h2 = harmonize(read_symbols(args.tesv_input), pc_genes,
                              id_to_symbol, id_to_type, args.allow_network,
                              "TESV")
    tesv_in = sorted(set(tesv_full) & background)
    tesv_out = sorted(set(tesv_full) - background)

    hypoxia, h3 = harmonize(read_symbols(args.hypoxia_input), pc_genes,
                            id_to_symbol, id_to_type, args.allow_network,
                            "HYPOXIA")
    hypoxia_in = sorted(set(hypoxia) & background)

    union = sorted(set(selection) | set(tesv_in))
    shared = sorted(set(selection) & set(tesv_in))
    composition = pd.DataFrame({
        "GENE": union,
        "IN_SELECTION": [int(g in set(selection)) for g in union],
        "IN_TESV": [int(g in set(tesv_in)) for g in union],
        "IN_SV5KB_BACKGROUND": [int(g in background) for g in union],
        "SOURCE": ["shared" if (g in set(selection) and g in set(tesv_in))
                   else ("selection_only" if g in set(selection) else "tesv_only")
                   for g in union]})

    write_list(os.path.join(args.outdir, "selection_148.txt"), selection)
    write_list(os.path.join(args.outdir, "tesv_full_proteinCoding.txt"), tesv_full)
    write_list(os.path.join(args.outdir, "tesv_in_SV5kb_background.txt"), tesv_in)
    write_list(os.path.join(args.outdir, "union297.txt"), union)
    write_list(os.path.join(args.outdir, "hypoxia_harmonized.txt"), hypoxia)
    write_list(os.path.join(args.outdir, "hypoxia_in_SV5kb_background.txt"),
               hypoxia_in)
    pd.DataFrame({"GENE": tesv_out, "REASON": "OUTSIDE_SV5KB_BACKGROUND"}).to_csv(
        os.path.join(args.outdir, "tesv_excluded_outside_background.tsv"),
        sep="\t", index=False)
    composition.to_csv(os.path.join(args.outdir, "union297_composition.tsv"),
                       sep="\t", index=False)
    pd.concat([h1, h2, h3]).to_csv(
        os.path.join(args.outdir, "symbol_harmonization.tsv"), sep="\t",
        index=False)

    print(f"Selection: {len(selection)} | TE-SV: {len(tesv_full)} total, "
          f"{len(tesv_in)} within background ({len(tesv_out)} excluded) | "
          f"Union: {len(union)} (shared {len(shared)}) | "
          f"Hypoxia: {len(hypoxia)} total, {len(hypoxia_in)} within background")


if __name__ == "__main__":
    main()
