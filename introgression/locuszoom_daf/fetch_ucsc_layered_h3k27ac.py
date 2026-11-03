#!/usr/bin/env python3
"""Fetch the UCSC hg38 seven-cell layered H3K27ac track for one region."""

from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import pandas as pd
import requests

TRACKS = {
    "GM12878": "wgEncodeRegMarkH3k27acGm12878",
    "H1-hESC": "wgEncodeRegMarkH3k27acH1hesc",
    "HSMM": "wgEncodeRegMarkH3k27acHsmm",
    "HUVEC": "wgEncodeRegMarkH3k27acHuvec",
    "K562": "wgEncodeRegMarkH3k27acK562",
    "NHEK": "wgEncodeRegMarkH3k27acNhek",
    "NHLF": "wgEncodeRegMarkH3k27acNhlf",
}


def get_json(url: str) -> dict:
    error = None
    for attempt in range(5):
        try:
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed to fetch {url}: {error}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    args = parser.parse_args()
    chrom = args.chrom if args.chrom.startswith("chr") else f"chr{args.chrom}"
    out = args.run_root / "data/UCSC_layered_H3K27ac_7cell.region.tsv.gz"
    raw_dir = args.run_root / "metadata/UCSC_layered_H3K27ac_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    rows, track_meta = [], {}
    for cell, track in TRACKS.items():
        url = ("https://api.genome.ucsc.edu/getData/track?genome=hg38;"
               f"track={track};chrom={chrom};start={args.start - 1};end={args.end}")
        payload = get_json(url)
        (raw_dir / f"{cell}.json.gz").write_bytes(gzip.compress(json.dumps(payload).encode()))
        records = payload.get(track, [])
        for record in records:
            rows.append({"chrom": record["chrom"], "start0": record["start"], "end0": record["end"],
                         "value": record["value"], "cell_line": cell, "track": track})
        track_meta[cell] = {"track": track, "sqlTable": payload.get("sqlTable"),
                            "bigDataUrl": payload.get("bigDataUrl"), "dataTime": payload.get("dataTime"),
                            "records": len(records), "api_url": url}
        print(f"{cell}: {len(records)} intervals")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False, compression="gzip")
    summary = {"source": "UCSC Layered H3K27ac on 7 ENCODE cell lines", "assembly": "hg38",
               "region": f"{chrom}:{args.start}-{args.end}", "tracks": track_meta}
    (args.run_root / "metadata/UCSC_layered_H3K27ac_7cell.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
