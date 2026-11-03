#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


ap = argparse.ArgumentParser()
ap.add_argument("--typing-dir", required=True)
ap.add_argument("--out-manifest", required=True)
ap.add_argument("--out-ids", required=True)
a = ap.parse_args()

rows_out = []
for group in ("den", "nea", "ghost"):
    original = read_rows(Path(a.typing_dir) / f"326.{group}_typing.tsv")
    hg38 = read_rows(Path(a.typing_dir) / f"326.{group}_typing.hg38.tsv")
    unmapped = {
        r["ID"] for r in read_rows(Path(a.typing_dir) / f"326.{group}_typing.hg38.tsv.unmapped.tsv")
    }
    kept = [r for r in original if r["ID"] not in unmapped]
    if len(kept) != len(hg38):
        raise SystemExit(f"{group}: mapped row count mismatch: {len(kept)} vs {len(hg38)}")
    for old, new in zip(kept, hg38):
        for col in ("type", "AF_TIB", "AF_HAN", "AF_DIFF", "Denisovan_Altai",
                    "Neanderthal_Altai", "Neanderthal_Cha", "Neanderthal_Vindija"):
            if old[col] != new[col]:
                raise SystemExit(f"{group}: order/content mismatch at {old['ID']} / {new['ID']} ({col})")
        rows_out.append({
            "group": group.upper(), "case_id": new["ID"], "original_sv_id": old["ID"],
            "hg38_sv_id": new["ID"], "chm13_chrom": old["chr"], "chm13_start": old["start"],
            "chm13_end": old["end"], "hg38_chrom": new["chr"], "hg38_start": new["start"],
            "hg38_end": new["end"], "type": new["type"], "AF_TIB": new["AF_TIB"],
            "AF_HAN": new["AF_HAN"], "AF_DIFF": new["AF_DIFF"],
        })

fields = list(rows_out[0])
Path(a.out_manifest).parent.mkdir(parents=True, exist_ok=True)
with open(a.out_manifest, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows_out)

seen = set()
with open(a.out_ids, "w") as fh:
    for row in rows_out:
        sid = row["original_sv_id"]
        if sid not in seen:
            print(sid, file=fh)
            seen.add(sid)

print(f"group_rows={len(rows_out)} unique_original_SVs={len(seen)}")

