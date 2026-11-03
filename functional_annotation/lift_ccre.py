#!/usr/bin/env python3
"""Lift BED3+3 cCRE annotations without interpreting accessions as BED scores."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed", required=True, type=Path)
    p.add_argument("--chain", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--liftover", default="liftOver")
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    output = a.outdir / "ccre_t2t.bed"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    attrs = []
    with tempfile.TemporaryDirectory(dir=a.outdir) as tmp:
        tmp = Path(tmp)
        source, lifted = tmp / "source.bed4", tmp / "lifted.bed4"
        with a.bed.open(encoding="utf-8") as src, source.open("w", encoding="utf-8") as dst:
            for fields in csv.reader(src, delimiter="\t"):
                if not fields or fields[0].startswith("#"):
                    continue
                if len(fields) < 6:
                    raise ValueError("Expected chromosome,start,end,rDHS,cCRE,class")
                if int(fields[1]) < 0 or int(fields[2]) <= int(fields[1]):
                    raise ValueError("Invalid source cCRE BED interval")
                attrs.append(fields[3:6])
                dst.write("\t".join(fields[:3] + [str(len(attrs)-1)]) + "\n")
        cmd = [a.liftover, "-minMatch=0.95", str(source), str(a.chain),
               str(lifted), str(a.outdir / "unmapped.bed4")]
        subprocess.run(cmd, check=True)
        seen = set()
        pending = tmp / "restored.bed"
        with lifted.open() as src, pending.open("w", encoding="utf-8") as dst:
            for fields in csv.reader(src, delimiter="\t"):
                i = int(fields[3])
                if i in seen or not 0 <= i < len(attrs):
                    raise ValueError("Unexpected multiple/unknown liftOver output record")
                seen.add(i)
                dst.write("\t".join(fields[:3] + attrs[i]) + "\n")
        with (a.outdir / "input_record_ids.tsv").open("w", encoding="utf-8") as dst:
            dst.write("record_id\trDHS\tcCRE\tclass\n")
            for i, values in enumerate(attrs):
                dst.write("\t".join([str(i)] + values) + "\n")
        pending.replace(output)
    qc = dict(input_records=len(attrs), mapped_records=len(seen),
              unmapped_records=len(attrs)-len(seen), minMatch=0.95,
              multiple=False, chain=str(a.chain), input_bed=str(a.bed))
    (a.outdir / "liftover_qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
