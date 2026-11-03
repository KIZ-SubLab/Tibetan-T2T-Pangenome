#!/usr/bin/env python3
"""
Convert RepeatMasker .out tables into a BED6 satellite-class track.

Kept classes:
  alpha_satellite   ALR/alpha
  beta_satellite    BSR/beta
  gamma_satellite   GSAT/gamma
  HSat1             HSATI / HSat1A / HSat1B
  HSat2             HSATII
  HSat3             HSATIII
  other_satellite   remaining satellite entries (SAR etc.)

The mapping mirrors the class strings used by the project's
downstream HOR/monomer scripts (BSat, GSat, HSat1A, HSat1B,
HSat2, HSat3, alpha).
"""
import argparse
import glob
import os

# RepeatMasker class/family -> project satellite class
CLASS_MAP = {
    "ALR/Alpha": "alpha_satellite",
    "BSR/Beta": "beta_satellite",
    "GSAT": "gamma_satellite",
    "HSATI": "HSat1",
    "HSAT1A": "HSat1",
    "HSAT1B": "HSat1",
    "HSATII": "HSat2",
    "HSATIII": "HSat3",
}


def classify(rep_class, rep_family):
    key = f"{rep_class}/{rep_family}" if rep_family else rep_class
    for k, v in CLASS_MAP.items():
        if key.upper().startswith(k.upper()):
            return v
    if "SAT" in rep_class.upper():
        return "other_satellite"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    n = 0
    with open(a.output, "w") as out:
        for f in sorted(glob.glob(os.path.join(a.outdir, "*.out"))):
            sample = os.path.basename(f).replace(".out", "")
            with open(f) as fh:
                for line in fh:
                    p = line.split()
                    if len(p) < 15 or not p[0].isdigit():
                        continue
                    cls = classify(p[10], p[11] if len(p) > 11 else "")
                    if cls is None:
                        continue
                    # .out is 1-based inclusive; convert to 0-based half-open
                    start = int(p[5]) - 1
                    end = int(p[6])
                    contig = p[4]
                    out.write(f"{contig}\t{start}\t{end}\t{cls}\t{sample}\t{p[8]}\n")
                    n += 1
    print(f"wrote {n} satellite intervals -> {a.output}")


if __name__ == "__main__":
    main()
