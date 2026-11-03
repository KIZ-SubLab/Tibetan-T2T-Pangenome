#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


p = argparse.ArgumentParser()
p.add_argument("--plotcoords", required=True)
p.add_argument("--mutation-branches", required=True)
p.add_argument("--mut", required=True)
p.add_argument("--chrom", required=True)
p.add_argument("--position", required=True, type=int)
p.add_argument("--years-per-generation", default=28.0, type=float)
p.add_argument("--sample-label", default="target population")
p.add_argument("--out", required=True)
a = p.parse_args()

mpl.rcParams.update(
    {"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"}
)

row = None
with open(a.mut, encoding="utf-8") as fh:
    fh.readline()
    for line in fh:
        fields = line.rstrip("\n").split(";")
        if len(fields) > 10 and int(fields[1]) == a.position:
            if row is not None:
                raise ValueError(f"multiple mutation records at {a.position}")
            row = fields
if row is None:
    raise ValueError(f"position {a.position} absent from {a.mut}")

source_id = row[3]
focal_branch = int(row[5].replace(",", " ").split()[0])
age_lo, age_hi = float(row[8]), float(row[9])
anc, der = row[10].split("/")

coords = {"h": [], "v": [], "t": [], "m": []}
with open(a.plotcoords, encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        if not line.strip():
            continue
        xb, xe, yb, ye, branch, style = line.split()
        coords.setdefault(style, []).append(
            (float(xb), float(xe), float(yb), float(ye), int(branch))
        )

derived_branches = {focal_branch}
with open(a.mutation_branches, encoding="utf-8") as fh:
    next(fh, None)
    for line in fh:
        fields = line.split()
        if fields:
            try:
                derived_branches.add(int(fields[-1]))
            except ValueError:
                pass

g = a.years_per_generation
fig, ax = plt.subplots(figsize=(14, 10))
for style in ("h", "v"):
    for xb, xe, yb, ye, branch in coords.get(style, []):
        derived = branch in derived_branches
        ax.plot(
            [xb, xe],
            [yb * g, ye * g],
            color="#D73027" if derived else "#333333",
            lw=1.25 if derived else 0.45,
            zorder=2 if derived else 1,
        )

mx = [(xb + xe) / 2 for xb, xe, _, _, _ in coords.get("m", [])]
my = [(yb + ye) / 2 * g for _, _, yb, ye, _ in coords.get("m", [])]
if mx:
    ax.scatter(mx, my, s=4, color="#4575B4", alpha=0.55, zorder=3)

age_lo_y, age_hi_y = age_lo * g, age_hi * g
ax.axhspan(age_lo_y, age_hi_y, color="#74ADD1", alpha=0.12, zorder=0)
ax.axhline(age_lo_y, color="#4575B4", ls="--", lw=0.8)
ax.axhline(age_hi_y, color="#4575B4", ls="--", lw=0.8)
ax.set_ylabel("Years before present")
ax.set_xlabel(f"{a.sample_label} haplotypes")
ax.set_title(
    f"RELATE marginal genealogy at {a.chrom}:{a.position}\n"
    f"VCF ID {source_id} | ancestral {anc}, derived {der} | "
    f"mutation branch interval {age_lo_y/1000:.1f}-{age_hi_y/1000:.1f} kya"
)
ax.tick_params(axis="x", bottom=False, labelbottom=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(
    handles=[
        Line2D([0], [0], color="#D73027", lw=2, label=f"Derived {der} lineage"),
        Line2D([0], [0], color="#333333", lw=0.7, label=f"Ancestral {anc} lineage"),
        Line2D(
            [0],
            [0],
            color="#4575B4",
            ls="--",
            lw=1,
            label="RELATE mutation-branch interval",
        ),
    ],
    loc="upper right",
    frameon=False,
    fontsize=8,
)
fig.tight_layout()

out = Path(a.out)
out.parent.mkdir(parents=True, exist_ok=True)
for ext, kwargs in [("png", {"dpi": 450}), ("svg", {}), ("pdf", {})]:
    fig.savefig(out.with_suffix("." + ext), bbox_inches="tight", **kwargs)
plt.close(fig)

