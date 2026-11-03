#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d


p = argparse.ArgumentParser()
p.add_argument("--freqs", required=True)
p.add_argument("--post", required=True)
p.add_argument("--out", required=True)
p.add_argument("--title", required=True)
p.add_argument("--max-years", required=True, type=float)
p.add_argument("--generation-years", default=28.0, type=float)
p.add_argument("--derived-label", default="derived")
p.add_argument("--sample-label", default="target population")
p.add_argument("--n-samples", required=True, type=int)
p.add_argument("--n-haplotypes", required=True, type=int)
a = p.parse_args()

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

grid = np.loadtxt(a.freqs, delimiter=",", ndmin=1)
post = np.loadtxt(a.post, delimiter=",", ndmin=2)
if post.shape[0] != grid.size and post.shape[1] == grid.size:
    post = post.T
if post.shape[0] != grid.size:
    raise ValueError(f"posterior rows {post.shape[0]} != frequency grid {grid.size}")
if not np.all(np.isfinite(post)) or not np.all(np.isfinite(grid)):
    raise ValueError("non-finite CLUES frequency/posterior value")
post = np.clip(post, 0, None)
column_sums = post.sum(axis=0, keepdims=True)
if np.any(column_sums <= 0):
    raise ValueError("CLUES posterior contains an empty time column")
post /= column_sums
cdf = np.cumsum(post, axis=0)


def quantile(q):
    return grid[np.argmax(cdf >= q, axis=0)]


low, median, high = quantile(0.025), quantile(0.5), quantile(0.975)
generations = np.arange(post.shape[1], dtype=float)
years = generations * a.generation_years
keep = years <= a.max_years
years, generations = years[keep], generations[keep]
low, median, high = low[keep], median[keep], high[keep]

out = Path(a.out)
out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(
    {
        "generations_before_present": generations,
        "years_before_present": years,
        "q2_5": low,
        "median": median,
        "q97_5": high,
    }
).to_csv(out.with_suffix(".source_data.tsv"), sep="\t", index=False)

sigma = 25
low_s = gaussian_filter1d(low, sigma=sigma, mode="nearest")
med_s = gaussian_filter1d(median, sigma=sigma, mode="nearest")
high_s = gaussian_filter1d(high, sigma=sigma, mode="nearest")
low_s = np.minimum(low_s, med_s)
high_s = np.maximum(high_s, med_s)
display_max = min(a.max_years, float(years.max()))
display_years = np.linspace(0, display_max, int(display_max / 50) + 1)
low_d = np.clip(PchipInterpolator(years, low_s)(display_years), 0, 1)
med_d = np.clip(PchipInterpolator(years, med_s)(display_years), 0, 1)
high_d = np.clip(PchipInterpolator(years, high_s)(display_years), 0, 1)
pd.DataFrame(
    {
        "years_before_present": display_years,
        "q2_5_display": low_d,
        "median_display": med_d,
        "q97_5_display": high_d,
        "display_smoothing_sigma_generations": sigma,
    }
).to_csv(out.with_suffix(".display_smooth.tsv"), sep="\t", index=False)

fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.fill_between(
    display_years,
    low_d,
    high_d,
    color="#F2C94C",
    alpha=0.72,
    lw=0,
    label="95% posterior interval",
)
ax.plot(display_years, med_d, color="#7A3E00", lw=1.6, label="Posterior median")
ax.set(
    xlim=(0, a.max_years),
    ylim=(0, 1),
    xlabel="Years before present",
    ylabel=f"Derived {a.derived_label} allele frequency",
    title=a.title,
)
ax.set_xticks(np.linspace(0, a.max_years, 5))
ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper right", frameon=False, fontsize=7)
ax.text(
    0.01,
    0.04,
    f"{a.sample_label}: {a.n_samples} diploid individuals ({a.n_haplotypes} haplotypes)",
    transform=ax.transAxes,
    fontsize=6.6,
    color="#333333",
)
fig.tight_layout()
for ext, kwargs in [
    ("png", {"dpi": 450}),
    ("svg", {}),
    ("pdf", {}),
    ("tiff", {"dpi": 600}),
]:
    fig.savefig(out.with_suffix("." + ext), bbox_inches="tight", **kwargs)
plt.close(fig)

