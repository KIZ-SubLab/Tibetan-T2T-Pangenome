#!/usr/bin/env python3
"""Plot a Science Fig. 5B-style local TMRCA posterior heatmap."""

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Mandatory editable-text and publication export settings.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({
    "font.size": 7,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
})


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--posterior", required=True)
    p.add_argument("--parameters", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--physical-mutation-rate", type=float, default=1.25e-8)
    p.add_argument("--generation-years", type=float, default=29.0)
    p.add_argument("--y-max-mya", type=float, default=2.0)
    p.add_argument("--title", default="TBN2 hap2 vs hap1")
    return p.parse_args()


def read_decoder(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if not header or header[0] != "#time_boundaries:":
        raise ValueError("decoder output lacks #time_boundaries header")
    boundaries = np.array([float(x) for x in header[1:] if x != ""], dtype=float)
    arr = np.loadtxt(path, comments="#", ndmin=2)
    if arr.shape[1] != len(boundaries) + 1:
        raise ValueError(
            f"posterior has {arr.shape[1]-1} states but {len(boundaries)} boundaries"
        )
    positions = arr[:, 0]
    posterior = arr[:, 1:]
    return positions, posterior, boundaries


def read_parameters(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(df["parameter"], df["value"]))


def state_conditional_means(left, lambdas):
    """Conditional mean time in each piecewise-exponential coalescence state."""
    means = np.empty_like(left)
    for i, (a, lam) in enumerate(zip(left, lambdas)):
        if lam <= 0:
            raise ValueError("all lambda values must be positive")
        if i + 1 == len(left):
            means[i] = a + 1.0 / lam
            continue
        width = left[i + 1] - a
        z = lam * width
        if z > 50.0:
            # width / expm1(z) is numerically zero in this limit.
            means[i] = a + 1.0 / lam
        elif abs(z) < 1e-6:
            means[i] = a + width / 2.0
        else:
            means[i] = a + 1.0 / lam - width / math.expm1(z)
    return means


def coordinate_edges(x):
    if len(x) == 1:
        return np.array([x[0] - 500.0, x[0] + 500.0])
    mids = (x[:-1] + x[1:]) / 2.0
    return np.r_[x[0] - (mids[0] - x[0]), mids, x[-1] + (x[-1] - mids[-1])]


def main():
    args = parse_args()
    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    pos, posterior, raw_left = read_decoder(args.posterior)
    params = read_parameters(args.parameters)
    scaled_mut = float(params["MSMC2_scaled_mutation_rate"])
    lambdas = np.array([float(x) for x in params["decoder_lambda_vector"].split(",")])
    if len(lambdas) != posterior.shape[1]:
        raise ValueError("lambda vector and posterior state counts differ")
    row_sums = posterior.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-5, rtol=1e-5):
        raise ValueError(f"posterior rows do not sum to one: range={row_sums.min()}..{row_sums.max()}")

    raw_means = state_conditional_means(raw_left, lambdas)
    time_factor = scaled_mut / args.physical_mutation_rate * args.generation_years / 1e6
    left_mya = raw_left * time_factor
    state_mya = raw_means * time_factor
    posterior_mean = posterior @ state_mya
    map_index = np.argmax(posterior, axis=1)
    map_mya = state_mya[map_index]

    # Construct finite heatmap edges. The last state is open-ended; its display
    # edge is extended beyond the requested y limit without changing posterior data.
    y_edges = np.r_[left_mya, max(args.y_max_mya * 1.05, state_mya[-1] * 1.05)]
    x_mbp = pos / 1e6
    x_edges = coordinate_edges(x_mbp)

    fig = plt.figure(figsize=(89 / 25.4, 67 / 25.4))
    ax = fig.add_axes([0.14, 0.16, 0.68, 0.68])
    mesh = ax.pcolormesh(
        x_edges, y_edges, posterior.T, cmap="Greys", shading="flat",
        vmin=0.0, vmax=float(np.nanmax(posterior)), rasterized=True,
    )
    ax.plot(x_mbp, posterior_mean, color="#C4433A", lw=0.9, zorder=3)
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(0, args.y_max_mya)
    ax.set_xlabel("chr1 (Mbp)")
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Time to most recent common ancestor\n(TMRCA; million years ago)", labelpad=7)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="y", left=False, right=True)
    ax.text(-0.20, 1.13, "B", transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="left", va="bottom")
    if args.title:
        ax.text(0.0, 1.04, args.title, transform=ax.transAxes, fontsize=7,
                ha="left", va="bottom")

    cax = fig.add_axes([0.42, 0.89, 0.30, 0.045])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.ax.xaxis.set_ticks_position("bottom")
    cb.ax.tick_params(labelsize=6, length=2, pad=1)
    cb.outline.set_visible(False)
    fig.text(0.40, 0.915, "Posterior probability", ha="right", va="center", fontsize=6.5)

    source = pd.DataFrame({
        "position_GRCh38": pos.astype(int),
        "position_Mbp": x_mbp,
        "posterior_mean_TMRCA_Mya": posterior_mean,
        "MAP_state": map_index,
        "MAP_TMRCA_Mya": map_mya,
        "max_posterior_probability": posterior[np.arange(len(pos)), map_index],
    })
    for i in range(posterior.shape[1]):
        source[f"posterior_state_{i}"] = posterior[:, i]
    source.to_csv(str(out_prefix) + ".source_data.tsv", sep="\t", index=False)

    state_table = pd.DataFrame({
        "state": np.arange(len(raw_left)),
        "left_boundary_raw": raw_left,
        "lambda": lambdas,
        "conditional_mean_raw": raw_means,
        "left_boundary_Mya": left_mya,
        "conditional_mean_Mya": state_mya,
    })
    state_table.to_csv(str(out_prefix) + ".time_states.tsv", sep="\t", index=False)

    qa = {
        "decoded_positions": int(len(pos)),
        "position_min_GRCh38": int(pos.min()),
        "position_max_GRCh38": int(pos.max()),
        "posterior_row_sum_min": float(row_sums.min()),
        "posterior_row_sum_max": float(row_sums.max()),
        "posterior_probability_max": float(posterior.max()),
        "posterior_mean_TMRCA_Mya_median": float(np.median(posterior_mean)),
        "posterior_mean_TMRCA_Mya_min": float(np.min(posterior_mean)),
        "posterior_mean_TMRCA_Mya_max": float(np.max(posterior_mean)),
        "MAP_TMRCA_Mya_median": float(np.median(map_mya)),
        "physical_mutation_rate_per_site_per_generation": args.physical_mutation_rate,
        "generation_years": args.generation_years,
        "MSMC2_scaled_mutation_rate": scaled_mut,
        "time_conversion": "raw_state_time * MSMC2_scaled_mutation_rate / physical_mutation_rate * generation_years / 1e6",
        "red_line": "posterior mean using conditional mean within each fitted piecewise-exponential state",
        "smoothing": "none",
    }
    with open(str(out_prefix) + ".qa.json", "w") as fh:
        json.dump(qa, fh, indent=2, sort_keys=True)
        fh.write("\n")

    for ext, kwargs in {
        "svg": {},
        "pdf": {},
        "tiff": {"dpi": 600},
        "png": {"dpi": 300},
    }.items():
        fig.savefig(str(out_prefix) + "." + ext, bbox_inches="tight", pad_inches=0.02, **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
