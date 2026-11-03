#!/usr/bin/env python3
"""Plot A+B+C D and fd window traces for DISC1."""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROWS = [
    ("A_relative", "A: Den, Nea, TBN, YRI"),
    ("B_Denisovan", "B: YRI, TBN, Denisovan, ancestral"),
    ("C_Neanderthal", "C: YRI, TBN, Neanderthal, ancestral"),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--windows", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--region-label", default="DISC1 ghost candidate")
    args = p.parse_args()
    w = pd.read_csv(args.windows, sep="\t")
    s = pd.read_csv(args.summary, sep="\t").set_index("comparison")
    x = w.midpoint / 1e6
    fig, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True,
                             gridspec_kw={"hspace": 0.24, "wspace": 0.18})
    for row_i, (code, label) in enumerate(ROWS):
        d_ax, fd_ax = axes[row_i]
        d = w[f"{code}_D_estimate"]
        fd = w[f"{code}_fd_estimate"]
        sig = w[f"{code}_D_q_BH"] < 0.05
        d_ax.axhline(0, color="#666666", lw=0.8)
        d_ax.plot(x, d, color="#D55E00", lw=1.2, drawstyle="steps-mid")
        d_ax.scatter(x[sig], d[sig], s=12, color="#A50026", zorder=3,
                     label="D BH-FDR q<0.05")
        d_ax.set_ylabel("D")
        d_ax.set_title(
            f"{label}\nregion D={s.loc[code,'D_estimate']:.4f}, "
            f"Z={s.loc[code,'D_z']:.2f}, P={s.loc[code,'D_p_two_sided']:.2g}",
            loc="left", fontsize=9)
        if row_i == 0:
            d_ax.legend(frameon=False, fontsize=7, loc="best")
        fd_ax.axhline(0, color="#666666", lw=0.8)
        fd_ax.plot(x, fd, color="#2878B5", lw=1.3, drawstyle="steps-mid")
        fd_ax.fill_between(x, 0, fd, step="mid", color="#5DA5DA", alpha=0.38)
        fd_ax.set_ylabel("$f_d$")
        fd_ax.set_title(
            f"region fd={s.loc[code,'fd_estimate']:.4f}; "
            "window significance uses D q",
            loc="left", fontsize=9)
        for ax in (d_ax, fd_ax):
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="x", color="#E0E0E0", lw=0.5)
    axes[-1, 0].set_xlabel("GRCh38 position on chr1 (Mb)")
    axes[-1, 1].set_xlabel("GRCh38 position on chr1 (Mb)")
    fig.suptitle(args.region_label + ": common ancestral-polarized TBN SNPs",
                 x=0.07, ha="left", fontsize=12)
    fig.savefig(args.out_prefix + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(args.out_prefix + ".pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
