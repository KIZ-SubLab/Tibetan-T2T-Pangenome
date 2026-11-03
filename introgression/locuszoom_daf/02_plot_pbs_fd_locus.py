#!/usr/bin/env python3
"""Plot a publication-grade PBS/fD locus with regulatory and gene tracks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd
import pyBigWig

PBS_CMAP = LinearSegmentedColormap.from_list("pbs_gray_orange", ["#B9B9B9", "#F28E2B"])
FD_COLOR = "#355EAA"
CCRE_COLORS = {"PLS": "#D73027", "pELS": "#F4A620", "dELS": "#F6D32D", "other": "#37A6C6"}
H3_COLORS = {
    "GM12878": "#6A3D9A", "H1-hESC": "#B159B2", "HSMM": "#33A02C",
    "HUVEC": "#1F9E89", "K562": "#E31A1C", "NHEK": "#FDBF6F", "NHLF": "#1F78B4",
}


def ccre_group(value: str) -> str:
    low = str(value).lower()
    if low.startswith("pls") or "promoter" in low:
        return "PLS"
    if low.startswith("pels") or "proximal enhancer" in low:
        return "pELS"
    if low.startswith("dels") or "distal enhancer" in low:
        return "dELS"
    return "other"


def load_ccres(root: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    path = root / "refs/encodeCcreCombined_ENCODE3.bb"
    with pyBigWig.open(str(path)) as bigbed:
        entries = bigbed.entries(chrom, start - 1, end) or []
    rows = []
    for start0, end0, payload in entries:
        fields = payload.split("\t")
        rows.append({
            "chrom": chrom, "start0": int(start0), "end0": int(end0),
            "start": int(start0) + 1, "end": int(end0), "accession": fields[0],
            "score": int(fields[1]), "ccre_class": fields[6], "encode_label": fields[7],
            "dnase_max_z": float(fields[8]),
        })
    columns = ["chrom", "start0", "end0", "start", "end", "accession", "score", "ccre_class", "encode_label", "dnase_max_z"]
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_csv(root / "data/ENCODE3_cCREs.region.tsv", sep="\t", index=False)
    return frame


def add_ccre_track(axis, ccres: pd.DataFrame):
    for row in ccres.itertuples():
        group = ccre_group(row.ccre_class)
        axis.add_patch(Rectangle((row.start, 0.10), row.end - row.start + 1, 0.80,
                                 facecolor=CCRE_COLORS[group], edgecolor="none"))
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.spines[:].set_visible(False)
    axis.text(-0.025, 0.50, "ENCODE3 cCREs", transform=axis.transAxes,
              ha="right", va="center", fontsize=7.0)


def add_h3_track(axis, h3: pd.DataFrame, start: int, end: int):
    n_bins = max(600, int(np.ceil((end - start + 1) / 250)))
    edges = np.linspace(start - 1, end, n_bins + 1)
    widths = np.diff(edges)
    x = (edges[:-1] + edges[1:]) / 2 + 0.5
    parameters = {}
    for cell_line, color in H3_COLORS.items():
        signal_sum = np.zeros(n_bins, dtype=float)
        for row in h3[h3.cell_line == cell_line].itertuples():
            interval_start = max(float(row.start0), edges[0])
            interval_end = min(float(row.end0), edges[-1])
            if interval_end <= interval_start:
                continue
            left = max(0, np.searchsorted(edges, interval_start, side="right") - 1)
            right = min(n_bins - 1, np.searchsorted(edges, interval_end, side="left"))
            for index in range(left, right + 1):
                overlap = max(0.0, min(interval_end, edges[index + 1]) - max(interval_start, edges[index]))
                signal_sum[index] += overlap * max(0.0, float(row.value))
        signal = signal_sum / widths
        baseline = float(np.median(signal))
        corrected = np.clip(signal - baseline, 0, None)
        positive = corrected[corrected > 0]
        scale = float(np.percentile(positive, 99.5)) if positive.size else 1.0
        scale = max(scale, np.finfo(float).eps)
        relative = np.clip(corrected / scale, 0, 1)
        parameters[cell_line] = {"median_background": baseline, "p99_5_positive_residual": scale}
        axis.fill_between(x, 0, relative, step="mid", color=color, alpha=0.16, linewidth=0)
        axis.plot(x, relative, color=color, alpha=0.72, lw=0.32)
    axis.set_ylim(0, 1.04)
    axis.set_yticks([0, 1])
    axis.tick_params(axis="y", labelsize=5.6, length=2, width=0.5, pad=1.5)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_linewidth(0.55)
    axis.text(-0.025, 0.52, "H3K27ac\nrelative signal", transform=axis.transAxes,
              ha="right", va="center", fontsize=6.7)
    axis.text(0.995, 0.94, "7 ENCODE cell lines; median-subtracted, P99.5 = 1",
              transform=axis.transAxes, ha="right", va="top", fontsize=5.3, color="#555555")
    return parameters


def add_gene_track(axis, genes: pd.DataFrame, exons: pd.DataFrame, start: int, end: int):
    if genes.empty:
        axis.text(0.5, 0.5, "No GENCODE genes in interval", transform=axis.transAxes,
                  ha="center", va="center")
        axis.set_yticks([])
        return
    lane_count = int(genes.lane.max()) + 1
    for gene in genes.itertuples():
        y = lane_count - 1 - int(gene.lane)
        x1, x2 = max(start, int(gene.start)), min(end, int(gene.end))
        color = "#18256B" if gene.gene_type == "protein_coding" else "#667089"
        axis.plot([x1, x2], [y, y], color=color, lw=0.8, solid_capstyle="round", zorder=1)
        for exon in exons[exons.transcript_id == gene.transcript_id].itertuples():
            axis.add_patch(Rectangle((exon.start, y - 0.11), exon.end - exon.start + 1, 0.22,
                                     facecolor=color, edgecolor=color, lw=0.3, zorder=2))
        if x2 > x1:
            arrow_x = np.linspace(x1 + 0.08 * (x2 - x1), x2 - 0.08 * (x2 - x1), 4)
            axis.scatter(arrow_x, np.full_like(arrow_x, y), marker=">" if gene.strand == "+" else "<",
                         s=7, color=color, linewidths=0, zorder=3)
        if not str(gene.gene_name).startswith("ENSG"):
            label_x = min(end - 10_000, max(start + 10_000, (x1 + x2) / 2))
            axis.text(label_x, y + 0.20, gene.gene_name, ha="center", va="bottom", fontsize=5.8,
                      color="#111111", fontstyle="italic" if gene.gene_type == "protein_coding" else "normal")
    axis.set_ylim(-0.55, lane_count - 0.15 + 0.65)
    axis.set_yticks([])
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.spines["bottom"].set_linewidth(0.6)
    axis.text(-0.025, 0.50, "GENCODE v50", transform=axis.transAxes,
              ha="right", va="center", fontsize=7.0)


def tick_positions(start: int, end: int):
    span = end - start
    rough_step = span / 6
    magnitude = 10 ** np.floor(np.log10(max(rough_step, 1)))
    step = min([1, 2, 5, 10], key=lambda value: abs(value * magnitude - rough_step)) * magnitude
    return np.arange(np.ceil(start / step) * step, end + 1, step, dtype=int)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--point-size", type=float, default=15.0)
    parser.add_argument("--archaic-code", choices=["DEN", "NEA", "GHOST"], required=True)
    parser.add_argument("--archaic-label", required=True)
    parser.add_argument("--pbs-ymax", type=float)
    args = parser.parse_args()
    chrom = args.chrom if args.chrom.startswith("chr") else f"chr{args.chrom}"
    root, start, end = args.run_root, args.start, args.end
    if args.point_size <= 0:
        raise ValueError("--point-size must be positive")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    pbs = pd.read_csv(root / "data/PBS_sites.region.tsv", sep="\t")
    fd_name = "GHOST_ABC_fD_windows.region.tsv" if args.archaic_code == "GHOST" else f"{args.archaic_code}_fD_windows.region.tsv"
    fd = pd.read_csv(root / "data" / fd_name, sep="\t")
    genes = pd.read_csv(root / "data/gencode_v50_basic.representative_genes.tsv", sep="\t")
    exons = pd.read_csv(root / "data/gencode_v50_basic.representative_exons.tsv", sep="\t")
    h3 = pd.read_csv(root / "data/UCSC_layered_H3K27ac_7cell.region.tsv.gz", sep="\t")
    ccres = load_ccres(root, chrom, start, end)
    summary = json.loads((root / "metadata/analysis_summary.json").read_text())

    figure = plt.figure(figsize=(7.10, 5.90))
    grid = figure.add_gridspec(4, 1, height_ratios=[5.6, 0.42, 1.12, 2.35], hspace=0.10)
    pbs_axis = figure.add_subplot(grid[0])
    ccre_axis = figure.add_subplot(grid[1], sharex=pbs_axis)
    h3_axis = figure.add_subplot(grid[2], sharex=pbs_axis)
    gene_axis = figure.add_subplot(grid[3], sharex=pbs_axis)
    fd_axis = pbs_axis.twinx()

    if args.archaic_code == "GHOST":
        fd_axis.bar(fd.midpoint, fd.A_relative_fd_estimate,
                    width=fd.display_width, color=FD_COLOR, alpha=0.25,
                    edgecolor=FD_COLOR, linewidth=0.25, align="center", zorder=1)
        fd_values = fd.A_relative_fd_estimate
    else:
        fd_axis.bar(fd.midpoint, fd.fD, width=fd.display_width, color=FD_COLOR, alpha=0.22,
                    edgecolor=FD_COLOR, linewidth=0.25, align="center", zorder=1)
        fd_values = fd.fD
    fd_axis.axhline(0, color="#D73027" if args.archaic_code == "GHOST" else FD_COLOR,
                    lw=0.60, alpha=0.80 if args.archaic_code == "GHOST" else 0.65, zorder=1)
    fd_low = min(-0.05, float(fd_values.min()) * 1.18)
    fd_high = max(0.10, float(fd_values.max()) * 1.12)
    fd_axis.set_ylim(fd_low, fd_high)
    fd_axis.set_ylabel(
        r"Test A relative $f_D$" if args.archaic_code == "GHOST" else rf"{args.archaic_label} $f_D$",
        fontsize=8.3, color=FD_COLOR, labelpad=8)
    fd_axis.tick_params(axis="y", colors=FD_COLOR, labelsize=6.8, width=0.55, length=3)
    fd_axis.spines["right"].set_color(FD_COLOR)
    fd_axis.spines["right"].set_linewidth(0.65)
    fd_axis.spines["top"].set_visible(False)
    fd_axis.patch.set_visible(False)
    fd_axis.set_zorder(1)

    norm = Normalize(vmin=0, vmax=float(pbs.PBS.max()))
    scatter = pbs_axis.scatter(pbs.POS, pbs.PBS, c=pbs.PBS, cmap=PBS_CMAP, norm=norm,
                               s=args.point_size, edgecolors="none", alpha=0.90, zorder=4)
    pbs_axis.set_zorder(2)
    pbs_axis.patch.set_alpha(0)
    pbs_axis.set_xlim(start, end)
    data_max = float(pbs.PBS.max())
    if args.pbs_ymax is not None and args.pbs_ymax < data_max:
        raise ValueError(
            f"--pbs-ymax ({args.pbs_ymax:g}) is smaller than the regional PBS maximum ({data_max:g}); "
            "increase it to avoid clipping points")
    pbs_ymax = args.pbs_ymax if args.pbs_ymax is not None else data_max * 1.14
    pbs_axis.set_ylim(0, pbs_ymax)
    pbs_axis.set_ylabel("Population branch statistic (PBS)", fontsize=8.3, color="#D97800")
    pbs_axis.tick_params(axis="y", colors="#D97800", labelsize=6.8, width=0.55, length=3)
    pbs_axis.tick_params(axis="x", labelbottom=False)
    pbs_axis.spines[["top", "right"]].set_visible(False)
    pbs_axis.spines["left"].set_color("#D97800")
    pbs_axis.spines["left"].set_linewidth(0.65)
    pbs_axis.spines["bottom"].set_linewidth(0.65)
    pbs_axis.grid(axis="y", color="#D9D9D9", lw=0.42, alpha=0.48, zorder=0)

    color_axis = pbs_axis.inset_axes([0.035, 0.835, 0.26, 0.045])
    colorbar = figure.colorbar(scatter, cax=color_axis, orientation="horizontal")
    colorbar.set_ticks([0, round(float(pbs.PBS.max()) / 2, 3), round(float(pbs.PBS.max()), 3)])
    colorbar.ax.tick_params(labelsize=5.8, length=2, width=0.45, pad=1)
    colorbar.outline.set_linewidth(0.4)
    colorbar.ax.set_title("PBS", fontsize=6.2, pad=2, loc="left")

    if args.archaic_code == "GHOST":
        legend_handles = [Patch(facecolor=FD_COLOR, edgecolor=FD_COLOR, alpha=0.25,
                                label=r"Test A: relative $f_D$")] 
    else:
        legend_handles = [Patch(facecolor=FD_COLOR, edgecolor=FD_COLOR, alpha=0.22,
                                label=rf"{args.archaic_code} $f_D$ (100-SNP windows; 50-SNP step)")]
    pbs_axis.legend(handles=legend_handles,
                    loc="upper right", bbox_to_anchor=(0.995, 0.995), fontsize=6.2,
                    handlelength=1.4, borderaxespad=0.2, ncol=1)
    title_subject = "ghost-candidate test A" if args.archaic_code == "GHOST" else f"{args.archaic_code} introgression"
    pbs_axis.text(0, 1.22, f"PBS and {title_subject} across {chrom}:{start:,}-{end:,}",
                  transform=pbs_axis.transAxes, fontsize=10.2, weight="bold", ha="left", va="bottom")
    pbs_axis.text(0, 1.09,
                  (f"n={len(pbs):,} PBS-positive SNPs; {len(fd)} test A windows; "
                   "background=test A, D/fD(Denisovan, Neanderthal, TBN, YRI); GRCh38")
                  if args.archaic_code == "GHOST" else
                  (f"n={len(pbs):,} PBS-positive SNPs; {len(fd)} overlapping {args.archaic_code} $f_D$ windows; "
                   f"source={args.archaic_label}; GRCh38"),
                  transform=pbs_axis.transAxes, fontsize=6.5, color="#555555", ha="left", va="bottom")

    add_ccre_track(ccre_axis, ccres)
    h3_parameters = add_h3_track(h3_axis, h3, start, end)
    add_gene_track(gene_axis, genes, exons, start, end)
    gene_axis.set_xlabel(f"GRCh38 {chrom} position (Mb)", fontsize=8.2)
    ticks = tick_positions(start, end)
    gene_axis.set_xticks(ticks)
    tick_step = float(np.min(np.diff(ticks))) if len(ticks) > 1 else float(end - start)
    tick_decimals = 1 if tick_step >= 100_000 else 2 if tick_step >= 10_000 else 3
    gene_axis.set_xticklabels([f"{tick / 1e6:.{tick_decimals}f}" for tick in ticks], fontsize=6.8)
    for linked_axis in [ccre_axis, h3_axis]:
        linked_axis.tick_params(axis="x", labelbottom=False, length=0)

    figure.subplots_adjust(left=0.155, right=0.855, top=0.805, bottom=0.090)
    prefix = root / "figures" / f"{chrom}_{start}_{end}_PBS_{args.archaic_code}_fD"
    figure.savefig(prefix.with_suffix(".svg"), format="svg")
    figure.savefig(prefix.with_suffix(".pdf"), format="pdf")
    figure.savefig(prefix.with_suffix(".png"), format="png", dpi=300)
    figure.savefig(prefix.with_suffix(".tiff"), format="tiff", dpi=600,
                   pil_kwargs={"compression": "tiff_lzw"})
    plt.close(figure)

    spec = {
        "figure_contract": {
            "core_conclusion": (
                "Regional PBS peaks can be localized relative to test A Denisovan-vs-Neanderthal affinity in a ghost-candidate region."
                if args.archaic_code == "GHOST" else
                f"Regional PBS peaks can be localized relative to {args.archaic_label} fD signal and regulatory/gene context."
            ),
            "archetype": "quantitative genomic track with one hero dual-axis panel and three supporting annotation tracks",
            "backend": "Python matplotlib only",
            "final_size_mm": [180.3, 149.9],
            "reviewer_risk": "PBS and fD are different estimators and scales; dual axes and source/window definitions are explicitly labeled.",
        },
        "point_y": "PBS from the supplied PBS-positive SNP table",
        "point_color": "PBS, fixed linear normalization from 0 to the regional maximum; low gray to high orange",
        "point_marker_area_points_squared": args.point_size,
        "pbs_y_axis": {"minimum": 0.0, "maximum": float(pbs_ymax),
                       "mode": "user-defined" if args.pbs_ymax is not None else "automatic"},
        "background_bars": (
            "Ghost mode displays only test A (A_relative_fd_estimate): D/fD(Denisovan,Neanderthal,TBN,YRI), matching panel A of the Ghost_DISC1 diagnostic figure. It is a relative sampled-archaic affinity statistic, not a directly observed ghost fD; tests B and C remain in the source results but are not drawn in the integrated figure."
            if args.archaic_code == "GHOST" else
            f"fD(EAS,TBN;{args.archaic_label},YRI), 100-SNP windows with 50-SNP step; bars centered on window midpoint"
        ),
        "fd_data_coverage_note": f"The first source fD window starts at {summary['fd_first_window_start']}; no fD value is imputed before it.",
        "ccre_track": "ENCODE3 combined cell-type-agnostic cCRE registry; all overlapping elements",
        "encode3_ccres_in_region": int(len(ccres)),
        "h3k27ac_track": "UCSC hg38 Layered H3K27ac, seven ENCODE cell lines; 250-bp base-overlap-weighted means",
        "h3k27ac_display_normalization": "Per cell line: subtract median bin mean, divide positive residuals by P99.5, clip to [0,1]",
        "h3k27ac_display_parameters_raw_units": h3_parameters,
        "gene_track": "one representative transcript per gene from GENCODE v50 basic; MANE Select/canonical preferred",
        "exports": ["editable SVG", "PDF with embedded TrueType fonts", "PNG 300 dpi", "TIFF 600 dpi LZW"],
    }
    (root / "metadata/figure_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    print(json.dumps(spec, indent=2))


if __name__ == "__main__":
    main()
