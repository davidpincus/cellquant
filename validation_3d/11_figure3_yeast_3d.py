"""validation_3d/11_figure3_yeast_3d.py

Rebuild manuscript Figure 3B (yeast temperature-series quantification) on the
NATIVE 3D metrics, in the published superplot/trend style: per-cell points
coloured by replicate, replicate medians as large diamonds, condition medians
joined by a line, significance brackets vs the 25 C reference.

Emits the MAIN figure only: Figure3B_yeast_3D, 6 panels (2x3) -- Sis1/Nsr1/Tif6
mean intensity, then the three pairwise Pearson R coefficients (Sis1-Tif6,
Nsr1-Sis1, Nsr1-Tif6), all computed natively on the z-stacks.

The three compartment-volume panels this script used to emit as a standalone
supplement are now row 1 of Figure S2 (16_figureS2_yeast_supp.py), so that one
script owns the whole supplement and the two cannot drift apart. SUPP_PANELS is
kept below only as the record of which columns those were.

NOTE ON NUCLEAR VOLUME. It is not in the supplement because it does not exist
in this dataset: the series was acquired as Tif6 / Nsr1 / Sis1 with Nsr1 in the
`nucleolus` role and no `nucleus` channel, so nuclear segmentation never ran and
`nucleus_volume_um3` is 0 for every cell. See SUPP_PANELS for why every other
volume column is degenerate or redundant. Any panel whose column is absent or
identically zero renders as an explicit "not measured" placeholder rather than
a flat line at zero, which would read as a real null result.

Metrics dropped relative to the published 2D figure: nucleolar circularity,
solidity and eccentricity are properties of a projection and have no 3D
counterpart; sphericity is not measurable at this axial resolution and was
removed from cellquant.

Statistics: two-sided Mann-Whitney U on replicate medians, each
temperature vs 25 C, Bonferroni-corrected across the four comparisons. At
n = 5-8 replicates the exact test is appropriate; the normal approximation
(scipy.stats.ranksums) gives materially different values at this n.

Usage:
    python 11_figure3_yeast_3d.py [--outputs-dir DIR] [--outdir DIR]

Needs sklearn-free deps only: pandas, scipy, matplotlib.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory

# One size for every text element in both figures -- titles, axis labels, tick
# labels and significance labels. Previously these ran 7 / 8 / 9 / 9.5 pt.
FONT_SIZE = 12

matplotlib.rcParams.update({
    "svg.fonttype": "none",     # editable text in Illustrator
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.linewidth": 0.8,
    "axes.unicode_minus": False,
})

TEMP_ORDER = ["25C", "30C", "32C", "36C", "40C"]
REF = "25C"
N_COMPARISONS = 4          # 30/32/36/40 vs 25 -> Bonferroni denominator

# (column, panel title, y-axis label)
# Main figure: the three marker intensities, then the three pairwise
# colocalization coefficients computed natively on the z-stacks.
MAIN_PANELS = [
    ("Sis1_cell_mean",          "Sis1 intensity",           "mean intensity (AU)"),
    ("Nsr1_cell_mean",          "Nsr1 intensity",           "mean intensity (AU)"),
    ("Tif6_cell_mean",          "Tif6 intensity",           "mean intensity (AU)"),
    ("pearson_r_Tif6_vs_Sis1",  "Sis1–Tif6 colocalization", "Pearson $R$"),
    ("pearson_r_Sis1_vs_Nsr1",  "Nsr1–Sis1 colocalization", "Pearson $R$"),
    ("pearson_r_Tif6_vs_Nsr1",  "Nsr1–Tif6 colocalization", "Pearson $R$"),
]

# Supplementary figure: compartment volumes.
#
# Nuclear volume is NOT here, and cannot be: this series was acquired as
# Tif6 / Nsr1(nucleolus) / Sis1 with no nucleus channel, so nuclear
# segmentation never ran and `nucleus_volume_um3` is 0 for all 8309 cells.
# Every other volume column is degenerate for the same reason or redundant with
# a panel already shown -- `cytosol_volume_um3` is exactly equal to
# `cell_volume_um3` (cytosol = cell minus a nucleus that was never segmented),
# `nucleolar_eq_diameter_um` is the cube root of panel 3, and `n_nucleoli` is 1
# for every cell. The nucleolar volume fraction is the one non-redundant
# quantity left, and it is the one that bears on compaction: it falls mostly
# because the cell grows (37 -> 65 um^3), not because the nucleolus shrinks
# (2.88 -> 3.75 um^3). It is non-monotonic, spiking at 30 C like the raw series.
SUPP_PANELS = [
    ("cell_volume_um3",           "Cell volume",              "volume (µm$^3$)"),
    ("nucleolar_volume_um3",      "Nucleolar volume",         "volume (µm$^3$)"),
    ("nucleolar_volume_fraction", "Nucleolar volume fraction", "nucleolus / cell"),
]

REP_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
              "#937860", "#DA8BC3", "#8C8C8C"]

# Per-cell point size. The published version used s=1.6, which rendered the
# cloud as a faint smear that the replicate-median diamonds swamped entirely --
# the reader could not see the distribution the medians were drawn from. Larger
# points need a wider jitter and a wider per-replicate offset to stay separable,
# so these four move together; changing DOT_SIZE alone will blob the strips.
DOT_SIZE = 6.0
DOT_ALPHA = 0.30
DOT_JITTER = 0.028         # half-width of the random x-spread within a replicate
REP_OFFSET = 0.085         # x-gap between adjacent replicate strips
DIAMOND_SIZE = 52          # replicate medians -- unchanged; they are the stat unit

# Significance-bracket type size and vertical stacking pitch, in axes fractions.
# Defaults suit the 6-panel main figure. A dense supplement (15 small panels)
# must lower the size and raise the pitch or the labels overlap each other --
# 16_figureS2_yeast_supp.py overrides both. Pitch has to exceed the label height
# as a fraction of panel height, so the two move together.
BRACKET_FONT_SIZE = FONT_SIZE
BRACKET_PITCH = 0.075

# Panels pinned to a shared y-range so members of a family can be read against
# each other directly. Requested explicitly, and it does clip: on the intensity
# panels a handful of bright cells fall above 350 (Tif6 reaches 410), and on the
# Pearson panels the small negative-R tail (min -0.21) falls below 0. The
# replicate medians -- the statistical unit -- are well inside both ranges.
# A fixed axis has no headroom to expand into, so significance brackets on
# these panels are stacked downward from the top edge instead; see draw_panel.
YLIM = {
    "Sis1_cell_mean":         (100.0, 350.0),
    "Nsr1_cell_mean":         (100.0, 350.0),
    "Tif6_cell_mean":         (100.0, 350.0),
    "pearson_r_Tif6_vs_Sis1": (0.0, 1.0),
    "pearson_r_Sis1_vs_Nsr1": (0.0, 1.0),
    "pearson_r_Tif6_vs_Nsr1": (0.0, 1.0),
}


def save_all(fig, stem: Path):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix("." + ext), bbox_inches="tight", dpi=300)
    print(f"  wrote {stem.name}.svg / .pdf / .png")


def load_cells(outputs_dir: Path) -> pd.DataFrame:
    files = sorted(outputs_dir.glob("*/cells.csv"))
    if not files:
        raise FileNotFoundError(f"no */cells.csv under {outputs_dir}")
    frames = []
    for f in files:
        df = pd.read_csv(f)
        rep_dir = f.parent.name
        if "condition" not in df.columns or df["condition"].isna().all():
            cond, rep = rep_dir.split("_series1_rep", 1)
            df["condition"], df["replicate"] = cond, rep
        df["rep_id"] = rep_dir
        frames.append(df)
    cells = pd.concat(frames, ignore_index=True)
    if "keep" in cells.columns:
        cells = cells[cells["keep"] == True].copy()   # noqa: E712
    cells["replicate"] = cells["replicate"].astype(str)

    # Derived: nucleolus as a share of cell volume. Guard the divisor rather
    # than trusting it -- a zero-volume cell would otherwise yield inf and drag
    # the panel's y-limits to infinity.
    if {"nucleolar_volume_um3", "cell_volume_um3"} <= set(cells.columns):
        denom = cells["cell_volume_um3"].where(cells["cell_volume_um3"] > 0)
        cells["nucleolar_volume_fraction"] = cells["nucleolar_volume_um3"] / denom

    print(f"loaded {len(files)} replicates, {len(cells)} cells (keep==True)")
    return cells.reset_index(drop=True)


def replicate_medians(cells: pd.DataFrame, col: str) -> pd.DataFrame:
    """One median per replicate image -- the unit of statistical analysis."""
    return (cells.groupby(["condition", "rep_id"], as_index=False)[col]
                 .median().rename(columns={col: "value"}))


def stats_vs_ref(rm: pd.DataFrame) -> dict[str, dict]:
    ref = rm.loc[rm["condition"] == REF, "value"].dropna().values
    out = {}
    for t in TEMP_ORDER:
        if t == REF:
            continue
        v = rm.loc[rm["condition"] == t, "value"].dropna().values
        if len(v) < 2 or len(ref) < 2:
            out[t] = {"p_raw": np.nan, "p_adj": np.nan, "sig": False}
            continue
        p = mannwhitneyu(v, ref, alternative="two-sided").pvalue
        p_adj = min(p * N_COMPARISONS, 1.0)
        out[t] = {"p_raw": float(p), "p_adj": float(p_adj), "sig": p_adj < 0.05}
    return out


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "n.s."
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}".rstrip("0").rstrip(".") if p < 0.01 else f"p = {p:.2f}"


def draw_panel(ax, cells, rm, st, col, title, ylab, reps_by_cond):
    x_of = {t: i for i, t in enumerate(TEMP_ORDER)}
    rng = np.random.default_rng(0)

    # per-cell points, coloured by replicate within each condition
    for t in TEMP_ORDER:
        sub = cells[cells["condition"] == t]
        for j, rid in enumerate(reps_by_cond[t]):
            v = sub.loc[sub["rep_id"] == rid, col].dropna().values
            if not len(v):
                continue
            off = (j - (len(reps_by_cond[t]) - 1) / 2) * REP_OFFSET
            x = x_of[t] + off + rng.uniform(-DOT_JITTER, DOT_JITTER, len(v))
            ax.scatter(x, v, s=DOT_SIZE, alpha=DOT_ALPHA, linewidths=0,
                       color=REP_COLORS[j % len(REP_COLORS)], rasterized=True)

    # replicate medians as diamonds + condition median line
    cond_med = []
    for t in TEMP_ORDER:
        vals = rm.loc[rm["condition"] == t]
        for j, rid in enumerate(reps_by_cond[t]):
            r = vals.loc[vals["rep_id"] == rid, "value"]
            if r.empty or not np.isfinite(r.iloc[0]):
                continue
            off = (j - (len(reps_by_cond[t]) - 1) / 2) * REP_OFFSET
            ax.scatter(x_of[t] + off, r.iloc[0], s=DIAMOND_SIZE, marker="D",
                       color=REP_COLORS[j % len(REP_COLORS)],
                       edgecolor="k", linewidths=0.7, zorder=5)
        cond_med.append(float(np.nanmedian(vals["value"].values)) if len(vals) else np.nan)
    ax.plot(range(len(TEMP_ORDER)), cond_med, "-", color="k", lw=1.3, zorder=6)
    ax.scatter(range(len(TEMP_ORDER)), cond_med, s=26, color="k", zorder=7)

    ax.set_xticks(range(len(TEMP_ORDER)))
    ax.set_xticklabels([t.replace("C", "") for t in TEMP_ORDER])
    ax.set_xlabel("Temperature (°C)", fontsize=FONT_SIZE)
    ax.set_ylabel(ylab, fontsize=FONT_SIZE)
    ax.set_title(title, fontsize=FONT_SIZE, fontweight="bold")
    ax.tick_params(labelsize=FONT_SIZE)

    fixed = YLIM.get(col)
    if fixed:
        ax.set_ylim(*fixed)
    else:
        # Scale to the bulk of the per-cell distribution, not to its outliers.
        # Auto-limits driven by a handful of extreme cells squash every trend
        # into the bottom of the panel; clip at the 0.5/99.5 percentiles and
        # make sure all replicate medians stay in view.
        vals = cells[col].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            lo, hi = np.nanpercentile(vals, [0.5, 99.5])
            meds = rm["value"].to_numpy(dtype=float)
            meds = meds[np.isfinite(meds)]
            if meds.size:
                lo, hi = min(lo, meds.min()), max(hi, meds.max())
            if hi <= lo:
                hi = lo + 1.0
            pad = (hi - lo) * 0.08
            lo = lo - pad if lo < 0 else max(0.0, lo - pad)
            ax.set_ylim(lo, hi + pad)

    sig = [(t, st[t]) for t in TEMP_ORDER if t != REF and st[t]["sig"]]
    if sig and not fixed:
        # Free axis: grow headroom so the bracket stack sits ABOVE the data,
        # then fall through to the same compact drawing the pinned axes use.
        # The old free-axis path stacked tall brackets at 0.13 of the span each,
        # which ran into the panel title as soon as a panel had three or four
        # of them -- fine in a 6-panel main figure, unreadable in a 15-panel
        # supplement. Reserve is capped so a panel with four brackets still
        # gives the data half its height.
        reserve = min(0.10 + BRACKET_PITCH * len(sig), 0.45)
        lo_, hi_ = ax.get_ylim()
        ax.set_ylim(lo_, lo_ + (hi_ - lo_) / (1.0 - reserve))
    if sig:
        # Brackets stacked downward from the top of the panel, label sitting ON
        # the bar with its white box masking the line behind.
        # Order is reversed against k so the widest span still sits highest --
        # nested, as when they were drawn above the data. Labels get a white
        # backing box because at this size they can cross the point cloud.
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        # The label sits ON the bar rather than above it, its white box masking
        # the line behind. That halves each bracket's vertical footprint to
        # roughly one line of text, which is what lets three of them stack
        # inside a pinned panel without reaching the title or the diamonds.
        # Ticks point down, toward the groups being compared.
        for k, (t, s) in enumerate(sig):
            y = 0.96 - BRACKET_PITCH * (len(sig) - 1 - k)
            x0, x1 = x_of[REF], x_of[t]
            ax.plot([x0, x0, x1, x1], [y - 0.016, y, y, y - 0.016],
                    lw=0.9, color="k", transform=tr, zorder=8,
                    path_effects=[pe.withStroke(linewidth=2.8,
                                                foreground="white")])
            ax.text((x0 + x1) / 2, y, fmt_p(s["p_adj"]),
                    ha="center", va="center", fontsize=BRACKET_FONT_SIZE,
                    transform=tr, zorder=9,
                    bbox=dict(boxstyle="square,pad=0.10", fc="white",
                              ec="none", alpha=0.85))


def build_figure(cells, reps_by_cond, panels, nrows, ncols, figsize, stem, rows):
    """Render one grid of superplot panels and append its stats to `rows`."""
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(nrows, ncols, figure=fig, hspace=0.42, wspace=0.34)

    for i, (col, title, ylab) in enumerate(panels):
        ax = fig.add_subplot(gs[i // ncols, i % ncols])

        # A column that is absent, or present but identically zero, is not a
        # flat trend -- it is a quantity this run never measured. Drawing it as
        # a line at zero would read as a real null result, so say so instead.
        note = None
        if col not in cells.columns:
            note = f"{col}\nnot in 3D output"
        else:
            v = cells[col].to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            if not v.size or np.all(v == 0):
                note = (f"{title}\nnot measured in this run\n"
                        f"({col} is 0 for all {len(cells)} cells)")
        if note:
            ax.text(0.5, 0.5, note, ha="center", va="center", fontsize=FONT_SIZE,
                    color="0.45", transform=ax.transAxes)
            ax.set_title(title, fontsize=FONT_SIZE, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_linestyle((0, (3, 3))); s.set_color("0.7")
            print(f"  !! {col}: {note.splitlines()[-1]}")
            continue

        rm = replicate_medians(cells, col)
        st = stats_vs_ref(rm)
        draw_panel(ax, cells, rm, st, col, title, ylab, reps_by_cond)

        med = {t: float(np.nanmedian(rm.loc[rm["condition"] == t, "value"].values))
               for t in TEMP_ORDER}
        print(f"\n{title}  [{col}]")
        print("   medians: " + "  ".join(f"{t}={med[t]:.4g}" for t in TEMP_ORDER))
        print("   vs 25C : " + "  ".join(
            f"{t}: p={st[t]['p_raw']:.4f} p_adj={st[t]['p_adj']:.4f}"
            f"{' *' if st[t]['sig'] else ''}" for t in TEMP_ORDER if t != REF))
        for t in TEMP_ORDER:
            rows.append({"figure": stem.name, "panel": title, "metric": col,
                         "temperature": t,
                         "replicate_median_of_medians": med[t],
                         "n_replicates": len(reps_by_cond[t]),
                         "p_raw_vs_25C": st[t]["p_raw"] if t != REF else np.nan,
                         "p_bonferroni_vs_25C": st[t]["p_adj"] if t != REF else np.nan})

    save_all(fig, stem)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--outputs-dir", type=Path,
                    default=(here.parent / "jcb_revision_midway3_cellquant_analysis"
                             / "validation_3d" / "outputs_yeast" / "3d"))
    ap.add_argument("--outdir", type=Path, default=here / "figure3_3d")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(args.outputs_dir)
    reps_by_cond = {t: sorted(cells.loc[cells["condition"] == t, "rep_id"].unique())
                    for t in TEMP_ORDER}
    print("replicates per condition:", {t: len(v) for t, v in reps_by_cond.items()})

    rows = []
    print("\n" + "=" * 70 + "\nMAIN FIGURE (Fig 3B)\n" + "=" * 70)
    build_figure(cells, reps_by_cond, MAIN_PANELS, 2, 3, (11.5, 7),
                 args.outdir / "Figure3B_yeast_3D", rows)

    pd.DataFrame(rows).to_csv(args.outdir / "figure3b_3d_stats.csv", index=False)
    print(f"\n  wrote figure3b_3d_stats.csv")


if __name__ == "__main__":
    np.random.seed(0)
    main()
