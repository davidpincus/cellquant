"""validation_3d/13_figure2_mammalian_2d.py

Rebuild manuscript Figure 2B (HCT116 stress granule quantification) and Table S1
from the reproducible 2D re-run of the published MIPs.

WHY 2D, NOT 3D. The mammalian 3D arm has 3 control vs 3 arsenite replicates
(the published MIP set has 3 vs 4), and at n=3 vs 3 the smallest achievable
two-sided exact Mann-Whitney p is 0.10 -- no panel could reach significance.
Figure 2 also reports no colocalization, so the projection argument that forces
the yeast colocalization into 3D does not apply here. Projected puncta counts
and intensity-ratio condensation measures are legitimate quantities.

SOURCE DATA. outputs_mammalian/published_2d/ -- cellquant re-run of the 7
published MIPs (control rep1/4/5, arsenite rep1-4). The per-cell output of the
original March 2026 run no longer exists on disk, and Cellpose is not
deterministic across versions, so values here differ from the published Table
S1. This run is the one the deposited code reproduces.

STATISTICS. Two-sided exact Mann-Whitney U on replicate medians, arsenite vs
control. A single comparison per metric, so no multiple-comparison correction
(contrast the yeast series, which is Bonferroni-corrected across its four
temperature comparisons). Aggregation is the median of replicate medians
throughout, matching the yeast figures and tables.

Usage:
    python 13_figure2_mammalian_2d.py [--outputs-dir DIR] [--outdir DIR]
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# One size for every text element, matching 11_figure3_yeast_3d.py so the two
# figures set the same type in print. Previously these ran 7.5 / 8 / 9 / 9.5 pt.
FONT_SIZE = 12

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import SymmetricalLogLocator
from matplotlib.transforms import blended_transform_factory

matplotlib.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": FONT_SIZE, "axes.linewidth": 0.8, "axes.unicode_minus": False,
})

ORDER = ["control", "arsenite"]
COND_COLOR = {"control": "#4C72B0", "arsenite": "#C44E52"}
REP_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]

# Panels mirror the published Fig 2B set.
PANELS = [
    ("G3BP1_puncta_n",                 "G3BP1 puncta per cell",   "puncta / cell"),
    ("G3BP1_puncta_area_px",           "G3BP1 puncta area",       "area (px)"),
    ("G3BP1_frac_intensity_in_puncta", "G3BP1 fraction condensed", "fraction of signal"),
    ("PABPC1_frac_intensity_in_puncta","PABPC1 fraction condensed","fraction of signal"),
]
# Symlog for the panels whose per-cell values span ~100x. On a linear axis the
# whole distribution -- and every replicate median with it -- collapses onto the
# floor, which defeats the point of the superplot: the reader is meant to judge
# the effect size by eye. Log alone will not do, because 10-12 of the 79 cells
# contain exactly zero puncta; symlog keeps a linear zone below `linthresh` so
# those cells stay on the plot instead of being silently dropped. `linthresh`
# sits below the smallest nonzero value in each column (area 3 px; fractions
# 2.3e-4), so every real measurement lands in the log region.
YSCALE = {
    "G3BP1_puncta_area_px":            ("symlog", 1.0),
    "G3BP1_frac_intensity_in_puncta":  ("symlog", 1e-4),
    "PABPC1_frac_intensity_in_puncta": ("symlog", 1e-4),
}
# Extra metrics reported in Table S1 but not plotted.
TABLE_EXTRA = [
    ("PABPC1_puncta_n",             "PABPC1 puncta per cell"),
    ("G3BP1_condensate_index_cell", "G3BP1 condensate index"),
    ("PABPC1_condensate_index_cell","PABPC1 condensate index"),
]

# Per-cell point size. Only ~11 cells land in each replicate strip here (79
# across 7 replicates), so the points can be far larger than in the yeast
# figure without merging. At the published s=5 they read as specks next to the
# s=70 diamonds; the reader could not see the distribution the medians came
# from. Bigger points need a wider jitter to stay separable, so the two move
# together.
DOT_SIZE = 13.0
DOT_ALPHA = 0.45
DOT_JITTER = 0.045         # half-width of the random x-spread within a replicate
REP_OFFSET = 0.13          # x-gap between adjacent replicate strips
DIAMOND_SIZE = 70          # replicate medians -- unchanged; they are the stat unit


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
        df["rep_id"] = f.parent.name
        frames.append(df)
    cells = pd.concat(frames, ignore_index=True)
    if "keep" in cells.columns:
        cells = cells[cells["keep"] == True].copy()   # noqa: E712
    print(f"loaded {len(files)} replicates, {len(cells)} cells (keep==True)")
    return cells.reset_index(drop=True)


def rep_medians(cells: pd.DataFrame, col: str) -> pd.DataFrame:
    return (cells.groupby(["condition", "rep_id"], as_index=False)[col]
                 .median().rename(columns={col: "value"}))


def test(rm: pd.DataFrame):
    a = rm.loc[rm["condition"] == "arsenite", "value"].dropna().values
    b = rm.loc[rm["condition"] == "control", "value"].dropna().values
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan
    p = mannwhitneyu(a, b, alternative="two-sided").pvalue
    return float(p), float(np.median(b)), float(np.median(a))


def effect_size(rm: pd.DataFrame) -> dict:
    """Effect sizes for arsenite vs control, on replicate medians.

    The figure reports these instead of p-values. At n=3 vs 4 the smallest
    attainable two-sided p is 0.057, so no panel can reach 0.05 no matter how
    large the effect: a p-value here measures the replicate count, not the
    biology. Two complements are computed --

      fold   ratio of the median-of-replicate-medians. Directly readable, and
             the quantity the Results text already speaks in. Every metric here
             is on a ratio scale with a meaningful zero, so a ratio is defined.
      delta  Cliff's delta, the rank effect size belonging to the Mann-Whitney
             test actually run: the probability an arsenite replicate exceeds a
             control one, minus the reverse. Bounded [-1, 1]; +1 means complete
             separation. Coarse at this n (steps of 1/6), which is why it goes
             in the stats CSV rather than on the panel.
    """
    a = rm.loc[rm["condition"] == "arsenite", "value"].dropna().values
    b = rm.loc[rm["condition"] == "control", "value"].dropna().values
    if not len(a) or not len(b):
        return {"fold": np.nan, "cliffs_delta": np.nan}
    mb = float(np.median(b))
    fold = float(np.median(a)) / mb if mb else np.nan
    diff = np.sign(a[:, None] - b[None, :])
    return {"fold": fold, "cliffs_delta": float(diff.sum() / diff.size)}


def fmt_fold(fold: float) -> str:
    if not np.isfinite(fold):
        return "n/a"
    if fold >= 10:
        return f"{fold:.0f}×"
    # Two decimals below 1 so a small decrease does not round to "1.0x".
    return f"{fold:.1f}×" if fold >= 1 else f"{fold:.2f}×"


def floor_p(n1: int, n2: int) -> float:
    """Smallest two-sided exact p attainable at these sample sizes."""
    a, b = np.arange(n1), np.arange(n1, n1 + n2)
    return float(mannwhitneyu(b, a, alternative="two-sided", method="exact").pvalue)


def draw(ax, cells, rm, col, title, ylab, reps_by_cond, fold):
    rng = np.random.default_rng(0)
    x_of = {c: i for i, c in enumerate(ORDER)}
    for c in ORDER:
        sub = cells[cells["condition"] == c]
        for j, rid in enumerate(reps_by_cond[c]):
            v = sub.loc[sub["rep_id"] == rid, col].dropna().values
            if not len(v):
                continue
            off = (j - (len(reps_by_cond[c]) - 1) / 2) * REP_OFFSET
            ax.scatter(x_of[c] + off + rng.uniform(-DOT_JITTER, DOT_JITTER, len(v)),
                       v, s=DOT_SIZE, alpha=DOT_ALPHA, linewidths=0,
                       color=REP_COLORS[j % len(REP_COLORS)], rasterized=True)
        vals = rm.loc[rm["condition"] == c]
        for j, rid in enumerate(reps_by_cond[c]):
            r = vals.loc[vals["rep_id"] == rid, "value"]
            if r.empty or not np.isfinite(r.iloc[0]):
                continue
            off = (j - (len(reps_by_cond[c]) - 1) / 2) * REP_OFFSET
            ax.scatter(x_of[c] + off, r.iloc[0], s=DIAMOND_SIZE, marker="D",
                       color=REP_COLORS[j % len(REP_COLORS)],
                       edgecolor="k", linewidths=0.8, zorder=5)
        med = float(np.nanmedian(vals["value"].values)) if len(vals) else np.nan
        ax.plot([x_of[c] - 0.28, x_of[c] + 0.28], [med, med], "-",
                color="k", lw=1.8, zorder=6)

    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(["control", "arsenite"])
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.set_ylabel(ylab, fontsize=FONT_SIZE)
    ax.set_title(title, fontsize=FONT_SIZE, fontweight="bold")
    ax.tick_params(labelsize=FONT_SIZE)

    vals = cells[col].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    meds = rm["value"].to_numpy(dtype=float)
    meds = meds[np.isfinite(meds)]

    scale, linthresh = YSCALE.get(col, ("linear", None))
    if scale == "symlog" and vals.size:
        # No percentile clipping here: on a log axis the tail costs a decade of
        # height, not most of the panel, so every cell can stay visible.
        ax.set_yscale("symlog", linthresh=linthresh, linscale=0.6)
        top = max(float(np.nanmax(vals)), float(meds.max()) if meds.size else 0.0)
        # Drop the floor just below zero so the 10-12 cells with no puncta are
        # drawn as whole markers instead of half-clipped on the spine. Ticks are
        # pinned at and above zero, so no negative label implies negative data.
        ax.set_ylim(-linthresh * 0.45, top * 4.0)
        ax.yaxis.set_minor_locator(SymmetricalLogLocator(
            base=10.0, linthresh=linthresh, subs=list(range(2, 10))))
        ax.tick_params(axis="y", which="minor", length=1.6, width=0.6)
        ax.set_yticks([t for t in ax.get_yticks() if t >= 0])
    else:
        if vals.size:
            lo, hi = np.nanpercentile(vals, [0.5, 99.5])
            if meds.size:
                lo, hi = min(lo, meds.min()), max(hi, meds.max())
            if hi <= lo:
                hi = lo + 1.0
            pad = (hi - lo) * 0.08
            ax.set_ylim(max(0.0, lo - pad) if lo >= 0 else lo - pad, hi + pad)
        lo, hi = ax.get_ylim(); span = hi - lo
        ax.set_ylim(lo, hi + span * 0.20)

    # Bracket height in axes fractions, x in data units. The old arithmetic
    # added a fraction of the linear span, which places the bracket off-panel
    # once the axis is symlog.
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([0, 0, 1, 1], [0.90, 0.935, 0.935, 0.90], lw=0.8, color="k",
            transform=tr, clip_on=False)
    ax.text(0.5, 0.945, fmt_fold(fold), ha="center", va="bottom",
            fontsize=FONT_SIZE, fontweight="bold", transform=tr)


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--outputs-dir", type=Path,
                    default=(here.parent / "jcb_revision_midway3_cellquant_analysis"
                             / "validation_3d" / "outputs_mammalian" / "published_2d"))
    ap.add_argument("--outdir", type=Path, default=here / "figure2_2d")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(args.outputs_dir)
    cells["condition"] = cells["condition"].astype(str)
    reps_by_cond = {c: sorted(cells.loc[cells["condition"] == c, "rep_id"].unique())
                    for c in ORDER}
    n1, n2 = len(reps_by_cond["control"]), len(reps_by_cond["arsenite"])
    print("replicates:", {c: len(v) for c, v in reps_by_cond.items()})
    print(f"smallest attainable two-sided exact p at n={n1} vs {n2}: {floor_p(n1, n2):.4f}")

    fig = plt.figure(figsize=(12.5, 3.9))
    gs = GridSpec(1, 4, figure=fig, wspace=0.46)

    rows = []
    for i, (col, title, ylab) in enumerate(PANELS):
        ax = fig.add_subplot(gs[0, i])
        if col not in cells.columns:
            ax.text(0.5, 0.5, f"{col}\nabsent", ha="center", va="center", fontsize=FONT_SIZE)
            ax.set_axis_off(); print(f"  !! {col} absent"); continue
        rm = rep_medians(cells, col)
        p, mc, ma = test(rm)
        es = effect_size(rm)
        draw(ax, cells, rm, col, title, ylab, reps_by_cond, es["fold"])
        print(f"  {title:30s} control={mc:.4g}  arsenite={ma:.4g}  "
              f"fold={es['fold']:.2f}  delta={es['cliffs_delta']:+.2f}  (p={p:.4f})")
        rows.append({"metric": col, "panel": title, "control_median": mc,
                     "arsenite_median": ma, "fold_change": es["fold"],
                     "cliffs_delta": es["cliffs_delta"], "p_exact_two_sided": p,
                     "n_control": n1, "n_arsenite": n2})

    save_all(fig, args.outdir / "Figure2B_mammalian_2D")

    for col, title in TABLE_EXTRA:
        if col not in cells.columns:
            continue
        rm = rep_medians(cells, col)
        p, mc, ma = test(rm)
        es = effect_size(rm)
        print(f"  {title:30s} control={mc:.4g}  arsenite={ma:.4g}  "
              f"fold={es['fold']:.2f}  delta={es['cliffs_delta']:+.2f}  "
              f"(p={p:.4f})  (table only)")
        rows.append({"metric": col, "panel": title, "control_median": mc,
                     "arsenite_median": ma, "fold_change": es["fold"],
                     "cliffs_delta": es["cliffs_delta"], "p_exact_two_sided": p,
                     "n_control": n1, "n_arsenite": n2})

    df = pd.DataFrame(rows)
    df.to_csv(args.outdir / "table_s1_2d_stats.csv", index=False)

    # LaTeX for Table S1, same convention as the yeast tables. Effect sizes come
    # before the P-value column: at n=3 vs 4 the smallest attainable two-sided p
    # is 0.057, so the P column reports the replicate count and the effect-size
    # columns report the biology. The figure shows fold change alone; the table
    # keeps P as well, because a supplementary table hiding its test invites the
    # question more loudly than printing it does.
    lines = [r"\begin{tabular}{lccccc}", r"\toprule",
             r"Metric & Control & Arsenite & Fold change & Cliff's $\delta$ "
             r"& $P$-value \\", r"\midrule"]
    for r in rows:
        def fmt(v):
            return f"{v:.4g}" if abs(v) >= 0.001 or v == 0 else f"{v:.2e}"
        fold = r["fold_change"]
        fold_s = (f"{fold:.1f}$\\times$" if fold >= 1 else f"{fold:.2f}$\\times$") \
            if np.isfinite(fold) else "--"
        lines.append(f"{r['panel']} & {fmt(r['control_median'])} & "
                     f"{fmt(r['arsenite_median'])} & {fold_s} & "
                     f"${r['cliffs_delta']:+.2f}$ & {r['p_exact_two_sided']:.2f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (args.outdir / "table_s1_body.tex").write_text("\n".join(lines) + "\n")
    print(f"\n  wrote table_s1_2d_stats.csv and table_s1_body.tex")


if __name__ == "__main__":
    np.random.seed(0)
    main()
