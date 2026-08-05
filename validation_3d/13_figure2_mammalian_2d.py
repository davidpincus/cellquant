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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

matplotlib.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8, "axes.unicode_minus": False,
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
# Extra metrics reported in Table S1 but not plotted.
TABLE_EXTRA = [
    ("PABPC1_puncta_n",             "PABPC1 puncta per cell"),
    ("G3BP1_condensate_index_cell", "G3BP1 condensate index"),
    ("PABPC1_condensate_index_cell","PABPC1 condensate index"),
]


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


def floor_p(n1: int, n2: int) -> float:
    """Smallest two-sided exact p attainable at these sample sizes."""
    a, b = np.arange(n1), np.arange(n1, n1 + n2)
    return float(mannwhitneyu(b, a, alternative="two-sided", method="exact").pvalue)


def draw(ax, cells, rm, col, title, ylab, reps_by_cond, p):
    rng = np.random.default_rng(0)
    x_of = {c: i for i, c in enumerate(ORDER)}
    for c in ORDER:
        sub = cells[cells["condition"] == c]
        for j, rid in enumerate(reps_by_cond[c]):
            v = sub.loc[sub["rep_id"] == rid, col].dropna().values
            if not len(v):
                continue
            off = (j - (len(reps_by_cond[c]) - 1) / 2) * 0.13
            ax.scatter(x_of[c] + off + rng.uniform(-0.035, 0.035, len(v)), v,
                       s=5, alpha=0.35, linewidths=0,
                       color=REP_COLORS[j % len(REP_COLORS)], rasterized=True)
        vals = rm.loc[rm["condition"] == c]
        for j, rid in enumerate(reps_by_cond[c]):
            r = vals.loc[vals["rep_id"] == rid, "value"]
            if r.empty or not np.isfinite(r.iloc[0]):
                continue
            off = (j - (len(reps_by_cond[c]) - 1) / 2) * 0.13
            ax.scatter(x_of[c] + off, r.iloc[0], s=70, marker="D",
                       color=REP_COLORS[j % len(REP_COLORS)],
                       edgecolor="k", linewidths=0.8, zorder=5)
        med = float(np.nanmedian(vals["value"].values)) if len(vals) else np.nan
        ax.plot([x_of[c] - 0.28, x_of[c] + 0.28], [med, med], "-",
                color="k", lw=1.8, zorder=6)

    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(["control", "arsenite"])
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.set_ylabel(ylab, fontsize=8)
    ax.set_title(title, fontsize=9.5, fontweight="bold")
    ax.tick_params(labelsize=8)

    vals = cells[col].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size:
        lo, hi = np.nanpercentile(vals, [0.5, 99.5])
        m = rm["value"].to_numpy(dtype=float); m = m[np.isfinite(m)]
        if m.size:
            lo, hi = min(lo, m.min()), max(hi, m.max())
        if hi <= lo:
            hi = lo + 1.0
        pad = (hi - lo) * 0.08
        ax.set_ylim(max(0.0, lo - pad) if lo >= 0 else lo - pad, hi + pad)

    lo, hi = ax.get_ylim(); span = hi - lo
    ax.set_ylim(lo, hi + span * 0.20)
    y = hi + span * 0.05
    ax.plot([0, 0, 1, 1], [y, y + span * 0.03, y + span * 0.03, y], lw=0.8, color="k")
    lbl = "n.s." if (not np.isfinite(p) or p >= 0.05) else f"p = {p:.3f}"
    if np.isfinite(p):
        lbl = f"p = {p:.2f}" if p >= 0.01 else f"p = {p:.3f}"
    ax.text(0.5, y + span * 0.045, lbl, ha="center", va="bottom", fontsize=7.5)


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

    fig = plt.figure(figsize=(11, 3.2))
    gs = GridSpec(1, 4, figure=fig, wspace=0.42)

    rows = []
    for i, (col, title, ylab) in enumerate(PANELS):
        ax = fig.add_subplot(gs[0, i])
        if col not in cells.columns:
            ax.text(0.5, 0.5, f"{col}\nabsent", ha="center", va="center", fontsize=8)
            ax.set_axis_off(); print(f"  !! {col} absent"); continue
        rm = rep_medians(cells, col)
        p, mc, ma = test(rm)
        draw(ax, cells, rm, col, title, ylab, reps_by_cond, p)
        print(f"  {title:30s} control={mc:.4g}  arsenite={ma:.4g}  p={p:.4f}")
        rows.append({"metric": col, "panel": title, "control_median": mc,
                     "arsenite_median": ma, "p_exact_two_sided": p,
                     "n_control": n1, "n_arsenite": n2})

    save_all(fig, args.outdir / "Figure2B_mammalian_2D")

    for col, title in TABLE_EXTRA:
        if col not in cells.columns:
            continue
        rm = rep_medians(cells, col)
        p, mc, ma = test(rm)
        print(f"  {title:30s} control={mc:.4g}  arsenite={ma:.4g}  p={p:.4f}  (table only)")
        rows.append({"metric": col, "panel": title, "control_median": mc,
                     "arsenite_median": ma, "p_exact_two_sided": p,
                     "n_control": n1, "n_arsenite": n2})

    df = pd.DataFrame(rows)
    df.to_csv(args.outdir / "table_s1_2d_stats.csv", index=False)

    # LaTeX for Table S1, same convention as the yeast tables.
    lines = [r"\begin{tabular}{lccc}", r"\toprule",
             r"Metric & Control & Arsenite & $P$-value \\", r"\midrule"]
    for r in rows:
        def fmt(v):
            return f"{v:.4g}" if abs(v) >= 0.001 or v == 0 else f"{v:.2e}"
        lines.append(f"{r['panel']} & {fmt(r['control_median'])} & "
                     f"{fmt(r['arsenite_median'])} & {r['p_exact_two_sided']:.2f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (args.outdir / "table_s1_body.tex").write_text("\n".join(lines) + "\n")
    print(f"\n  wrote table_s1_2d_stats.csv and table_s1_body.tex")


if __name__ == "__main__":
    np.random.seed(0)
    main()
