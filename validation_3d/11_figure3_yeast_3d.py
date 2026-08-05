"""validation_3d/11_figure3_yeast_3d.py

Rebuild manuscript Figure 3B (yeast temperature-series quantification) on the
NATIVE 3D metrics, in the published superplot/trend style: per-cell points
coloured by replicate, replicate medians as large diamonds, condition medians
joined by a line, significance brackets vs the 25 C reference.

Panel mapping, published 2D -> 3D:
    cell area (px)              -> cell volume (um^3)
    nucleolar area (px)         -> nucleolar volume (um^3)
    nucleolar circularity       -> NO 3D ANALOG.  Circularity, solidity and
                                   eccentricity are properties of a projection;
                                   sphericity is not measurable at this axial
                                   resolution and was removed from cellquant.
                                   Substituted: Sis1 fraction of puncta within
                                   0.5 um of the nucleolus (3D-native, and a
                                   natural companion to the distance panel).
    Sis1 mean distance (px)     -> Sis1 distance to nucleolus (um)
    Sis1/Nsr1/Tif6 cell mean    -> unchanged
    Pearson R Tif6 vs Nsr1      -> computed natively on the z-stack

Statistics: two-sided exact Mann-Whitney U on replicate medians, each
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
from matplotlib.gridspec import GridSpec

matplotlib.rcParams.update({
    "svg.fonttype": "none",     # editable text in Illustrator
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.unicode_minus": False,
})

TEMP_ORDER = ["25C", "30C", "32C", "36C", "40C"]
REF = "25C"
N_COMPARISONS = 4          # 30/32/36/40 vs 25 -> Bonferroni denominator

# (column, panel title, y-axis label)
PANELS = [
    ("cell_volume_um3",         "Cell volume",              "volume (µm$^3$)"),
    ("nucleolar_volume_um3",    "Nucleolar volume",         "volume (µm$^3$)"),
    ("Sis1_mean_distance",      "Sis1 puncta to nucleolus", "distance (µm)"),
    ("Sis1_fraction_proximal",  "Sis1 puncta proximal",     "fraction $\\leq$0.5 µm"),
    ("Sis1_cell_mean",          "Sis1 intensity",           "mean intensity (AU)"),
    ("Nsr1_cell_mean",          "Nsr1 intensity",           "mean intensity (AU)"),
    ("Tif6_cell_mean",          "Tif6 intensity",           "mean intensity (AU)"),
    ("pearson_r_Tif6_vs_Nsr1",  "Nsr1–Tif6 colocalization", "Pearson $R$"),
]

REP_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
              "#937860", "#DA8BC3", "#8C8C8C"]


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
        p = mannwhitneyu(v, ref, alternative="two-sided", method="exact").pvalue
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
            off = (j - (len(reps_by_cond[t]) - 1) / 2) * 0.075
            x = x_of[t] + off + rng.uniform(-0.018, 0.018, len(v))
            ax.scatter(x, v, s=1.6, alpha=0.28, linewidths=0,
                       color=REP_COLORS[j % len(REP_COLORS)], rasterized=True)

    # replicate medians as diamonds + condition median line
    cond_med = []
    for t in TEMP_ORDER:
        vals = rm.loc[rm["condition"] == t]
        for j, rid in enumerate(reps_by_cond[t]):
            r = vals.loc[vals["rep_id"] == rid, "value"]
            if r.empty or not np.isfinite(r.iloc[0]):
                continue
            off = (j - (len(reps_by_cond[t]) - 1) / 2) * 0.075
            ax.scatter(x_of[t] + off, r.iloc[0], s=52, marker="D",
                       color=REP_COLORS[j % len(REP_COLORS)],
                       edgecolor="k", linewidths=0.7, zorder=5)
        cond_med.append(float(np.nanmedian(vals["value"].values)) if len(vals) else np.nan)
    ax.plot(range(len(TEMP_ORDER)), cond_med, "-", color="k", lw=1.3, zorder=6)
    ax.scatter(range(len(TEMP_ORDER)), cond_med, s=26, color="k", zorder=7)

    ax.set_xticks(range(len(TEMP_ORDER)))
    ax.set_xticklabels([t.replace("C", "") for t in TEMP_ORDER])
    ax.set_xlabel("Temperature (°C)", fontsize=8)
    ax.set_ylabel(ylab, fontsize=8)
    ax.set_title(title, fontsize=9.5, fontweight="bold")
    ax.tick_params(labelsize=8)

    # Scale to the bulk of the per-cell distribution, not to its outliers.
    # Auto-limits driven by a handful of extreme cells squash every trend into
    # the bottom of the panel; clip at the 0.5/99.5 percentiles and make sure
    # all replicate medians stay in view.
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

    # significance brackets, stacked in headroom added above the data
    sig = [(t, st[t]) for t in TEMP_ORDER if t != REF and st[t]["sig"]]
    if sig:
        lo, hi = ax.get_ylim()
        span = hi - lo
        ax.set_ylim(lo, hi + span * (0.13 * len(sig) + 0.06))
        for k, (t, s) in enumerate(sig):
            y = hi + span * (0.05 + 0.13 * k)
            x0, x1 = x_of[REF], x_of[t]
            ax.plot([x0, x0, x1, x1],
                    [y, y + span * 0.025, y + span * 0.025, y],
                    lw=0.8, color="k")
            ax.text((x0 + x1) / 2, y + span * 0.035, fmt_p(s["p_adj"]),
                    ha="center", va="bottom", fontsize=7)


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

    fig = plt.figure(figsize=(14, 7))
    gs = GridSpec(2, 4, figure=fig, hspace=0.55, wspace=0.34)

    rows = []
    for i, (col, title, ylab) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        if col not in cells.columns:
            ax.text(0.5, 0.5, f"{col}\nnot in 3D output", ha="center", va="center",
                    fontsize=8, color="0.5")
            ax.set_axis_off()
            print(f"  !! {col} absent")
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
            rows.append({"panel": title, "metric": col, "temperature": t,
                         "replicate_median_of_medians": med[t],
                         "n_replicates": len(reps_by_cond[t]),
                         "p_raw_vs_25C": st[t]["p_raw"] if t != REF else np.nan,
                         "p_bonferroni_vs_25C": st[t]["p_adj"] if t != REF else np.nan})

    save_all(fig, args.outdir / "Figure3B_yeast_3D")
    pd.DataFrame(rows).to_csv(args.outdir / "figure3b_3d_stats.csv", index=False)
    print(f"\n  wrote figure3b_3d_stats.csv")


if __name__ == "__main__":
    np.random.seed(0)
    main()
