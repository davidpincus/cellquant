"""validation_3d/05_compare_mammalian_coloc.py

Focused colocalization comparison for the mammalian HCT116 stress-granule data
(R3 rebuttal). Pair = G3BP1_vs_PABPC1. Compartment = whole-cell (matches the
published mammalian run; see config_used.yml).

Two axes from identical source voxels / same cellquant version:
  - 3d          : native z-stack voxel-wise Pearson/Manders per cell.
  - matched_2d  : SAME 3D cell masks projected to XY, channels MIP'd, 2D coloc.
  - published_2d: secondary axis (independent 2D segmentation on the MIP); its
                  replicates are NOT the same cells, so it appears in the
                  per-condition table only, not in the paired MIP-3D delta.

Reuses the existing colocalization.csv outputs (per-cell coloc is written to its
own file by cellquant). Nothing is recomputed here.

Stats are replicate-level only (honest n):
  - Mann-Whitney U: arsenite vs control on per-replicate medians, per axis.
  - Wilcoxon signed-rank: per-replicate median (MIP - 3D) deltas vs 0.

Outputs:
  outputs_mammalian/comparisons/coloc_replicate_table.csv
  outputs_mammalian/comparisons/coloc_paired_deltas.csv
  05_mammalian_coloc_comparison.png
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

PAIR = "G3BP1_vs_PABPC1"
METRICS = ["pearson_r", "manders_m1", "manders_m2"]
CONDITIONS = ["control", "arsenite"]
AXES = {
    "3d": setup_paths.OUT_MAMM_3D,
    "matched_2d": setup_paths.OUT_MAMM_MATCHED2D,
    "published_2d": setup_paths.OUT_MAMM_PUBLISHED2D,
}
# tab10 by replicate number (mirrors cellquant superplot styling)
REP_CMAP = plt.get_cmap("tab10")


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def load_axis(axis_dir: Path) -> pd.DataFrame:
    """Concatenate every rep's colocalization.csv under an axis dir."""
    frames = []
    for sub in sorted(axis_dir.glob("*/")):
        f = sub / "colocalization.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        df = df[df["pair"] == PAIR].copy()
        if df.empty:
            continue
        df["rep_id"] = sub.name  # e.g. arsenite_rep1
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Replicate-level summary table
# ---------------------------------------------------------------------------
def replicate_table(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for axis, df in data.items():
        if df.empty:
            continue
        for cond in CONDITIONS:
            sub = df[df["condition"] == cond]
            if sub.empty:
                continue
            for metric in METRICS:
                # median over cells within each rep, then median + IQR across reps
                rep_med = sub.groupby("rep_id")[metric].median()
                q1, q3 = np.percentile(rep_med.values, [25, 75])
                rows.append({
                    "axis": axis,
                    "condition": cond,
                    "metric": metric,
                    "rep_median_of_medians": float(np.median(rep_med.values)),
                    "iqr_low": float(q1),
                    "iqr_high": float(q3),
                    "n_reps": int(rep_med.size),
                    "n_cells": int(sub.shape[0]),
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Paired MIP - 3D deltas (same cells)
# ---------------------------------------------------------------------------
def paired_deltas(df3d: pd.DataFrame, dfm: pd.DataFrame) -> pd.DataFrame:
    """Per cell (MIP - 3D), then per-replicate median delta. Inner-join cells."""
    keys = ["rep_id", "condition", "cell_id"]
    a = df3d[keys + METRICS].rename(columns={m: f"{m}_3d" for m in METRICS})
    b = dfm[keys + METRICS].rename(columns={m: f"{m}_2d" for m in METRICS})
    merged = a.merge(b, on=keys, how="inner")
    for m in METRICS:
        merged[f"{m}_delta"] = merged[f"{m}_2d"] - merged[f"{m}_3d"]
    return merged


def per_rep_median_delta(paired: pd.DataFrame) -> pd.DataFrame:
    agg = {f"{m}_delta": "median" for m in METRICS}
    g = paired.groupby(["rep_id", "condition"], as_index=False).agg(agg)
    return g


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------
def stats_block(data: dict[str, pd.DataFrame],
                rep_delta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # Mann-Whitney arsenite vs control on per-replicate medians, per axis/metric
    for axis, df in data.items():
        if df.empty:
            continue
        for metric in METRICS:
            ctrl = (df[df["condition"] == "control"]
                    .groupby("rep_id")[metric].median().values)
            ars = (df[df["condition"] == "arsenite"]
                   .groupby("rep_id")[metric].median().values)
            if len(ctrl) < 1 or len(ars) < 1:
                continue
            try:
                u, p = mannwhitneyu(ars, ctrl, alternative="two-sided")
            except ValueError:
                u, p = np.nan, np.nan
            rows.append({
                "test": "MannWhitney_arsenite_vs_control",
                "axis": axis, "metric": metric,
                "stat": float(u) if u == u else np.nan,
                "p_value": float(p) if p == p else np.nan,
                "n_arsenite_reps": len(ars), "n_control_reps": len(ctrl),
                "median_arsenite": float(np.median(ars)),
                "median_control": float(np.median(ctrl)),
                "arsenite_gt_control": bool(np.median(ars) > np.median(ctrl)),
            })
    # Wilcoxon signed-rank: per-replicate median (MIP - 3D) delta vs 0 (all reps)
    for metric in METRICS:
        d = rep_delta[f"{metric}_delta"].dropna().values
        try:
            w, p = wilcoxon(d, alternative="two-sided")
        except ValueError:
            w, p = np.nan, np.nan
        rows.append({
            "test": "Wilcoxon_MIPminus3D_vs_0",
            "axis": "matched_2d_vs_3d", "metric": metric,
            "stat": float(w) if w == w else np.nan,
            "p_value": float(p) if p == p else np.nan,
            "n_reps": len(d),
            "median_delta_MIP_minus_3D": float(np.median(d)),
            "mean_delta_MIP_minus_3D": float(np.mean(d)),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure: 1x3 superplot (Pearson | M1 | M2), matched + 3d paired per condition
# ---------------------------------------------------------------------------
def make_figure(data: dict[str, pd.DataFrame], rep_delta: pd.DataFrame,
                out_png: Path) -> None:
    df3d, dfm = data["3d"], data["matched_2d"]
    # x layout: control_3D=0, control_MIP=0.55, arsenite_3D=1.6, arsenite_MIP=2.15
    layout = {
        ("control", "3d"): 0.0, ("control", "matched_2d"): 0.55,
        ("arsenite", "3d"): 1.6, ("arsenite", "matched_2d"): 2.15,
    }
    axis_df = {"3d": df3d, "matched_2d": dfm}
    rng = np.random.default_rng(0)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
    for ax, metric in zip(axes, METRICS):
        for (cond, axis), xc in layout.items():
            sub = axis_df[axis]
            sub = sub[sub["condition"] == cond]
            if sub.empty:
                continue
            # small points colored by replicate
            for rid, grp in sub.groupby("rep_id"):
                repnum = int(str(rid).split("rep")[-1])
                vals = grp[metric].dropna().astype(float).values
                xj = xc + rng.normal(0, 0.04, size=len(vals))
                ax.scatter(xj, vals, s=10, alpha=0.55,
                           color=REP_CMAP((repnum - 1) % 10), edgecolors="none")
            # large diamonds = replicate medians
            rep_med = sub.groupby("rep_id")[metric].median()
            offs = (np.linspace(-0.08, 0.08, len(rep_med))
                    if len(rep_med) > 1 else np.array([0.0]))
            for (rid, mv), o in zip(rep_med.items(), offs):
                repnum = int(str(rid).split("rep")[-1])
                ax.scatter(xc + o, mv, s=140, marker="D", linewidths=1.2,
                           edgecolors="black",
                           facecolors=REP_CMAP((repnum - 1) % 10), zorder=5)
        ax.set_xticks([0.275, 1.875])
        ax.set_xticklabels(["control", "arsenite"])
        for xc in (0.0, 0.55, 1.6, 2.15):
            ax.text(xc, -0.06, "3D" if xc in (0.0, 1.6) else "MIP",
                    ha="center", va="top", fontsize=7, color="gray",
                    transform=ax.get_xaxis_transform())
        ax.set_ylabel(metric)
        ax.set_title(metric)
        # headroom so the annotation box clears the top diamonds
        y0, y1 = ax.get_ylim()
        ax.set_ylim(y0, y1 + 0.20 * (y1 - y0))
        # annotate MIP - 3D shift per condition (per-rep median delta, mean over reps)
        ann = []
        for cond in CONDITIONS:
            d = rep_delta.loc[rep_delta["condition"] == cond,
                              f"{metric}_delta"].values
            if len(d):
                ann.append(f"{cond}: MIP-3D={np.median(d):+.3f}")
        ax.annotate("\n".join(ann), xy=(0.5, 0.98), xycoords="axes fraction",
                    ha="center", va="top", fontsize=8,
                    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    fig.suptitle(f"Mammalian HCT116 coloc ({PAIR}), whole-cell: MIP vs 3D",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main() -> None:
    data = {axis: load_axis(d) for axis, d in AXES.items()}
    out_dir = setup_paths.OUT_MAMM_COMPARISONS
    out_dir.mkdir(parents=True, exist_ok=True)

    table = replicate_table(data)
    paired = paired_deltas(data["3d"], data["matched_2d"])
    rep_delta = per_rep_median_delta(paired)
    stats = stats_block(data, rep_delta)

    table_path = out_dir / "coloc_replicate_table.csv"
    delta_path = out_dir / "coloc_paired_deltas.csv"
    stats_path = out_dir / "coloc_stats.csv"
    fig_path = setup_paths.VALIDATION_ROOT / "05_mammalian_coloc_comparison.png"
    table.to_csv(table_path, index=False)
    rep_delta.to_csv(delta_path, index=False)
    stats.to_csv(stats_path, index=False)
    make_figure(data, rep_delta, fig_path)

    # --- printed output: table, p-values, paths only ---
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    print("=== COMPARISON TABLE (axis x condition x metric; replicate medians) ===")
    show = table.copy()
    show["median (IQR)"] = show.apply(
        lambda r: f"{r['rep_median_of_medians']:.3f} "
                  f"({r['iqr_low']:.3f}-{r['iqr_high']:.3f})", axis=1)
    print(show[["axis", "condition", "metric", "median (IQR)",
                "n_reps", "n_cells"]].to_string(index=False))

    print("\n=== PER-REPLICATE MEDIAN DELTA (MIP - 3D) ===")
    dd = rep_delta.copy()
    print(dd.to_string(index=False))

    print("\n=== P-VALUES ===")
    print(stats.to_string(index=False))

    # Direction-flip flag
    print("\n=== DIRECTION CHECK (arsenite vs control, 3d vs matched_2d) ===")
    for metric in METRICS:
        d3 = table[(table.axis == "3d") & (table.metric == metric)]
        dm = table[(table.axis == "matched_2d") & (table.metric == metric)]
        def _gt(t):
            a = t[t.condition == "arsenite"]["rep_median_of_medians"].values
            c = t[t.condition == "control"]["rep_median_of_medians"].values
            return (a[0] > c[0]) if len(a) and len(c) else None
        g3, gm = _gt(d3), _gt(dm)
        flag = "  <-- FLIP" if (g3 is not None and gm is not None and g3 != gm) else ""
        print(f"  {metric}: 3d arsenite>control={g3}; matched_2d arsenite>control={gm}{flag}")

    print("\n=== OUTPUT PATHS ===")
    for p in (table_path, delta_path, stats_path, fig_path):
        print(f"  {p}")


if __name__ == "__main__":
    main()
