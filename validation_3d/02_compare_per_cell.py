"""validation_3d/02_compare_per_cell.py

Per-cell paired comparison: matched-2D-from-3D vs 3D-native for every cell
in every image, plus scatters, Bland-Altman plots, and summary stats.

For colocalization metrics, also produces a panel coloring points by the
Z-extent of the cell — the predicted bias direction is "2D-from-MIP is
inflated for cells that span more Z planes."
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
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

# Metrics common to 2D and 3D output. 3D-only metrics (cell_volume_um3,
# nucleolar_sphericity) are excluded since there's no 2D counterpart to
# compare against; they are reported separately in the report.
PAIRED_METRICS = [
    # Cell-level basics
    "cell_area_px",  # 2D-only key; 3D adds cell_area_px=None — handle gracefully
    # Channel-level intensities
    "{ch}_cell_mean", "{ch}_nucleus_mean", "{ch}_cytosol_mean",
    "{ch}_condensate_index_cell", "{ch}_condensate_index_cytosol",
    # Puncta metrics
    "{ch}_puncta_n", "{ch}_puncta_mean_intensity",
    "{ch}_diffuse_mean_intensity", "{ch}_frac_intensity_in_puncta",
    "{ch}_puncta_integrated_intensity",
    "{ch}_fragmentation_index_simple", "{ch}_fragmentation_index_persistence",
    # Colocalization (per pair)
    "pearson_r_{pair}", "manders_m1_{pair}", "manders_m2_{pair}",
]

# Channels for the mammalian dataset (matches MAMMALIAN_CHANNELS in setup)
MAMMALIAN_QUANT_CHANNELS = ["G3BP1", "PABPC1"]
MAMMALIAN_COLOC_PAIRS = ["G3BP1_vs_PABPC1"]


def _expand_metric_templates(channels: list[str],
                             pairs: list[str]) -> list[str]:
    cols: list[str] = []
    for tmpl in PAIRED_METRICS:
        if "{ch}" in tmpl:
            for ch in channels:
                cols.append(tmpl.replace("{ch}", ch))
        elif "{pair}" in tmpl:
            for pair in pairs:
                cols.append(tmpl.replace("{pair}", pair))
        else:
            cols.append(tmpl)
    return cols


def _load_paired(stem: str, dataset: str) -> pd.DataFrame:
    """Load 3D and matched-2D cells.csv for one image, inner-join on cell_id."""
    if dataset == "mammalian":
        d3 = setup_paths.OUT_MAMM_3D / stem / "cells.csv"
        d2 = setup_paths.OUT_MAMM_MATCHED2D / stem / "cells.csv"
    else:
        raise NotImplementedError(f"dataset {dataset}")
    if not d3.exists() or not d2.exists():
        return pd.DataFrame()
    df3 = pd.read_csv(d3).add_suffix("_3d")
    df3 = df3.rename(columns={"cell_id_3d": "cell_id"})
    df2 = pd.read_csv(d2).add_suffix("_2d")
    df2 = df2.rename(columns={"cell_id_2d": "cell_id"})
    merged = df3.merge(df2, on="cell_id", how="inner")
    merged["image"] = stem
    return merged


def _scatter_with_y_eq_x(ax, x, y, color_by=None, label=""):
    ax.scatter(x, y, c=color_by if color_by is not None else "C0",
               s=20, alpha=0.6, edgecolors="none")
    lo = min(np.nanmin(x), np.nanmin(y))
    hi = max(np.nanmax(x), np.nanmax(y))
    ax.plot([lo, hi], [lo, hi], "--", color="gray", linewidth=0.8)
    ax.set_xlabel(f"2D-from-3D-mask ({label})")
    ax.set_ylabel(f"3D-native ({label})")


def _bland_altman(ax, x, y, label=""):
    mean = 0.5 * (x + y)
    diff = y - x  # 3D - 2D-from-MIP
    md = np.nanmean(diff)
    sd = np.nanstd(diff)
    ax.scatter(mean, diff, s=18, alpha=0.6, edgecolors="none")
    ax.axhline(md, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(md + 1.96 * sd, color="gray", linestyle=":", linewidth=0.8)
    ax.axhline(md - 1.96 * sd, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel(f"Mean of 2D-from-MIP and 3D ({label})")
    ax.set_ylabel("3D - 2D-from-MIP")


def per_cell_summary(dataset: str = "mammalian") -> pd.DataFrame:
    """Build paired summary across all images for one dataset.

    Output: DataFrame with one row per metric, plus per-image scatter+BA panels
    written to outputs_{dataset}/comparisons/.
    """
    if dataset == "mammalian":
        zstacks = setup_paths.list_mammalian_zstacks()
        out_dir = setup_paths.OUT_MAMM_COMPARISONS
        channels = MAMMALIAN_QUANT_CHANNELS
        coloc_pairs = MAMMALIAN_COLOC_PAIRS
    else:
        raise NotImplementedError(f"dataset {dataset}")
    out_dir.mkdir(parents=True, exist_ok=True)

    paired_dfs = []
    for p in zstacks:
        m = _load_paired(p.stem, dataset)
        if not m.empty:
            paired_dfs.append(m)
    if not paired_dfs:
        return pd.DataFrame()
    all_paired = pd.concat(paired_dfs, ignore_index=True)
    all_paired.to_csv(out_dir / "paired_all_cells.csv", index=False)

    metric_cols = _expand_metric_templates(channels, coloc_pairs)

    summary_rows = []
    for metric in metric_cols:
        col_3d = f"{metric}_3d"
        col_2d = f"{metric}_2d"
        if col_3d not in all_paired.columns or col_2d not in all_paired.columns:
            continue
        x = all_paired[col_2d].astype(float).values
        y = all_paired[col_3d].astype(float).values
        mask = ~(np.isnan(x) | np.isnan(y))
        x, y = x[mask], y[mask]
        if len(x) < 5:
            continue
        rho, _ = spearmanr(x, y)
        rel_err = np.abs(y - x) / (np.abs(y) + np.abs(x) + 1e-9) * 2.0
        summary_rows.append({
            "metric": metric,
            "n_cells": len(x),
            "spearman_r": float(rho),
            "mean_bias_2d_minus_3d": float(np.mean(x - y)),
            "mean_3d": float(np.mean(y)),
            "mean_2d": float(np.mean(x)),
            "median_abs_relative_error": float(np.median(rel_err)),
            "p90_abs_relative_error": float(np.percentile(rel_err, 90)),
        })
        # Scatter + Bland-Altman per metric
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        _scatter_with_y_eq_x(axes[0], x, y, label=metric)
        axes[0].set_title(
            f"{metric}\nSpearman={rho:.3f}, n={len(x)}")
        _bland_altman(axes[1], x, y, label=metric)
        axes[1].set_title(
            f"Bland-Altman: 3D - 2D-from-MIP\n"
            f"bias={np.nanmean(y-x):.3g}")
        fig.tight_layout()
        # Sanitize metric name for filename
        safe = metric.replace("/", "_").replace(" ", "_")
        fig.savefig(out_dir / f"per_cell_{safe}.png", dpi=120)
        plt.close(fig)

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        # Agreement classification
        def _classify(row):
            if row["spearman_r"] >= 0.9 and row["median_abs_relative_error"] <= 0.10:
                return "high agreement"
            if row["spearman_r"] >= 0.7 and row["median_abs_relative_error"] <= 0.25:
                return "moderate agreement"
            return "discrepant"
        summary["agreement"] = summary.apply(_classify, axis=1)
    summary.to_csv(out_dir / "per_cell_summary.csv", index=False)
    return summary


def main() -> None:
    print("=== Mammalian per-cell paired comparison ===")
    summary = per_cell_summary("mammalian")
    if summary.empty:
        print("(no paired data — did the 3D and matched-2D passes complete?)")
        return
    print(summary.to_string(index=False))
    print(f"\nWrote {len(summary)} per-metric panels + summary to "
          f"{setup_paths.OUT_MAMM_COMPARISONS}/")


if __name__ == "__main__":
    main()
