"""validation_3d/03_compare_aggregates.py

Per-condition aggregate recapitulation: does the 3D pipeline recover the
published biological conclusion?

Mammalian: replicate-level Wilcoxon rank-sum on each metric for
arsenite vs control. Produces a comparison table across {3D, matched-2D,
published-2D}, with effect direction and corrected p-value per axis.

Output:
  outputs_mammalian/comparisons/aggregate_recapitulation.csv
  outputs_mammalian/comparisons/aggregate_bars.png
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
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

MAMMALIAN_METRICS_OF_INTEREST = [
    "G3BP1_puncta_n", "G3BP1_frac_intensity_in_puncta",
    "G3BP1_condensate_index_cell",
    "PABPC1_puncta_n", "PABPC1_frac_intensity_in_puncta",
    "PABPC1_condensate_index_cell",
    "pearson_r_G3BP1_vs_PABPC1",
    "manders_m1_G3BP1_vs_PABPC1",
    "manders_m2_G3BP1_vs_PABPC1",
    "G3BP1_fragmentation_index_persistence",
    "PABPC1_fragmentation_index_persistence",
]


def _load_axis_cells(axis: str) -> pd.DataFrame:
    """Concatenate cells.csv from every image of one axis.

    axis ∈ {"3d", "matched_2d", "published_2d"}
    """
    if axis == "3d":
        root = setup_paths.OUT_MAMM_3D
    elif axis == "matched_2d":
        root = setup_paths.OUT_MAMM_MATCHED2D
    elif axis == "published_2d":
        root = setup_paths.OUT_MAMM_PUBLISHED2D
    else:
        raise ValueError(axis)

    frames = []
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        cf = sub / "cells.csv"
        if not cf.exists():
            continue
        df = pd.read_csv(cf)
        # Ensure condition + replicate columns present (cellquant writes
        # them via filename pattern). If missing, infer from directory name.
        if "condition" not in df.columns or df["condition"].isna().all():
            cond, _, rep = sub.name.rpartition("_rep")
            df["condition"] = cond
            df["replicate"] = rep
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _wilcoxon_rank_sum_replicate(cells: pd.DataFrame,
                                 metric: str,
                                 cond_a: str,
                                 cond_b: str) -> dict[str, float]:
    """Replicate-level Wilcoxon rank-sum on the per-image medians of `metric`.

    Mirrors the manuscript's statistical convention (Lord et al. superplots).
    Returns {p_value, n_a, n_b, median_a, median_b, sign}.
    """
    if metric not in cells.columns:
        return {"p_value": np.nan, "n_a": 0, "n_b": 0,
                "median_a": np.nan, "median_b": np.nan, "sign": ""}
    # keep gate from cellquant: include only "kept" cells if column present
    if "keep" in cells.columns:
        df = cells[cells["keep"]].copy()
    else:
        df = cells.copy()
    # Coerce metric to numeric (strings get NaN; lots of NaNs is OK)
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    # Per-image median
    rep_meds = (
        df.groupby(["condition", "image"], as_index=False)[metric]
          .median()
          .rename(columns={metric: "rep_median"})
    )
    a = rep_meds.loc[rep_meds["condition"] == cond_a, "rep_median"].dropna().values
    b = rep_meds.loc[rep_meds["condition"] == cond_b, "rep_median"].dropna().values
    if len(a) < 2 or len(b) < 2:
        return {"p_value": np.nan, "n_a": len(a), "n_b": len(b),
                "median_a": np.nanmedian(a) if len(a) else np.nan,
                "median_b": np.nanmedian(b) if len(b) else np.nan,
                "sign": "ns (small N)"}
    try:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
    except Exception:
        p = np.nan
    med_a = float(np.median(a))
    med_b = float(np.median(b))
    sign = "+" if med_b > med_a else ("-" if med_b < med_a else "0")
    return {"p_value": float(p), "n_a": len(a), "n_b": len(b),
            "median_a": med_a, "median_b": med_b, "sign": sign}


def mammalian_aggregate_table() -> pd.DataFrame:
    axes = {"3d": _load_axis_cells("3d"),
            "matched_2d": _load_axis_cells("matched_2d"),
            "published_2d": _load_axis_cells("published_2d")}
    rows = []
    for metric in MAMMALIAN_METRICS_OF_INTEREST:
        row = {"metric": metric}
        for axis, df in axes.items():
            if df.empty:
                row[f"{axis}_n_a"] = 0
                row[f"{axis}_n_b"] = 0
                row[f"{axis}_median_control"] = np.nan
                row[f"{axis}_median_arsenite"] = np.nan
                row[f"{axis}_sign"] = "no data"
                row[f"{axis}_p"] = np.nan
                continue
            stat = _wilcoxon_rank_sum_replicate(df, metric, "control", "arsenite")
            row[f"{axis}_n_a"] = stat["n_a"]
            row[f"{axis}_n_b"] = stat["n_b"]
            row[f"{axis}_median_control"] = stat["median_a"]
            row[f"{axis}_median_arsenite"] = stat["median_b"]
            row[f"{axis}_sign"] = stat["sign"]
            row[f"{axis}_p"] = stat["p_value"]
        # Concordance: are the three axes pointing the same direction?
        signs = [row.get(f"{a}_sign") for a in ("3d", "matched_2d", "published_2d")]
        signs = [s for s in signs if s in ("+", "-")]
        row["all_axes_agree_direction"] = len(set(signs)) <= 1 and len(signs) > 0
        rows.append(row)
    return pd.DataFrame(rows)


def _bar_chart(table: pd.DataFrame, out_path: Path) -> None:
    """Side-by-side bar of median differences (arsenite - control) per axis."""
    fig, ax = plt.subplots(figsize=(11, 5))
    metrics = table["metric"].tolist()
    width = 0.27
    x = np.arange(len(metrics))
    for i, axis in enumerate(("3d", "matched_2d", "published_2d")):
        diffs = (table[f"{axis}_median_arsenite"]
                 - table[f"{axis}_median_control"])
        ax.bar(x + (i - 1) * width, diffs, width, label=axis)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=45, ha="right", fontsize=8)
    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_ylabel("Median (arsenite) - median (control)")
    ax.set_title("Mammalian: arsenite effect direction by axis")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def main() -> None:
    out_dir = setup_paths.OUT_MAMM_COMPARISONS
    out_dir.mkdir(parents=True, exist_ok=True)
    table = mammalian_aggregate_table()
    table.to_csv(out_dir / "aggregate_recapitulation.csv", index=False)
    if not table.empty:
        _bar_chart(table, out_dir / "aggregate_bars.png")
    print(table.to_string(index=False))
    print(f"\nWrote aggregate table + bars to {out_dir}/")


if __name__ == "__main__":
    main()
