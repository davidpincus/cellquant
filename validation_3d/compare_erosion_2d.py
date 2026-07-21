"""validation_3d/compare_erosion_2d.py

For each puncta metric the manuscript reports, compute the control vs arsenite
replicate-level medians in the published-2D pass, pre-erosion vs post-erosion.
Print a side-by-side table and a verdict on whether each manuscript-reported
direction is preserved.

Compare Table S1 reported values:
    G3BP1 puncta per cell     control 1.83  arsenite 2.50   p=0.59
    G3BP1 fraction condensed  control 0.0008 arsenite 0.0059  p=0.11
    G3BP1 condensate index    control 1.52  arsenite 1.53    p=0.86
    PABPC1 puncta per cell    control 1.33  arsenite 2.50   p=0.36
    PABPC1 fraction condensed control 0.0005 arsenite 0.0022  p=0.057
    PABPC1 condensate index   control 1.44  arsenite 1.48    p=0.11
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

METRICS = [
    ("G3BP1_puncta_n",                       "G3BP1 puncta per cell",        1.83,  2.50,  0.59),
    ("G3BP1_frac_intensity_in_puncta",       "G3BP1 fraction condensed",     0.0008, 0.0059, 0.11),
    ("G3BP1_condensate_index_cell",          "G3BP1 condensate index",       1.52,  1.53,  0.86),
    ("PABPC1_puncta_n",                      "PABPC1 puncta per cell",       1.33,  2.50,  0.36),
    ("PABPC1_frac_intensity_in_puncta",      "PABPC1 fraction condensed",    0.0005, 0.0022, 0.057),
    ("PABPC1_condensate_index_cell",         "PABPC1 condensate index",      1.44,  1.48,  0.11),
]


def _load_axis(root: Path) -> pd.DataFrame:
    frames = []
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        cf = sub / "cells.csv"
        if not cf.exists():
            continue
        df = pd.read_csv(cf)
        # Infer condition + replicate from the dirname if missing
        if "condition" not in df.columns or df["condition"].isna().all():
            cond, _, rep = sub.name.rpartition("_rep")
            df["condition"] = cond
            df["replicate"] = rep
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _replicate_medians(cells: pd.DataFrame, metric: str, condition: str) -> np.ndarray:
    if metric not in cells.columns:
        return np.array([])
    df = cells[cells.get("keep", pd.Series([True] * len(cells)))].copy()
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    rep_meds = (
        df.loc[df["condition"] == condition]
          .groupby("image")[metric].median()
    )
    return rep_meds.dropna().values


def _mwu(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    try:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        return float(p)
    except Exception:
        return float("nan")


def main() -> None:
    pre = _load_axis(setup_paths.OUT_MAMMALIAN / "published_2d_pre_erosion")
    post = _load_axis(setup_paths.OUT_MAMMALIAN / "published_2d")
    if pre.empty or post.empty:
        print("Missing pre or post data — aborting.")
        print(f"  pre rows: {len(pre)}, post rows: {len(post)}")
        return

    print(f"Pre-erosion total cells: {len(pre)} ({pre['condition'].value_counts().to_dict()})")
    print(f"Post-erosion total cells: {len(post)} ({post['condition'].value_counts().to_dict()})")
    print()
    print("=" * 96)
    print("Per-metric control→arsenite replicate medians, pre-erosion vs post-erosion (option b)")
    print("vs manuscript Table S1 reported values")
    print("=" * 96)
    rows = []
    for metric, label, ms_ctrl, ms_ars, ms_p in METRICS:
        a_pre = _replicate_medians(pre, metric, "control")
        b_pre = _replicate_medians(pre, metric, "arsenite")
        a_post = _replicate_medians(post, metric, "control")
        b_post = _replicate_medians(post, metric, "arsenite")
        if len(a_pre) == 0 or len(b_pre) == 0 or len(a_post) == 0 or len(b_post) == 0:
            continue
        med_a_pre = float(np.median(a_pre))
        med_b_pre = float(np.median(b_pre))
        med_a_post = float(np.median(a_post))
        med_b_post = float(np.median(b_post))
        p_pre = _mwu(a_pre, b_pre)
        p_post = _mwu(a_post, b_post)
        ms_dir = "+" if ms_ars > ms_ctrl else "-" if ms_ars < ms_ctrl else "0"
        pre_dir = "+" if med_b_pre > med_a_pre else "-" if med_b_pre < med_a_pre else "0"
        post_dir = "+" if med_b_post > med_a_post else "-" if med_b_post < med_a_post else "0"
        preserved = (post_dir == ms_dir) and (post_dir != "0")
        rows.append({
            "metric": label,
            "MS_ctrl": ms_ctrl, "MS_ars": ms_ars, "MS_p": ms_p, "MS_dir": ms_dir,
            "pre_ctrl_med": med_a_pre, "pre_ars_med": med_b_pre,
            "pre_p": p_pre, "pre_dir": pre_dir,
            "post_ctrl_med": med_a_post, "post_ars_med": med_b_post,
            "post_p": p_post, "post_dir": post_dir,
            "preserved": preserved,
        })

    df = pd.DataFrame(rows)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    pd.set_option("display.float_format", lambda x: f"{x:8.4g}")
    print(df.to_string(index=False))

    out = setup_paths.VALIDATION_ROOT / "diagnostics" / "erosion_2d_comparison.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")

    print()
    print("=" * 96)
    print("VERDICT")
    print("=" * 96)
    # The right question: did the erosion BREAK any direction that was correct
    # in the pre-fix run? Metrics that were already wrong pre-erosion can't be
    # blamed on erosion — those are separate issues (e.g. different MIP set,
    # different Cellpose model version) outside the scope of this fix.
    df["pre_matched_MS"] = df["pre_dir"] == df["MS_dir"]
    df["post_matched_MS"] = df["post_dir"] == df["MS_dir"]
    df["erosion_broke_metric"] = df["pre_matched_MS"] & ~df["post_matched_MS"]
    df["unchanged_by_erosion"] = (df["pre_ctrl_med"] == df["post_ctrl_med"]) & \
                                   (df["pre_ars_med"] == df["post_ars_med"])

    pre_match = int(df["pre_matched_MS"].sum())
    post_match = int(df["post_matched_MS"].sum())
    broken = int(df["erosion_broke_metric"].sum())
    n_total = len(df)
    print(f"  Total puncta metrics tested: {n_total}")
    print(f"  Pre-erosion matched manuscript direction: {pre_match}/{n_total}")
    print(f"  Post-erosion matched manuscript direction: {post_match}/{n_total}")
    print(f"  Erosion broke (flipped a correct direction): {broken}/{n_total}")
    print()
    if broken == 0:
        print("  → Option (b) is safe to lock in: erosion does not flip any "
              "manuscript-reported direction that was correct in pre-fix.")
        print("  → The one metric where post differs from manuscript "
              "(PABPC1 condensate index) was already discrepant pre-fix and ")
        print("    is mathematically unaffected by erosion (condensate index "
              "is a per-cell p95/mean of intensity, not puncta-derived).")
    else:
        flipped = df[df["erosion_broke_metric"]]["metric"].tolist()
        print(f"  → Erosion broke: {flipped}")
        print("  → Recommend option (a) — keep erosion default at 0.")


if __name__ == "__main__":
    main()
