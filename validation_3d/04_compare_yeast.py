"""validation_3d/04_compare_yeast.py

Per-temperature aggregate recapitulation for the yeast series. For each
metric and each axis (3D, matched-2D-from-3D, published-2D), compute
replicate-level Wilcoxon rank-sum tests vs the 25C reference for the four
elevated temperatures (30, 32, 36, 40C). Compare direction and significance
against published Tables S2 and S3 and flag disagreements.

Output:
  outputs_yeast/comparisons/per_temperature_medians.csv
  outputs_yeast/comparisons/wilcoxon_vs_25C.csv
  outputs_yeast/comparisons/directional_agreement.csv
  outputs_yeast/comparisons/disagreements.txt
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

REF_TEMP = "25C"
ELEVATED_TEMPS = ["30C", "32C", "36C", "40C"]
ALL_TEMPS = [REF_TEMP] + ELEVATED_TEMPS

# Metrics that exist with the same column name in 3D and 2D outputs.
COMMON_METRICS = [
    "Sis1_condensate_index_cell",
    "Sis1_frac_intensity_in_puncta",
    "Sis1_puncta_n",
    "Tif6_condensate_index_cell",
    "Tif6_frac_intensity_in_puncta",
    "Tif6_puncta_n",
    "Nsr1_condensate_index_cell",
    "pearson_r_Tif6_vs_Sis1",
    "pearson_r_Tif6_vs_Nsr1",
    "pearson_r_Sis1_vs_Nsr1",
    "manders_m1_Tif6_vs_Sis1",
    "manders_m1_Tif6_vs_Nsr1",
    "manders_m1_Sis1_vs_Nsr1",
    # Distance units differ across modes (3D µm vs 2D px) but direction is
    # comparable.
    "Sis1_mean_distance",
    "Sis1_fraction_proximal",
    "Tif6_mean_distance",
    "Tif6_fraction_proximal",
]

# Mode-specific nucleolar morphology metrics; direction interpretation matches
# the published 2D metric in the same row.
# nucleolar_sphericity and nucleolar_surface_um2 were removed from the pipeline
# (unreliable on a voxelized mask at anisotropic resolution). The 3D axis now
# carries only the volume-based size metrics; there is no 3D shape metric to
# line up against published circularity/solidity/eccentricity, so those rows
# simply have no 3D counterpart.
MORPHOLOGY_METRICS_BY_AXIS = {
    "3d": ["nucleolar_volume_um3", "nucleolar_eq_diameter_um"],
    "matched_2d": ["nucleolar_area", "nucleolar_circularity",
                   "nucleolar_solidity", "nucleolar_eccentricity"],
    "published_2d": ["nucleolar_area", "nucleolar_circularity",
                     "nucleolar_solidity", "nucleolar_eccentricity"],
}

# Published Table S2 means (rows = metric, cols = 25C/30C/32C/36C/40C). Used
# to derive an expected direction vs 25C per elevated temperature.
PUBLISHED_S2 = {
    "Sis1_condensate_index_cell": [1.68, 1.61, 1.71, 1.88, 1.67],
    "Sis1_frac_intensity_in_puncta": [0.012, 0.023, 0.018, 0.006, 0.017],
    "Sis1_puncta_n": [0.41, 0.82, 1.19, 0.20, 0.80],
    "Tif6_condensate_index_cell": [1.85, 1.95, 1.89, 1.87, 1.70],
    "Tif6_frac_intensity_in_puncta": [0.010, 0.018, 0.034, 0.004, 0.028],
    "Tif6_puncta_n": [0.28, 0.39, 1.39, 0.08, 2.58],
    "pearson_r_Tif6_vs_Sis1": [0.87, 0.85, 0.84, 0.79, 0.71],
    "pearson_r_Tif6_vs_Nsr1": [0.84, 0.83, 0.82, 0.72, 0.54],
    "pearson_r_Sis1_vs_Nsr1": [0.77, 0.74, 0.73, 0.76, 0.59],
    "manders_m1_Tif6_vs_Sis1": [1.00, 0.79, 1.00, 0.99, 0.66],
    "manders_m1_Tif6_vs_Nsr1": [0.85, 0.85, 0.78, 0.86, 0.49],
    "manders_m1_Sis1_vs_Nsr1": [0.82, 0.81, 0.79, 0.86, 0.68],
    "Sis1_fraction_proximal": [0.41, 0.44, 0.23, 0.96, 0.71],
    "Sis1_mean_distance": [8.32, 8.55, 10.82, 1.06, 4.54],
    "Tif6_fraction_proximal": [0.56, 0.59, 0.30, 0.85, 0.35],
    "Tif6_mean_distance": [6.28, 6.91, 10.31, 1.69, 9.94],
}

# Published Table S3 Bonferroni-corrected p-values vs 25C (cols 30/32/36/40C).
# None = not reported in S3 (e.g. Nsr1_condensate_index_cell is in S3 but not
# S2). Used for significance-pattern comparison.
PUBLISHED_S3 = {
    "Sis1_condensate_index_cell": [1.00, 0.50, 0.009, 0.03],
    "Sis1_puncta_n": [1.00, 0.009, 1.00, 0.02],
    "Sis1_frac_intensity_in_puncta": [1.00, 0.02, 1.00, 0.02],
    "Tif6_condensate_index_cell": [0.02, 0.02, 1.00, 1.00],
    "Tif6_puncta_n": [1.00, 0.009, 1.00, 0.006],
    "Tif6_frac_intensity_in_puncta": [1.00, 0.02, 1.00, 0.006],
    "Nsr1_condensate_index_cell": [1.00, 0.33, 0.10, 0.005],
    "pearson_r_Tif6_vs_Sis1": [0.02, 0.02, 0.009, 0.003],
    "pearson_r_Tif6_vs_Nsr1": [1.00, 1.00, 0.009, 0.003],
    "pearson_r_Sis1_vs_Nsr1": [0.07, 0.07, 0.10, 0.08],
    "manders_m1_Tif6_vs_Sis1": [1.00, 0.60, 1.00, 1.00],
    "manders_m1_Tif6_vs_Nsr1": [1.00, 1.00, 0.53, 0.003],
    "manders_m1_Sis1_vs_Nsr1": [1.00, 1.00, 1.00, 0.003],
    "Sis1_mean_distance": [1.00, 0.02, 0.01, 0.91],
    "Sis1_fraction_proximal": [1.00, 0.18, 0.01, 1.00],
    "Tif6_mean_distance": [1.00, 0.03, 0.03, 0.009],
    "Tif6_fraction_proximal": [0.55, 0.01, 1.00, 0.005],
    "nucleolar_area": [0.05, 1.00, 0.02, 0.009],
    "nucleolar_solidity": [1.00, 1.00, 0.03, 0.12],
    "nucleolar_circularity": [0.50, 1.00, 0.03, 1.00],
    "nucleolar_eccentricity": [0.07, 0.03, 0.009, 0.03],
}

# Direction of the published morphology metrics derived from text/Table S2
# context (volume↑ at high T from compaction text? actually compaction reduces
# volume — leave neutral and compute from data).
ALPHA = 0.05  # for "significant" classification


def _sign(a: float, b: float) -> str:
    if not (np.isfinite(a) and np.isfinite(b)):
        return ""
    if b > a:
        return "+"
    if b < a:
        return "-"
    return "0"


def _load_axis_cells(axis: str) -> pd.DataFrame:
    """Concatenate cells.csv across all per-image runs of one axis."""
    if axis == "3d":
        root = setup_paths.OUT_YEAST_3D
    elif axis == "matched_2d":
        root = setup_paths.OUT_YEAST_MATCHED2D
    elif axis == "published_2d":
        root = setup_paths.OUT_YEAST_PUBLISHED2D
    else:
        raise ValueError(axis)
    frames = []
    if not root.exists():
        return pd.DataFrame()
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        cf = sub / "cells.csv"
        if not cf.exists():
            continue
        df = pd.read_csv(cf)
        if "condition" not in df.columns or df["condition"].isna().all():
            # Yeast stem "25C_series1_rep1": condition before _series1_rep
            stem = sub.name
            if "_series1_rep" in stem:
                cond, rep = stem.split("_series1_rep", 1)
                df["condition"] = cond
                df["replicate"] = rep
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _per_image_medians(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Group to (condition, image) and take median of `metric`."""
    if metric not in df.columns:
        return pd.DataFrame(columns=["condition", "image", "rep_median"])
    work = df.copy()
    if "keep" in work.columns:
        work = work[work["keep"]]
    work[metric] = pd.to_numeric(work[metric], errors="coerce")
    return (work.groupby(["condition", "image"], as_index=False)[metric]
                .median()
                .rename(columns={metric: "rep_median"}))


def _wilcoxon_vs_ref(reps_a: np.ndarray, reps_b: np.ndarray) -> float:
    a = reps_a[np.isfinite(reps_a)]
    b = reps_b[np.isfinite(reps_b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    try:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        return float(p)
    except Exception:
        return float("nan")


def _per_temp_summary_for_axis(cells: pd.DataFrame, metric: str
                               ) -> dict[str, dict]:
    """For one (axis, metric), return {temp: {median, n_reps, p_vs_25C, sign_vs_25C}}."""
    out: dict[str, dict] = {t: {"median": float("nan"), "n_reps": 0,
                                "p_vs_25C": float("nan"), "sign_vs_25C": ""}
                            for t in ALL_TEMPS}
    if cells.empty or metric not in cells.columns:
        return out
    reps = _per_image_medians(cells, metric)
    ref_vals = reps.loc[reps["condition"] == REF_TEMP, "rep_median"].values
    ref_med = float(np.nanmedian(ref_vals)) if len(ref_vals) else float("nan")
    out[REF_TEMP]["median"] = ref_med
    out[REF_TEMP]["n_reps"] = int(len(ref_vals))
    for t in ELEVATED_TEMPS:
        vals = reps.loc[reps["condition"] == t, "rep_median"].values
        med = float(np.nanmedian(vals)) if len(vals) else float("nan")
        out[t]["median"] = med
        out[t]["n_reps"] = int(len(vals))
        out[t]["p_vs_25C"] = _wilcoxon_vs_ref(ref_vals, vals)
        out[t]["sign_vs_25C"] = _sign(ref_med, med)
    return out


def _published_sign(metric: str, temp: str) -> str:
    """Sign of published Table S2 mean vs 25C at this temperature."""
    if metric not in PUBLISHED_S2:
        return ""
    ref = PUBLISHED_S2[metric][0]
    idx = ALL_TEMPS.index(temp)
    val = PUBLISHED_S2[metric][idx]
    return _sign(ref, val)


def _published_sig(metric: str, temp: str) -> bool | None:
    """True if Table S3 p < 0.05 at this temp; None if metric not reported."""
    if metric not in PUBLISHED_S3:
        return None
    idx = ELEVATED_TEMPS.index(temp)
    p = PUBLISHED_S3[metric][idx]
    return p < ALPHA


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    """Build the three CSV tables + the disagreement list."""
    axes = ("3d", "matched_2d", "published_2d")
    cells_by_axis = {a: _load_axis_cells(a) for a in axes}

    # Per-axis metrics: common + axis-specific morphology
    metrics_by_axis = {
        a: COMMON_METRICS + MORPHOLOGY_METRICS_BY_AXIS[a] for a in axes
    }

    medians_rows = []
    wilcoxon_rows = []
    agree_rows = []
    disagreements: list[str] = []

    all_metrics = sorted({m for ms in metrics_by_axis.values() for m in ms})
    for metric in all_metrics:
        per_axis_summary = {}
        for a in axes:
            if metric not in metrics_by_axis[a]:
                continue
            per_axis_summary[a] = _per_temp_summary_for_axis(
                cells_by_axis[a], metric)

        # Medians table — long form
        for a, summary in per_axis_summary.items():
            for t in ALL_TEMPS:
                medians_rows.append({
                    "metric": metric, "axis": a, "temperature": t,
                    "median": summary[t]["median"],
                    "n_reps": summary[t]["n_reps"],
                })

        # Wilcoxon table — long form, only elevated temps
        for a, summary in per_axis_summary.items():
            for t in ELEVATED_TEMPS:
                wilcoxon_rows.append({
                    "metric": metric, "axis": a, "temperature": t,
                    "p_vs_25C": summary[t]["p_vs_25C"],
                    "sign_vs_25C": summary[t]["sign_vs_25C"],
                    "significant": (
                        summary[t]["p_vs_25C"] < ALPHA
                        if np.isfinite(summary[t]["p_vs_25C"]) else None),
                })

        # Directional agreement table — one row per (metric, temp)
        for t in ELEVATED_TEMPS:
            pub_sign = _published_sign(metric, t)
            pub_sig = _published_sig(metric, t)
            row = {"metric": metric, "temperature": t,
                   "published_sign": pub_sign,
                   "published_significant": pub_sig}
            for a in axes:
                if a in per_axis_summary:
                    row[f"{a}_sign"] = per_axis_summary[a][t]["sign_vs_25C"]
                    row[f"{a}_p"] = per_axis_summary[a][t]["p_vs_25C"]
                else:
                    row[f"{a}_sign"] = ""
                    row[f"{a}_p"] = float("nan")
            # Flag direction disagreement: only flag where published direction
            # is non-empty AND axis direction is non-empty AND they differ.
            if pub_sign in ("+", "-"):
                for a in axes:
                    s = row.get(f"{a}_sign", "")
                    if s in ("+", "-") and s != pub_sign:
                        msg = (f"{metric} @ {t}: published {pub_sign}, "
                               f"{a} reports {s} (p={row[f'{a}_p']:.3g})")
                        disagreements.append(msg)
            row["any_direction_disagreement"] = any(
                row.get(f"{a}_sign", "") in ("+", "-")
                and pub_sign in ("+", "-")
                and row[f"{a}_sign"] != pub_sign for a in axes)
            agree_rows.append(row)

    return (pd.DataFrame(medians_rows),
            pd.DataFrame(wilcoxon_rows),
            pd.DataFrame(agree_rows),
            disagreements)


def main() -> None:
    out_dir = setup_paths.OUT_YEAST_COMPARISONS
    out_dir.mkdir(parents=True, exist_ok=True)

    medians, wilcoxon, agree, disagreements = build_tables()

    medians.to_csv(out_dir / "per_temperature_medians.csv", index=False)
    wilcoxon.to_csv(out_dir / "wilcoxon_vs_25C.csv", index=False)
    agree.to_csv(out_dir / "directional_agreement.csv", index=False)

    disagree_path = out_dir / "disagreements.txt"
    if disagreements:
        with open(disagree_path, "w") as fh:
            fh.write("DIRECTION DISAGREEMENTS vs published Table S2\n")
            fh.write("=" * 60 + "\n")
            for d in disagreements:
                fh.write(d + "\n")
        print(f"\n{len(disagreements)} direction disagreement(s) flagged:")
        for d in disagreements:
            print(" ", d)
    else:
        with open(disagree_path, "w") as fh:
            fh.write("No direction disagreements vs Table S2 detected.\n")
        print("\nNo direction disagreements vs Table S2 detected.")

    print(f"\nWrote {medians.shape[0]} median rows, {wilcoxon.shape[0]} "
          f"Wilcoxon rows, {agree.shape[0]} agreement rows to {out_dir}/")


# ===========================================================================
# FOCUSED COLOCALIZATION ANALYSIS
# Q1 PRESERVATION: does 3D reproduce the published Pearson decline + S3 sig?
# Q2 ARTIFACT: MIP vs 3D paired per-rep deltas (masks held fixed).
# ===========================================================================
from scipy.stats import wilcoxon as _wilcoxon_signed

COLOC_PAIRS = ["Tif6_vs_Sis1", "Tif6_vs_Nsr1", "Sis1_vs_Nsr1"]
COLOC_METRICS = ["pearson_r", "manders_m1"]


def _load_coloc(axis: str) -> pd.DataFrame:
    """Concatenate colocalization.csv across reps of one axis.
    Returns long df: image, condition, replicate, cell_id, pair, metric cols."""
    root = {"3d": setup_paths.OUT_YEAST_3D,
            "matched_2d": setup_paths.OUT_YEAST_MATCHED2D,
            "published_2d": setup_paths.OUT_YEAST_PUBLISHED2D}[axis]
    frames = []
    if not root.exists():
        return pd.DataFrame()
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        cf = sub / "colocalization.csv"
        if not cf.exists() or cf.stat().st_size == 0:
            continue
        df = pd.read_csv(cf)
        stem = sub.name
        if "_series1_rep" in stem:
            cond, rep = stem.split("_series1_rep", 1)
            df["condition"] = cond
            df["replicate"] = rep
            df["rep_id"] = stem
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _rep_medians(df: pd.DataFrame, pair: str, metric: str) -> pd.DataFrame:
    """Per-rep median of metric for one pair. Returns condition, rep_id, value, n_cells."""
    sub = df[df["pair"] == pair].copy()
    sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
    g = (sub.groupby(["condition", "rep_id"], as_index=False)
            .agg(value=(metric, "median"), n_cells=(metric, "count")))
    return g


def focused_coloc() -> None:
    out_dir = setup_paths.OUT_YEAST_COMPARISONS
    out_dir.mkdir(parents=True, exist_ok=True)
    d3 = _load_coloc("3d")
    d2 = _load_coloc("matched_2d")

    # ---- Q1: 3D preservation vs published S2/S3 ----
    q1_rows = []
    nbonf = len(ELEVATED_TEMPS)  # per-pair Bonferroni factor = 4
    for pair in COLOC_PAIRS:
        for metric in COLOC_METRICS:
            key = f"{metric}_{pair}"
            pub = PUBLISHED_S2.get(key)
            pubsig = PUBLISHED_S3.get(key)  # p-values [30,32,36,40]
            rm = _rep_medians(d3, pair, metric)
            ref = rm[rm["condition"] == REF_TEMP]["value"].values
            ref_med = float(np.nanmedian(ref)) if len(ref) else float("nan")
            for ti, t in enumerate(ALL_TEMPS):
                vals = rm[rm["condition"] == t]["value"].values
                ncells = int(rm[rm["condition"] == t]["n_cells"].sum())
                med = float(np.nanmedian(vals)) if len(vals) else float("nan")
                q25, q75 = (float(np.nanpercentile(vals, 25)),
                            float(np.nanpercentile(vals, 75))) if len(vals) else (np.nan, np.nan)
                pub_val = pub[ti] if pub else np.nan
                if t == REF_TEMP:
                    p_raw = np.nan; p_bonf = np.nan; sign = ""
                    sig3d = None
                else:
                    ei = ELEVATED_TEMPS.index(t)
                    p_raw = _wilcoxon_vs_ref(ref, vals)
                    p_bonf = min(1.0, p_raw * nbonf) if np.isfinite(p_raw) else np.nan
                    sign = _sign(ref_med, med)
                    sig3d = (p_bonf < ALPHA) if np.isfinite(p_bonf) else None
                # published direction vs 25C and S3 significance
                pub_sign = _sign(pub[0], pub[ti]) if pub else ""
                pub_sig = (pubsig[ei] < ALPHA) if (pubsig and t != REF_TEMP) else None
                match = None
                if t != REF_TEMP and sign in ("+", "-") and pub_sign in ("+", "-"):
                    dir_ok = (sign == pub_sign)
                    sig_ok = (sig3d == pub_sig) if (sig3d is not None and pub_sig is not None) else None
                    match = bool(dir_ok and (sig_ok if sig_ok is not None else True))
                q1_rows.append({
                    "pair": pair, "metric": metric, "temperature": t,
                    "n_reps": int(len(vals)), "n_cells": ncells,
                    "median_3d": round(med, 4),
                    "iqr_lo": round(q25, 4), "iqr_hi": round(q75, 4),
                    "published": pub_val,
                    "sign_3d_vs25": sign, "p_raw": p_raw, "p_bonf": p_bonf,
                    "sig_3d": sig3d,
                    "pub_sign": pub_sign, "pub_sig": pub_sig,
                    "dir_AND_sig_match": match,
                })
    q1 = pd.DataFrame(q1_rows)
    q1.to_csv(out_dir / "focused_coloc_Q1.csv", index=False)

    # ---- Q2: MIP - 3D paired per-rep delta ----
    q2_rows = []
    for pair in COLOC_PAIRS:
        for metric in COLOC_METRICS:
            r3 = _rep_medians(d3, pair, metric).set_index("rep_id")["value"]
            r2 = _rep_medians(d2, pair, metric).set_index("rep_id")["value"]
            for t in ALL_TEMPS:
                ids3 = _rep_medians(d3, pair, metric)
                reps_t = ids3[ids3["condition"] == t]["rep_id"].tolist()
                paired = [(r2.get(rid, np.nan), r3.get(rid, np.nan)) for rid in reps_t]
                deltas = [a - b for a, b in paired if np.isfinite(a) and np.isfinite(b)]
                deltas = np.array(deltas, dtype=float)
                if len(deltas) >= 1:
                    md = float(np.median(deltas))
                else:
                    md = np.nan
                if len(deltas) >= 2 and np.any(deltas != 0):
                    try:
                        _, p = _wilcoxon_signed(deltas)
                    except Exception:
                        p = np.nan
                else:
                    p = np.nan
                q2_rows.append({
                    "pair": pair, "metric": metric, "temperature": t,
                    "n_pairs": int(len(deltas)),
                    "median_delta_MIP_minus_3D": round(md, 4) if np.isfinite(md) else np.nan,
                    "sign": ("+" if md > 0 else "-" if md < 0 else "0") if np.isfinite(md) else "",
                    "p_wilcoxon_vs0": p,
                })
    q2 = pd.DataFrame(q2_rows)
    q2.to_csv(out_dir / "focused_coloc_Q2.csv", index=False)

    # ---- print ----
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    print("\n" + "=" * 70)
    print("Q1 PEARSON R — 3D replicate median (IQR) vs published, Bonf Wilcoxon vs 25C")
    print("=" * 70)
    pe = q1[q1["metric"] == "pearson_r"]
    for pair in COLOC_PAIRS:
        print(f"\n  {pair}")
        for _, r in pe[pe["pair"] == pair].iterrows():
            tag = "" if r["temperature"] == REF_TEMP else (
                f" sign={r['sign_3d_vs25']} pBonf={r['p_bonf']:.3g} "
                f"sig3d={r['sig_3d']} | pub_sign={r['pub_sign']} pub_sig={r['pub_sig']} "
                f"MATCH={r['dir_AND_sig_match']}")
            print(f"    {r['temperature']:>4} n={r['n_reps']}r/{r['n_cells']:>4}c  "
                  f"3D={r['median_3d']:.3f} [{r['iqr_lo']:.3f},{r['iqr_hi']:.3f}]  "
                  f"pub={r['published']}{tag}")
    print("\n" + "=" * 70)
    print("Q1 MANDERS M1 [MARK: half-saturated in S2, projection-fragile]")
    print("=" * 70)
    mm = q1[q1["metric"] == "manders_m1"]
    for pair in COLOC_PAIRS:
        print(f"\n  {pair}")
        for _, r in mm[mm["pair"] == pair].iterrows():
            tag = "" if r["temperature"] == REF_TEMP else (
                f" sign={r['sign_3d_vs25']} pBonf={r['p_bonf']:.3g} "
                f"sig3d={r['sig_3d']} | pub_sign={r['pub_sign']} pub_sig={r['pub_sig']} "
                f"MATCH={r['dir_AND_sig_match']}")
            print(f"    {r['temperature']:>4} n={r['n_reps']}r/{r['n_cells']:>4}c  "
                  f"3D={r['median_3d']:.3f} [{r['iqr_lo']:.3f},{r['iqr_hi']:.3f}]  "
                  f"pub={r['published']}{tag}")
    print("\n" + "=" * 70)
    print("Q2 MIP - 3D paired per-rep delta (Wilcoxon signed-rank vs 0)")
    print("=" * 70)
    for metric in COLOC_METRICS:
        print(f"\n  -- {metric} --")
        sub = q2[q2["metric"] == metric]
        for pair in COLOC_PAIRS:
            row = "    " + f"{pair:14s}"
            for t in ALL_TEMPS:
                r = sub[(sub["pair"] == pair) & (sub["temperature"] == t)].iloc[0]
                row += f" {t}:{r['median_delta_MIP_minus_3D']:+.3f}(p{r['p_wilcoxon_vs0']:.2g})"
            print(row)
    print(f"\nWrote focused_coloc_Q1.csv, focused_coloc_Q2.csv to {out_dir}/")


if __name__ == "__main__":
    if "--coloc" in sys.argv:
        focused_coloc()
    else:
        main()
