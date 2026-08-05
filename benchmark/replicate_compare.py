#!/usr/bin/env python
"""Replicate-level agreement between cellquant and CellProfiler (HCT116 stress granules).

Per-cell puncta count/area-fraction do not correlate cell-by-cell (the tools define a
punctum differently and per-cell counts are low). The paper makes its claims at the
REPLICATE/POPULATION level (medians per image, Table S1 convention), so this asks:
does each tool INDEPENDENTLY recover the same replicate-level pattern?

For EACH tool's OWN full per-image cell population (not the matched subset), per image:
median per-cell puncta count and median per-cell puncta area fraction -> 7 paired
points per measure (3 control + 4 arsenite). Then:
  * cross-tool tracking: Spearman rho (foregrounded — rank agreement on 7 points is the
    right test and is robust to the scale/definition difference), Pearson r + OLS slope.
  * condition-level per tool: control vs arsenite median (of replicate medians),
    direction, and Wilcoxon rank-sum on replicate medians (3v4 — a direction/verdict
    check, NOT strong evidence; near-powerless).
  * key hypothesis: does area fraction track BETTER across tools than count?

cellquant population = its analysis set (keep filter, as in Table S1); CellProfiler
population = all secondary objects (no keep filter) — i.e. each tool's standalone output.
Nothing is tuned. Fails loud on empty per-image population or non-overlapping image sets.
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, ranksums

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bench_io import COLMAP, load_cellquant, load_cellprofiler

MEASURES = [("count", "median per-cell puncta count"),
            ("area_fraction", "median per-cell puncta area fraction")]
COND_COLOR = {"control": "#4C78A8", "arsenite": "#E45756"}


def condition_of(base):
    m = re.match(r"MAX_([A-Za-z]+)_rep", str(base))
    if not m:
        raise SystemExit(f"ERROR: cannot parse condition from image base '{base}'")
    return m.group(1)


def per_image_medians(df, measure_col):
    """median of measure_col over each image's cells; returns {base: median}, {base: n}."""
    med = df.groupby("base")[measure_col].median().to_dict()
    n = df.groupby("base").size().to_dict()
    return med, n


def ols(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or np.allclose(x, x[0]):
        return float("nan"), float("nan")
    s, b = np.polyfit(x, y, 1)
    return float(s), float(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cellquant-cells", required=True)
    ap.add_argument("--cellquant-masks", required=True)
    ap.add_argument("--cp-cells", required=True)
    ap.add_argument("--cp-puncta", required=True)
    ap.add_argument("--cp-image", default=None)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    cp_image = args.cp_image or os.path.join(os.path.dirname(os.path.abspath(args.cp_cells)), "cp_Image.csv")

    print("=" * 78)
    print("Replicate-level agreement — HCT116 G3BP1 stress granules (medians per image)")
    print("=" * 78)
    cq = load_cellquant(args.cellquant_cells, args.cellquant_masks, COLMAP["cellquant"])
    cp = load_cellprofiler(args.cp_cells, args.cp_puncta, cp_image,
                           COLMAP["cp_cells"], COLMAP["cp_puncta"], COLMAP["cp_image"])
    print(f"[pop] cellquant analysis cells: {len(cq)} (keep filter) | "
          f"CellProfiler cells: {len(cp)} (all secondary objects)")

    # ---- fail loud: overlap + non-empty per-image populations ----
    cq_bases, cp_bases = set(cq["base"]), set(cp["base"])
    shared = sorted(cq_bases & cp_bases)
    if not shared:
        raise SystemExit(f"ERROR: image sets do not overlap.\n  cq={sorted(cq_bases)}\n  cp={sorted(cp_bases)}")
    if cq_bases ^ cp_bases:
        print(f"[align] WARNING: non-shared images ignored: "
              f"cq-only={sorted(cq_bases - cp_bases)} cp-only={sorted(cp_bases - cq_bases)}")
    cq_n = cq.groupby("base").size().to_dict(); cp_n = cp.groupby("base").size().to_dict()
    for b in shared:
        if cq_n.get(b, 0) == 0 or cp_n.get(b, 0) == 0:
            raise SystemExit(f"ERROR: image {b} has an empty cell population "
                             f"(cq={cq_n.get(b,0)}, cp={cp_n.get(b,0)}) — fail loud.")

    # ---- per-image (replicate) medians, both tools, both measures ----
    rows = []
    for b in shared:
        cqi, cpi = cq[cq["base"] == b], cp[cp["base"] == b]
        rows.append({
            "image": b, "condition": condition_of(b),
            "cq_median_count": float(cqi["measure_count"].median()),
            "cp_median_count": float(cpi["measure_count"].median()),
            "cq_median_areafrac": float(cqi["measure_area_fraction"].median()),
            "cp_median_areafrac": float(cpi["measure_area_fraction"].median()),
            "cq_n_cells": int(len(cqi)), "cp_n_cells": int(len(cpi)),
        })
    rep = pd.DataFrame(rows).sort_values(["condition", "image"]).reset_index(drop=True)
    rep.to_csv(os.path.join(args.out_dir, "replicate_agreement.csv"), index=False)

    # ---- cross-tool tracking (7 points) per measure ----
    tracking = {}
    for key, label in MEASURES:
        x = rep[f"cp_median_{'areafrac' if key=='area_fraction' else 'count'}"].values.astype(float)
        y = rep[f"cq_median_{'areafrac' if key=='area_fraction' else 'count'}"].values.astype(float)
        sr, sp = spearmanr(x, y)
        pr, pp = pearsonr(x, y)
        slope, intercept = ols(x, y)
        tracking[key] = {"label": label, "n": int(len(x)),
                         "spearman_rho": float(sr), "spearman_p": float(sp),
                         "pearson_r": float(pr), "pearson_p": float(pp),
                         "ols_slope": slope, "ols_intercept": intercept}

    # ---- condition-level per tool: medians of replicate medians + Wilcoxon (3v4) ----
    condition = {}
    for tool in ("cq", "cp"):
        condition[tool] = {}
        for key in ("count", "areafrac"):
            col = f"{tool}_median_{key}"
            ctrl = rep.loc[rep.condition == "control", col].values
            ars = rep.loc[rep.condition == "arsenite", col].values
            stat, p = ranksums(ars, ctrl)   # arsenite vs control
            condition[tool][key] = {
                "control_median": float(np.median(ctrl)), "arsenite_median": float(np.median(ars)),
                "direction": "arsenite>control" if np.median(ars) > np.median(ctrl)
                             else ("arsenite<control" if np.median(ars) < np.median(ctrl) else "equal"),
                "wilcoxon_stat": float(stat), "wilcoxon_p": float(p),
                "n_control": int(len(ctrl)), "n_arsenite": int(len(ars)),
            }

    # ---- key hypothesis: area fraction tracks better than count across tools? ----
    hyp = {
        "spearman_count": tracking["count"]["spearman_rho"],
        "spearman_area_fraction": tracking["area_fraction"]["spearman_rho"],
        "area_fraction_tracks_better": bool(
            tracking["area_fraction"]["spearman_rho"] > tracking["count"]["spearman_rho"]),
    }

    with open(os.path.join(args.out_dir, "replicate_summary.json"), "w") as fh:
        json.dump({"population_note": "cq = keep-filtered analysis set (Table S1); cp = all secondary objects",
                   "shared_images": shared, "tracking": tracking,
                   "condition": condition, "hypothesis": hyp}, fh, indent=2)

    _figures(rep, tracking, condition, args.out_dir)

    # ---- stdout ----
    print("\n--- replicate_agreement (median per-cell, per image) ---")
    with pd.option_context("display.width", 160, "display.max_columns", None):
        print(rep.to_string(index=False))
    print("\n--- cross-tool tracking (7 replicate medians; Spearman foregrounded) ---")
    print(f"  {'measure':34s} {'Spearman rho (p)':>22} {'Pearson r (p)':>20} {'OLS slope':>10}")
    for key, label in MEASURES:
        t = tracking[key]
        print(f"  {label:34s} {t['spearman_rho']:>8.3f} (p={t['spearman_p']:.3f})  "
              f"{t['pearson_r']:>8.3f} (p={t['pearson_p']:.3f})  {t['ols_slope']:>10.3g}")
    print("\n--- condition medians (of replicate medians) + Wilcoxon rank-sum (3v4; direction check) ---")
    for tool, nm in (("cq", "cellquant"), ("cp", "CellProfiler")):
        for key, kl in (("count", "count"), ("areafrac", "area fraction")):
            c = condition[tool][key]
            print(f"  {nm:12s} {kl:12s}: control={c['control_median']:.4g}  arsenite={c['arsenite_median']:.4g}  "
                  f"[{c['direction']}]  Wilcoxon p={c['wilcoxon_p']:.3f}")
    print("  (3v4 Wilcoxon is near-powerless; read direction/verdict, not significance.)")
    print("\n--- key hypothesis: area fraction tracks better across tools than count? ---")
    print(f"  Spearman count={hyp['spearman_count']:.3f}  vs  area fraction={hyp['spearman_area_fraction']:.3f}"
          f"  -> area fraction better: {hyp['area_fraction_tracks_better']}")
    if hyp["area_fraction_tracks_better"]:
        print("  Interpretation: the tools agree on HOW MUCH signal is condensed (area fraction)")
        print("  but partition it into different NUMBERS of discrete puncta (count) — i.e. the")
        print("  disagreement is object splitting, not detection of the condensed signal.")
    print(f"\nwrote: {args.out_dir}/replicate_agreement.csv, replicate_summary.json, "
          f"replicate_scatter_count.{{png,pdf}}, replicate_scatter_areafrac.{{png,pdf}}, "
          f"condition_medians.{{png,pdf}}, replicate_combined.{{png,pdf}}")


def _scatter(ax, rep, key, tracking):
    kk = "areafrac" if key == "area_fraction" else "count"
    x = rep[f"cp_median_{kk}"].values.astype(float)
    y = rep[f"cq_median_{kk}"].values.astype(float)
    for cond in ("control", "arsenite"):
        mask = (rep["condition"] == cond).values
        ax.scatter(x[mask], y[mask], s=70, c=COND_COLOR[cond], edgecolor="k", lw=0.4, label=cond, zorder=3)
    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    pad = 0.08 * (hi - lo + 1e-9)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1, label="identity", zorder=1)
    t = tracking[key]
    if np.isfinite(t["ols_slope"]):
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, t["ols_slope"] * xs + t["ols_intercept"], "b-", lw=1.2, label="OLS", zorder=2)
    ax.set_xlabel("CellProfiler"); ax.set_ylabel("cellquant")
    ax.set_title(f"{t['label']}\nSpearman rho={t['spearman_rho']:.3f} (p={t['spearman_p']:.3f}); "
                 f"Pearson r={t['pearson_r']:.3f}", fontsize=9)
    ax.legend(fontsize=7)


def _condbars(ax, condition, key, title):
    conds = ["control", "arsenite"]
    cqv = [condition["cq"][key]["control_median"], condition["cq"][key]["arsenite_median"]]
    cpv = [condition["cp"][key]["control_median"], condition["cp"][key]["arsenite_median"]]
    x = np.arange(2); w = 0.35
    ax.bar(x - w / 2, cqv, w, label="cellquant", color="#54A24B", edgecolor="k", lw=0.4)
    ax.bar(x + w / 2, cpv, w, label="CellProfiler", color="#B279A2", edgecolor="k", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(conds)
    ax.set_ylabel("median of replicate medians"); ax.set_title(title, fontsize=9); ax.legend(fontsize=7)


def _figures(rep, tracking, condition, out_dir):
    def save(fig, name):
        fig.savefig(os.path.join(out_dir, f"{name}.png"), dpi=150)
        fig.savefig(os.path.join(out_dir, f"{name}.pdf")); plt.close(fig)
    # individual
    for key, name in (("count", "replicate_scatter_count"), ("area_fraction", "replicate_scatter_areafrac")):
        f, a = plt.subplots(figsize=(5, 4.6)); _scatter(a, rep, key, tracking); f.tight_layout(); save(f, name)
    f, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    _condbars(ax[0], condition, "count", "puncta count (per-cell median)")
    _condbars(ax[1], condition, "areafrac", "puncta area fraction (per-cell median)")
    f.suptitle("Condition medians (of replicate medians), both tools", fontsize=11)
    f.tight_layout(rect=[0, 0, 1, 0.94]); save(f, "condition_medians")
    # combined
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    _scatter(axes[0], rep, "count", tracking)
    _scatter(axes[1], rep, "area_fraction", tracking)
    _condbars(axes[2], condition, "areafrac", "condition medians — area fraction")
    fig.suptitle("Replicate-level agreement — HCT116 G3BP1 (n=7 images: 3 control + 4 arsenite)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93]); save(fig, "replicate_combined")


if __name__ == "__main__":
    main()
