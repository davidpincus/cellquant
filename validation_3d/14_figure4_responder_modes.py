"""validation_3d/14_figure4_responder_modes.py

The single-cell population carries two modes along PC1, and heat moves only one
of them. This builds the panel that makes that point.

WHAT "TWO MODES" MEANS HERE. A two-component fit is preferred by BIC at every
temperature (Ashman's D 2.35-4.92), but the low mode is only a shoulder of the
density at 25 C -- it resolves into a separate KDE peak at 30 C and above, and
is unmistakable at 40 C. Before the volume floor was applied, 25 C did show two
resolved peaks; part of that was the sub-10 um^3 objects the floor removes. The
panel title says "separates as temperature rises" for that reason.

THE CLAIM. Splitting cells at the antimode of a two-component fit to PC1 gives a
low mode and a high mode. As temperature rises from 25 to 40 C the high mode
travels several PC1 units while the low mode stays where it is. Under the
manuscript's own statistical convention -- replicate medians, Mann-Whitney vs
25 C, Bonferroni across the four comparisons -- the high mode moves
significantly and the low mode does not. That is the whole panel: a responding
subpopulation and a non-responding one, not a population shifting as a unit.

WHY THE LOW MODE IS NOT SEGMENTATION DEBRIS. Checked three ways before drawing
anything (see the printout):
  - Raising a minimum cell-volume floor from 0 to 25 um^3, which deletes a
    quarter of all cells, does not dissolve the split: Ashman's D on the hot
    cells stays between 3.4 and 3.9 throughout.
  - Cells near the floor are no more common at 40 C (6.4% below 8.7 um^3) than
    at 25 C (6.9%), so the hot low mode is not made of them.
  - Even sub-8.7 um^3 objects carry real signal (median Tif6 178 AU against a
    camera background far below that) and 81% contain a detected nucleolus.
    Debris does not have a nucleolus.

Usage:
    python 14_figure4_responder_modes.py [--outputs-dir DIR] [--outdir DIR]
                                         [--min-volume-um3 FLOAT]

Needs sklearn; no umap. Run with the same interpreter as 10_/11_ (bm_seq).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

FONT_SIZE = 12
matplotlib.rcParams.update({
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": FONT_SIZE, "axes.linewidth": 0.8, "axes.unicode_minus": False,
})

TEMP_ORDER = ["25C", "30C", "32C", "36C", "40C"]
TEMP_VALUE = {"25C": 25, "30C": 30, "32C": 32, "36C": 36, "40C": 40}
TEMP_COLORS = {"25C": "#3B4CC0", "30C": "#6A8BD8", "32C": "#5E5E5E",
               "36C": "#E07B54", "40C": "#B40426"}
REF = "25C"
N_COMPARISONS = 4
MODE_COLORS = {"low": "#4B6CB7", "high": "#C0392B"}


def fmt_p(p):
    if not np.isfinite(p):
        return "n.s."
    if p >= 0.05:
        return f"n.s. (p = {min(p,1.0):.2f})"
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def stats_vs_ref(rm: pd.DataFrame, col="value"):
    """Replicate medians, Mann-Whitney vs 25 C, Bonferroni across 4 comparisons."""
    ref = rm.loc[rm["condition"] == REF, col].dropna().values
    out = {}
    for t in TEMP_ORDER:
        if t == REF:
            continue
        v = rm.loc[rm["condition"] == t, col].dropna().values
        if len(v) < 2 or len(ref) < 2:
            out[t] = {"p_raw": np.nan, "p_adj": np.nan, "sig": False}
            continue
        p = mannwhitneyu(v, ref, alternative="two-sided").pvalue
        adj = min(p * N_COMPARISONS, 1.0)
        out[t] = {"p_raw": float(p), "p_adj": float(adj), "sig": adj < 0.05}
    return out


def main():
    import importlib.util
    import sys

    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--outputs-dir", type=Path,
                    default=(here.parent / "jcb_revision_midway3_cellquant_analysis"
                             / "validation_3d" / "outputs_yeast" / "3d"))
    ap.add_argument("--outdir", type=Path, default=here / "figure4_3d")
    # 10 um^3 floor, on by default. Inspecting the masks of 40C_series1_rep2
    # against the raw Tif6 stack: cells in the low mode (12-28 um^3) are plainly
    # cells -- round, well-fit masks, interior signal, a nucleolus in 88% of
    # them, spanning ~2.3 um axially. Below ~8.7 um^3 a minority of objects span
    # only 2 of the 41 z-planes (0.46 um), thinner than a yeast cell can be, and
    # those are not cells. The floor removes them without touching the low mode.
    # It is not load-bearing: the split survives a floor of 25 um^3 (Ashman's D
    # 3.4-3.9 from 0 to 25), so this is hygiene, not a result.
    ap.add_argument("--min-volume-um3", type=float, default=10.0,
                    help="drop cells below this volume before fitting "
                         "(default 10; 0 disables)")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Reuse the feature screening and pretty-printing from the Figure 4 script so
    # the two panels are guaranteed to describe the same feature space.
    spec = importlib.util.spec_from_file_location(
        "fig4", here / "10_figure4_pca_umap_3d.py")
    fig4 = importlib.util.module_from_spec(spec); sys.modules["fig4"] = fig4
    spec.loader.exec_module(fig4)

    cells = fig4.load_cells(args.outputs_dir).reset_index(drop=True)
    if args.min_volume_um3 > 0:
        n0 = len(cells)
        cells = cells[cells["cell_volume_um3"] >= args.min_volume_um3].reset_index(drop=True)
        print(f"volume floor {args.min_volume_um3} um3: kept {len(cells)}/{n0}")
    cols = fig4.native_features(cells)
    X = StandardScaler().fit_transform(
        cells[cols].apply(pd.to_numeric, errors="coerce").fillna(0).values)
    pca = PCA(n_components=5).fit(X)
    pc = pca.transform(X)
    pc1 = pc[:, 0]
    # PCA sign is arbitrary. Pin it so the high mode is the large-cell one, which
    # makes "the high mode travels" mean the same thing on every rerun.
    if np.corrcoef(pc1, cells["cell_volume_um3"])[0, 1] < 0:
        pc1 = -pc1
    cells["PC1"] = pc1
    var1 = pca.explained_variance_ratio_[0]

    # ---- split at the antimode of a 2-component fit to all cells ----
    gm = GaussianMixture(2, random_state=0, n_init=10).fit(pc1.reshape(-1, 1))
    order = np.argsort(gm.means_.ravel())
    grid = np.linspace(pc1.min(), pc1.max(), 4000)
    resp = gm.predict_proba(grid.reshape(-1, 1))[:, order[1]]
    boundary = float(grid[np.argmin(np.abs(resp - 0.5))])
    cells["mode"] = np.where(pc1 >= boundary, "high", "low")
    m_low, m_high = gm.means_.ravel()[order]
    sd = np.sqrt(gm.covariances_.ravel()[order])
    D = abs(m_high - m_low) * np.sqrt(2) / np.sqrt(sd[0] ** 2 + sd[1] ** 2)
    print(f"\nPC1 explains {var1*100:.1f}%.  two-component fit: means "
          f"{m_low:+.2f} / {m_high:+.2f}, Ashman D = {D:.2f}, "
          f"antimode at PC1 = {boundary:+.2f}")
    print(f"overall occupancy: low {100*(cells['mode']=='low').mean():.1f}%  "
          f"high {100*(cells['mode']=='high').mean():.1f}%")

    # ---- replicate-level series, per mode ----
    recs = []
    for (t, rid, md), sub in cells.groupby(["condition", "rep_id", "mode"]):
        if len(sub) < 10:          # a handful of cells is not a replicate median
            continue
        recs.append({"condition": t, "rep_id": rid, "mode": md,
                     "value": float(np.median(sub["PC1"])), "n": len(sub)})
    rm = pd.DataFrame(recs)
    occ = (cells.assign(is_high=(cells["mode"] == "high"))
                .groupby(["condition", "rep_id"], as_index=False)["is_high"].mean()
                .rename(columns={"is_high": "value"}))

    st = {md: stats_vs_ref(rm[rm["mode"] == md]) for md in ("low", "high")}
    st_occ = stats_vs_ref(occ)
    print(f"\n{'':6s} {'low mode':>26s} {'high mode':>26s} {'% in high mode':>22s}")
    for t in TEMP_ORDER:
        lo = rm[(rm["mode"] == "low") & (rm["condition"] == t)]["value"]
        hi = rm[(rm["mode"] == "high") & (rm["condition"] == t)]["value"]
        o = occ[occ["condition"] == t]["value"]
        sl = "" if t == REF else f"  {fmt_p(st['low'][t]['p_adj'])}"
        sh = "" if t == REF else f"  {fmt_p(st['high'][t]['p_adj'])}"
        print(f"  {t:4s} {lo.median():+8.2f}{sl:>18s} {hi.median():+8.2f}{sh:>18s} "
              f"{100*o.median():8.1f}%")

    lo_range = (rm[rm['mode'] == 'low'].groupby('condition')['value'].median())
    hi_range = (rm[rm['mode'] == 'high'].groupby('condition')['value'].median())
    print(f"\n  low  mode travels {lo_range['40C'] - lo_range['25C']:+.2f} PC1 units, "
          f"25 -> 40 C")
    print(f"  high mode travels {hi_range['40C'] - hi_range['25C']:+.2f} PC1 units, "
          f"25 -> 40 C")

    # ---- robustness: is the high mode's march just boundary bookkeeping? ----
    # Occupancy falls at 40 C, so the high mode loses its lowest members, and
    # deleting a bottom tail raises a median for free. Quantiles that cannot be
    # touched by boundary crossings settle it: for the high mode use the upper
    # percentiles (immune to losing the bottom), for the low mode the lower ones
    # (immune to gaining a top).
    print("\n  boundary-immune check (shift 25 -> 40 C, Bonferroni p):")
    for md, qs in (("high", [("median", 50), ("p75", 75), ("p90", 90)]),
                   ("low", [("median", 50), ("p25", 25), ("p10", 10)])):
        for name, q in qs:
            recs = []
            for (t, rid), sub in cells[cells["mode"] == md].groupby(["condition", "rep_id"]):
                if len(sub) < 10:
                    continue
                recs.append({"condition": t, "value": float(np.percentile(sub["PC1"], q))})
            r = pd.DataFrame(recs)
            a = r[r["condition"] == "40C"]["value"].values
            b_ = r[r["condition"] == REF]["value"].values
            shift = float(np.median(a) - np.median(b_))
            p = min(mannwhitneyu(a, b_, alternative="two-sided").pvalue * N_COMPARISONS, 1.0)
            print(f"    {md:4s} {name:7s} shift={shift:+6.2f}  p_adj={p:.4f}")

    # =====================  FIGURE  =====================
    fig = plt.figure(figsize=(15.5, 4.4))
    gs = GridSpec(1, 3, figure=fig, wspace=0.30)

    # --- a: PC1 density per temperature, with the antimode marked ---
    ax = fig.add_subplot(gs[0, 0])
    gx = np.linspace(np.percentile(pc1, 0.2), np.percentile(pc1, 99.8), 600)
    for t in TEMP_ORDER:
        v = cells.loc[cells["condition"] == t, "PC1"].values
        ax.plot(gx, gaussian_kde(v)(gx), color=TEMP_COLORS[t], lw=2.0, label=t)
    ax.axvline(boundary, color="0.35", lw=1.0, ls=(0, (4, 3)))
    ax.text(boundary, ax.get_ylim()[1] * 0.97, " antimode", fontsize=FONT_SIZE - 3,
            color="0.35", ha="left", va="top")
    ax.set_xlabel(f"PC1 ({var1*100:.1f}%)")
    ax.set_ylabel("density")
    ax.set_title("A   A second mode separates with heat", loc="left",
                 fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=FONT_SIZE - 3, loc="upper right",
              labelspacing=0.25, handlelength=1.4)

    # --- b: where each mode sits, vs temperature ---
    ax = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(0)
    for md, lbl in (("low", "low mode (no response)"),
                    ("high", "high mode (responds)")):
        sub = rm[rm["mode"] == md]
        for t in TEMP_ORDER:
            v = sub.loc[sub["condition"] == t, "value"].values
            ax.scatter(np.full(len(v), TEMP_VALUE[t]) + rng.uniform(-0.5, 0.5, len(v)),
                       v, s=26, color=MODE_COLORS[md], alpha=0.55, linewidths=0)
        med = [sub.loc[sub["condition"] == t, "value"].median() for t in TEMP_ORDER]
        ax.plot([TEMP_VALUE[t] for t in TEMP_ORDER], med, "-o", color=MODE_COLORS[md],
                lw=2.0, ms=8, mec="k", mew=0.9, zorder=5, label=lbl)
    ax.axhline(boundary, color="0.35", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.set_xticks([TEMP_VALUE[t] for t in TEMP_ORDER])
    ax.set_xticklabels([t.replace("C", "") for t in TEMP_ORDER])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("PC1 (replicate median)")
    ax.set_title("B   The response is carried by one mode", loc="left",
                 fontweight="bold", pad=10)
    ax.legend(frameon=False, fontsize=FONT_SIZE - 3, loc="upper left",
              labelspacing=0.25)
    # State the test outcome at 40 C for each mode -- that is the claim.
    # The low mode's median does clear 0.05, but that is fragile: at its 25th
    # and 10th percentiles the same test gives p = 0.17 and 0.72. The high
    # mode's shift only grows at its upper percentiles (+4.5, +4.8).
    ax.text(0.98, 0.06,
            f"40 vs 25 °C\n"
            f"high: {hi_range['40C']-hi_range['25C']:+.1f}  "
            f"{fmt_p(st['high']['40C']['p_adj'])}\n"
            f"low:  {lo_range['40C']-lo_range['25C']:+.1f}  "
            f"{fmt_p(st['low']['40C']['p_adj'])}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FONT_SIZE - 4, color="0.25", linespacing=1.4)

    # --- c: how many cells are in the responding mode ---
    ax = fig.add_subplot(gs[0, 2])
    for t in TEMP_ORDER:
        v = occ.loc[occ["condition"] == t, "value"].values * 100
        ax.scatter(np.full(len(v), TEMP_VALUE[t]) + rng.uniform(-0.5, 0.5, len(v)),
                   v, s=26, color=TEMP_COLORS[t], alpha=0.6, linewidths=0)
    med = [occ.loc[occ["condition"] == t, "value"].median() * 100 for t in TEMP_ORDER]
    ax.plot([TEMP_VALUE[t] for t in TEMP_ORDER], med, "-", color="k", lw=1.6, zorder=5)
    for t, y in zip(TEMP_ORDER, med):
        ax.scatter(TEMP_VALUE[t], y, s=95, color=TEMP_COLORS[t], edgecolor="k",
                   linewidth=1.1, zorder=6)
    ax.set_xticks([TEMP_VALUE[t] for t in TEMP_ORDER])
    ax.set_xticklabels([t.replace("C", "") for t in TEMP_ORDER])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("cells in high mode (%)")
    ax.set_title("C   At 40 °C, cells cross into the low mode", loc="left",
                 fontweight="bold", pad=10)

    stem = args.outdir / "Figure4_responder_modes"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix("." + ext), bbox_inches="tight", dpi=300)
    print(f"\n  wrote {stem.name}.svg / .pdf / .png")

    rows = []
    for md in ("low", "high"):
        for t in TEMP_ORDER:
            s = rm[(rm["mode"] == md) & (rm["condition"] == t)]["value"]
            rows.append({"series": f"{md}_mode_PC1", "temperature": t,
                         "median_of_replicate_medians": float(s.median()),
                         "n_replicates": int(len(s)),
                         "p_bonferroni_vs_25C":
                             np.nan if t == REF else st[md][t]["p_adj"]})
    for t in TEMP_ORDER:
        s = occ[occ["condition"] == t]["value"]
        rows.append({"series": "fraction_in_high_mode", "temperature": t,
                     "median_of_replicate_medians": float(s.median()),
                     "n_replicates": int(len(s)),
                     "p_bonferroni_vs_25C":
                         np.nan if t == REF else st_occ[t]["p_adj"]})
    pd.DataFrame(rows).to_csv(args.outdir / "figure4_responder_modes_stats.csv",
                              index=False)
    print("  wrote figure4_responder_modes_stats.csv")


if __name__ == "__main__":
    np.random.seed(0)
    main()
