#!/usr/bin/env python
"""cellquant vs CellProfiler per-cell agreement harness (HCT116 G3BP1 stress granules).

Runs AFTER both tools have processed the *same* images. Matches cells between the
tools by centroid (mutual nearest neighbour, per image), then reports per-cell
agreement on three measures: G3BP1 puncta count, puncta area fraction, and mean
G3BP1 intensity. The claim under test is EQUIVALENCE (no accuracy advantage), so
this succeeds by showing agreement; disagreement is a reported result, not a bug.

Nothing here is tuned to improve agreement. Column names are fixed in COLMAP below
and can be overridden on the CLI. The harness fails loud (raises) if an output is
missing/empty, if the two tools share no images, or if zero cells match.

Only numpy + scipy + pandas + matplotlib are used (no scikit-learn).

Design notes:
  * cellquant's cells.csv has no centroid column, so cellquant centroids are derived
    from its saved integer cell masks ({image}_cellmask.tif; label == cell_id) via
    --cellquant-masks. This is why the CLI has one arg beyond the minimal spec.
  * CellProfiler intensity (Intensity_MeanIntensity_G3BP1) is rescaled to 0-1;
    cellquant reports raw 16-bit counts. Pearson r and the OLS slope are
    scale-invariant and are the primary intensity-agreement metrics. The intensity
    Bland-Altman bias mostly reflects the unit difference and is labelled as such.
    Count and area fraction are on shared units, so their Bland-Altman is direct.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Shared IO (loaders + verified column mapping) lives in bench_io so the two
# diagnostic scripts reuse the exact same parsing. Behaviour here is unchanged.
from bench_io import COLMAP, load_cellquant, load_cellprofiler

# Measures compared; shared_units controls whether Bland-Altman is a direct
# quantity or dominated by a unit difference (intensity).
MEASURES = [
    ("count",         "G3BP1 puncta count",       True),
    ("area_fraction", "puncta area fraction",     True),
    ("intensity",     "mean cell G3BP1 intensity", False),
]


def mutual_nn(A, B, tol):
    """Mutual nearest-neighbour pairs within tol. A,B: (n,2) arrays. -> [(i,j,dist)]."""
    if len(A) == 0 or len(B) == 0:
        return []
    ta, tb = cKDTree(A), cKDTree(B)
    da, ia = tb.query(A)      # nearest B for each A
    _, ib = ta.query(B)       # nearest A for each B
    out = []
    for i, (j, d) in enumerate(zip(ia, da)):
        if d <= tol and ib[j] == i:
            out.append((i, int(j), float(d)))
    return out


def ols(x, y):
    """OLS y ~ x -> (slope, intercept). numpy only."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2 or np.allclose(x, x[0]):
        return np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def measure_stats(cq, cp, shared_units):
    """cq, cp: 1-D arrays (cellquant, cellprofiler) on matched pairs."""
    x, y = np.asarray(cp, float), np.asarray(cq, float)   # OLS cellquant ~ cellprofiler
    n = len(x)
    pr, pp = pearsonr(x, y) if n > 2 else (np.nan, np.nan)
    sr, sp = spearmanr(x, y) if n > 2 else (np.nan, np.nan)
    slope, intercept = ols(x, y)
    diff = y - x                                          # cellquant - cellprofiler
    bias = float(np.mean(diff)); sd = float(np.std(diff, ddof=1)) if n > 1 else np.nan
    return {
        "n": int(n),
        "pearson_r": float(pr), "pearson_p": float(pp),
        "spearman_r": float(sr), "spearman_p": float(sp),
        "ols_slope": slope, "ols_intercept": intercept,
        "ba_bias": bias, "ba_sd": sd,
        "ba_loa_low": bias - 1.96 * sd, "ba_loa_high": bias + 1.96 * sd,
        "ba_units": "direct (shared units)" if shared_units
                    else "reflects unit difference (CP 0-1 vs cellquant raw counts); use Pearson r / OLS slope",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cellquant-cells", required=True, help="cellquant cells.csv")
    ap.add_argument("--cellquant-masks", required=True,
                    help="cellquant masks/ dir ({image}_cellmask.tif) — centroids are derived here")
    ap.add_argument("--cp-cells", required=True, help="CellProfiler cp_Cells.csv")
    ap.add_argument("--cp-puncta", required=True, help="CellProfiler cp_Puncta.csv")
    ap.add_argument("--cp-image", default=None,
                    help="CellProfiler cp_Image.csv (default: cp_Image.csv beside --cp-cells)")
    ap.add_argument("--out-dir", required=True, help="output dir for summary/matches/figures")
    ap.add_argument("--match-tol-px", type=float, default=None,
                    help="centroid match tolerance in px (default: 0.25 * median cell equiv diameter)")
    # column-mapping overrides (name=value), e.g. --cq-intensity-col G3BP1_cell_mean
    ap.add_argument("--cq-count-col"); ap.add_argument("--cq-intensity-col")
    ap.add_argument("--cp-count-col"); ap.add_argument("--cp-intensity-col")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    cm_cq = dict(COLMAP["cellquant"]); cm_c = dict(COLMAP["cp_cells"])
    cm_p = dict(COLMAP["cp_puncta"]); cm_i = dict(COLMAP["cp_image"])
    if args.cq_count_col: cm_cq["count"] = args.cq_count_col
    if args.cq_intensity_col: cm_cq["intensity"] = args.cq_intensity_col
    if args.cp_count_col: cm_c["count"] = args.cp_count_col
    if args.cp_intensity_col: cm_c["intensity"] = args.cp_intensity_col
    cp_image = args.cp_image or os.path.join(os.path.dirname(os.path.abspath(args.cp_cells)), "cp_Image.csv")

    print("=" * 78)
    print("cellquant vs CellProfiler agreement — HCT116 G3BP1 stress granules")
    print("=" * 78)
    cq = load_cellquant(args.cellquant_cells, args.cellquant_masks, cm_cq)
    cp = load_cellprofiler(args.cp_cells, args.cp_puncta, cp_image, cm_c, cm_p, cm_i)

    # ---- align images ----
    cq_bases, cp_bases = set(cq["base"]), set(cp["base"])
    shared = sorted(cq_bases & cp_bases)
    if not shared:
        raise SystemExit(
            f"ERROR: the two tools share NO images.\n"
            f"       cellquant images: {sorted(cq_bases)}\n"
            f"       cellprofiler images: {sorted(cp_bases)}"
        )
    only_cq = sorted(cq_bases - cp_bases); only_cp = sorted(cp_bases - cq_bases)
    if only_cq: print(f"[align] WARNING: {len(only_cq)} image(s) only in cellquant: {only_cq}")
    if only_cp: print(f"[align] WARNING: {len(only_cp)} image(s) only in cellprofiler: {only_cp}")
    print(f"[align] {len(shared)} shared image(s): {shared}")

    # ---- match tolerance ----
    if args.match_tol_px is not None:
        tol = float(args.match_tol_px); tol_src = "user (--match-tol-px)"
    else:
        med_d = float(np.median(np.concatenate([cq["equiv_d"].values, cp["equiv_d"].values])))
        tol = 0.25 * med_d; tol_src = f"0.25 * median equiv diameter ({med_d:.1f} px)"
    print(f"[match] tolerance = {tol:.1f} px  [{tol_src}]")

    # ---- match cells within each image (mutual NN) ----
    rows = []
    per_image = []
    all_dist = []
    n_cq_tot = n_cp_tot = n_match_tot = 0
    for b in shared:
        A = cq[cq["base"] == b]; B = cp[cp["base"] == b]
        Axy = A[["cx", "cy"]].values; Bxy = B[["cx", "cy"]].values
        pairs = mutual_nn(Axy, Bxy, tol)
        n_cq_tot += len(A); n_cp_tot += len(B); n_match_tot += len(pairs)
        per_image.append({"image": b, "n_cellquant": len(A), "n_cellprofiler": len(B),
                          "n_matched": len(pairs),
                          "cq_match_rate": len(pairs) / max(len(A), 1),
                          "cp_match_rate": len(pairs) / max(len(B), 1)})
        if len(pairs) == 0:
            print(f"[match] WARNING: image {b} has ZERO matched cells "
                  f"(cq={len(A)}, cp={len(B)}).")
        Ai = A.reset_index(drop=True); Bi = B.reset_index(drop=True)
        for i, j, d in pairs:
            all_dist.append(d)
            rows.append({
                "image": b, "dist_px": d,
                "cq_cell_id": Ai.iloc[i].get(cm_cq["cell_id"]),
                "cp_object": Bi.iloc[j].get(cm_c["object"]),
                "cq_count": Ai.iloc[i]["measure_count"], "cp_count": Bi.iloc[j]["measure_count"],
                "cq_area_fraction": Ai.iloc[i]["measure_area_fraction"], "cp_area_fraction": Bi.iloc[j]["measure_area_fraction"],
                "cq_intensity": Ai.iloc[i]["measure_intensity"], "cp_intensity": Bi.iloc[j]["measure_intensity"],
            })

    matches = pd.DataFrame(rows)
    if len(matches) == 0:
        raise SystemExit("ERROR: ZERO cells matched across the whole dataset — refusing to report agreement.")
    matches.to_csv(os.path.join(args.out_dir, "matches.csv"), index=False)

    all_dist = np.array(all_dist)
    match_summary = {
        "n_cellquant": int(n_cq_tot), "n_cellprofiler": int(n_cp_tot), "n_matched": int(n_match_tot),
        "cq_match_rate": n_match_tot / max(n_cq_tot, 1), "cp_match_rate": n_match_tot / max(n_cp_tot, 1),
        "tolerance_px": tol, "tolerance_source": tol_src,
        "match_dist_median": float(np.median(all_dist)), "match_dist_mean": float(np.mean(all_dist)),
        "match_dist_max": float(np.max(all_dist)),
        "images_only_cellquant": only_cq, "images_only_cellprofiler": only_cp,
        "unmatched_cellquant": int(n_cq_tot - n_match_tot),
        "unmatched_cellprofiler": int(n_cp_tot - n_match_tot),
    }

    # ---- per-measure stats ----
    stats = {}
    for key, label, shared_units in MEASURES:
        s = measure_stats(matches[f"cq_{key}"].values, matches[f"cp_{key}"].values, shared_units)
        s["label"] = label; s["shared_units"] = shared_units
        stats[key] = s

    # sanity: Pearson(intensity) >= Pearson(area_fraction) >= Pearson(count), roughly
    pr = {k: stats[k]["pearson_r"] for k in ("intensity", "area_fraction", "count")}
    sane = (pr["intensity"] >= pr["area_fraction"] - 0.05) and (pr["area_fraction"] >= pr["count"] - 0.05)
    sanity = {"expected_order": "pearson(intensity) >= pearson(area_fraction) >= pearson(count)",
              "pearson_intensity": pr["intensity"], "pearson_area_fraction": pr["area_fraction"],
              "pearson_count": pr["count"], "holds": bool(sane)}

    # ---- write summary ----
    summ_rows = []
    for key, label, _ in MEASURES:
        s = stats[key]
        summ_rows.append({"measure": key, "label": label, **{k: s[k] for k in
            ("n", "pearson_r", "pearson_p", "spearman_r", "spearman_p",
             "ols_slope", "ols_intercept", "ba_bias", "ba_loa_low", "ba_loa_high", "ba_units")}})
    pd.DataFrame(summ_rows).to_csv(os.path.join(args.out_dir, "benchmark_summary.csv"), index=False)
    with open(os.path.join(args.out_dir, "benchmark_summary.json"), "w") as fh:
        json.dump({"match": match_summary, "per_image": per_image,
                   "measures": stats, "sanity": sanity}, fh, indent=2)

    # ---- figures ----
    _make_figures(matches, stats, match_summary, args.out_dir)

    # ---- stdout summary ----
    print("\n" + "-" * 78)
    print(f"MATCHING: cellquant={n_cq_tot} cells, cellprofiler={n_cp_tot} cells, matched={n_match_tot}")
    print(f"  match rate: cellquant {match_summary['cq_match_rate']:.2f} | "
          f"cellprofiler {match_summary['cp_match_rate']:.2f}  "
          f"(unmatched cq={match_summary['unmatched_cellquant']}, cp={match_summary['unmatched_cellprofiler']})")
    print(f"  match distance px: median {match_summary['match_dist_median']:.1f}, "
          f"mean {match_summary['match_dist_mean']:.1f}, max {match_summary['match_dist_max']:.1f} "
          f"(tol {tol:.1f})")
    print("  per image:")
    for pi in per_image:
        print(f"    {pi['image']:22s} cq={pi['n_cellquant']:3d} cp={pi['n_cellprofiler']:3d} "
              f"matched={pi['n_matched']:3d} (cq {pi['cq_match_rate']:.2f} / cp {pi['cp_match_rate']:.2f})")
    print("\nAGREEMENT (on matched pairs):")
    print(f"  {'measure':22s} {'n':>4} {'Pearson r':>12} {'Spearman rho':>13} {'OLS slope':>10} {'BA bias':>12}")
    for key, label, _ in MEASURES:
        s = stats[key]
        print(f"  {label:22s} {s['n']:>4} {s['pearson_r']:>7.3f}(p={s['pearson_p']:.1e}) "
              f"{s['spearman_r']:>13.3f} {s['ols_slope']:>10.3f} {s['ba_bias']:>12.4g}")
        if not s["shared_units"]:
            print(f"    ^ intensity BA bias {s['ba_bias']:.4g}: {s['ba_units']}")
    print(f"\nSANITY (pearson intensity>=area_frac>=count, ~): {'OK' if sane else 'VIOLATED'} "
          f"[int={pr['intensity']:.3f}, area={pr['area_fraction']:.3f}, count={pr['count']:.3f}]")
    if not sane:
        print("    ^ a violation is more likely a harness bug than biology — inspect before trusting.")
    print(f"\nwrote: {args.out_dir}/benchmark_summary.csv, .json, matches.csv, "
          f"agreement_combined.png/.pdf, agreement_scatter.png, agreement_blandaltman.png")
    print("-" * 78)


def _make_figures(matches, stats, match_summary, out_dir):
    n_match = match_summary["n_matched"]
    # combined 2x3: scatter (top) + Bland-Altman (bottom), one column per measure
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for col, (key, label, shared_units) in enumerate(MEASURES):
        x = matches[f"cp_{key}"].values.astype(float)   # cellprofiler
        y = matches[f"cq_{key}"].values.astype(float)   # cellquant
        s = stats[key]
        # scatter
        ax = axes[0, col]
        ax.scatter(x, y, s=16, alpha=0.6, edgecolor="none")
        if shared_units:
            lo = float(min(x.min(), y.min())); hi = float(max(x.max(), y.max()))
            ax.plot([lo, hi], [lo, hi], "r--", lw=1, label="identity")
        if np.isfinite(s["ols_slope"]):
            xs = np.linspace(x.min(), x.max(), 50)
            ax.plot(xs, s["ols_slope"] * xs + s["ols_intercept"], "b-", lw=1.2,
                    label=f"OLS (slope={s['ols_slope']:.2f})")
        ax.set_xlabel(f"CellProfiler\n{label}"); ax.set_ylabel(f"cellquant\n{label}")
        ttl = f"{label}\nPearson r={s['pearson_r']:.3f} (p={s['pearson_p']:.1e})\nSpearman rho={s['spearman_r']:.3f}"
        if not shared_units:
            ttl += "  [scales differ]"
        ax.set_title(ttl, fontsize=9); ax.legend(fontsize=7)
        # Bland-Altman
        axb = axes[1, col]
        mean = (x + y) / 2 if shared_units else None
        diff = y - x
        if shared_units:
            axb.scatter(mean, diff, s=16, alpha=0.6, edgecolor="none")
            axb.set_xlabel(f"mean of tools\n{label}")
        else:
            axb.scatter(range(len(diff)), diff, s=16, alpha=0.6, edgecolor="none")
            axb.set_xlabel("matched cell index")
        axb.axhline(s["ba_bias"], color="b", lw=1, label=f"bias={s['ba_bias']:.3g}")
        axb.axhline(s["ba_loa_low"], color="gray", ls="--", lw=1, label="±1.96 SD")
        axb.axhline(s["ba_loa_high"], color="gray", ls="--", lw=1)
        axb.set_ylabel("cellquant − CellProfiler")
        bt = f"Bland-Altman  ({'direct' if shared_units else 'unit-diff — see r/slope'})"
        axb.set_title(bt, fontsize=9); axb.legend(fontsize=7)
    fig.suptitle(f"cellquant vs CellProfiler per-cell agreement — HCT116 G3BP1 stress granules "
                 f"(n={n_match} matched cells)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(out_dir, "agreement_combined.png"), dpi=150)
    fig.savefig(os.path.join(out_dir, "agreement_combined.pdf"))
    # standalone scatter row and BA row (for flexibility)
    for which, r in (("agreement_scatter", 0), ("agreement_blandaltman", 1)):
        f2, ax2 = plt.subplots(1, 3, figsize=(15, 4.6))
        for col in range(3):
            src = axes[r, col]
            for line in src.get_lines():
                ax2[col].plot(line.get_xdata(), line.get_ydata(), line.get_linestyle(),
                              color=line.get_color(), lw=line.get_linewidth(), label=line.get_label())
            for coll in src.collections:
                off = coll.get_offsets()
                if len(off):
                    ax2[col].scatter(off[:, 0], off[:, 1], s=16, alpha=0.6, edgecolor="none")
            ax2[col].set_xlabel(src.get_xlabel()); ax2[col].set_ylabel(src.get_ylabel())
            ax2[col].set_title(src.get_title(), fontsize=9); ax2[col].legend(fontsize=7)
        f2.tight_layout(); f2.savefig(os.path.join(out_dir, f"{which}.png"), dpi=150)
        plt.close(f2)
    plt.close(fig)


if __name__ == "__main__":
    main()
