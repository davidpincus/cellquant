"""validation_3d/09_pca_3d_yeast.py

Rebuild the manuscript yeast PCA on 3D-derived per-cell metrics, plotted like
the published Figure 4 (PC1 vs PC2 colored by temperature). Mirrors published
recipe: per-cell, keep==True, fillna(0), StandardScaler z-score, sklearn PCA.

Published 39 features -> 3D mapping: 30 map 1:1; 7 are 2D->3D analogs
(area_px->volume_vox, circularity->sphericity); 2 (nucleolar_solidity,
nucleolar_eccentricity) have no 3D analog and are dropped. -> 37 features.

Outputs (validation_3d/):
  pca_3d_loadings.csv, pca_3d_variance_explained.csv, pca_3d_scores.csv
  pca_3d_yeast.png
"""
from __future__ import annotations
import importlib, os, sys
from pathlib import Path
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
setup_paths = importlib.import_module("00_setup_paths")
PUB_DIR = (HERE.parent / "AVCGTIA_manuscript" / "workups")

TEMP_ORDER = ["25C", "30C", "32C", "36C", "40C"]

# Published 39 features in published-PCA order.
PUB_FEATURES = [l.split(",")[0] for l in
                (PUB_DIR / "pca_loadings.csv").read_text().splitlines()[1:]]

# Published feature -> 3D column. Identity unless remapped. None = drop (no analog).
FEATURE_MAP = {
    "nucleolar_circularity": "nucleolar_sphericity",
    "cytosol_area_px": "cytosol_volume_vox",
    "cell_area_px": "cell_volume_vox",
    "nucleolar_area": "nucleolar_volume_vox",
    "Tif6_puncta_area_px": "Tif6_puncta_volume_vox",
    "Sis1_puncta_area_px": "Sis1_puncta_volume_vox",
    "nucleus_area_px": "nucleus_volume_vox",       # degenerate (no nucleus ch) in both
    "nucleolar_solidity": None,                     # no 3D analog -> drop
    "nucleolar_eccentricity": None,                 # no 3D analog -> drop
}

# Fragile puncta/condensate-derived features (STEP 4b robustness drop).
FRAGILE = {
    "Sis1_condensate_index_cell", "Sis1_condensate_index_cytosol",
    "Tif6_condensate_index_cell", "Tif6_condensate_index_cytosol",
    "Nsr1_condensate_index_cell", "Nsr1_condensate_index_cytosol",
    "Sis1_puncta_n", "Tif6_puncta_n",
    "Sis1_frac_intensity_in_puncta", "Tif6_frac_intensity_in_puncta",
    "Sis1_puncta_integrated_intensity", "Tif6_puncta_integrated_intensity",
    "Sis1_puncta_volume_vox", "Tif6_puncta_volume_vox",
}


def build_matrix():
    """Return (df_meta, X df with 3D columns, list of (pub_name, 3d_col))."""
    frames = []
    for sub in sorted(setup_paths.OUT_YEAST_3D.iterdir()):
        cf = sub / "cells.csv"
        if not cf.exists():
            continue
        df = pd.read_csv(cf)
        if "condition" not in df or df["condition"].isna().all():
            cond, rep = sub.name.split("_series1_rep", 1)
            df["condition"] = cond
            df["replicate"] = rep
        df["rep_id"] = sub.name
        frames.append(df)
    cells = pd.concat(frames, ignore_index=True)
    if "keep" in cells:
        cells = cells[cells["keep"] == True].copy()

    mapping = []   # (published_name, 3d_col)
    dropped = []
    for f in PUB_FEATURES:
        col = FEATURE_MAP.get(f, f)
        if col is None:
            dropped.append((f, "no 3D analog"))
            continue
        if col not in cells.columns:
            dropped.append((f, f"3D col '{col}' absent"))
            continue
        if cells[col].notna().mean() <= 0.5:   # mirror published >50% non-NaN QC
            dropped.append((f, f"{col} >50% NaN"))
            continue
        mapping.append((f, col))
    feat_cols = [c for _, c in mapping]
    Xdf = cells[feat_cols].copy()
    meta = cells[["rep_id", "image", "condition", "replicate", "cell_id"]].reset_index(drop=True)
    return meta, Xdf.reset_index(drop=True), mapping, dropped


def run_pca(Xdf, drop_cols=None):
    cols = [c for c in Xdf.columns if not (drop_cols and c in drop_cols)]
    X = Xdf[cols].fillna(0).values
    Xs = StandardScaler().fit_transform(X)
    pca = PCA()
    scores = pca.fit_transform(Xs)
    return pca, scores, cols


def temp_separation(scores, meta, pc=0):
    """Spearman-like monotonic check: mean PC score per temp, ordered."""
    means = [float(np.mean(scores[meta["condition"].values == t, pc]))
             for t in TEMP_ORDER]
    return means


def main():
    meta, Xdf, mapping, dropped = build_matrix()
    print(f"3D matrix: {Xdf.shape[0]} cells x {Xdf.shape[1]} features "
          f"(from {len(PUB_FEATURES)} published; dropped {len(dropped)})")
    print("MAPPING (published -> 3D):")
    for p, c in mapping:
        tag = "" if p == c else "   <-- analog"
        print(f"  {p:38s} -> {c}{tag}")
    print("DROPPED:")
    for p, why in dropped:
        print(f"  {p:38s} ({why})")

    # ---- full PCA (37 features) ----
    pca, scores, cols = run_pca(Xdf)
    pub_name = {c: p for p, c in mapping}
    var = pca.explained_variance_ratio_

    # save loadings + variance
    load = pd.DataFrame({"feature_3d": cols,
                         "feature_published": [pub_name[c] for c in cols]})
    for i in range(min(10, pca.n_components_)):
        load[f"PC{i+1}_loading"] = pca.components_[i]
    load["PC1_abs"] = load["PC1_loading"].abs()
    load = load.sort_values("PC1_abs", ascending=False).drop(columns="PC1_abs")
    load.to_csv(HERE / "pca_3d_loadings.csv", index=False)
    var_df = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(var))],
                           "variance_explained": var,
                           "cumulative_variance": np.cumsum(var)})
    var_df.to_csv(HERE / "pca_3d_variance_explained.csv", index=False)
    sc = meta.copy()
    sc["PC1"], sc["PC2"] = scores[:, 0], scores[:, 1]
    sc.to_csv(HERE / "pca_3d_scores.csv", index=False)

    # ---- robustness: drop fragile puncta features ----
    drop_cols = [c for c in cols if c in FRAGILE]
    pca_r, scores_r, cols_r = run_pca(Xdf, drop_cols=set(drop_cols))
    var_r = pca_r.explained_variance_ratio_

    # ---- UMAP on the same standardized 37-feature matrix (published params) ----
    import umap
    Xs_full = StandardScaler().fit_transform(Xdf[cols].fillna(0).values)
    emb = umap.UMAP(n_components=2, random_state=42, n_neighbors=30,
                    min_dist=0.3).fit_transform(Xs_full)
    sc["UMAP1"], sc["UMAP2"] = emb[:, 0], emb[:, 1]
    sc.to_csv(HERE / "pca_3d_scores.csv", index=False)

    # ---- figure: mirror published A (PCA) + C (loadings) + D (heatmap) + B (UMAP) ----
    cmap = plt.cm.coolwarm
    colors = {t: cmap(i / (len(TEMP_ORDER) - 1)) for i, t in enumerate(TEMP_ORDER)}
    fig = plt.figure(figsize=(18, 10))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.82], hspace=0.30)
    gs = outer[0].subgridspec(1, 3, width_ratios=[1.15, 1.0, 1.25], wspace=0.32)

    # Panel A: PC1 vs PC2 scatter
    axA = fig.add_subplot(gs[0, 0])
    for t in TEMP_ORDER:
        m = meta["condition"].values == t
        axA.scatter(scores[m, 0], scores[m, 1], c=[colors[t]], label=t,
                    s=4, alpha=0.4, rasterized=True, linewidths=0)
    axA.set_xlabel(f"PC1 ({var[0]*100:.1f}%)")
    axA.set_ylabel(f"PC2 ({var[1]*100:.1f}%)")
    axA.set_title("3D PCA — cells by temperature")
    axA.legend(markerscale=4, frameon=False, fontsize=10)
    # temp-mean PC1 trajectory overlay
    mp1 = temp_separation(scores, meta, 0)
    mp2 = temp_separation(scores, meta, 1)
    axA.plot(mp1, mp2, "-", color="k", lw=1.0, zorder=5)
    for t, x, y in zip(TEMP_ORDER, mp1, mp2):
        axA.scatter(x, y, c=[colors[t]], s=120, edgecolor="k", linewidth=1.2, zorder=6)

    # Panel C: top PC1 loadings bar
    axC = fig.add_subplot(gs[0, 1])
    top = load.reindex(load["PC1_loading"].abs().sort_values(ascending=False).index).head(12)
    yp = np.arange(len(top))[::-1]
    axC.barh(yp, top["PC1_loading"], color=["#b2182b" if v > 0 else "#2166ac"
                                            for v in top["PC1_loading"]])
    axC.set_yticks(yp)
    axC.set_yticklabels(top["feature_published"], fontsize=8)
    axC.set_xlabel("PC1 loading")
    axC.set_title("Top PC1 loadings (3D)")
    axC.axvline(0, color="grey", lw=0.6)

    # Panel D: loadings heatmap PC1-PC5 across top-PC1 features
    axD = fig.add_subplot(gs[0, 2])
    topf = load.reindex(load["PC1_loading"].abs().sort_values(ascending=False).index).head(15)
    H = topf[[f"PC{i+1}_loading" for i in range(5)]].values
    im = axD.imshow(H, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    axD.set_yticks(range(len(topf)))
    axD.set_yticklabels(topf["feature_published"], fontsize=7)
    axD.set_xticks(range(5))
    axD.set_xticklabels([f"PC{i+1}" for i in range(5)])
    axD.set_title("Loadings PC1-5 (3D)")
    fig.colorbar(im, ax=axD, fraction=0.046, pad=0.04)

    # Panel B: UMAP — merge + per-temperature small multiples (mirror published B)
    gb = outer[1].subgridspec(1, 6, wspace=0.12)
    cond = meta["condition"].values
    # merge panel (all temps overlaid)
    axM = fig.add_subplot(gb[0, 0])
    for t in TEMP_ORDER:
        m = cond == t
        axM.scatter(emb[m, 0], emb[m, 1], c=[colors[t]], label=t, s=3,
                    alpha=0.4, rasterized=True, linewidths=0)
    axM.set_title("UMAP (merge)", fontsize=11)
    axM.set_xlabel("UMAP1"); axM.set_ylabel("UMAP2")
    axM.legend(markerscale=4, frameon=False, fontsize=8, loc="best")
    xlim, ylim = axM.get_xlim(), axM.get_ylim()
    # per-temperature panels: focal temp colored over grey background
    for i, t in enumerate(TEMP_ORDER):
        ax = fig.add_subplot(gb[0, i + 1])
        ax.scatter(emb[:, 0], emb[:, 1], c="lightgrey", s=2, alpha=0.3,
                   rasterized=True, linewidths=0)
        m = cond == t
        ax.scatter(emb[m, 0], emb[m, 1], c=[colors[t]], s=3, alpha=0.6,
                   rasterized=True, linewidths=0)
        ax.set_title(t, fontsize=11)
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle("3D yeast PCA + UMAP (37 comparable features) — mirrors published Fig 4A,B,C,D",
                 fontsize=14)
    fig.savefig(HERE / "pca_3d_yeast.png", dpi=150, bbox_inches="tight")

    # ---- print STEP 4/5 summary ----
    print("\n=== VARIANCE EXPLAINED (3D vs published) ===")
    pub_var = pd.read_csv(PUB_DIR / "pca_variance_explained.csv")
    for i in range(5):
        print(f"  PC{i+1}: 3D={var[i]*100:5.1f}%   published={pub_var['variance_explained'][i]*100:5.1f}%")
    print(f"  cum5: 3D={np.cumsum(var)[4]*100:.1f}%  published={pub_var['cumulative_variance'][4]*100:.1f}%")

    print("\n=== PC1 TOP LOADINGS (3D) ===")
    for _, r in load.head(10).iterrows():
        print(f"  {r['feature_published']:38s} {r['PC1_loading']:+.3f}")
    print("=== PC1 TOP LOADINGS (published) ===")
    pub_load = pd.read_csv(PUB_DIR / "pca_loadings.csv")
    pub_load["a"] = pub_load["PC1_loading"].abs()
    for _, r in pub_load.sort_values("a", ascending=False).head(10).iterrows():
        print(f"  {r['feature']:38s} {r['PC1_loading']:+.3f}")

    print("\n=== TEMPERATURE SEPARATION along PC1 (mean score per temp) ===")
    print("  3D full   :", "  ".join(f"{t}={v:+.2f}" for t, v in zip(TEMP_ORDER, mp1)))
    mono = all(np.diff(mp1) > 0) or all(np.diff(mp1) < 0)
    print(f"  monotonic in temperature? {mono}  (range {max(mp1)-min(mp1):.2f})")

    print("\n=== ROBUSTNESS: drop fragile puncta/condensate features ===")
    print(f"  dropped {len(drop_cols)} fragile -> {len(cols_r)} robust features")
    print(f"  robust PCs: PC1={var_r[0]*100:.1f}%  PC2={var_r[1]*100:.1f}%")
    mpr = temp_separation(scores_r, meta, 0)
    print("  robust PC1:", "  ".join(f"{t}={v:+.2f}" for t, v in zip(TEMP_ORDER, mpr)))
    mono_r = all(np.diff(mpr) > 0) or all(np.diff(mpr) < 0)
    print(f"  monotonic in temperature (robust)? {mono_r}  (range {max(mpr)-min(mpr):.2f})")
    # robust top PC1 loadings
    rl = pd.DataFrame({"feat": [pub_name[c] for c in cols_r],
                       "PC1": pca_r.components_[0]})
    rl["a"] = rl["PC1"].abs()
    print("  robust PC1 top loadings:")
    for _, r in rl.sort_values("a", ascending=False).head(8).iterrows():
        print(f"    {r['feat']:36s} {r['PC1']:+.3f}")

    print("\nwrote pca_3d_yeast.png, pca_3d_loadings.csv, "
          "pca_3d_variance_explained.csv, pca_3d_scores.csv")


if __name__ == "__main__":
    np.random.seed(0)
    main()
