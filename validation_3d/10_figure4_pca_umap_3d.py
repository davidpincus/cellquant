"""validation_3d/10_figure4_pca_umap_3d.py

Rebuild manuscript Figure 4 (PCA + UMAP of single-cell cellquant features) on the
NATIVE 3D yeast metrics.

Supersedes 09_pca_3d_yeast.py, which mapped nucleolar_circularity ->
nucleolar_sphericity. Sphericity has since been removed from cellquant (it is not
measurable at the axial resolution of these data), so that mapping no longer
resolves and the comparable-feature count drops from 37 to 36.

Two feature sets are computed:

  NATIVE  (default, used for the figure) -- every non-degenerate, non-redundant
          numeric per-cell column the 3D pipeline emits. This is the honest
          "multi-parameter output of cellquant" set for a 3D analysis. It ADDS
          nucleolar proximity and fragmentation indices, which the published 2D
          feature list did not contain, and uses physical volumes (um^3).

  MIRROR  (comparability check) -- the published 39-feature list mapped onto 3D
          columns; 3 drop out because all three nucleolar SHAPE descriptors
          (circularity, solidity, eccentricity) are projection-only quantities
          with no 3D counterpart. 36 remain.

Redundancy rules applied to NATIVE (a PCA double-weights collinear columns):
  - keep *_um3, drop the perfectly collinear *_vox twin
  - drop nucleolar_eq_diameter_um (monotone transform of nucleolar_volume_um3)
  - drop fraction_proximal_0.4um / _0.6um (sensitivity band around the 0.5 um
    default; keeping all three triple-weights one quantity)
  - drop zero-variance and all-NaN columns (n_nuclei, nucleus_volume_*,
    *_nucleus_mean -- there is no nucleus channel in this experiment)

Usage:
    python 10_figure4_pca_umap_3d.py [--outputs-dir DIR] [--outdir DIR]

Needs sklearn + umap-learn (e.g. miniforge3/envs/bm_seq).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# Illustrator-friendly SVG/PDF: keep text as TEXT (not outlined paths) so panel
# labels and tick labels stay editable, and use TrueType rather than Type 3.
# Matches 11_figure3_yeast_3d.py and 13_figure2_mammalian_2d.py so all three
# figures set the same type in print.
FONT_SIZE = 12
# The two dense label stacks -- panel C's 12 feature names and panel D's 15
# heatmap rows -- would collide at 12 pt in any sane panel width, so they get
# their own smaller size. Everything a reader actually reads is FONT_SIZE.
DENSE_LABEL_SIZE = 9

matplotlib.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.linewidth": 0.8,
    # U+2212 MINUS SIGN often renders as a missing-glyph box in Illustrator;
    # fall back to ASCII hyphen so axis labels survive the round trip.
    "axes.unicode_minus": False,
})


def save_all(fig, stem: Path):
    """Write .svg (editable text) + .pdf + .png. Dense scatters are rasterized
    inside an otherwise-vector figure, so Illustrator stays responsive."""
    for ext in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix("." + ext), bbox_inches="tight", dpi=300)
    print(f"  wrote {stem.name}.svg / .pdf / .png")

TEMP_ORDER = ["25C", "30C", "32C", "36C", "40C"]
TEMP_VALUE = {"25C": 25, "30C": 30, "32C": 32, "36C": 36, "40C": 40}

# Explicit temperature palette, replacing a straight sample of `coolwarm`.
# A diverging map puts its lightest colour at the midpoint, so 32 C came out
# near-white (0.865 grey) and vanished -- both against the page and against the
# 0.90 grey context cloud in the UMAP small multiples. The two flanking steps
# were nearly as pale. These five hold the blue->red temperature reading but
# keep luminance roughly constant, so every condition is visible at 2-3 px.
TEMP_COLORS = {
    "25C": "#3B4CC0",   # coolwarm endpoint, already dark
    "30C": "#6A8BD8",   # darkened from #8DB0FE
    "32C": "#5E5E5E",   # was #DCDCDC -- the invisible one
    "36C": "#E07B54",   # darkened from #F6A385
    "40C": "#B40426",   # coolwarm endpoint, already dark
}

METADATA = {"image", "condition", "replicate", "cell_id", "n_nuclei", "keep", "rep_id"}

# Redundant twins: drop the left column, keep the right one.
VOX_TWINS = {
    "cell_volume_vox": "cell_volume_um3",
    "cytosol_volume_vox": "cytosol_volume_um3",
    "nucleus_volume_vox": "nucleus_volume_um3",
    "nucleolar_volume_vox": "nucleolar_volume_um3",
    "Tif6_puncta_volume_vox": "Tif6_puncta_volume_um3",
    "Sis1_puncta_volume_vox": "Sis1_puncta_volume_um3",
}
EXPLICIT_DROP = {
    "nucleolar_eq_diameter_um",          # monotone in nucleolar_volume_um3
    "Tif6_fraction_proximal_0.4um", "Tif6_fraction_proximal_0.6um",
    "Sis1_fraction_proximal_0.4um", "Sis1_fraction_proximal_0.6um",
}

# Published 39-feature name -> 3D column, for the MIRROR set. None = no analog.
MIRROR_MAP = {
    "cell_area_px": "cell_volume_vox",
    "cytosol_area_px": "cytosol_volume_vox",
    "nucleus_area_px": "nucleus_volume_vox",
    "nucleolar_area": "nucleolar_volume_vox",
    "Tif6_puncta_area_px": "Tif6_puncta_volume_vox",
    "Sis1_puncta_area_px": "Sis1_puncta_volume_vox",
    "nucleolar_circularity": None,   # was -> sphericity; sphericity removed from cellquant
    "nucleolar_solidity": None,      # projection-only
    "nucleolar_eccentricity": None,  # projection-only
}

# Pretty labels for the figure.
PRETTY = {
    "cell_volume_um3": "cell volume", "cytosol_volume_um3": "cytosol volume",
    "nucleolar_volume_um3": "nucleolar volume", "n_nucleoli": "nucleoli per cell",
    "frac_intensity_in_puncta": "frac. signal in puncta",
    "puncta_integrated_intensity": "puncta integ. intensity",
    "puncta_mean_intensity": "puncta mean int.",
    "diffuse_mean_intensity": "diffuse mean int.",
    "puncta_volume_um3": "puncta volume", "puncta_n": "puncta count",
    "condensate_index_cell": "condensate index (cell)",
    "condensate_index_cytosol": "condensate index (cyto)",
    "fragmentation_index_simple": "fragmentation (simple)",
    "fragmentation_index_persistence": "fragmentation (persist.)",
    "cell_mean": "mean int. (cell)", "cytosol_mean": "mean int. (cyto)",
    "mean_distance": "dist. to nucleolus", "fraction_proximal": "frac. proximal",
}


def prettify(col: str) -> str:
    for ch in ("Tif6", "Nsr1", "Sis1"):
        if col.startswith(ch + "_"):
            rest = col[len(ch) + 1:]
            return f"{ch} {PRETTY.get(rest, rest.replace('_', ' '))}"
    if col.startswith("pearson_r_"):
        return "Pearson R " + col[len("pearson_r_"):].replace("_vs_", "/")
    if col.startswith("manders_m1_"):
        return "Manders M1 " + col[len("manders_m1_"):].replace("_vs_", "/")
    if col.startswith("manders_m2_"):
        return "Manders M2 " + col[len("manders_m2_"):].replace("_vs_", "/")
    return PRETTY.get(col, col.replace("_", " "))


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
    print(f"loaded {len(files)} replicates, {len(cells)} cells (keep==True)")
    return cells.reset_index(drop=True)


def usable(cells: pd.DataFrame, cols: list[str]) -> list[str]:
    """Drop all-NaN / mostly-NaN / zero-variance columns, preserving order."""
    out = []
    for c in cols:
        if c not in cells.columns:
            print(f"    drop {c:42s} (absent)")
            continue
        s = pd.to_numeric(cells[c], errors="coerce")
        if s.notna().mean() <= 0.5:
            print(f"    drop {c:42s} (>50% NaN)")
            continue
        if not np.isfinite(s.std()) or s.std() == 0:
            print(f"    drop {c:42s} (zero variance)")
            continue
        out.append(c)
    return out


def native_features(cells: pd.DataFrame) -> list[str]:
    num = cells.select_dtypes("number").columns
    cand = [c for c in num
            if c not in METADATA and c not in VOX_TWINS and c not in EXPLICIT_DROP]
    print("  NATIVE feature screening:")
    return usable(cells, cand)


def mirror_features(cells: pd.DataFrame, pub_loadings: Path) -> tuple[list[str], list[str]]:
    pub = [l.split(",")[0] for l in pub_loadings.read_text().splitlines()[1:]]
    cols, dropped = [], []
    print("  MIRROR feature screening:")
    for f in pub:
        c = MIRROR_MAP.get(f, f)
        if c is None:
            dropped.append(f)
            print(f"    drop {f:42s} (no 3D analog)")
            continue
        cols.append(c)
    return usable(cells, cols), dropped


def run_pca(cells: pd.DataFrame, cols: list[str]):
    X = cells[cols].apply(pd.to_numeric, errors="coerce").fillna(0).values
    Xs = StandardScaler().fit_transform(X)
    pca = PCA()
    scores = pca.fit_transform(Xs)
    return pca, scores, Xs


def temp_corr(scores: np.ndarray, cells: pd.DataFrame, n_pc: int = 5):
    t = cells["condition"].map(TEMP_VALUE).values.astype(float)
    return [float(spearmanr(t, scores[:, i]).statistic) for i in range(n_pc)]


def centroids(scores: np.ndarray, cells: pd.DataFrame, pcx: int, pcy: int):
    xs, ys = [], []
    for t in TEMP_ORDER:
        m = cells["condition"].values == t
        xs.append(float(np.mean(scores[m, pcx])))
        ys.append(float(np.mean(scores[m, pcy])))
    return np.array(xs), np.array(ys)


def summarize(tag, pca, scores, cells, cols):
    var = pca.explained_variance_ratio_
    rho = temp_corr(scores, cells)
    print(f"\n=== {tag}: {len(cols)} features, {scores.shape[0]} cells ===")
    print("  variance: " + "  ".join(f"PC{i+1}={var[i]*100:.1f}%" for i in range(5))
          + f"  cum5={var[:5].sum()*100:.1f}%")
    print("  Spearman(temperature, PC): "
          + "  ".join(f"PC{i+1}={rho[i]:+.2f}" for i in range(5)))
    lead = int(np.argmax(np.abs(rho)))
    print(f"  --> temperature is carried by PC{lead+1} (rho={rho[lead]:+.2f})")
    return var, rho, lead


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    default_out = (here.parent / "jcb_revision_midway3_cellquant_analysis"
                   / "validation_3d" / "outputs_yeast" / "3d")
    ap.add_argument("--outputs-dir", type=Path, default=default_out)
    ap.add_argument("--pub-loadings", type=Path,
                    default=here.parent / "AVCGTIA_manuscript" / "workups" / "pca_loadings.csv")
    ap.add_argument("--outdir", type=Path, default=here / "figure4_3d")
    ap.add_argument("--umap-neighbors", type=int, default=30)
    ap.add_argument("--umap-min-dist", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(args.outputs_dir)
    print("cells per temperature:",
          {t: int((cells["condition"] == t).sum()) for t in TEMP_ORDER})

    nat_cols = native_features(cells)
    mir_cols, mir_dropped = mirror_features(cells, args.pub_loadings)

    pca_n, sc_n, Xs_n = run_pca(cells, nat_cols)
    pca_m, sc_m, _ = run_pca(cells, mir_cols)
    var_n, rho_n, lead_n = summarize("NATIVE", pca_n, sc_n, cells, nat_cols)
    var_m, rho_m, lead_m = summarize("MIRROR (published 39 -> 3D)", pca_m, sc_m, cells, mir_cols)
    print(f"\n  MIRROR dropped {len(mir_dropped)}: {', '.join(mir_dropped)}")

    # ---- UMAP on the NATIVE standardized matrix ----
    import umap
    emb = umap.UMAP(n_components=2, random_state=args.seed,
                    n_neighbors=args.umap_neighbors,
                    min_dist=args.umap_min_dist).fit_transform(Xs_n)

    # ---- persist ----
    load = pd.DataFrame({"feature": nat_cols, "label": [prettify(c) for c in nat_cols]})
    for i in range(min(10, pca_n.n_components_)):
        load[f"PC{i+1}_loading"] = pca_n.components_[i]
    load.to_csv(args.outdir / "figure4_3d_loadings.csv", index=False)
    pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(var_n))],
                  "variance_explained": var_n,
                  "cumulative_variance": np.cumsum(var_n)}
                 ).to_csv(args.outdir / "figure4_3d_variance.csv", index=False)
    sc = cells[["rep_id", "image", "condition", "replicate", "cell_id"]].copy()
    for i in range(5):
        sc[f"PC{i+1}"] = sc_n[:, i]
    sc["UMAP1"], sc["UMAP2"] = emb[:, 0], emb[:, 1]
    sc.to_csv(args.outdir / "figure4_3d_scores.csv", index=False)

    # =====================  MAIN FIGURE (A-E)  =====================
    # Plot the PC that actually carries temperature on the y-axis.
    pcx, pcy = 0, (lead_n if lead_n != 0 else 1)
    colors = {t: matplotlib.colors.to_rgba(TEMP_COLORS[t]) for t in TEMP_ORDER}
    cond = cells["condition"].values
    n_per_t = {t: int((cond == t).sum()) for t in TEMP_ORDER}

    # Portrait. The old 15x11 landscape had to squeeze A, C and D into one row,
    # which is what drove panel C's feature names under its own title. Here each
    # of those gets half the page width, and the six UMAP panels take two rows of
    # three instead of one row of six.
    fig = plt.figure(figsize=(10, 15))
    gs = GridSpec(5, 6, figure=fig,
                  height_ratios=[1.10, 1.10, 1.00, 1.00, 0.80],
                  hspace=0.45, wspace=0.55)

    # --- A: PCA scatter + centroid trajectory ---
    axA = fig.add_subplot(gs[0, 0:3])
    for t in TEMP_ORDER:
        m = cond == t
        axA.scatter(sc_n[m, pcx], sc_n[m, pcy], c=[colors[t]], label=t, s=4,
                    alpha=0.35, linewidths=0, rasterized=True)
    cx, cy = centroids(sc_n, cells, pcx, pcy)
    axA.plot(cx, cy, "-", color="k", lw=1.2, zorder=5)
    for t, x, y in zip(TEMP_ORDER, cx, cy):
        axA.scatter(x, y, c=[colors[t]], s=150, edgecolor="k", linewidth=1.3, zorder=6)
    axA.set_xlabel(f"PC{pcx+1} ({var_n[pcx]*100:.1f}%)")
    axA.set_ylabel(f"PC{pcy+1} ({var_n[pcy]*100:.1f}%)")
    axA.set_title(f"A   Single-cell PCA, {len(nat_cols)} 3D features", loc="left",
                  fontweight="bold", pad=10)
    # Hand-built handles: reusing the scatter artists inherits their alpha=0.35,
    # which washed the swatches out to roughly the colours this palette was
    # introduced to fix. Pinned upper-left rather than "best" -- the cloud
    # shifts between reruns and the legend was landing on data.
    axA.legend(handles=[Line2D([0], [0], marker="o", linestyle="none",
                               markerfacecolor=colors[t], markeredgecolor="none",
                               markersize=7, label=t) for t in TEMP_ORDER],
               frameon=False, fontsize=FONT_SIZE - 2, loc="upper left",
               borderaxespad=0.3, labelspacing=0.35, handletextpad=0.4)

    # --- E: the temperature axis vs temperature ---
    axE = fig.add_subplot(gs[0, 3:6])

    # --- C: top PC loadings for the temperature-carrying PC ---
    axC = fig.add_subplot(gs[1, 0:3])
    key = f"PC{pcy+1}_loading"
    top = load.reindex(load[key].abs().sort_values(ascending=False).index).head(12)
    yp = np.arange(len(top))[::-1]
    axC.barh(yp, top[key], color=["#b2182b" if v > 0 else "#2166ac" for v in top[key]])
    axC.set_yticks(yp)
    axC.set_yticklabels(top["label"], fontsize=DENSE_LABEL_SIZE)
    axC.axvline(0, color="grey", lw=0.6)
    axC.set_xlabel(f"PC{pcy+1} loading")
    # Titles anchored to the FIGURE-left of the panel rather than the axes-left,
    # so a long feature name cannot run underneath its own title.
    axC.set_title(f"C   Top PC{pcy+1} loadings (temperature axis)", loc="left",
                  fontweight="bold", pad=10)

    # --- D: loadings heatmap PC1-PC5 ---
    axD = fig.add_subplot(gs[1, 3:6])
    topf = load.reindex(load[key].abs().sort_values(ascending=False).index).head(15)
    H = topf[[f"PC{i+1}_loading" for i in range(5)]].values
    im = axD.imshow(H, cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    axD.set_yticks(range(len(topf)))
    axD.set_yticklabels(topf["label"], fontsize=DENSE_LABEL_SIZE)
    # Row names on the RIGHT. On the left they are ~1.5 in long and ran straight
    # across the gutter into panel C's bars; there is no column gap wide enough
    # to hold them that does not waste the page.
    axD.yaxis.tick_right()
    axD.set_xticks(range(5)); axD.set_xticklabels([f"PC{i+1}" for i in range(5)])
    axD.tick_params(axis="x", labelsize=FONT_SIZE - 2)
    axD.set_title("D   Feature loadings, PC1-PC5", loc="left", fontweight="bold",
                  pad=10)
    # Colorbar moved below the heatmap; on the right it would sit between the
    # map and its own row labels.
    cb = fig.colorbar(im, ax=axD, orientation="horizontal",
                      fraction=0.055, pad=0.16, aspect=32)
    cb.ax.tick_params(labelsize=FONT_SIZE - 3)
    cb.set_label("loading", fontsize=FONT_SIZE - 2)

    # --- B: UMAP merge + per-temperature small multiples ---
    # Point size scales with 1/sqrt(n) so the low-n conditions (32 C has 916
    # cells vs 2577 at 36 C) stay legible instead of vanishing into the grey.
    n_ref = max(n_per_t.values())
    # Six panels over two rows of three: (all cells), then the five conditions.
    slots = [(2, 0), (2, 2), (2, 4), (3, 0), (3, 2), (3, 4)]
    axM = fig.add_subplot(gs[slots[0][0], slots[0][1]:slots[0][1] + 2])
    for t in TEMP_ORDER:
        m = cond == t
        axM.scatter(emb[m, 0], emb[m, 1], c=[colors[t]], s=3, alpha=0.45,
                    linewidths=0, rasterized=True)
    axM.set_title("B   UMAP (all cells)", loc="left", fontweight="bold",
                  fontsize=FONT_SIZE, pad=8)
    axM.set_xlabel("UMAP1", fontsize=FONT_SIZE - 2)
    axM.set_ylabel("UMAP2", fontsize=FONT_SIZE - 2)
    axM.set_xticks([]); axM.set_yticks([])
    xlim, ylim = axM.get_xlim(), axM.get_ylim()

    def darken(c, f=0.55):
        return tuple(np.clip(np.asarray(c[:3]) * f, 0, 1))

    for i, t in enumerate(TEMP_ORDER):
        r, c0 = slots[i + 1]
        ax = fig.add_subplot(gs[r, c0:c0 + 2])
        ax.scatter(emb[:, 0], emb[:, 1], c="0.90", s=2, alpha=0.55,
                   linewidths=0, rasterized=True)
        m = cond == t
        s_t = float(np.clip(3.0 * np.sqrt(n_ref / max(n_per_t[t], 1)), 3.5, 7.0))
        ax.scatter(emb[m, 0], emb[m, 1], c=[colors[t]], s=s_t, alpha=0.85,
                   edgecolors=[darken(colors[t])], linewidths=0.28,
                   rasterized=True)
        ax.set_title(f"{t}  (n={n_per_t[t]})", fontsize=FONT_SIZE - 1, pad=6)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xticks([]); ax.set_yticks([])

    # --- E: the temperature axis vs temperature (axes created with row 0) ---
    rng = np.random.default_rng(0)
    for t in TEMP_ORDER:
        m = cond == t
        axE.scatter(np.full(m.sum(), TEMP_VALUE[t]) + rng.uniform(-0.55, 0.55, m.sum()),
                    sc_n[m, pcy], c=[colors[t]], s=3, alpha=0.28, linewidths=0,
                    rasterized=True)
    axE.plot([TEMP_VALUE[t] for t in TEMP_ORDER], cy, "-", color="k", lw=1.5, zorder=5)
    for t, y in zip(TEMP_ORDER, cy):
        axE.scatter(TEMP_VALUE[t], y, c=[colors[t]], s=110, edgecolor="k",
                    linewidth=1.2, zorder=6)
    axE.set_xticks([TEMP_VALUE[t] for t in TEMP_ORDER])
    axE.set_xticklabels([t.replace("C", "") for t in TEMP_ORDER])
    axE.set_xlabel("Temperature (°C)")
    axE.set_ylabel(f"PC{pcy+1} score")
    axE.set_title(f"E   PC{pcy+1} tracks temperature "
                  f"(Spearman $\\rho$ = {rho_n[pcy]:+.2f})", loc="left",
                  fontweight="bold", pad=10)

    # Reserve the whole bottom row for the hand-drawn model panel.
    # Delete this placeholder in Illustrator once the schematic is dropped in.
    axPH = fig.add_subplot(gs[4, 0:6])
    axPH.set_xticks([]); axPH.set_yticks([])
    for s in axPH.spines.values():
        s.set_linestyle((0, (6, 6))); s.set_color("0.75"); s.set_linewidth(0.9)
    axPH.text(0.5, 0.5, "F   model schematic\n(placeholder — delete in Illustrator)",
              ha="center", va="center", color="0.6", fontsize=FONT_SIZE)

    save_all(fig, args.outdir / "Figure4_PCA_UMAP_3D")

    # =====================  SUPPLEMENT (scree + which-PC) =====================
    figS = plt.figure(figsize=(11, 4.4))
    gsS = GridSpec(1, 2, figure=figS, wspace=0.28)

    axS1 = figS.add_subplot(gsS[0, 0])
    xs = np.arange(1, 11)
    axS1.bar(xs, var_n[:10] * 100, color="#4C72B0")
    axS1.plot(xs, np.cumsum(var_n[:10]) * 100, "o-", color="#B40426", ms=5,
              label="cumulative")
    axS1.set_xticks(xs); axS1.set_xlabel("Principal component")
    axS1.set_ylabel("Variance explained (%)")
    axS1.legend(frameon=False, fontsize=FONT_SIZE - 2)
    axS1.set_title("A   Variance explained", loc="left", fontweight="bold", pad=10)

    axS2 = figS.add_subplot(gsS[0, 1])
    axS2.bar(np.arange(1, 6), rho_n[:5],
             color=["#B40426" if i == pcy else "#bbbbbb" for i in range(5)])
    axS2.axhline(0, color="k", lw=0.6)
    axS2.set_xticks(np.arange(1, 6))
    axS2.set_xticklabels([f"PC{i+1}" for i in range(5)])
    axS2.set_ylabel("Spearman $\\rho$ vs temperature")
    axS2.set_ylim(-1, 1)
    axS2.set_title("B   Which component carries temperature", loc="left",
                   fontweight="bold", pad=10)

    save_all(figS, args.outdir / "FigureS_PCA_diagnostics_3D")

    # ---- comparability note ----
    print("\n=== NATIVE vs MIRROR agreement ===")
    print(f"  temperature-carrying PC: NATIVE=PC{lead_n+1} ({rho_n[lead_n]:+.2f})   "
          f"MIRROR=PC{lead_m+1} ({rho_m[lead_m]:+.2f})")
    print(f"  PC1 variance:            NATIVE={var_n[0]*100:.1f}%   MIRROR={var_m[0]*100:.1f}%")
    print("\n  top 8 loadings on the temperature axis (NATIVE):")
    for _, r in load.reindex(load[key].abs().sort_values(ascending=False).index).head(8).iterrows():
        print(f"    {r['label']:38s} {r[key]:+.3f}")


if __name__ == "__main__":
    np.random.seed(0)
    main()
