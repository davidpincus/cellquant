"""validation_3d/05_algorithmic_sanity.py

Six by-eye QC panels per the validation plan:
  (a) 3D cell segmentation overlay — 2 mammalian (control + arsenite)
  (b) 3D cell segmentation overlay — 2 yeast (25C + 40C)  [skipped without sld reader]
  (c) 3D LoG puncta detection — mammalian arsenite
  (d) 3D LoG puncta detection — yeast 36C  [skipped]
  (e) ROI colocalization check — Pearson/Manders three ways for one mammalian cell
  (f) Marching-cubes nucleolar mesh viz  [skipped — yeast]

Each panel is a multi-panel PNG with a short caption.
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
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")
sys.path.insert(0, str(setup_paths.CELLQUANT_SCRIPT.parent))

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "cellquant", setup_paths.CELLQUANT_SCRIPT)
cq = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cq)


def _robust_norm(im: np.ndarray, p_lo=1.0, p_hi=99.5) -> np.ndarray:
    lo, hi = np.percentile(im, [p_lo, p_hi])
    if hi <= lo:
        return np.zeros_like(im, dtype=np.float32)
    return np.clip((im - lo) / (hi - lo), 0, 1).astype(np.float32)


def panel_segmentation_overlay_mammalian() -> Path | None:
    """3D segmentation overlay on 1 control + 1 arsenite mammalian image."""
    out_dir = setup_paths.OUT_MAMM_COMPARISONS / "sanity"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for row, stem in enumerate(["control_rep1", "arsenite_rep1"]):
        zstack_path = setup_paths.MAMMALIAN_ZSTACK_DIR / f"{stem}.tif"
        run_dir = setup_paths.OUT_MAMM_3D / stem
        cell_mask_path = run_dir / "masks" / f"{stem}_cellmask.tif"
        if not zstack_path.exists() or not cell_mask_path.exists():
            for ax in axes[row]:
                ax.text(0.5, 0.5, f"missing\n{stem}", ha="center", va="center")
                ax.set_axis_off()
            continue
        arr = tifffile.imread(str(zstack_path))   # (Z, C, Y, X)
        cell_mask = tifffile.imread(str(cell_mask_path)).astype(np.int32)
        if cell_mask.ndim == 4:
            cell_mask = cell_mask[:, 0]  # drop singleton
        # DAPI MIP
        dapi_mip = arr[:, 0].max(axis=0).astype(np.float32)
        # cell mask MIP (any-Z presence)
        mask_mip = (cell_mask > 0).any(axis=0).astype(np.uint8)
        # middle Z slice
        z_mid = arr.shape[0] // 2
        dapi_mid = arr[z_mid, 0].astype(np.float32)
        mask_mid = (cell_mask[z_mid] > 0).astype(np.uint8)
        # XZ projection (middle Y)
        y_mid = arr.shape[2] // 2
        dapi_xz = arr[:, 0, y_mid, :].astype(np.float32)
        mask_xz = (cell_mask[:, y_mid, :] > 0).astype(np.uint8)

        a0, a1, a2, a3 = axes[row]
        a0.imshow(_robust_norm(dapi_mip), cmap="gray")
        a0.contour(mask_mip, levels=[0.5], colors="cyan", linewidths=0.6)
        a0.set_title(f"{stem}: DAPI MIP + cell mask MIP")
        a0.set_axis_off()
        a1.imshow(_robust_norm(dapi_mid), cmap="gray")
        a1.contour(mask_mid, levels=[0.5], colors="cyan", linewidths=0.6)
        a1.set_title(f"middle Z = {z_mid}")
        a1.set_axis_off()
        a2.imshow(_robust_norm(dapi_xz), cmap="gray", aspect="auto")
        a2.contour(mask_xz, levels=[0.5], colors="cyan", linewidths=0.6)
        a2.set_title("XZ orthogonal view")
        a2.set_axis_off()
        # n_cells panel
        n_cells = int(cell_mask.max())
        keep_df_path = run_dir / "cells.csv"
        if keep_df_path.exists():
            df = pd.read_csv(keep_df_path)
            n_keep = int(df["keep"].sum()) if "keep" in df.columns else len(df)
        else:
            n_keep = -1
        a3.bar(["raw cells", "kept"], [n_cells, n_keep], color=["C0", "C2"])
        a3.set_title(f"cells: {n_cells} raw, {n_keep} kept")
    fig.suptitle("Mammalian 3D segmentation overlay")
    fig.tight_layout()
    out_path = out_dir / "sanity_a_mammalian_segmentation.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


def panel_puncta_detection_mammalian() -> Path | None:
    """3D LoG puncta overlay on a single arsenite-treated image."""
    out_dir = setup_paths.OUT_MAMM_COMPARISONS / "sanity"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "arsenite_rep1"
    zstack_path = setup_paths.MAMMALIAN_ZSTACK_DIR / f"{stem}.tif"
    run_dir = setup_paths.OUT_MAMM_3D / stem
    if not zstack_path.exists() or not run_dir.exists():
        return None
    arr = tifffile.imread(str(zstack_path))  # (Z, C, Y, X)
    g3bp1 = arr[:, 1]
    g3bp1_mip = g3bp1.max(axis=0).astype(np.float32)
    z_mid = arr.shape[0] // 2
    g3bp1_mid = g3bp1[z_mid].astype(np.float32)

    puncta_path = run_dir / "masks" / f"{stem}_G3BP1_punctamask.tif"
    if not puncta_path.exists():
        return None
    puncta_mask = tifffile.imread(str(puncta_path))
    # 3D vs 2D-on-MIP puncta detection comparison
    puncta_mip = (puncta_mask > 0).any(axis=0)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(_robust_norm(g3bp1_mip), cmap="gray")
    axes[0].contour(puncta_mip, levels=[0.5], colors="magenta", linewidths=0.5)
    axes[0].set_title(f"{stem}: G3BP1 MIP + 3D puncta any-Z")
    axes[0].set_axis_off()

    axes[1].imshow(_robust_norm(g3bp1_mid), cmap="gray")
    puncta_z = (puncta_mask[z_mid] > 0)
    axes[1].contour(puncta_z, levels=[0.5], colors="magenta", linewidths=0.5)
    axes[1].set_title(f"middle Z={z_mid}: G3BP1 slice + 3D puncta in slice")
    axes[1].set_axis_off()

    # Per-cell puncta count distribution
    cells_csv = run_dir / "cells.csv"
    if cells_csv.exists():
        df = pd.read_csv(cells_csv)
        if "G3BP1_puncta_n" in df.columns:
            axes[2].hist(df["G3BP1_puncta_n"].dropna(), bins=20, color="C1")
            axes[2].set_xlabel("G3BP1 puncta per cell (3D)")
            axes[2].set_ylabel("# cells")
            axes[2].set_title(
                f"median = {df['G3BP1_puncta_n'].median():.1f}")
    fig.suptitle("Mammalian 3D puncta detection")
    fig.tight_layout()
    out_path = out_dir / "sanity_c_mammalian_puncta.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


def panel_roi_colocalization_three_ways() -> Path | None:
    """Compute Pearson/Manders three ways on a single arsenite cell."""
    out_dir = setup_paths.OUT_MAMM_COMPARISONS / "sanity"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "arsenite_rep1"
    zstack_path = setup_paths.MAMMALIAN_ZSTACK_DIR / f"{stem}.tif"
    run_dir = setup_paths.OUT_MAMM_3D / stem
    if not (run_dir / "cells.csv").exists():
        return None
    arr = tifffile.imread(str(zstack_path))  # (Z, C, Y, X)
    cell_mask = tifffile.imread(str(run_dir / "masks" / f"{stem}_cellmask.tif")).astype(np.int32)
    df = pd.read_csv(run_dir / "cells.csv")
    if df.empty:
        return None
    # Pick a single kept cell with the most puncta
    df_keep = df[df.get("keep", pd.Series([True] * len(df)))].copy()
    if df_keep.empty:
        return None
    df_keep = df_keep.sort_values("G3BP1_puncta_n", ascending=False)
    target_cid = int(df_keep["cell_id"].iloc[0])
    g3bp1 = arr[:, 1].astype(np.float64)
    pabpc1 = arr[:, 2].astype(np.float64)
    cell_3d = (cell_mask == target_cid)
    g3bp1_v = g3bp1[cell_3d]
    pabpc1_v = pabpc1[cell_3d]

    # Hand-computed Pearson
    def _pearson(a, b):
        if a.size < 2:
            return float("nan")
        am, bm = a.mean(), b.mean()
        astd, bstd = a.std(), b.std()
        if astd == 0 or bstd == 0:
            return float("nan")
        return float(((a - am) * (b - bm)).mean() / (astd * bstd))

    r_hand_3d = _pearson(g3bp1_v, pabpc1_v)

    # MIP equivalent (collapse Z, then restrict to cell footprint)
    cell_2d = cell_3d.any(axis=0)
    g3bp1_mip = arr[:, 1].max(axis=0).astype(np.float64)
    pabpc1_mip = arr[:, 2].max(axis=0).astype(np.float64)
    r_hand_mip = _pearson(g3bp1_mip[cell_2d], pabpc1_mip[cell_2d])

    # cellquant's own 3D value from the run's colocalization.csv
    coloc_csv = run_dir / "colocalization.csv"
    r_cq_3d = float("nan")
    m1_cq_3d = float("nan")
    if coloc_csv.exists():
        coloc = pd.read_csv(coloc_csv)
        mask = (coloc["cell_id"] == target_cid) & (coloc["pair"] == "G3BP1_vs_PABPC1")
        if mask.any():
            r_cq_3d = float(coloc.loc[mask, "pearson_r"].iloc[0])
            m1_cq_3d = float(coloc.loc[mask, "manders_m1"].iloc[0])

    fig, ax = plt.subplots(figsize=(8, 4))
    labels = ["Hand 3D", "Hand 2D-on-MIP", "cellquant 3D"]
    vals = [r_hand_3d, r_hand_mip, r_cq_3d]
    bars = ax.bar(labels, vals, color=["C0", "C3", "C2"])
    ax.set_ylabel("Pearson R (G3BP1 vs PABPC1)")
    ax.set_title(f"{stem} cell {target_cid}: three-way coloc check")
    ax.axhline(0, color="black", linewidth=0.5)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.3f}",
                ha="center", va="bottom" if v >= 0 else "top")
    fig.tight_layout()
    out_path = out_dir / "sanity_e_mammalian_coloc_threeway.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)

    # Also persist the raw numbers
    summary_csv = out_dir / "sanity_e_summary.csv"
    pd.DataFrame([{
        "image": stem, "cell_id": target_cid,
        "pearson_hand_3d": r_hand_3d,
        "pearson_hand_mip": r_hand_mip,
        "pearson_cellquant_3d": r_cq_3d,
        "manders_m1_cellquant_3d": m1_cq_3d,
        "abs_diff_hand_vs_cq_3d": abs(r_hand_3d - r_cq_3d),
        "mip_minus_3d_pearson": r_hand_mip - r_hand_3d,
    }]).to_csv(summary_csv, index=False)
    return out_path


def main() -> None:
    print("=== Algorithmic sanity panels ===")
    panels = {
        "a_mammalian_segmentation": panel_segmentation_overlay_mammalian(),
        "c_mammalian_puncta": panel_puncta_detection_mammalian(),
        "e_mammalian_coloc_threeway": panel_roi_colocalization_three_ways(),
    }
    for k, v in panels.items():
        print(f"  {k}: {v if v else '(skipped — input missing)'}")


if __name__ == "__main__":
    main()
