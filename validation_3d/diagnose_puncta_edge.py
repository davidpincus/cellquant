"""validation_3d/diagnose_puncta_edge.py

Diagnoses whether the LoG puncta detector lights up along cell-membrane
gradients. For each detected punctum (3D and 2D-on-MIP), compute the minimum
distance from its centroid to the cell-mask boundary in physical units (µm).
If the membrane-edge hypothesis is right, puncta will cluster within ~0.5 µm
of the boundary in 3D; the 2D path may or may not show the same artifact
depending on whether seg_downsample=3 suppresses the edge gradient.

Outputs:
  - validation_3d/diagnostics/edge_artifact_summary.csv
      per-image, per-axis (2D vs 3D), per-channel summary stats
  - validation_3d/diagnostics/edge_artifact_<image>_<channel>.png
      boundary-distance histogram, 3D vs 2D side-by-side
  - validation_3d/diagnostics/edge_artifact_cell_view_<image>.png
      single-cell visualization: boundary contour + puncta + distance heatmap

Run WITHOUT touching cellquant.py — pre-fix baseline diagnostic.
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
from scipy.ndimage import distance_transform_edt
from skimage import measure

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

# Physical band width to call "membrane-adjacent" (microns)
MEMBRANE_BAND_UM = 0.5

# Histogram bin edges in microns
HIST_EDGES_UM = np.linspace(0, 5.0, 26)


def _boundary_distance_field_3d(cell_mask: np.ndarray,
                                voxel_xy_um: float,
                                voxel_z_um: float) -> np.ndarray:
    """In-cell distance to nearest non-same-cell voxel, in microns.

    For each in-cell voxel, returns the distance to the nearest non-cell
    (background OR different-cell) voxel. Uses anisotropic EDT spacing.
    """
    spacing = (voxel_z_um, voxel_xy_um, voxel_xy_um)
    # Treat any background voxel as the boundary. For two touching cells, the
    # boundary between them is also background-adjacent in the typical case
    # (Cellpose leaves a 1-voxel gap). Where it doesn't, we slightly
    # underestimate distance — fine for diagnostic.
    bg = (cell_mask == 0)
    # distance from each voxel to nearest bg voxel; inside cell this is
    # distance to boundary.
    return distance_transform_edt(~bg, sampling=spacing)


def _boundary_distance_field_2d(cell_mask: np.ndarray,
                                voxel_xy_um: float) -> np.ndarray:
    spacing = (voxel_xy_um, voxel_xy_um)
    bg = (cell_mask == 0)
    return distance_transform_edt(~bg, sampling=spacing)


def _per_punctum_distance(puncta_mask: np.ndarray,
                          dist_field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For each labeled punctum, return centroid distance + total voxel count."""
    if puncta_mask.max() == 0:
        return np.array([]), np.array([])
    centroid_distances = []
    voxel_counts = []
    for prop in measure.regionprops(puncta_mask):
        # nearest-voxel centroid index
        centroid = tuple(int(round(c)) for c in prop.centroid)
        centroid_distances.append(float(dist_field[centroid]))
        voxel_counts.append(int(prop.area))
    return np.asarray(centroid_distances), np.asarray(voxel_counts)


def _voxel_band_fraction(cell_mask: np.ndarray,
                         puncta_mask: np.ndarray,
                         dist_field: np.ndarray) -> tuple[float, float, float]:
    """Voxel-level membrane-band statistics.

    Returns:
        cell_band_fraction      — fraction of in-cell voxels that are in band
        puncta_band_fraction    — fraction of puncta voxels in band
        enrichment_ratio        — puncta_band_fraction / cell_band_fraction
                                   (>1 means puncta over-enriched at edge)
    """
    in_cell = (cell_mask > 0)
    in_band = (dist_field <= MEMBRANE_BAND_UM) & in_cell
    in_punct = (puncta_mask > 0)
    n_cell = int(in_cell.sum())
    n_band = int(in_band.sum())
    n_punct = int(in_punct.sum())
    n_punct_in_band = int((in_punct & in_band).sum())
    if n_cell == 0:
        return 0.0, 0.0, float("nan")
    cell_band_frac = n_band / n_cell
    if n_punct == 0:
        return cell_band_frac, 0.0, float("nan")
    punct_band_frac = n_punct_in_band / n_punct
    enrichment = (punct_band_frac / cell_band_frac
                  if cell_band_frac > 0 else float("nan"))
    return cell_band_frac, punct_band_frac, enrichment


def analyze_image(stem: str,
                  voxel_xy_um: float = 0.094,
                  voxel_z_um: float = 0.25,
                  ) -> list[dict]:
    """Run the diagnostic on one image for both 3D and 2D-on-MIP."""
    rows = []
    run_3d = setup_paths.OUT_MAMM_3D / stem
    run_2d = setup_paths.OUT_MAMM_PUBLISHED2D / stem
    if not (run_3d / "masks").exists():
        print(f"[skip] {stem}: no 3D masks")
        return rows

    # ---- 3D ----
    cell_3d = tifffile.imread(run_3d / "masks" / f"{stem}_cellmask.tif").astype(np.int32)
    if cell_3d.ndim == 4:
        cell_3d = cell_3d[:, 0]
    print(f"  {stem} 3D cell mask: shape={cell_3d.shape}, {cell_3d.max()} cells")
    dist_3d = _boundary_distance_field_3d(cell_3d, voxel_xy_um, voxel_z_um)
    for ch in ("G3BP1", "PABPC1"):
        pp = run_3d / "masks" / f"{stem}_{ch}_punctamask.tif"
        if not pp.exists():
            continue
        p3 = tifffile.imread(pp).astype(np.int32)
        d_3d, v_3d = _per_punctum_distance(p3, dist_3d)
        cell_band_3d, punct_band_3d, enrich_3d = _voxel_band_fraction(
            cell_3d, p3, dist_3d)
        rows.append({
            "image": stem, "axis": "3d", "channel": ch,
            "n_puncta": int(len(d_3d)),
            "median_centroid_dist_um": float(np.median(d_3d)) if len(d_3d) else float("nan"),
            "frac_centroids_in_band_<0.5um": float((d_3d < MEMBRANE_BAND_UM).mean()) if len(d_3d) else float("nan"),
            "cell_band_voxel_fraction": float(cell_band_3d),
            "puncta_band_voxel_fraction": float(punct_band_3d),
            "band_enrichment_ratio": float(enrich_3d),
        })

    # ---- 2D-on-MIP ----
    if not (run_2d / "masks").exists():
        print(f"  [skip 2D] {stem}: no 2D masks")
        return rows
    cell_2d = tifffile.imread(run_2d / "masks" / f"MAX_{stem}_cellmask.tif").astype(np.int32)
    if cell_2d.ndim == 3:
        cell_2d = cell_2d[0]
    print(f"  {stem} 2D cell mask: shape={cell_2d.shape}, {cell_2d.max()} cells")
    dist_2d = _boundary_distance_field_2d(cell_2d, voxel_xy_um)
    for ch in ("G3BP1", "PABPC1"):
        pp = run_2d / "masks" / f"MAX_{stem}_{ch}_punctamask.tif"
        if not pp.exists():
            continue
        p2 = tifffile.imread(pp).astype(np.int32)
        d_2d, v_2d = _per_punctum_distance(p2, dist_2d)
        cell_band_2d, punct_band_2d, enrich_2d = _voxel_band_fraction(
            cell_2d, p2, dist_2d)
        rows.append({
            "image": stem, "axis": "2d", "channel": ch,
            "n_puncta": int(len(d_2d)),
            "median_centroid_dist_um": float(np.median(d_2d)) if len(d_2d) else float("nan"),
            "frac_centroids_in_band_<0.5um": float((d_2d < MEMBRANE_BAND_UM).mean()) if len(d_2d) else float("nan"),
            "cell_band_voxel_fraction": float(cell_band_2d),
            "puncta_band_voxel_fraction": float(punct_band_2d),
            "band_enrichment_ratio": float(enrich_2d),
        })

    return rows


def _hist_panel(stem: str, channel: str, out_path: Path) -> None:
    """3D vs 2D centroid-distance histogram + band-fraction annotation."""
    rows_3d = []
    rows_2d = []
    # Re-load masks (small enough)
    cell_3d = tifffile.imread(
        setup_paths.OUT_MAMM_3D / stem / "masks" / f"{stem}_cellmask.tif"
    ).astype(np.int32)
    if cell_3d.ndim == 4:
        cell_3d = cell_3d[:, 0]
    dist_3d = _boundary_distance_field_3d(cell_3d, 0.094, 0.25)
    pp3 = setup_paths.OUT_MAMM_3D / stem / "masks" / f"{stem}_{channel}_punctamask.tif"
    if pp3.exists():
        p3 = tifffile.imread(pp3).astype(np.int32)
        d_3d, _ = _per_punctum_distance(p3, dist_3d)
    else:
        d_3d = np.array([])

    cell_2d = tifffile.imread(
        setup_paths.OUT_MAMM_PUBLISHED2D / stem / "masks" / f"MAX_{stem}_cellmask.tif"
    ).astype(np.int32)
    if cell_2d.ndim == 3:
        cell_2d = cell_2d[0]
    dist_2d = _boundary_distance_field_2d(cell_2d, 0.094)
    pp2 = setup_paths.OUT_MAMM_PUBLISHED2D / stem / "masks" / f"MAX_{stem}_{channel}_punctamask.tif"
    if pp2.exists():
        p2 = tifffile.imread(pp2).astype(np.int32)
        d_2d, _ = _per_punctum_distance(p2, dist_2d)
    else:
        d_2d = np.array([])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, dists, label, color in [
        (axes[0], d_3d, f"3D (n={len(d_3d)})", "C0"),
        (axes[1], d_2d, f"2D-on-MIP (n={len(d_2d)})", "C3"),
    ]:
        ax.hist(dists, bins=HIST_EDGES_UM, color=color, alpha=0.85,
                edgecolor="black", linewidth=0.4)
        ax.axvline(MEMBRANE_BAND_UM, color="red", linestyle="--", linewidth=1.0,
                   label=f"membrane band ({MEMBRANE_BAND_UM} µm)")
        ax.set_xlabel("Punctum centroid distance to cell boundary (µm)")
        ax.set_ylabel("# puncta")
        ax.set_title(f"{stem} — {channel} — {label}")
        ax.legend(fontsize=8)
        if len(dists):
            frac = (dists < MEMBRANE_BAND_UM).mean()
            ax.text(0.97, 0.97,
                    f"median = {np.median(dists):.2f} µm\n"
                    f"frac < {MEMBRANE_BAND_UM} µm = {frac:.2%}",
                    transform=ax.transAxes, va="top", ha="right",
                    fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="gray", alpha=0.85))
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _cell_view_panel(stem: str, out_path: Path) -> None:
    """Single-cell visualization: boundary, puncta, distance heatmap (MIP + mid-Z)."""
    cell_3d = tifffile.imread(
        setup_paths.OUT_MAMM_3D / stem / "masks" / f"{stem}_cellmask.tif"
    ).astype(np.int32)
    if cell_3d.ndim == 4:
        cell_3d = cell_3d[:, 0]
    # Pick a cell with the most G3BP1 puncta — biggest signal
    cells_csv = setup_paths.OUT_MAMM_3D / stem / "cells.csv"
    if not cells_csv.exists():
        return
    df = pd.read_csv(cells_csv)
    if df.empty or "G3BP1_puncta_n" not in df.columns:
        return
    df_keep = df[df.get("keep", pd.Series([True] * len(df)))].copy()
    if df_keep.empty:
        return
    df_keep = df_keep.sort_values("G3BP1_puncta_n", ascending=False)
    target_cid = int(df_keep["cell_id"].iloc[0])

    cell_bin = (cell_3d == target_cid)
    if cell_bin.sum() == 0:
        return
    dist_3d = _boundary_distance_field_3d(cell_3d, 0.094, 0.25)
    p3 = tifffile.imread(
        setup_paths.OUT_MAMM_3D / stem / "masks" / f"{stem}_G3BP1_punctamask.tif"
    ).astype(np.int32)
    p3_in_cell = p3 * cell_bin

    # Slim to cell bounding box for clarity
    zz, yy, xx = np.where(cell_bin)
    z0, z1 = zz.min(), zz.max() + 1
    y0, y1 = yy.min(), yy.max() + 1
    x0, x1 = xx.min(), xx.max() + 1
    sub_cell = cell_bin[z0:z1, y0:y1, x0:x1]
    sub_dist = dist_3d[z0:z1, y0:y1, x0:x1]
    sub_dist[~sub_cell] = np.nan  # mask out background
    sub_punct = p3_in_cell[z0:z1, y0:y1, x0:x1]
    z_mid = sub_cell.shape[0] // 2

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # Top-left: MIP boundary distance heatmap
    mip_dist = np.nanmax(sub_dist, axis=0)
    im0 = axes[0, 0].imshow(mip_dist, cmap="viridis")
    axes[0, 0].set_title(f"MIP: max(boundary distance) within cell {target_cid}")
    axes[0, 0].set_axis_off()
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, label="µm")

    # Top-right: middle-Z slice boundary-distance heatmap
    mid_dist = sub_dist[z_mid]
    im1 = axes[0, 1].imshow(mid_dist, cmap="viridis")
    axes[0, 1].set_title(f"Z={z_mid + z0}: boundary-distance heatmap")
    axes[0, 1].set_axis_off()
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046, label="µm")

    # Bottom-left: MIP cell boundary + puncta MIP overlay
    mip_cell = sub_cell.any(axis=0)
    mip_punct = (sub_punct > 0).any(axis=0)
    axes[1, 0].imshow(mip_cell, cmap="gray", alpha=0.5)
    axes[1, 0].contour(mip_cell, levels=[0.5], colors="cyan", linewidths=1.0)
    axes[1, 0].contour(mip_punct, levels=[0.5], colors="magenta", linewidths=0.8)
    axes[1, 0].set_title("MIP: cell boundary (cyan) + puncta any-Z (magenta)")
    axes[1, 0].set_axis_off()

    # Bottom-right: middle-Z slice with same overlays
    mid_cell = sub_cell[z_mid]
    mid_punct = (sub_punct[z_mid] > 0)
    axes[1, 1].imshow(mid_cell, cmap="gray", alpha=0.5)
    axes[1, 1].contour(mid_cell, levels=[0.5], colors="cyan", linewidths=1.0)
    axes[1, 1].contour(mid_punct, levels=[0.5], colors="magenta", linewidths=0.8)
    axes[1, 1].set_title(f"Z={z_mid + z0}: cell + puncta-in-slice")
    axes[1, 1].set_axis_off()

    fig.suptitle(f"{stem}: cell {target_cid} ({int(cell_bin.sum())} voxels)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def main() -> None:
    out_dir = setup_paths.VALIDATION_ROOT / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)

    target_images = ["control_rep1", "arsenite_rep1"]
    all_rows: list[dict] = []
    for stem in target_images:
        print(f"\n=== {stem} ===")
        rows = analyze_image(stem)
        all_rows.extend(rows)
        for ch in ("G3BP1", "PABPC1"):
            _hist_panel(stem, ch, out_dir / f"edge_artifact_{stem}_{ch}.png")
        _cell_view_panel(stem, out_dir / f"edge_artifact_cell_view_{stem}.png")

    df = pd.DataFrame(all_rows)
    df.to_csv(out_dir / "edge_artifact_summary.csv", index=False)
    print()
    print("=" * 72)
    print("EDGE-ARTIFACT DIAGNOSTIC SUMMARY")
    print("=" * 72)
    print(df.to_string(index=False))
    print()

    # Verdict
    print("VERDICT:")
    if df.empty:
        print("  (no data — diagnostic could not run)")
        return
    for stem in target_images:
        for ch in ("G3BP1", "PABPC1"):
            sub = df[(df["image"] == stem) & (df["channel"] == ch)]
            if sub.empty:
                continue
            line_parts = [f"  {stem} / {ch}:"]
            for axis in ("3d", "2d"):
                r = sub[sub["axis"] == axis]
                if r.empty:
                    continue
                row = r.iloc[0]
                enrich = row["band_enrichment_ratio"]
                frac = row["frac_centroids_in_band_<0.5um"]
                verdict = ("EDGE-BIASED" if enrich >= 1.5
                           else "uniform-ish" if enrich >= 0.7
                           else "interior-biased")
                line_parts.append(
                    f"{axis}={verdict} (enrich={enrich:.2f}, "
                    f"frac<{MEMBRANE_BAND_UM}µm={frac:.0%})")
            print(" | ".join(line_parts))
    print()
    print(f"Wrote {out_dir}/edge_artifact_summary.csv + panels")


if __name__ == "__main__":
    main()
