"""validation_3d/diagnostics/yeast_axial_profile.py

Characterize the axial (Z) cell distribution across the 5 yeast one-rep-per-temp
z-stacks, then derive a CPU-feasible crop window.

Two passes:
  (a) For each of the 5 z-stacks: per-Z 95th-percentile intensity per channel +
      fraction of integrated signal per channel.
  (b) From the 25C smoke 3D cell mask: per-Z unique cell count, per-cell
      z-extent (max_z - min_z + 1), and per-cell centroid z. The mask was
      written by the full-fidelity smoke run with seg_3d_method=full at full
      XY+Z resolution.

Outputs (all under validation_3d/diagnostics/):
  yeast_axial_profiles.png       — channel-wise per-Z intensity profiles (5x3)
  yeast_axial_profiles.csv       — raw profile values
  yeast_cells_per_z.png          — 25C per-Z cell count + z-extent histogram
  yeast_cells_per_z.csv          — 25C per-Z cell count + per-cell z-extent
  AXIAL_SUMMARY.txt              — proposed crop window + cell-height rationale
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
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
setup_paths = importlib.import_module("00_setup_paths")

OUT_DIR = Path(__file__).resolve().parent
SMOKE_DIR = setup_paths.VALIDATION_ROOT / "_smoke_yeast_3d" / "out"


def _load_zstack_czyx(path: Path) -> np.ndarray:
    """Return (C, Z, Y, X) float32. The yeast TIFFs are ZCYX uint16."""
    arr = np.asarray(tifffile.imread(str(path)))
    if arr.ndim != 4:
        raise ValueError(f"{path.name}: expected 4D, got {arr.shape}")
    # Detect axis order by smallest non-(Y,X) dim being C.
    if arr.shape[1] < arr.shape[0] and arr.shape[1] <= 10:
        return np.moveaxis(arr, 1, 0).astype(np.float32, copy=False)
    return arr.astype(np.float32, copy=False)


def per_z_signal_profiles(paths: list[Path]) -> pd.DataFrame:
    """For each image and each channel, per-Z p95 intensity + fraction of total."""
    channels = ["Tif6", "Nsr1", "Sis1"]
    rows: list[dict] = []
    for p in paths:
        print(f"  loading {p.name} for signal profile...")
        czyx = _load_zstack_czyx(p)  # (C, Z, Y, X)
        C, Z, _, _ = czyx.shape
        for c_idx, name in enumerate(channels):
            stack_c = czyx[c_idx]
            integrated_per_z = stack_c.sum(axis=(1, 2))
            total = float(integrated_per_z.sum()) or 1.0
            for z in range(Z):
                rows.append({
                    "image": p.stem,
                    "channel": name,
                    "z": z,
                    "p95_intensity": float(
                        np.percentile(stack_c[z].ravel(), 95)),
                    "integrated_intensity": float(integrated_per_z[z]),
                    "fraction_of_total": float(integrated_per_z[z] / total),
                })
        del czyx  # free 4 GB before next image
    return pd.DataFrame(rows)


def plot_signal_profiles(df: pd.DataFrame, out_png: Path) -> None:
    images = df["image"].unique().tolist()
    channels = df["channel"].unique().tolist()
    fig, axes = plt.subplots(len(images), len(channels),
                             figsize=(4.0 * len(channels), 2.0 * len(images)),
                             sharex=True, sharey="col")
    if len(images) == 1:
        axes = np.array([axes])
    if len(channels) == 1:
        axes = axes[:, None]
    for i, img in enumerate(images):
        for j, ch in enumerate(channels):
            ax = axes[i, j]
            sub = df[(df["image"] == img) & (df["channel"] == ch)]
            ax.plot(sub["z"], sub["fraction_of_total"], color="C0",
                    label="frac of total")
            ax2 = ax.twinx()
            ax2.plot(sub["z"], sub["p95_intensity"], color="C1", linestyle="--",
                     label="p95")
            ax.set_title(f"{img.split('_')[0]}  {ch}", fontsize=9)
            if i == len(images) - 1:
                ax.set_xlabel("Z slice")
            if j == 0:
                ax.set_ylabel("fraction", fontsize=8)
            if j == len(channels) - 1:
                ax2.set_ylabel("p95 (au)", color="C1", fontsize=8)
    fig.suptitle("Per-Z signal profile (yeast 1 rep/temp)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def per_z_cell_distribution(mask_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """From a 3D labelled cell mask, return:
       cells_per_z: (z, n_cells_present, frac_voxels_labelled)
       per_cell:    (cell_id, z_min, z_max, z_extent, z_centroid, vox)
    """
    print(f"  loading 3D cell mask {mask_path.name} ({mask_path.stat().st_size/1e6:.0f} MB)...")
    cellmask = np.asarray(tifffile.imread(str(mask_path)))  # (Z, Y, X) int
    Z = cellmask.shape[0]
    # Per-Z: how many distinct labels are present
    cpz: list[dict] = []
    for z in range(Z):
        ids = np.unique(cellmask[z])
        ids = ids[ids != 0]
        labelled_frac = float((cellmask[z] > 0).mean())
        cpz.append({"z": z, "n_cells_present": int(len(ids)),
                    "labelled_fraction": labelled_frac})
    # Per-cell z-extent + centroid via ndi
    # objects = list of slice tuples (z_slice, y_slice, x_slice) per label
    objects = ndi.find_objects(cellmask)
    rows: list[dict] = []
    for label_idx, sl in enumerate(objects, start=1):
        if sl is None:
            continue
        z_slice = sl[0]
        z_min = int(z_slice.start)
        z_max = int(z_slice.stop - 1)
        z_extent = z_max - z_min + 1
        # Voxel count for that label (sum of mask == label)
        # Restrict to the bounding box for speed
        sub = cellmask[sl]
        vox = int((sub == label_idx).sum())
        # Centroid Z (mean over voxels)
        if vox > 0:
            zs = np.where(sub == label_idx)[0] + z_min
            z_centroid = float(zs.mean())
        else:
            z_centroid = float("nan")
        rows.append({
            "cell_id": label_idx,
            "z_min": z_min,
            "z_max": z_max,
            "z_extent": z_extent,
            "z_centroid": z_centroid,
            "voxel_count": vox,
        })
    return pd.DataFrame(cpz), pd.DataFrame(rows)


def plot_cells_per_z(cpz: pd.DataFrame, per_cell: pd.DataFrame,
                     out_png: Path, voxel_z_um: float) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    ax = axes[0]
    ax.plot(cpz["z"], cpz["n_cells_present"], color="C2")
    ax.set_xlabel("Z slice")
    ax.set_ylabel("# distinct cell labels present")
    ax.set_title("Cells present per Z slice  (25C smoke, full-fidelity)")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.hist(per_cell["z_extent"], bins=range(1, int(per_cell["z_extent"].max()) + 2),
            color="C3", alpha=0.7)
    med = per_cell["z_extent"].median()
    p95 = per_cell["z_extent"].quantile(0.95)
    ax.axvline(med, color="k", linestyle="--", label=f"median={med:.0f} sl")
    ax.axvline(p95, color="k", linestyle=":", label=f"p95={p95:.0f} sl")
    ax.set_xlabel("Per-cell z-extent (slices)")
    ax.set_ylabel("# cells")
    ax.set_title(f"Cell z-extent  (Z step {voxel_z_um:.2f} µm)")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.hist(per_cell["z_centroid"], bins=40, color="C4", alpha=0.7)
    ax.set_xlabel("Per-cell centroid Z (slice)")
    ax.set_ylabel("# cells")
    ax.set_title("Cell centroid Z distribution")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def propose_window(signal_df: pd.DataFrame, cpz: pd.DataFrame,
                   per_cell: pd.DataFrame, total_z: int,
                   voxel_z_um: float) -> dict:
    """Choose a central crop window sized to contain whole middle-layer cells.

    Strategy:
      window_half = p95 cell z-extent / 2 (so a typical big cell fully fits)
      window_size >= p95 z-extent (so cell can fit anywhere within window
                                    without straddling its boundary)
      center on the densest centroid region (the median centroid z)
    """
    z_extent_p95 = float(per_cell["z_extent"].quantile(0.95))
    z_extent_med = float(per_cell["z_extent"].median())
    z_extent_max = float(per_cell["z_extent"].max())
    centroid_med = float(per_cell["z_centroid"].median())

    # Find centroid range covering most cells
    centroid_p05 = float(per_cell["z_centroid"].quantile(0.05))
    centroid_p95 = float(per_cell["z_centroid"].quantile(0.95))

    # Window must accommodate (centroid spread) + (cell half-height on each end)
    # so cells whose centroids are at the centroid_p95 don't poke out the top.
    half_cell_p95 = z_extent_p95 / 2.0
    half_cell_max = z_extent_max / 2.0
    # The window needed to contain centroids p05..p95 + 1 cell-half on each side
    proposed_lo = max(0, int(np.floor(centroid_p05 - half_cell_p95)))
    proposed_hi = min(total_z - 1, int(np.ceil(centroid_p95 + half_cell_p95)))
    proposed_size = proposed_hi - proposed_lo + 1
    # Center the window around the median centroid
    proposed_center = int(round(centroid_med))

    return {
        "total_z": total_z,
        "z_extent_median_slices": z_extent_med,
        "z_extent_p95_slices": z_extent_p95,
        "z_extent_max_slices": z_extent_max,
        "z_extent_median_um": z_extent_med * voxel_z_um,
        "z_extent_p95_um": z_extent_p95 * voxel_z_um,
        "centroid_median": centroid_med,
        "centroid_p05": centroid_p05,
        "centroid_p95": centroid_p95,
        "proposed_window_lo": proposed_lo,
        "proposed_window_hi": proposed_hi,
        "proposed_window_size_slices": proposed_size,
        "proposed_window_size_um": proposed_size * voxel_z_um,
        "proposed_window_center": proposed_center,
    }


def write_summary(prop: dict, summary_path: Path,
                  signal_df: pd.DataFrame, voxel_z_um: float) -> None:
    lines = []
    lines.append("Axial cell distribution + proposed crop window")
    lines.append("=" * 60)
    lines.append("")
    lines.append("Cell z-extent (from 25C_series1_rep1 smoke 3D mask, n=678):")
    lines.append(f"  median z-extent : {prop['z_extent_median_slices']:.1f} slices "
                 f"= {prop['z_extent_median_um']:.2f} µm")
    lines.append(f"  p95    z-extent : {prop['z_extent_p95_slices']:.1f} slices "
                 f"= {prop['z_extent_p95_um']:.2f} µm")
    lines.append(f"  max    z-extent : {prop['z_extent_max_slices']:.1f} slices "
                 f"= {prop['z_extent_max_slices'] * voxel_z_um:.2f} µm")
    lines.append("")
    lines.append("Cell centroid-Z distribution (where the cells live):")
    lines.append(f"  p05 / median / p95 centroid Z: "
                 f"{prop['centroid_p05']:.1f} / {prop['centroid_median']:.1f} / "
                 f"{prop['centroid_p95']:.1f}  (slice indices, 0..{prop['total_z']-1})")
    lines.append("")
    lines.append("Per-channel signal envelope (across all 5 images, summary):")
    for ch in signal_df["channel"].unique():
        ch_df = signal_df[signal_df["channel"] == ch]
        # Z range that contains the central 90% of integrated signal,
        # averaged across images
        cum = (ch_df.groupby("image")
                    .apply(lambda d: d.sort_values("z")["fraction_of_total"].cumsum())
                    .reset_index())
        cum.columns = ["image", "z", "cum_frac"]
        lo_each = cum[cum["cum_frac"] >= 0.05].groupby("image")["z"].min()
        hi_each = cum[cum["cum_frac"] >= 0.95].groupby("image")["z"].min()
        lines.append(f"  {ch}: 5%..95% signal between Z="
                     f"{int(lo_each.median())}..{int(hi_each.median())} (across images)")
    lines.append("")
    lines.append("Proposed central crop window:")
    lines.append(f"  Z range : [{prop['proposed_window_lo']}, "
                 f"{prop['proposed_window_hi']}]  "
                 f"({prop['proposed_window_size_slices']} slices = "
                 f"{prop['proposed_window_size_um']:.2f} µm of the "
                 f"{prop['total_z'] * voxel_z_um:.1f} µm stack)")
    lines.append(f"  centered around Z={prop['proposed_window_center']} "
                 f"(median cell centroid)")
    lines.append(f"  rationale: spans the centroid-p05..p95 of cells +"
                 f" half a p95-tall cell on each end, so a typical cell fits "
                 f"entirely inside the window without straddling. Cells whose "
                 f"centroids fall OUTSIDE this range — i.e., the top and bottom "
                 f"cell layers — are below 5% / above 95% centroid quantiles.")
    lines.append("")
    lines.append("Speedup vs full-stack baseline:")
    speedup = prop['total_z'] / prop['proposed_window_size_slices']
    lines.append(f"  Slices: {prop['total_z']} -> "
                 f"{prop['proposed_window_size_slices']} = {speedup:.2f}x fewer "
                 f"per-slice Cellpose calls (stitch path).")
    lines.append(f"  Note: XY downsampling stacks multiplicatively on top.")

    summary_path.write_text("\n".join(lines))
    for line in lines:
        print(line)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    voxel_z_um = setup_paths.YEAST_VOXEL_Z_UM
    zstacks = setup_paths.list_yeast_zstacks_one_per_temp()
    print(f"signal profiles on {len(zstacks)} images:")
    sig_df = per_z_signal_profiles(zstacks)
    sig_df.to_csv(OUT_DIR / "yeast_axial_profiles.csv", index=False)
    plot_signal_profiles(sig_df, OUT_DIR / "yeast_axial_profiles.png")
    print(f"wrote {OUT_DIR / 'yeast_axial_profiles.png'}")

    mask_path = SMOKE_DIR / "masks" / "25C_series1_rep1_cellmask.tif"
    if not mask_path.exists():
        print(f"!! smoke cell mask missing at {mask_path}; skipping cell-dist pass")
        return
    print(f"per-Z cell distribution on smoke 25C mask:")
    cpz, per_cell = per_z_cell_distribution(mask_path)
    cpz.to_csv(OUT_DIR / "yeast_cells_per_z.csv", index=False)
    per_cell.to_csv(OUT_DIR / "yeast_cells_zextent.csv", index=False)
    plot_cells_per_z(cpz, per_cell, OUT_DIR / "yeast_cells_per_z.png",
                     voxel_z_um=voxel_z_um)
    print(f"wrote {OUT_DIR / 'yeast_cells_per_z.png'}")

    total_z = int(cpz["z"].max() + 1)
    prop = propose_window(sig_df, cpz, per_cell, total_z, voxel_z_um)
    write_summary(prop, OUT_DIR / "AXIAL_SUMMARY.txt", sig_df, voxel_z_um)


if __name__ == "__main__":
    main()
