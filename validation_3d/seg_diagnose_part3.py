"""validation_3d/seg_diagnose_part3.py

Run 2D Cellpose-SAM segmentation independently on each Z-slice of
control_rep3's cell-seg composite. Count cells found per slice. Plot.

If signal is concentrated in only a few Z-slices and dim slices return
zero/few masks, the stitch-by-IoU 3D linkage may fail because the per-slice
masks are too sparse to link into coherent objects — explaining why 3D
loses cells that the MIP (which collapses all signal into one bright plane)
recovers fine.
"""
from __future__ import annotations

import importlib
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

# Load cellquant as a module so we can reuse its model init + segment_cells
_spec = importlib.util.spec_from_file_location(
    "cellquant", setup_paths.CELLQUANT_SCRIPT)
cq = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cq)


def composite_seg_3d(channels_zyx: list[np.ndarray]) -> np.ndarray:
    """Replicate cellquant's build_composite_seg_image, on a 3D stack."""
    acc = np.zeros_like(channels_zyx[0], dtype=np.float32)
    for ch in channels_zyx:
        lo, hi = np.percentile(ch, [1.0, 99.8])
        if hi > lo:
            r = np.clip((ch - lo) / (hi - lo), 0, 1).astype(np.float32)
        else:
            r = np.zeros_like(ch, dtype=np.float32)
        acc += r
    lo, hi = acc.min(), acc.max()
    if hi > lo:
        acc = (acc - lo) / (hi - lo)
    return acc.astype(np.float32)


def main() -> None:
    stem = "control_rep3"
    zpath = setup_paths.MAMMALIAN_ZSTACK_DIR / f"{stem}.tif"
    arr = tifffile.imread(zpath)  # (Z, C, Y, X)
    print(f"loaded {stem}: shape={arr.shape}")
    z = arr.shape[0]
    chs = [arr[:, c].astype(np.float32) for c in range(arr.shape[1])]
    comp = composite_seg_3d(chs)
    print(f"composite shape: {comp.shape}, range [{comp.min():.3f}, {comp.max():.3f}]")

    # Match the mammalian preset's 2D segmentation params (seg_downsample=3)
    cfg = dict(cq.DEFAULTS)
    cfg.update(cq.CELL_TYPE_PRESETS["mammalian"])
    cfg["use_gpu"] = False  # CPU
    cfg.pop("_3d", None)
    ds = max(1, int(cfg.get("seg_downsample", 1)))
    print(f"using seg_downsample={ds}")

    print("initializing Cellpose-SAM model (one-time ~10-30s on CPU)…")
    model = cq.init_model(cfg)

    per_slice_counts = []
    per_slice_kept_areas = []
    middle_z = z // 2
    middle_mask = None
    mip_mask = None

    # Also run on the MIP and on the middle Z for reference
    mip = comp.max(axis=0)
    mip_ds = mip[::ds, ::ds] if ds > 1 else mip
    print(f"segmenting MIP ({mip_ds.shape}) for reference…")
    mip_mask_ds = cq.segment_cells(model, mip_ds, cfg)
    mip_mask = mip_mask_ds
    n_mip = int(mip_mask.max())
    print(f"  MIP: {n_mip} cells")

    for zi in range(z):
        slice_img = comp[zi]
        slice_ds = slice_img[::ds, ::ds] if ds > 1 else slice_img
        mask = cq.segment_cells(model, slice_ds, cfg)
        n_cells = int(mask.max())
        per_slice_counts.append(n_cells)
        per_slice_kept_areas.append(
            np.bincount(mask.ravel())[1:].mean() if n_cells > 0 else 0.0)
        print(f"  z={zi}: {n_cells} cells (mean area = {per_slice_kept_areas[-1]:.0f} px)")
        if zi == middle_z:
            middle_mask = mask

    # 3D stitch result already on disk
    cells_3d = setup_paths.OUT_MAMM_3D / stem / "cells.csv"
    if cells_3d.exists():
        import pandas as pd
        d3 = pd.read_csv(cells_3d)
        n3 = int(d3["keep"].sum())
        n3_total = len(d3)
    else:
        n3 = n3_total = -1

    out_dir = setup_paths.VALIDATION_ROOT / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].bar(range(z), per_slice_counts, color="C0", edgecolor="black", linewidth=0.4)
    axes[0].axhline(n_mip, color="C2", linestyle="--", linewidth=1.0,
                    label=f"MIP-based: {n_mip} cells")
    axes[0].axhline(n3, color="C3", linestyle="--", linewidth=1.0,
                    label=f"3D stitch kept: {n3} (total {n3_total})")
    axes[0].set_xlabel("Z slice")
    axes[0].set_ylabel("Cells found (2D Cellpose-SAM)")
    axes[0].set_title(f"{stem}: per-slice 2D segmentation cell counts")
    axes[0].set_xticks(range(z))
    axes[0].legend(fontsize=8)

    # Per-slice mean signal (rescaled) overlay
    axial = comp.mean(axis=(1, 2))
    axes[1].plot(range(z), axial, marker="o", color="C0", label="composite mean intensity")
    axes[1].set_xlabel("Z slice")
    axes[1].set_ylabel("Mean composite signal")
    axes[1].set_title(f"{stem}: axial signal profile")
    axes[1].set_xticks(range(z))
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_dir / "control_rep3_cells_per_z.png", dpi=130)
    plt.close(fig)

    # Write the numeric per-slice table
    import pandas as pd
    pd.DataFrame({
        "z": range(z),
        "cells_found_2d": per_slice_counts,
        "mean_kept_area_px": per_slice_kept_areas,
        "axial_mean_signal": axial,
    }).to_csv(out_dir / "control_rep3_per_slice.csv", index=False)

    print()
    print(f"=== {stem} summary ===")
    print(f"  cells found per Z slice (2D, downsample={ds}): {per_slice_counts}")
    print(f"  median per-slice count: {int(np.median(per_slice_counts))}")
    print(f"  MIP-based 2D segmentation: {n_mip} cells")
    print(f"  3D stitch kept: {n3} cells (total {n3_total})")
    print()
    print(f"  Wrote {out_dir}/control_rep3_cells_per_z.png")
    print(f"  Wrote {out_dir}/control_rep3_per_slice.csv")


if __name__ == "__main__":
    main()
