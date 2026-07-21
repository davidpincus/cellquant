"""validation_3d/seg_diagnose_part1.py — Part 1 of segmentation diagnostic.

Per-image characterization. Produces:
  - diagnostics/seg_characterization.csv
  - diagnostics/per_z_signal_profiles.png   (axial mean intensity per image)
  - diagnostics/rejected_cell_size_dist.png (volume hist for nuclei-rejected cells)
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
from skimage import filters, measure

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

# Same logic cellquant uses to build the composite cell-seg image:
# rescale every channel to [0,1] using p1/p99.8, sum, renormalize.
def composite_seg(channels_zyx: list[np.ndarray]) -> np.ndarray:
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


def quick_nuclei_estimate_2d(dapi_mip: np.ndarray) -> int:
    """Threshold + label + count for a ballpark nuclei count from the DAPI MIP."""
    try:
        thr = filters.threshold_otsu(dapi_mip)
    except Exception:
        return 0
    bw = dapi_mip > thr
    # remove small objects (< 50 px = below realistic nucleus area at 0.094 µm/px)
    lab, n = measure.label(bw, return_num=True)
    if n == 0:
        return 0
    props = measure.regionprops(lab)
    return int(sum(1 for p in props if p.area >= 100))


def characterize_image(zstack_path: Path,
                       voxel_xy_um: float = 0.094,
                       voxel_z_um: float = 0.25) -> dict:
    stem = zstack_path.stem
    print(f"  characterizing {stem}…")
    arr = tifffile.imread(zstack_path)  # (Z, C, Y, X)
    z, c, y, x = arr.shape
    dapi = arr[:, 0]
    g3bp1 = arr[:, 1]
    pabpc1 = arr[:, 2]

    # Composite the same way cellquant 3D does (sum-of-all-non-skip-channels)
    comp = composite_seg([dapi.astype(np.float32),
                          g3bp1.astype(np.float32),
                          pabpc1.astype(np.float32)])

    # SNR — defined as p95 / p5 on cell-seg input
    comp_mip = comp.max(axis=0)
    z_mid = z // 2
    comp_mid = comp[z_mid]

    def _snr(im):
        p5, p95 = np.percentile(im, [5, 95])
        return float(p95 / max(p5, 1e-6))

    snr_mip = _snr(comp_mip)
    snr_mid = _snr(comp_mid)
    # Per-channel SNR on MIP (for context)
    snr_dapi_mip = _snr(dapi.max(axis=0).astype(np.float32))
    snr_g3bp1_mip = _snr(g3bp1.max(axis=0).astype(np.float32))
    snr_pabpc1_mip = _snr(pabpc1.max(axis=0).astype(np.float32))

    # Axial signal profile: mean of composite per Z
    axial_profile = comp.mean(axis=(1, 2))
    axial_peak_z = int(np.argmax(axial_profile))
    axial_peak_val = float(axial_profile.max())
    # How "narrow" is the signal in Z? Count slices > 0.5 * peak.
    n_bright_slices = int((axial_profile > 0.5 * axial_peak_val).sum())

    # Cell counts from existing 3D run
    cells_3d_csv = setup_paths.OUT_MAMM_3D / stem / "cells.csv"
    if cells_3d_csv.exists():
        d3 = pd.read_csv(cells_3d_csv)
        n3_total = len(d3)
        n3_keep = int(d3["keep"].sum())
        n3_zero_nuc = int((d3["n_nuclei"] == 0).sum())
        n3_one_to_four = int((d3["n_nuclei"].between(1, 4)).sum())
        n3_over_four = int((d3["n_nuclei"] > 4).sum())
        rejected = d3[~d3["keep"]]
        if "cell_volume_vox" in d3.columns:
            rej_vol_med = float(rejected["cell_volume_vox"].median()) if len(rejected) else float("nan")
            rej_vol_p25 = float(rejected["cell_volume_vox"].quantile(0.25)) if len(rejected) else float("nan")
            rej_vol_p75 = float(rejected["cell_volume_vox"].quantile(0.75)) if len(rejected) else float("nan")
            kept_vol_med = float(d3.loc[d3["keep"], "cell_volume_vox"].median()) if n3_keep else float("nan")
        else:
            rej_vol_med = rej_vol_p25 = rej_vol_p75 = kept_vol_med = float("nan")
    else:
        n3_total = n3_keep = n3_zero_nuc = n3_one_to_four = n3_over_four = 0
        rej_vol_med = rej_vol_p25 = rej_vol_p75 = kept_vol_med = float("nan")

    # 2D-published cell counts on the matching MIP, when available
    pub_dir = setup_paths.OUT_MAMM_PUBLISHED2D / stem
    if (pub_dir / "cells.csv").exists():
        d2 = pd.read_csv(pub_dir / "cells.csv")
        n2_total = len(d2)
        n2_keep = int(d2["keep"].sum()) if "keep" in d2.columns else n2_total
    else:
        n2_total = -1
        n2_keep = -1

    # Quick DAPI-MIP nuclei estimate (ground truth ballpark)
    n_nuclei_estimate = quick_nuclei_estimate_2d(dapi.max(axis=0).astype(np.float32))

    return {
        "image": stem,
        "z_depth": z,
        "y": y, "x": x,
        # SNR
        "snr_composite_mip": snr_mip,
        "snr_composite_mid_z": snr_mid,
        "snr_dapi_mip": snr_dapi_mip,
        "snr_g3bp1_mip": snr_g3bp1_mip,
        "snr_pabpc1_mip": snr_pabpc1_mip,
        # axial
        "axial_peak_z": axial_peak_z,
        "n_bright_slices_gt_half_peak": n_bright_slices,
        "axial_profile": axial_profile.tolist(),
        # cell counts
        "n_3d_total": n3_total,
        "n_3d_keep": n3_keep,
        "n_3d_zero_nuclei": n3_zero_nuc,
        "n_3d_1to4_nuclei": n3_one_to_four,
        "n_3d_over_4_nuclei": n3_over_four,
        "n_2d_published_total": n2_total,
        "n_2d_published_keep": n2_keep,
        "n_dapi_mip_estimate": n_nuclei_estimate,
        # debris signature
        "rejected_volume_vox_median": rej_vol_med,
        "rejected_volume_vox_p25": rej_vol_p25,
        "rejected_volume_vox_p75": rej_vol_p75,
        "kept_volume_vox_median": kept_vol_med,
    }


def main() -> None:
    out_dir = setup_paths.VALIDATION_ROOT / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in setup_paths.list_mammalian_zstacks():
        rows.append(characterize_image(p))

    # Strip axial_profile into a separate plot, write the rest to CSV
    profiles = {r["image"]: np.array(r.pop("axial_profile")) for r in rows}
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "seg_characterization.csv", index=False)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 220)
    pd.set_option("display.float_format", lambda x: f"{x:8.3g}")
    print()
    print("=" * 80)
    print(df[[
        "image", "z_depth",
        "snr_composite_mip", "snr_composite_mid_z",
        "axial_peak_z", "n_bright_slices_gt_half_peak",
        "n_3d_total", "n_3d_keep",
        "n_3d_zero_nuclei", "n_3d_1to4_nuclei", "n_3d_over_4_nuclei",
        "n_2d_published_keep", "n_dapi_mip_estimate",
        "kept_volume_vox_median", "rejected_volume_vox_median",
    ]].to_string(index=False))

    # Axial profile plot
    fig, ax = plt.subplots(figsize=(9, 4))
    for stem, prof in profiles.items():
        zs = np.arange(len(prof))
        ax.plot(zs, prof, marker="o", label=stem)
    ax.set_xlabel("Z slice")
    ax.set_ylabel("Mean composite-seg intensity (normalized)")
    ax.set_title("Axial signal distribution per image")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "per_z_signal_profiles.png", dpi=130)
    plt.close(fig)

    # Rejected-vs-kept volume distribution
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    rej_meds = df["rejected_volume_vox_median"].values
    kept_meds = df["kept_volume_vox_median"].values
    axes[0].scatter(df["image"], df["n_3d_total"], label="total objects found", color="C0", s=60)
    axes[0].scatter(df["image"], df["n_3d_keep"], label="kept after nuclei gate", color="C2", s=60)
    axes[0].set_xticks(range(len(df)))
    axes[0].set_xticklabels(df["image"].tolist(), rotation=25, ha="right", fontsize=8)
    axes[0].set_ylabel("Cell count")
    axes[0].set_title("3D segmentation: total vs kept")
    axes[0].legend(fontsize=8)
    width = 0.4
    x = np.arange(len(df))
    axes[1].bar(x - width / 2, kept_meds, width, label="kept median", color="C2")
    axes[1].bar(x + width / 2, rej_meds, width, label="rejected median", color="C3")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(df["image"].tolist(), rotation=25, ha="right", fontsize=8)
    axes[1].set_ylabel("Median object volume (voxels)")
    axes[1].set_title("Object size: kept vs nuclei-rejected")
    axes[1].set_yscale("log")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "rejected_cell_size_dist.png", dpi=130)
    plt.close(fig)

    print()
    print(f"Wrote {out_dir}/seg_characterization.csv")
    print(f"Wrote {out_dir}/per_z_signal_profiles.png")
    print(f"Wrote {out_dir}/rejected_cell_size_dist.png")


if __name__ == "__main__":
    main()
