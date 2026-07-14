#!/usr/bin/env python3
"""
view_3d_segmentation.py

Visualize the masks produced by cellquant3d.py to verify the segmentation
is genuinely 3D (not just 2D slices stacked). Run this AFTER cellquant3d.py.

Generates four figures in <out_dir>/viz/:
    1. mip_overlay.png       — XY MIP of channels with cell + puncta outlines
    2. z_gallery.png         — same masks at six Z planes (cells should appear
                               and disappear with depth)
    3. orthogonal_views.png  — XY + XZ + YZ slices (proves cells aren't infinite
                               vertical bars)
    4. render_3d.png         — actual surface meshes via marching cubes,
                               with puncta centroids as points

Usage:
    python view_3d_segmentation.py <input_tif> <out_dir> [--voxel-size XY Z]

Example:
    python cellquant3d.py stacks/ "1:Sis1:quantify" "2:Hsp70:quantify" \\
        --out test_run/ --cell-type yeast --voxel-size 0.1 0.1 \\
        --seg-3d-method full --skip-plots --colocalization

    python view_3d_segmentation.py stacks/sample.tif test_run/ \\
        --voxel-size 0.1 0.1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import tifffile as tiff


def detect_axes(arr: np.ndarray, n_channels: int, declared: str | None = None) -> str:
    """Trust the file's declared axes if present; otherwise guess from shape."""
    if declared:
        d = declared.upper().replace("S", "C")
        kept = "".join(c for c in d if c in "ZCYX")
        if len(kept) == arr.ndim and set(kept) >= {"Y", "X"}:
            return kept
    if arr.ndim == 3:
        return "ZYX"
    shape = arr.shape
    cand_c = [i for i, s in enumerate(shape) if s == n_channels and s <= 10]
    if not cand_c:
        small = [(i, s) for i, s in enumerate(shape) if s <= 10]
        cand_c = [min(small, key=lambda x: x[1])[0]]
    c_axis = cand_c[0]
    spatial = sorted(
        [(i, s) for i, s in enumerate(shape) if i != c_axis],
        key=lambda x: x[1], reverse=True,
    )
    yx_axes = {spatial[0][0], spatial[1][0]}
    z_axis = next(i for i in range(4) if i != c_axis and i not in yx_axes)
    layout = [None, None, None, None]
    layout[z_axis] = "Z"; layout[c_axis] = "C"
    layout[spatial[0][0]] = "Y"; layout[spatial[1][0]] = "X"
    return "".join(layout)  # type: ignore[arg-type]


def load_zcyx(path: Path) -> tuple[np.ndarray, str]:
    """Load a 3D/4D TIFF and reshape to (C, Z, Y, X). Returns (arr, axes)."""
    declared = None
    try:
        with tiff.TiffFile(str(path)) as tf:
            if tf.series and tf.series[0].axes:
                declared = tf.series[0].axes
    except Exception:
        pass
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 3:
        return arr[np.newaxis], "CZYX"
    n_ch = 2 if arr.ndim == 4 and min(arr.shape) <= 10 else 1
    axes = detect_axes(arr, n_ch, declared=declared)
    perm = [axes.index(a) for a in "CZYX"]
    return np.transpose(arr, perm), axes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_tif", help="Original input TIFF (Z-stack)")
    ap.add_argument("out_dir", help="cellquant3d.py output dir (contains masks/)")
    ap.add_argument("--voxel-size", nargs=2, type=float, default=[0.1, 0.1],
                    metavar=("XY_UM", "Z_UM"))
    ap.add_argument("--channel-names", nargs="+", default=None,
                    help="Optional channel names for axis labels (e.g. Sis1 Hsp70)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage import measure

    in_path = Path(args.input_tif)
    out_root = Path(args.out_dir)
    masks_dir = out_root / "masks"
    viz_dir = out_root / "viz"
    viz_dir.mkdir(parents=True, exist_ok=True)
    stem = in_path.stem

    # Load original image
    czyx, src_axes = load_zcyx(in_path)
    print(f"Input: {in_path.name}, shape {czyx.shape} (CZYX), source axes {src_axes}")
    n_channels, Z, Y, X = czyx.shape
    voxel_xy, voxel_z = args.voxel_size

    if args.channel_names and len(args.channel_names) == n_channels:
        ch_names = args.channel_names
    else:
        ch_names = [f"ch{i+1}" for i in range(n_channels)]

    # Load cell mask
    cell_path = masks_dir / f"{stem}_cellmask.tif"
    if not cell_path.exists():
        print(f"ERROR: cell mask not found at {cell_path}")
        print("Run cellquant3d.py first.")
        return 1
    cell_mask = np.asarray(tiff.imread(str(cell_path)))
    n_cells = int(cell_mask.max())
    print(f"Cells: {n_cells}")

    # Load puncta masks per channel (if present)
    puncta_masks: dict[str, np.ndarray] = {}
    for name in ch_names:
        for candidate in masks_dir.glob(f"{stem}_*{name}*_punctamask.tif"):
            puncta_masks[name] = np.asarray(tiff.imread(str(candidate)))
            break
        else:
            # Try generic patterns
            for candidate in masks_dir.glob(f"{stem}_*_punctamask.tif"):
                if name.lower() in candidate.name.lower():
                    puncta_masks[name] = np.asarray(tiff.imread(str(candidate)))
                    break
    for name, pm in puncta_masks.items():
        print(f"  {name} puncta: {int(pm.max())}")

    # Pre-compute display ranges
    p99 = [float(np.percentile(czyx[c], 99.5)) for c in range(n_channels)]

    # ---- Figure 1: MIP overlay ----
    fig, ax = plt.subplots(figsize=(9, 8.5), dpi=120)
    rgb = np.zeros((Y, X, 3))
    if n_channels >= 1:
        rgb[..., 1] = np.clip(czyx[0].max(0) / max(p99[0], 1), 0, 1)
    if n_channels >= 2:
        rgb[..., 0] = np.clip(czyx[1].max(0) / max(p99[1], 1), 0, 1)
    ax.imshow(rgb); ax.set_axis_off()
    mip_cell = cell_mask.max(0)
    for cid in np.unique(mip_cell):
        if cid == 0: continue
        ax.contour(mip_cell == cid, levels=[0.5], colors="cyan", linewidths=0.9)
        ys, xs = np.where(mip_cell == cid)
        ax.text(xs.mean(), ys.mean(), str(int(cid)), color="yellow",
                fontsize=11, ha="center", weight="bold")
    colors_p = ["magenta", "white", "orange", "deepskyblue"]
    for i, (name, pm) in enumerate(puncta_masks.items()):
        mip_p = pm.max(0)
        if mip_p.max() == 0: continue
        ax.contour(mip_p > 0, levels=[0.5], colors=colors_p[i % len(colors_p)], linewidths=0.6)
    ax.set_title(f"MIP overlay — {n_cells} cells, "
                 + ", ".join(f"{n}: {int(p.max())} puncta" for n, p in puncta_masks.items()),
                 fontsize=11)
    plt.tight_layout()
    plt.savefig(viz_dir / "mip_overlay.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {viz_dir / 'mip_overlay.png'}")

    # ---- Figure 2: Z gallery ----
    z_planes = np.linspace(Z * 0.15, Z * 0.85, 6).astype(int)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=120)
    for ax, z in zip(axes.flat, z_planes):
        rgb = np.zeros((Y, X, 3))
        if n_channels >= 1:
            rgb[..., 1] = np.clip(czyx[0, z] / max(p99[0], 1), 0, 1)
        if n_channels >= 2:
            rgb[..., 0] = np.clip(czyx[1, z] / max(p99[1], 1), 0, 1)
        ax.imshow(rgb); ax.set_axis_off()
        ax.set_title(f"Z = {z}  ({z * voxel_z:.1f} µm)", fontsize=11)
        cell_z = cell_mask[z]
        for cid in np.unique(cell_z):
            if cid == 0: continue
            ax.contour(cell_z == cid, levels=[0.5], colors="cyan", linewidths=1.0)
            ys, xs = np.where(cell_z == cid)
            ax.text(xs.mean(), ys.mean(), str(int(cid)), color="yellow",
                    fontsize=10, ha="center", weight="bold")
        for i, (name, pm) in enumerate(puncta_masks.items()):
            if pm[z].max() > 0:
                ax.contour(pm[z] > 0, levels=[0.5],
                           colors=colors_p[i % len(colors_p)], linewidths=0.6)
    plt.suptitle("Same masks at six Z planes — cells should appear, peak, and disappear "
                 "with depth (proves segmentation is volumetric, not 2D-stacked)",
                 fontsize=11, y=0.995)
    plt.tight_layout()
    plt.savefig(viz_dir / "z_gallery.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {viz_dir / 'z_gallery.png'}")

    # ---- Figure 3: orthogonal views ----
    from scipy import ndimage as ndi
    if n_cells > 0:
        centroids = ndi.center_of_mass(
            cell_mask > 0, cell_mask, range(1, n_cells + 1)
        )
        centroids = np.array(centroids)
        y_slice = int(np.median(centroids[:, 1]))
        x_slice = int(np.median(centroids[:, 2]))
    else:
        y_slice, x_slice = Y // 2, X // 2

    def _rgb(s2d, h2d):
        rgb = np.zeros((*s2d.shape, 3))
        if n_channels >= 1:
            rgb[..., 1] = np.clip(s2d / max(p99[0], 1), 0, 1)
        if n_channels >= 2:
            rgb[..., 0] = np.clip(h2d / max(p99[1], 1), 0, 1)
        return rgb

    s = czyx[0]; h = czyx[1] if n_channels >= 2 else czyx[0]
    xy_rgb = _rgb(s.max(0), h.max(0))
    xz_rgb = _rgb(s[:, y_slice, :], h[:, y_slice, :])
    yz_rgb = _rgb(s[:, :, x_slice], h[:, :, x_slice])
    xy_mask = cell_mask.max(0)
    xz_mask = cell_mask[:, y_slice, :]
    yz_mask = cell_mask[:, :, x_slice]

    fig = plt.figure(figsize=(13, 12), dpi=120)
    gs = fig.add_gridspec(2, 2, width_ratios=[X, Z], height_ratios=[Y, Z],
                          wspace=0.05, hspace=0.05)
    ax_xy = fig.add_subplot(gs[0, 0])
    ax_yz = fig.add_subplot(gs[0, 1])
    ax_xz = fig.add_subplot(gs[1, 0])
    ax_xy.imshow(xy_rgb, aspect=1)
    ax_xy.set_title("XY (max projection through Z)", fontsize=11)
    for cid in np.unique(xy_mask):
        if cid == 0: continue
        ax_xy.contour(xy_mask == cid, levels=[0.5], colors="cyan", linewidths=0.9)
        ys, xs = np.where(xy_mask == cid)
        ax_xy.text(xs.mean(), ys.mean(), str(int(cid)),
                   color="yellow", fontsize=11, ha="center", weight="bold")
    ax_xy.axhline(y_slice, color="lime", linewidth=0.5, linestyle="--", alpha=0.6)
    ax_xy.axvline(x_slice, color="lime", linewidth=0.5, linestyle="--", alpha=0.6)
    ax_xy.set_xticks([]); ax_xy.set_yticks([])

    ax_xz.imshow(xz_rgb, aspect=1)
    ax_xz.set_title(f"XZ slice at Y={y_slice}", fontsize=11)
    ax_xz.set_ylabel("Z (slices)")
    for cid in np.unique(xz_mask):
        if cid == 0: continue
        ax_xz.contour(xz_mask == cid, levels=[0.5], colors="cyan", linewidths=0.9)
        zs, xs = np.where(xz_mask == cid)
        ax_xz.text(xs.mean(), zs.mean(), str(int(cid)),
                   color="yellow", fontsize=11, ha="center", weight="bold")
    ax_xz.axvline(x_slice, color="lime", linewidth=0.5, linestyle="--", alpha=0.6)
    ax_xz.set_xticks([])

    ax_yz.imshow(yz_rgb.transpose(1, 0, 2), aspect=1)
    ax_yz.set_title(f"YZ slice at X={x_slice}", fontsize=11)
    yz_mask_t = yz_mask.T
    for cid in np.unique(yz_mask_t):
        if cid == 0: continue
        ax_yz.contour(yz_mask_t == cid, levels=[0.5], colors="cyan", linewidths=0.9)
        ys, zs = np.where(yz_mask_t == cid)
        ax_yz.text(zs.mean(), ys.mean(), str(int(cid)),
                   color="yellow", fontsize=11, ha="center", weight="bold")
    ax_yz.axhline(y_slice, color="lime", linewidth=0.5, linestyle="--", alpha=0.6)
    ax_yz.set_xlabel("Z (slices)")
    ax_yz.set_yticks([])

    plt.suptitle("Orthogonal views — true 3D cells are round in XY AND finite in Z. "
                 "Green dashed lines mark slice positions.", fontsize=11, y=0.995)
    plt.savefig(viz_dir / "orthogonal_views.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {viz_dir / 'orthogonal_views.png'}")

    # ---- Figure 4: 3D surface render ----
    spacing = (voxel_z, voxel_xy, voxel_xy)
    fig = plt.figure(figsize=(16, 7), dpi=120)
    cmap = plt.cm.tab10

    # Compute puncta centroids in physical units once
    punct_pts: dict[str, np.ndarray] = {}
    for name, pm in puncta_masks.items():
        ids = np.unique(pm); ids = ids[ids != 0]
        cents = []
        for pid in ids:
            zs, ys, xs = np.where(pm == pid)
            cents.append((zs.mean() * voxel_z, ys.mean() * voxel_xy, xs.mean() * voxel_xy))
        if cents:
            punct_pts[name] = np.array(cents)

    for view_idx, (elev, azim) in enumerate([(20, -60), (15, 30)]):
        ax = fig.add_subplot(1, 2, view_idx + 1, projection="3d")
        for cid in range(1, n_cells + 1):
            sub = (cell_mask == cid).astype(np.uint8)
            sub_pad = np.pad(sub, 1, mode="constant", constant_values=0)
            try:
                verts, faces, _, _ = measure.marching_cubes(sub_pad, level=0.5, spacing=spacing)
                verts -= np.array([voxel_z, voxel_xy, voxel_xy])
            except Exception:
                continue
            mesh = Poly3DCollection(verts[faces], alpha=0.35,
                                    facecolor=cmap((cid - 1) % 10),
                                    edgecolor="none")
            ax.add_collection3d(mesh)
        markers = [("o", "lime"), ("^", "red"), ("s", "deepskyblue"), ("D", "orange")]
        for i, (name, pts) in enumerate(punct_pts.items()):
            m, c = markers[i % len(markers)]
            ax.scatter(pts[:, 2], pts[:, 1], pts[:, 0],
                       c=c, s=18, marker=m, alpha=0.95,
                       edgecolors="black", linewidth=0.3,
                       label=f"{name} ({len(pts)})")
        ax.set_xlim(0, X * voxel_xy); ax.set_ylim(0, Y * voxel_xy); ax.set_zlim(0, Z * voxel_z)
        ax.set_xlabel("X (µm)"); ax.set_ylabel("Y (µm)"); ax.set_zlabel("Z (µm)")
        ax.set_box_aspect((X, Y, Z))
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(f"View {view_idx + 1} (elev={elev}°, azim={azim}°)", fontsize=10)
        if punct_pts:
            ax.legend(loc="upper left", fontsize=9)
    plt.suptitle(f"3D surface meshes — {n_cells} cells, with detected puncta as points",
                 fontsize=12, y=0.98)
    plt.tight_layout()
    plt.savefig(viz_dir / "render_3d.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {viz_dir / 'render_3d.png'}")

    print(f"\nDone. Open {viz_dir}/ to inspect.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
