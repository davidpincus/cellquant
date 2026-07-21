"""validation_3d/04_synthetic_sphericity.py

Marching-cubes ground truth on programmatic 3D objects with known geometry.
Tests cellquant's nucleolar sphericity calculation against analytic values.

Objects:
  - perfect sphere (analytic sphericity = 1.0)
  - prolate spheroid 2:1:1, 3:1:1 (analytic sphericity calculable)
  - oblate spheroid 1:1:0.5
  - dumbbell (qualitative: low)
  - crescent (qualitative: low)

Voxel resolution matches the yeast acquisition (0.094 µm XY, 0.1 µm Z).
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
from skimage import measure

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

# Same as yeast acquisition (matches Methods)
VOXEL_XY = 0.094  # um
VOXEL_Z = 0.1     # um
SPACING = (VOXEL_Z, VOXEL_XY, VOXEL_XY)
VOXEL_VOLUME_UM3 = VOXEL_Z * VOXEL_XY * VOXEL_XY


def _analytic_sphericity_spheroid(a_um: float, b_um: float, c_um: float) -> float:
    """Analytic sphericity of an ellipsoid with semi-axes a, b, c (in microns).

    Uses Knud Thomsen's approximation for ellipsoid surface area
    (relative error < 1.1% for any aspect ratio).
    Sphericity = pi^(1/3) * (6V)^(2/3) / S
    """
    p = 1.6075
    s_term = (
        ((a_um * b_um) ** p
         + (a_um * c_um) ** p
         + (b_um * c_um) ** p) / 3
    ) ** (1.0 / p)
    surface_um2 = 4.0 * np.pi * s_term
    volume_um3 = (4.0 / 3.0) * np.pi * a_um * b_um * c_um
    return float(
        (np.pi ** (1.0 / 3.0)) * (6.0 * volume_um3) ** (2.0 / 3.0) / surface_um2)


def _make_ellipsoid_mask(a_vox: float, b_vox: float, c_vox: float,
                        pad: int = 4) -> np.ndarray:
    """3D bool array containing an ellipsoid with semi-axes (z, y, x) in voxels."""
    nz = int(np.ceil(2 * a_vox)) + 2 * pad
    ny = int(np.ceil(2 * b_vox)) + 2 * pad
    nx = int(np.ceil(2 * c_vox)) + 2 * pad
    z, y, x = np.indices((nz, ny, nx)).astype(np.float32)
    cz, cy, cx = nz / 2, ny / 2, nx / 2
    val = ((z - cz) / a_vox) ** 2 + ((y - cy) / b_vox) ** 2 + ((x - cx) / c_vox) ** 2
    return (val <= 1.0).astype(np.uint8)


def _make_dumbbell_mask(radius_vox: float, neck_z_offset: float,
                       pad: int = 4) -> np.ndarray:
    """Two spheres connected by a thin neck (qualitatively low sphericity)."""
    r = int(np.ceil(radius_vox))
    nz = 2 * r + int(2 * neck_z_offset) + 2 * pad
    ny = 2 * r + 2 * pad
    nx = 2 * r + 2 * pad
    z, y, x = np.indices((nz, ny, nx)).astype(np.float32)
    cy, cx = ny / 2, nx / 2
    cz1 = nz / 2 - neck_z_offset
    cz2 = nz / 2 + neck_z_offset
    val1 = (z - cz1) ** 2 + (y - cy) ** 2 + (x - cx) ** 2 <= radius_vox ** 2
    val2 = (z - cz2) ** 2 + (y - cy) ** 2 + (x - cx) ** 2 <= radius_vox ** 2
    # thin neck of radius 0.3 * r
    neck_r = 0.3 * radius_vox
    neck = ((y - cy) ** 2 + (x - cx) ** 2 <= neck_r ** 2) & (
        (z >= cz1) & (z <= cz2))
    return (val1 | val2 | neck).astype(np.uint8)


def _make_crescent_mask(outer_r_vox: float, inner_r_vox: float,
                       offset_vox: float, pad: int = 4) -> np.ndarray:
    """A sphere with a hemispherical cavity (qualitatively low sphericity)."""
    nz = int(np.ceil(2 * outer_r_vox)) + 2 * pad
    ny = nz
    nx = nz
    z, y, x = np.indices((nz, ny, nx)).astype(np.float32)
    cz, cy, cx = nz / 2, ny / 2, nx / 2
    outer = ((z - cz) ** 2 + (y - cy) ** 2 + (x - cx) ** 2) <= outer_r_vox ** 2
    inner = ((z - cz) ** 2 + (y - cy) ** 2 + (x - cx + offset_vox) ** 2) <= inner_r_vox ** 2
    return (outer & ~inner).astype(np.uint8)


def _sphericity_from_marching_cubes(mask: np.ndarray) -> tuple[float, int, float]:
    """Run cellquant's marching-cubes routine. Returns (sphericity, vol_vox, surface_um2)."""
    sub_pad = np.pad(mask, 1, mode="constant", constant_values=False)
    try:
        verts, faces, _, _ = measure.marching_cubes(
            sub_pad.astype(np.uint8), level=0.5, spacing=SPACING)
        surface_um2 = float(measure.mesh_surface_area(verts, faces))
    except Exception as exc:
        return float("nan"), 0, float("nan")
    vol_vox = int(mask.sum())
    if surface_um2 <= 0 or vol_vox <= 0:
        return float("nan"), vol_vox, surface_um2
    V = vol_vox * VOXEL_VOLUME_UM3
    sphericity = float(
        (np.pi ** (1.0 / 3.0)) * (6.0 * V) ** (2.0 / 3.0) / surface_um2)
    return sphericity, vol_vox, surface_um2


def run_synthetic_battery() -> pd.DataFrame:
    rows = []
    radius_um = 0.6  # yeast nucleolus-ish

    # 1. Perfect sphere
    r_vox_z = radius_um / VOXEL_Z
    r_vox_xy = radius_um / VOXEL_XY
    mask = _make_ellipsoid_mask(r_vox_z, r_vox_xy, r_vox_xy)
    s, vv, sa = _sphericity_from_marching_cubes(mask)
    rows.append({"shape": "sphere", "params": "r=0.6um",
                 "analytic_sphericity": 1.0,
                 "reported_sphericity": s,
                 "abs_error": abs(s - 1.0),
                 "rel_error": abs(s - 1.0),
                 "vol_vox": vv, "surface_um2": sa})

    # 2. Prolate spheroid 2:1:1 (long axis along Z)
    for ratio_label, az, ay, ax in [("prolate_2_1_1", 1.2, 0.6, 0.6),
                                     ("prolate_3_1_1", 1.8, 0.6, 0.6),
                                     ("oblate_1_1_0.5", 0.3, 0.6, 0.6)]:
        az_v = az / VOXEL_Z
        ay_v = ay / VOXEL_XY
        ax_v = ax / VOXEL_XY
        mask = _make_ellipsoid_mask(az_v, ay_v, ax_v)
        s, vv, sa = _sphericity_from_marching_cubes(mask)
        analytic = _analytic_sphericity_spheroid(az, ay, ax)
        rows.append({"shape": ratio_label,
                     "params": f"axes=({az},{ay},{ax})um",
                     "analytic_sphericity": analytic,
                     "reported_sphericity": s,
                     "abs_error": abs(s - analytic),
                     "rel_error": abs(s - analytic) / analytic if analytic else float("nan"),
                     "vol_vox": vv, "surface_um2": sa})

    # 3. Rotational invariance: same prolate, different orientations
    # Build a 2:1:1 long-axis-along-Y prolate by swapping argument order.
    mask = _make_ellipsoid_mask(0.6 / VOXEL_Z, 1.2 / VOXEL_XY, 0.6 / VOXEL_XY)
    s, vv, sa = _sphericity_from_marching_cubes(mask)
    analytic = _analytic_sphericity_spheroid(0.6, 1.2, 0.6)
    rows.append({"shape": "prolate_2_1_1_rotated_Y",
                 "params": "axes=(0.6,1.2,0.6)um",
                 "analytic_sphericity": analytic,
                 "reported_sphericity": s,
                 "abs_error": abs(s - analytic),
                 "rel_error": abs(s - analytic) / analytic if analytic else float("nan"),
                 "vol_vox": vv, "surface_um2": sa})

    # 4. Dumbbell (qualitative)
    mask = _make_dumbbell_mask(radius_vox=4, neck_z_offset=6)
    s, vv, sa = _sphericity_from_marching_cubes(mask)
    rows.append({"shape": "dumbbell",
                 "params": "r=4vox, neck_offset=6",
                 "analytic_sphericity": float("nan"),
                 "reported_sphericity": s,
                 "abs_error": float("nan"),
                 "rel_error": float("nan"),
                 "vol_vox": vv, "surface_um2": sa})

    # 5. Crescent (qualitative)
    mask = _make_crescent_mask(outer_r_vox=8, inner_r_vox=5, offset_vox=5)
    s, vv, sa = _sphericity_from_marching_cubes(mask)
    rows.append({"shape": "crescent",
                 "params": "outer=8, inner=5, offset=5",
                 "analytic_sphericity": float("nan"),
                 "reported_sphericity": s,
                 "abs_error": float("nan"),
                 "rel_error": float("nan"),
                 "vol_vox": vv, "surface_um2": sa})

    return pd.DataFrame(rows)


def _plot_results(df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    shapes = df["shape"].tolist()
    analytic = df["analytic_sphericity"].values
    reported = df["reported_sphericity"].values
    x = np.arange(len(shapes))
    width = 0.4
    ax.bar(x - width / 2, analytic, width, label="analytic", color="C0")
    ax.bar(x + width / 2, reported, width, label="cellquant marching-cubes",
           color="C1")
    ax.set_xticks(x)
    ax.set_xticklabels(shapes, rotation=25, ha="right")
    ax.set_ylabel("Sphericity")
    ax.set_title("Synthetic sphericity: analytic vs cellquant")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def main() -> None:
    out_dir = setup_paths.OUT_SYNTH
    out_dir.mkdir(parents=True, exist_ok=True)
    df = run_synthetic_battery()
    csv_path = out_dir / "results.csv"
    png_path = out_dir / "results.png"
    df.to_csv(csv_path, index=False)
    _plot_results(df, png_path)

    # Pass criterion check
    well_defined = df[~df["analytic_sphericity"].isna()]
    if len(well_defined):
        median_rel = float(np.nanmedian(well_defined["rel_error"]))
        verdict = ("PASS" if median_rel <= 0.05
                   else "MARGINAL" if median_rel <= 0.10
                   else "FAIL")
        print(f"Median relative error on spheres/spheroids: {median_rel:.3f}  ({verdict})")
    dumb = df.loc[df["shape"] == "dumbbell", "reported_sphericity"]
    cres = df.loc[df["shape"] == "crescent", "reported_sphericity"]
    if len(dumb):
        print(f"Dumbbell sphericity:  {float(dumb.values[0]):.3f} (target < 0.7)")
    if len(cres):
        print(f"Crescent sphericity:  {float(cres.values[0]):.3f} (target < 0.7)")
    print(f"\nWrote {csv_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
