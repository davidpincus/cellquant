"""validation_3d/diagnose_sphericity.py

Test cellquant's sphericity routine at varying sphere radii to discriminate
between (a) a hidden unit/spacing bug and (b) marching-cubes discretization
error at small radii.

Logic:
  - Generate perfect synthetic spheres at radii 0.3, 0.6, 1.0, 1.5, 2.0, 3.0,
    5.0 µm with the yeast voxel size (0.094 µm XY, 0.1 µm Z).
  - Compute reported vs analytic sphericity for each.
  - Plot relative error vs radius.

  - If the unit handling is wrong, relative error stays constant or grows with
    radius — it doesn't care about how many voxels resolve the surface.
  - If marching-cubes discretization is the cause, relative error scales
    approximately as 1/radius: more voxels per surface area give a finer mesh
    and the error converges toward 0 as resolution increases.

Outputs:
  - validation_3d/diagnostics/sphericity_vs_radius.csv
  - validation_3d/diagnostics/sphericity_vs_radius.png
  - Verdict printed to stdout.
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

VOXEL_XY = 0.094
VOXEL_Z = 0.1
SPACING = (VOXEL_Z, VOXEL_XY, VOXEL_XY)
VOXEL_UM3 = VOXEL_Z * VOXEL_XY * VOXEL_XY

RADII_UM = [0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]


def make_sphere_mask(radius_um: float, pad: int = 6) -> np.ndarray:
    rz = radius_um / VOXEL_Z
    rxy = radius_um / VOXEL_XY
    nz = int(np.ceil(2 * rz)) + 2 * pad
    ny = int(np.ceil(2 * rxy)) + 2 * pad
    nx = ny
    z, y, x = np.indices((nz, ny, nx)).astype(np.float32)
    cz, cy, cx = nz / 2, ny / 2, nx / 2
    val = ((z - cz) / rz) ** 2 + ((y - cy) / rxy) ** 2 + ((x - cx) / rxy) ** 2
    return (val <= 1.0).astype(np.uint8)


def sphericity_via_cellquant(mask: np.ndarray) -> tuple[float, int, float]:
    """Same algorithm cellquant uses inside compute_nucleolar_morphology."""
    sub_pad = np.pad(mask, 1, mode="constant", constant_values=False)
    verts, faces, _, _ = measure.marching_cubes(
        sub_pad.astype(np.uint8), level=0.5, spacing=SPACING)
    surface_um2 = float(measure.mesh_surface_area(verts, faces))
    vol_vox = int(mask.sum())
    V = vol_vox * VOXEL_UM3
    sphericity = (np.pi ** (1.0 / 3.0)) * (6.0 * V) ** (2.0 / 3.0) / surface_um2
    return float(sphericity), vol_vox, surface_um2


def main() -> None:
    out_dir = setup_paths.VALIDATION_ROOT / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for r in RADII_UM:
        mask = make_sphere_mask(r)
        s, vol_vox, surf_um2 = sphericity_via_cellquant(mask)
        rxy_vox = r / VOXEL_XY
        rz_vox = r / VOXEL_Z
        analytic_vol_um3 = (4.0 / 3.0) * np.pi * r ** 3
        analytic_surf_um2 = 4.0 * np.pi * r ** 2
        rel_err_sphericity = abs(s - 1.0)
        vol_err = (vol_vox * VOXEL_UM3 - analytic_vol_um3) / analytic_vol_um3
        surf_err = (surf_um2 - analytic_surf_um2) / analytic_surf_um2
        rows.append({
            "radius_um": r,
            "rxy_vox": rxy_vox,
            "rz_vox": rz_vox,
            "n_vox": vol_vox,
            "vol_reported_um3": vol_vox * VOXEL_UM3,
            "vol_analytic_um3": analytic_vol_um3,
            "vol_rel_err": vol_err,
            "surf_reported_um2": surf_um2,
            "surf_analytic_um2": analytic_surf_um2,
            "surf_rel_err": surf_err,
            "sphericity_reported": s,
            "sphericity_rel_err": rel_err_sphericity,
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "sphericity_vs_radius.csv", index=False)
    print("Sphericity diagnostic (perfect spheres, yeast voxel = 0.094 XY / 0.1 Z µm):")
    print()
    print(df[[
        "radius_um", "rxy_vox", "n_vox",
        "vol_rel_err", "surf_rel_err", "sphericity_reported", "sphericity_rel_err",
    ]].to_string(index=False, float_format=lambda x: f"{x:8.4g}"))
    print()

    # 1/r fit on relative error: log(err) ~ -log(r) + const
    log_r = np.log(df["radius_um"].values)
    log_err = np.log(np.clip(df["sphericity_rel_err"].values, 1e-12, None))
    slope, intercept = np.polyfit(log_r, log_err, 1)
    print(f"log(rel_err) vs log(radius) slope = {slope:+.3f}")
    print("  (a slope near -1 ⇒ 1/r scaling ⇒ discretization-dominated;")
    print("   a slope near 0 ⇒ constant error ⇒ unit bug)")
    print()

    if slope < -0.5:
        verdict = "DISCRETIZATION (marching-cubes resolution-limited at small radii)"
    elif slope > -0.2:
        verdict = "UNIT BUG (error doesn't scale with radius)"
    else:
        verdict = "AMBIGUOUS"
    print(f"VERDICT: {verdict}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(df["radius_um"], df["sphericity_reported"], "o-", color="C0")
    axes[0].axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("Sphere radius (µm)")
    axes[0].set_ylabel("Reported sphericity")
    axes[0].set_title("Sphericity convergence to 1.0")

    axes[1].loglog(df["radius_um"], df["sphericity_rel_err"], "o-", color="C1",
                   label="data")
    # Reference 1/r line through the first point
    r0, e0 = df["radius_um"].iloc[0], df["sphericity_rel_err"].iloc[0]
    ref_r = np.linspace(r0, df["radius_um"].max(), 50)
    ref_err = e0 * (r0 / ref_r)
    axes[1].loglog(ref_r, ref_err, "--", color="gray", linewidth=0.8,
                   label="reference 1/r")
    axes[1].set_xlabel("Sphere radius (µm)")
    axes[1].set_ylabel("|sphericity − 1|")
    axes[1].set_title(
        f"Relative sphericity error vs radius (log-log)\n"
        f"slope = {slope:+.3f}")
    axes[1].legend(fontsize=8)
    fig.suptitle(verdict)
    fig.tight_layout()
    fig.savefig(out_dir / "sphericity_vs_radius.png", dpi=130)
    plt.close(fig)
    print(f"Wrote {out_dir}/sphericity_vs_radius.{{csv,png}}")


if __name__ == "__main__":
    main()
