"""validation_3d/seg_diagnose_part2.py — seg_3d_method comparison.

Re-run cellquant 3D on control_rep3 and arsenite_rep3 with both:
  --seg-3d-method stitch  (default mammalian behavior, already on disk; we re-run
                          to a comparison directory for clean side-by-side)
  --seg-3d-method full    (volumetric Cellpose-SAM with do_3D=True)

Saves outputs to validation_3d/diagnostics/seg_mode_compare/{image}/{mode}/
"""
from __future__ import annotations

import importlib
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

IMAGES = ["control_rep3", "arsenite_rep3"]
MODES = ["stitch", "full"]


def run_one(stem: str, mode: str) -> dict:
    out_dir = setup_paths.VALIDATION_ROOT / "diagnostics" / "seg_mode_compare" / stem / mode
    out_dir.mkdir(parents=True, exist_ok=True)

    staging = setup_paths.VALIDATION_ROOT / "_staging_seg_compare" / f"{stem}_{mode}"
    staging.mkdir(parents=True, exist_ok=True)
    zstack = setup_paths.MAMMALIAN_ZSTACK_DIR / f"{stem}.tif"
    link = staging / zstack.name
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(zstack)

    argv = [
        sys.executable, str(setup_paths.CELLQUANT_SCRIPT),
        str(staging),
        *setup_paths.MAMMALIAN_CHANNELS,
        "--cell-type", setup_paths.MAMMALIAN_CELL_TYPE,
        "--mode", "3d",
        "--voxel-size", str(setup_paths.MAMMALIAN_VOXEL_XY_UM),
        str(setup_paths.MAMMALIAN_VOXEL_Z_UM),
        "--seg-3d-method", mode,
        "--condensate-index",
        "--colocalization",
        "--filename-pattern", setup_paths.MAMMALIAN_FILENAME_PATTERN_3D,
        "--skip-plots",
        "--no-gpu",
        "--out", str(out_dir),
    ]
    env = dict(os.environ)
    env.update(setup_paths.ENV_OVERRIDES)
    t0 = time.time()
    try:
        result = subprocess.run(
            argv, capture_output=True, text=True, timeout=setup_paths.PER_IMAGE_TIMEOUT_SEC,
            env=env,
        )
        dt = time.time() - t0
        status = "ok" if result.returncode == 0 else f"failed ({result.returncode})"
        stdout_tail = result.stdout[-1500:]
    except Exception as exc:
        dt = time.time() - t0
        status = f"crash: {type(exc).__name__}: {exc}"
        stdout_tail = ""
    print(f"  {stem} / {mode}: {status} ({dt:.0f}s)")
    return {"image": stem, "mode": mode, "duration_sec": dt, "status": status,
            "stdout_tail": stdout_tail}


def main() -> None:
    print("=== Part 2: seg_3d_method comparison ===")
    records = []
    for stem in IMAGES:
        for mode in MODES:
            records.append(run_one(stem, mode))

    # Summarize from each cells.csv
    import pandas as pd
    rows = []
    out_root = setup_paths.VALIDATION_ROOT / "diagnostics" / "seg_mode_compare"
    for stem in IMAGES:
        for mode in MODES:
            cf = out_root / stem / mode / "cells.csv"
            if not cf.exists():
                rows.append({"image": stem, "mode": mode, "n_total": 0, "n_keep": 0,
                             "median_kept_vol_vox": float("nan"),
                             "median_rejected_vol_vox": float("nan")})
                continue
            df = pd.read_csv(cf)
            kept = df[df["keep"]] if "keep" in df.columns else df
            rejected = df[~df["keep"]] if "keep" in df.columns else df.iloc[0:0]
            rows.append({
                "image": stem,
                "mode": mode,
                "n_total": len(df),
                "n_keep": len(kept),
                "n_zero_nuclei": int((df.get("n_nuclei", 0) == 0).sum())
                                  if "n_nuclei" in df.columns else 0,
                "median_kept_vol_vox": float(kept["cell_volume_vox"].median())
                                        if len(kept) and "cell_volume_vox" in kept.columns
                                        else float("nan"),
                "median_rejected_vol_vox": float(rejected["cell_volume_vox"].median())
                                            if len(rejected) and "cell_volume_vox" in rejected.columns
                                            else float("nan"),
            })
    summary = pd.DataFrame(rows)
    pd.set_option("display.width", 140)
    print()
    print("=" * 78)
    print(summary.to_string(index=False, float_format=lambda x: f"{x:8.4g}"))

    summary.to_csv(out_root / "summary.csv", index=False)
    print(f"\nWrote {out_root}/summary.csv")


if __name__ == "__main__":
    main()
