"""Reusable cellquant runner for Part B experiments. Spawns cellquant with
a per-experiment set of override flags and writes outputs to a labelled dir.
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


def run(label: str, stem: str, *, stitch_threshold: float | None = None,
        min_cell_volume_vox: int | None = None,
        cell_seg_channel: str | None = None,
        nuclei_diameter: float | None = None) -> dict:
    out = (setup_paths.VALIDATION_ROOT / "diagnostics" / "seg_experiments" /
           label / stem)
    out.mkdir(parents=True, exist_ok=True)
    staging = setup_paths.VALIDATION_ROOT / "_staging_seg_exp" / f"{label}_{stem}"
    staging.mkdir(parents=True, exist_ok=True)
    zsta = setup_paths.MAMMALIAN_ZSTACK_DIR / f"{stem}.tif"
    link = staging / zsta.name
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(zsta)

    argv = [
        sys.executable, str(setup_paths.CELLQUANT_SCRIPT),
        str(staging),
        *setup_paths.MAMMALIAN_CHANNELS,
        "--cell-type", setup_paths.MAMMALIAN_CELL_TYPE,
        "--mode", "3d",
        "--voxel-size", str(setup_paths.MAMMALIAN_VOXEL_XY_UM),
        str(setup_paths.MAMMALIAN_VOXEL_Z_UM),
        "--condensate-index",
        "--colocalization",
        "--filename-pattern", setup_paths.MAMMALIAN_FILENAME_PATTERN_3D,
        "--skip-plots",
        "--no-gpu",
        "--out", str(out),
    ]
    if stitch_threshold is not None:
        argv += ["--stitch-threshold", str(stitch_threshold)]
    if min_cell_volume_vox is not None:
        argv += ["--min-cell-volume-vox", str(min_cell_volume_vox)]
    if cell_seg_channel is not None:
        argv += ["--cell-seg-channel", cell_seg_channel]
    if nuclei_diameter is not None:
        argv += ["--nuclei-diameter", str(nuclei_diameter)]

    env = dict(os.environ)
    env.update(setup_paths.ENV_OVERRIDES)
    t0 = time.time()
    result = subprocess.run(
        argv, capture_output=True, text=True,
        timeout=setup_paths.PER_IMAGE_TIMEOUT_SEC, env=env,
    )
    dt = time.time() - t0
    return {"label": label, "stem": stem,
            "duration_sec": dt,
            "returncode": result.returncode,
            "stdout_tail": result.stdout[-1000:],
            "stderr_tail": result.stderr[-500:]}


def summarize_cells(label: str, stem: str) -> dict:
    """Read cells.csv from one experiment output and report headline stats."""
    import pandas as pd
    cf = (setup_paths.VALIDATION_ROOT / "diagnostics" / "seg_experiments" /
          label / stem / "cells.csv")
    if not cf.exists():
        return {"label": label, "stem": stem, "n_total": 0, "n_keep": 0,
                "n_0nuc": 0, "median_kept_vol": float("nan"),
                "median_rej_vol": float("nan")}
    df = pd.read_csv(cf)
    kept = df[df["keep"]]
    rej = df[~df["keep"]]
    return {
        "label": label, "stem": stem,
        "n_total": int(len(df)),
        "n_keep": int(len(kept)),
        "n_0nuc": int((df["n_nuclei"] == 0).sum()),
        "median_kept_vol": float(kept["cell_volume_vox"].median()) if len(kept) else float("nan"),
        "median_rej_vol": float(rej["cell_volume_vox"].median()) if len(rej) else float("nan"),
    }
