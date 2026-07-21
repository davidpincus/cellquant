"""validation_3d/01_run_pipelines.py

Runs the three input passes per dataset:
  (a) cellquant in 3D on each z-stack
  (b) matched-2D-from-3D quantification: project each 3D cell mask to its XY
      footprint, compute per-cell metrics on a fresh MIP using the SAME cells
      as the 3D run, save cells_2d_matched.csv per image
  (c) cellquant in 2D on each published MIP (independent cell identities)

Per-image timings and any failures are logged. Run with --image to do a single
image; otherwise iterates the whole dataset.

USAGE
    python 01_run_pipelines.py                       # both datasets, all images
    python 01_run_pipelines.py --dataset mammalian   # mammalian only
    python 01_run_pipelines.py --dataset yeast       # yeast only (needs sld reader)
    python 01_run_pipelines.py --dataset mammalian --image arsenite_rep1
    python 01_run_pipelines.py --skip-3d             # only do matched-2D + published-2D
    python 01_run_pipelines.py --skip-matched        # only 3D and published-2D
    python 01_run_pipelines.py --skip-published      # only 3D and matched-2D
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any

# Suppress duplicate libomp crashes on macOS (must happen before any C
# extension that links against OpenMP is imported)
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
import tifffile
import yaml

# Local paths
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
setup_paths = importlib.import_module("00_setup_paths")  # noqa: E402

# Load cellquant.py as a module so we can call its internal helpers
# (per_cell_metrics, detect_puncta, costes_threshold, compute_colocalization).
# This guarantees the matched-2D quantification uses byte-identical algorithms.
_spec = importlib.util.spec_from_file_location(
    "cellquant", setup_paths.CELLQUANT_SCRIPT)
cq = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cq)


# Resume support: when True, skip any pass whose out_dir already has a
# cells.csv (crash recovery / incremental runs over a long CPU batch).
RESUME = False


def _skip_if_done(out_dir: Path, label: str) -> dict[str, Any] | None:
    if RESUME and (out_dir / "cells.csv").exists():
        return {"label": label,
                "status": "skipped (resume: cells.csv exists)",
                "duration_sec": 0.0}
    return None


# ---------------------------------------------------------------------------
# Subprocess wrapper for cellquant invocations
# ---------------------------------------------------------------------------
def _run_cellquant(
    argv: list[str],
    out_dir: Path,
    label: str,
) -> dict[str, Any]:
    """Invoke cellquant.py as a subprocess. Returns timing + status."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(setup_paths.CELLQUANT_SCRIPT),
        *argv,
        "--out", str(out_dir),
    ]
    env = dict(os.environ)
    env.update(setup_paths.ENV_OVERRIDES)

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=setup_paths.PER_IMAGE_TIMEOUT_SEC,
            env=env,
        )
        dt = time.time() - t0
        status = "ok" if result.returncode == 0 else f"failed (exit {result.returncode})"
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.TimeoutExpired as exc:
        dt = time.time() - t0
        status = f"timeout after {dt:.0f}s"
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
    except Exception as exc:  # pragma: no cover
        dt = time.time() - t0
        status = f"crash: {type(exc).__name__}: {exc}"
        stdout = ""
        stderr = traceback.format_exc()

    record = {
        "label": label,
        "argv": cmd,
        "status": status,
        "duration_sec": dt,
        "stdout_tail": stdout[-2000:],
        "stderr_tail": stderr[-2000:],
    }
    with open(out_dir / "_run_record.json", "w") as fh:
        json.dump(record, fh, indent=2)
    return record


def _log_failure(label: str, record: dict[str, Any]) -> None:
    """Append a failure summary to the global failures.log."""
    with open(setup_paths.OUT_FAILURES, "a") as fh:
        fh.write(f"\n=== {label} : {record['status']} ===\n")
        fh.write(f"argv: {' '.join(record['argv'])}\n")
        fh.write("stderr (tail):\n")
        fh.write(record["stderr_tail"])
        fh.write("\n")


# ---------------------------------------------------------------------------
# 3D pass
# ---------------------------------------------------------------------------
def run_3d_mammalian(zstack_path: Path) -> dict[str, Any]:
    """Run cellquant 3D on one mammalian z-stack."""
    stem = zstack_path.stem
    out_dir = setup_paths.OUT_MAMM_3D / stem
    in_dir = setup_paths.VALIDATION_ROOT / "_staging_mammalian" / stem
    in_dir.mkdir(parents=True, exist_ok=True)
    # cellquant expects a directory of TIFFs. Symlink the single z-stack in.
    link_path = in_dir / zstack_path.name
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    link_path.symlink_to(zstack_path)
    argv = [
        str(in_dir),
        *setup_paths.MAMMALIAN_CHANNELS,
        "--cell-type", setup_paths.MAMMALIAN_CELL_TYPE,
        "--mode", "3d",
        "--voxel-size",
        str(setup_paths.MAMMALIAN_VOXEL_XY_UM),
        str(setup_paths.MAMMALIAN_VOXEL_Z_UM),
        "--condensate-index",
        "--colocalization",
        "--filename-pattern", setup_paths.MAMMALIAN_FILENAME_PATTERN_3D,
        "--skip-plots",
    ]
    if setup_paths.NO_GPU:
        argv.append("--no-gpu")
    label = f"mammalian/3d/{stem}"
    record = _run_cellquant(argv, out_dir, label)
    if record["status"] != "ok":
        _log_failure(label, record)
    return record


def run_published_2d_mammalian(mip_path: Path) -> dict[str, Any]:
    """Run cellquant 2D on one mammalian MIP (independent cell identities)."""
    stem = mip_path.stem.replace("MAX_", "", 1)
    out_dir = setup_paths.OUT_MAMM_PUBLISHED2D / stem
    in_dir = setup_paths.VALIDATION_ROOT / "_staging_published_mammalian" / stem
    in_dir.mkdir(parents=True, exist_ok=True)
    link_path = in_dir / mip_path.name
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    link_path.symlink_to(mip_path)
    argv = [
        str(in_dir),
        *setup_paths.MAMMALIAN_CHANNELS,
        "--cell-type", setup_paths.MAMMALIAN_CELL_TYPE,
        "--mode", "2d",
        "--voxel-size", str(setup_paths.MAMMALIAN_VOXEL_XY_UM),
        "--condensate-index",
        "--colocalization",
        "--allow-2d-colocalization",
        "--filename-pattern", setup_paths.MAMMALIAN_FILENAME_PATTERN_2D,
        "--skip-plots",
    ]
    if setup_paths.NO_GPU:
        argv.append("--no-gpu")
    label = f"mammalian/published_2d/{stem}"
    record = _run_cellquant(argv, out_dir, label)
    if record["status"] != "ok":
        _log_failure(label, record)
    return record


# ---------------------------------------------------------------------------
# Matched-2D-from-3D quantification
# ---------------------------------------------------------------------------
def _load_zstack_channels(path: Path, channels: list[dict]) -> dict[str, np.ndarray]:
    """Load a 3-channel z-stack as {channel_name: (Z, Y, X) float32}."""
    arr = np.asarray(tifffile.imread(str(path)))
    # SG_zstacks shape is (Z, C, Y, X)
    if arr.ndim != 4:
        raise ValueError(f"{path.name}: expected 4D ZCYX, got shape {arr.shape}")
    # Detect axis order by shape: smallest non-(Y,X) dim is C, the other is Z
    if arr.shape[1] < arr.shape[0] and arr.shape[1] <= 10:
        czyx = np.moveaxis(arr, 1, 0)  # (Z, C, Y, X) -> (C, Z, Y, X)
    elif arr.shape[0] <= arr.shape[1] and arr.shape[0] <= 10:
        czyx = arr  # already (C, Z, Y, X)
    else:
        raise ValueError(f"{path.name}: ambiguous axes for shape {arr.shape}")
    out = {}
    for ch in channels:
        idx = ch["position"] - 1
        out[ch["name"]] = czyx[idx].astype(np.float32, copy=False)
    return out


def _config_from_3d_run(run_3d_dir: Path) -> dict[str, Any]:
    """Load the config_used.yml from a finished 3D run, adapt for 2D usage.

    The matched-2D quantification uses the same puncta / coloc parameters
    as the 3D run, just operating on 2D arrays. We need to overlay the few
    keys that differ between modes (puncta_min_area_px vs *_volume_vox,
    nucleolar_min_area_px vs *_volume_vox).
    """
    cfg_path = run_3d_dir / "config_used.yml"
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    # Strip 3D-only keys; the puncta_min_area_px / max_area_px from the
    # mammalian preset are still in cfg because cellquant overlays both
    # before mode dispatch. We rely on those for the 2D path.
    cfg["mode"] = "2d"
    return cfg


def _channels_from_3d_run(run_3d_dir: Path) -> list[dict]:
    """Load channel definitions back from the config_used.yml."""
    with open(run_3d_dir / "config_used.yml") as fh:
        cfg = yaml.safe_load(fh)
    chans = []
    for s in cfg.get("_channels", []):
        # "1:DAPI:nucleus" -> {"position":1, "name":"DAPI", "role":"nucleus"}
        p, n, r = s.split(":")
        chans.append({"position": int(p), "name": n, "role": r})
    return chans


def matched_2d_mammalian(zstack_path: Path) -> dict[str, Any]:
    """Project the 3D cell mask, MIP the z-stack, compute per-cell 2D metrics."""
    stem = zstack_path.stem
    run_3d_dir = setup_paths.OUT_MAMM_3D / stem
    out_dir = setup_paths.OUT_MAMM_MATCHED2D / stem
    out_dir.mkdir(parents=True, exist_ok=True)
    cells_3d_csv = run_3d_dir / "cells.csv"
    if not cells_3d_csv.exists():
        record = {"label": f"mammalian/matched_2d/{stem}",
                  "status": f"3D cells.csv not found at {cells_3d_csv}",
                  "duration_sec": 0.0}
        _log_failure(record["label"], {**record, "argv": [], "stderr_tail": "",
                                       "stdout_tail": ""})
        return record

    t0 = time.time()
    try:
        cfg = _config_from_3d_run(run_3d_dir)
        channels = _channels_from_3d_run(run_3d_dir)

        # Load the cell + nuc masks from the 3D run
        mask_dir = run_3d_dir / "masks"
        cell_mask_3d = np.asarray(
            tifffile.imread(str(mask_dir / f"{stem}_cellmask.tif"))).astype(np.int32)
        nuc_mask_path = mask_dir / f"{stem}_nucmask.tif"
        if nuc_mask_path.exists():
            nuc_mask_3d = np.asarray(
                tifffile.imread(str(nuc_mask_path))).astype(np.int32)
        else:
            nuc_mask_3d = np.zeros_like(cell_mask_3d, dtype=np.int32)

        # Project to 2D footprints. For per-cell ID continuity we use the
        # maximum along Z so that any voxel labelled c stamps the XY footprint
        # with c. This may merge IDs in unlucky overlapping cases — for now
        # assume cells are reasonably XY-separated (true for the HCT116 data).
        cell_mask_2d = cell_mask_3d.max(axis=0).astype(np.int32, copy=False)
        nuc_mask_2d = nuc_mask_3d.max(axis=0).astype(np.int32, copy=False)

        # Load z-stack and MIP each channel
        zstack = _load_zstack_channels(zstack_path, channels)
        images_2d = {name: ch.max(axis=0) for name, ch in zstack.items()}

        # Apply cellquant's 2D per-cell metric machinery on the projected
        # masks + MIP'd channels. This reuses identical algorithms.
        compartment = cfg.get("puncta_compartment", "cytosol")
        compartment_mask_2d = cq.make_compartment_mask(
            cell_mask_2d, nuc_mask_2d, compartment, cfg)

        puncta_masks_2d: dict[str, np.ndarray] = {}
        puncta_chs = [ch for ch in channels
                      if ch["role"] == "quantify"]
        for ch in puncta_chs:
            name = ch["name"]
            ch_cfg = cq.resolve_per_channel_cfg(cfg, name)
            puncta_masks_2d[name] = cq.detect_puncta(
                images_2d[name], compartment_mask_2d, ch_cfg)

        # Fragmentation indices (mirrors cellquant main loop)
        frag = {ch["name"]: cq.compute_fragmentation_indices(
            images_2d[ch["name"]], cell_mask_2d,
            cq.resolve_per_channel_cfg(cfg, ch["name"])) for ch in puncta_chs}

        cell_to_nucs = cq.map_nuclei_to_cells(nuc_mask_2d, cell_mask_2d)
        meta = {"image": zstack_path.name, "condition": "", "replicate": ""}
        # Try to parse condition+replicate from stem like "control_rep1"
        if "_rep" in stem:
            cond, rep = stem.rsplit("_rep", 1)
            meta["condition"] = cond
            meta["replicate"] = rep

        cells_df = cq.per_cell_metrics(
            zstack_path.name, images_2d, cell_mask_2d, nuc_mask_2d,
            puncta_masks_2d, cell_to_nucs, meta, channels, cfg,
            fragmentation_per_channel=frag,
        )

        # Colocalization (matches the 3D run's coloc_channels = quantify + nucleolus)
        if cfg.get("colocalization", False):
            quant_chs = cq.get_quantify_channels(channels)
            nol_chs = cq.get_nucleolus_channels(channels)
            coloc_channels = []
            seen = set()
            for ch in quant_chs + nol_chs:
                if ch["name"] not in seen:
                    coloc_channels.append(ch)
                    seen.add(ch["name"])
            coloc_df = cq.compute_colocalization(
                images_2d, coloc_channels, cell_mask_2d, nuc_mask_2d,
                cfg.get("colocalization_compartment", "whole-cell"),
                cfg, meta, zstack_path.name)
            # Pivot into cells_df mirroring cellquant's main-loop logic
            if not coloc_df.empty:
                for pair_name in coloc_df["pair"].unique():
                    pair_df = coloc_df[coloc_df["pair"] == pair_name]
                    for metric_col in ["pearson_r", "manders_m1", "manders_m2"]:
                        pivot_col = f"{metric_col}_{pair_name}"
                        merge_df = pair_df[["image", "cell_id", metric_col]].rename(
                            columns={metric_col: pivot_col})
                        cells_df = cells_df.merge(
                            merge_df, on=["image", "cell_id"], how="left")
            coloc_df.to_csv(out_dir / "colocalization.csv", index=False)

        cells_df.to_csv(out_dir / "cells.csv", index=False)
        dt = time.time() - t0
        return {"label": f"mammalian/matched_2d/{stem}",
                "status": "ok",
                "duration_sec": dt,
                "n_cells": int(len(cells_df))}
    except Exception as exc:
        dt = time.time() - t0
        record = {"label": f"mammalian/matched_2d/{stem}",
                  "status": f"crash: {type(exc).__name__}: {exc}",
                  "duration_sec": dt}
        _log_failure(record["label"],
                     {**record, "argv": [], "stderr_tail": traceback.format_exc(),
                      "stdout_tail": ""})
        return record


# ---------------------------------------------------------------------------
# Yeast 3D / published-2D / matched-2D
# ---------------------------------------------------------------------------
def run_3d_yeast(zstack_path: Path) -> dict[str, Any]:
    """Run cellquant 3D on one yeast z-stack (Nsr1=nucleolus, no nucleus channel)."""
    stem = zstack_path.stem
    out_dir = setup_paths.OUT_YEAST_3D / stem
    _sk = _skip_if_done(out_dir, f"yeast/3d/{stem}")
    if _sk is not None:
        return _sk
    in_dir = setup_paths.VALIDATION_ROOT / "_staging_yeast" / stem
    in_dir.mkdir(parents=True, exist_ok=True)
    link_path = in_dir / zstack_path.name
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    link_path.symlink_to(zstack_path)
    argv = [
        str(in_dir),
        *setup_paths.YEAST_CHANNELS,
        "--cell-type", setup_paths.YEAST_CELL_TYPE,
        "--mode", "3d",
        "--seg-3d-method", setup_paths.YEAST_SEG_3D_METHOD,
        "--z-crop-center", str(setup_paths.YEAST_Z_CROP_CENTER),
        "--seg-downsample", str(setup_paths.YEAST_SEG_DOWNSAMPLE),
        "--voxel-size",
        str(setup_paths.YEAST_VOXEL_XY_UM),
        str(setup_paths.YEAST_VOXEL_Z_UM),
        "--colocalization",
        "--nucleolar-proximity", "Nsr1",
        "--condensate-index",
        "--filename-pattern", setup_paths.YEAST_FILENAME_PATTERN_3D,
        "--skip-plots",
    ]
    if setup_paths.NO_GPU:
        argv.append("--no-gpu")
    label = f"yeast/3d/{stem}"
    record = _run_cellquant(argv, out_dir, label)
    if record["status"] != "ok":
        _log_failure(label, record)
    return record


def run_published_2d_yeast(mip_path: Path) -> dict[str, Any]:
    """Run cellquant 2D on one published yeast MIP (independent cell identities)."""
    stem = mip_path.stem  # MIPs share the z-stack naming for yeast
    out_dir = setup_paths.OUT_YEAST_PUBLISHED2D / stem
    _sk = _skip_if_done(out_dir, f"yeast/published_2d/{stem}")
    if _sk is not None:
        return _sk
    in_dir = setup_paths.VALIDATION_ROOT / "_staging_published_yeast" / stem
    in_dir.mkdir(parents=True, exist_ok=True)
    link_path = in_dir / mip_path.name
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    link_path.symlink_to(mip_path)
    argv = [
        str(in_dir),
        *setup_paths.YEAST_CHANNELS,
        "--cell-type", setup_paths.YEAST_CELL_TYPE,
        "--mode", "2d",
        "--voxel-size", str(setup_paths.YEAST_VOXEL_XY_UM),
        "--colocalization",
        # This MIP-coloc pass exists precisely to quantify the projection
        # artifact vs 3D (compare script Q2); the Phase C guardrail otherwise
        # aborts 2D colocalization, so opt in explicitly here.
        "--allow-2d-colocalization",
        "--nucleolar-proximity", "Nsr1",
        "--condensate-index",
        "--filename-pattern", setup_paths.YEAST_FILENAME_PATTERN_2D,
        "--skip-plots",
    ]
    if setup_paths.NO_GPU:
        argv.append("--no-gpu")
    label = f"yeast/published_2d/{stem}"
    record = _run_cellquant(argv, out_dir, label)
    if record["status"] != "ok":
        _log_failure(label, record)
    return record


def matched_2d_yeast(zstack_path: Path) -> dict[str, Any]:
    """Project yeast 3D cell + nucleolar masks, MIP the z-stack, redo 2D quant.

    Differs from mammalian: no nucleus channel → nuc_mask is zeros; nucleolar
    proximity is computed on the projected nucleolar mask; whole-cell
    compartment (no nuclear erosion).
    """
    stem = zstack_path.stem
    run_3d_dir = setup_paths.OUT_YEAST_3D / stem
    out_dir = setup_paths.OUT_YEAST_MATCHED2D / stem
    out_dir.mkdir(parents=True, exist_ok=True)
    _sk = _skip_if_done(out_dir, f"yeast/matched_2d/{stem}")
    if _sk is not None:
        return _sk
    cells_3d_csv = run_3d_dir / "cells.csv"
    if not cells_3d_csv.exists():
        record = {"label": f"yeast/matched_2d/{stem}",
                  "status": f"3D cells.csv not found at {cells_3d_csv}",
                  "duration_sec": 0.0}
        _log_failure(record["label"], {**record, "argv": [], "stderr_tail": "",
                                       "stdout_tail": ""})
        return record

    t0 = time.time()
    try:
        cfg = _config_from_3d_run(run_3d_dir)
        channels = _channels_from_3d_run(run_3d_dir)

        mask_dir = run_3d_dir / "masks"
        cell_mask_3d = np.asarray(
            tifffile.imread(str(mask_dir / f"{stem}_cellmask.tif"))).astype(np.int32)
        # No nucleus channel for yeast
        nuc_mask_2d = np.zeros(cell_mask_3d.shape[1:], dtype=np.int32)

        nucleolar_path = mask_dir / f"{stem}_nucleolarmask.tif"
        if nucleolar_path.exists():
            nucleolar_mask_3d = np.asarray(
                tifffile.imread(str(nucleolar_path))).astype(np.uint8)
        else:
            nucleolar_mask_3d = np.zeros_like(cell_mask_3d, dtype=np.uint8)

        # Project Z → XY footprints
        cell_mask_2d = cell_mask_3d.max(axis=0).astype(np.int32, copy=False)
        nucleolar_mask_2d = (nucleolar_mask_3d.max(axis=0) > 0).astype(np.uint8)

        # Load z-stack and MIP each channel
        zstack = _load_zstack_channels(zstack_path, channels)
        images_2d = {name: ch.max(axis=0) for name, ch in zstack.items()}

        # Compartment for puncta detection (yeast preset → whole-cell)
        compartment = cfg.get("puncta_compartment", "whole-cell")
        compartment_mask_2d = cq.make_compartment_mask(
            cell_mask_2d, nuc_mask_2d, compartment, cfg)

        puncta_masks_2d: dict[str, np.ndarray] = {}
        puncta_chs = [ch for ch in channels if ch["role"] == "quantify"]
        for ch in puncta_chs:
            name = ch["name"]
            ch_cfg = cq.resolve_per_channel_cfg(cfg, name)
            # Per-channel compartment override (unlikely for yeast, mirror mammalian)
            if ch_cfg.get("puncta_compartment") != compartment:
                ch_comp = cq.make_compartment_mask(
                    cell_mask_2d, nuc_mask_2d,
                    ch_cfg.get("puncta_compartment", compartment), cfg)
            else:
                ch_comp = compartment_mask_2d
            puncta_masks_2d[name] = cq.detect_puncta(
                images_2d[name], ch_comp, ch_cfg)

        frag = {ch["name"]: cq.compute_fragmentation_indices(
            images_2d[ch["name"]], cell_mask_2d,
            cq.resolve_per_channel_cfg(cfg, ch["name"])) for ch in puncta_chs}

        cell_to_nucs = cq.map_nuclei_to_cells(nuc_mask_2d, cell_mask_2d)
        # cellquant's parser handles "{condition}_series1_rep{replicate}" cleanly
        meta = cq.parse_filename_metadata(stem, cfg)
        meta["image"] = zstack_path.name

        cells_df = cq.per_cell_metrics(
            zstack_path.name, images_2d, cell_mask_2d, nuc_mask_2d,
            puncta_masks_2d, cell_to_nucs, meta, channels, cfg,
            fragmentation_per_channel=frag,
        )

        # Colocalization (matches the 3D run's quantify + nucleolus pairs)
        if cfg.get("colocalization", False):
            quant_chs = cq.get_quantify_channels(channels)
            nol_chs = cq.get_nucleolus_channels(channels)
            coloc_channels = []
            seen = set()
            for ch in quant_chs + nol_chs:
                if ch["name"] not in seen:
                    coloc_channels.append(ch)
                    seen.add(ch["name"])
            coloc_df = cq.compute_colocalization(
                images_2d, coloc_channels, cell_mask_2d, nuc_mask_2d,
                cfg.get("colocalization_compartment", "whole-cell"),
                cfg, meta, zstack_path.name)
            if not coloc_df.empty:
                for pair_name in coloc_df["pair"].unique():
                    pair_df = coloc_df[coloc_df["pair"] == pair_name]
                    for metric_col in ["pearson_r", "manders_m1", "manders_m2"]:
                        pivot_col = f"{metric_col}_{pair_name}"
                        merge_df = pair_df[["image", "cell_id", metric_col]].rename(
                            columns={metric_col: pivot_col})
                        cells_df = cells_df.merge(
                            merge_df, on=["image", "cell_id"], how="left")
            coloc_df.to_csv(out_dir / "colocalization.csv", index=False)

        # Nucleolar proximity (2D path uses pixel threshold; the function
        # picks the path from puncta_mask.ndim, not cfg["mode"])
        nuc_prox_ch_name = cfg.get("nucleolar_proximity_channel")
        if nuc_prox_ch_name and nucleolar_mask_2d.max() > 0:
            prox_threshold_px = float(cfg.get("proximity_threshold_px", 5))
            prox_dfs = []
            for pch in puncta_chs:
                prox_df = cq.compute_nucleolar_proximity(
                    puncta_masks_2d[pch["name"]], nucleolar_mask_2d, cell_mask_2d,
                    prox_threshold_px, meta, zstack_path.name, pch["name"],
                    cfg=cfg)
                prox_dfs.append(prox_df)
            prox_out = pd.concat(prox_dfs, ignore_index=True) if prox_dfs else pd.DataFrame()
            if not prox_out.empty:
                prox_out.to_csv(out_dir / "nucleolar_proximity.csv", index=False)
                for ch_name in prox_out["channel"].unique():
                    ch_df = prox_out[prox_out["channel"] == ch_name]
                    for metric_col in ["mean_distance", "fraction_proximal"]:
                        pivot_col = f"{ch_name}_{metric_col}"
                        merge_df = ch_df[["image", "cell_id", metric_col]].rename(
                            columns={metric_col: pivot_col})
                        cells_df = cells_df.merge(
                            merge_df, on=["image", "cell_id"], how="left")

        cells_df.to_csv(out_dir / "cells.csv", index=False)
        dt = time.time() - t0
        return {"label": f"yeast/matched_2d/{stem}",
                "status": "ok",
                "duration_sec": dt,
                "n_cells": int(len(cells_df))}
    except Exception as exc:
        dt = time.time() - t0
        record = {"label": f"yeast/matched_2d/{stem}",
                  "status": f"crash: {type(exc).__name__}: {exc}",
                  "duration_sec": dt}
        _log_failure(record["label"],
                     {**record, "argv": [], "stderr_tail": traceback.format_exc(),
                      "stdout_tail": ""})
        return record


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run_mammalian(image_filter: str | None,
                  skip_3d: bool, skip_matched: bool, skip_published: bool
                  ) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    zstacks = setup_paths.list_mammalian_zstacks()
    mips = setup_paths.list_mammalian_mips()
    if image_filter:
        zstacks = [p for p in zstacks if image_filter in p.stem]
        mips = [p for p in mips if image_filter in p.stem]
        if not zstacks and not mips:
            raise SystemExit(
                f"--image {image_filter!r} matched none of the mammalian "
                f"inputs (0 z-stacks, 0 MIPs); refusing to run 0 images.")

    # 3D pass
    if not skip_3d:
        print(f"\n=== Mammalian 3D pass ({len(zstacks)} images) ===")
        for i, p in enumerate(zstacks, 1):
            print(f"  [{i}/{len(zstacks)}] 3D: {p.name}")
            rec = run_3d_mammalian(p)
            print(f"    -> {rec['status']} ({rec['duration_sec']:.1f}s)")
            records.append(rec)

    # Matched-2D pass (depends on 3D having completed)
    if not skip_matched:
        print(f"\n=== Mammalian matched-2D pass ({len(zstacks)} images) ===")
        for i, p in enumerate(zstacks, 1):
            print(f"  [{i}/{len(zstacks)}] matched-2D: {p.name}")
            rec = matched_2d_mammalian(p)
            print(f"    -> {rec['status']} ({rec['duration_sec']:.1f}s)")
            records.append(rec)

    # Published-2D pass on every MIP
    if not skip_published:
        print(f"\n=== Mammalian published-2D pass ({len(mips)} MIPs) ===")
        for i, p in enumerate(mips, 1):
            print(f"  [{i}/{len(mips)}] published-2D: {p.name}")
            rec = run_published_2d_mammalian(p)
            print(f"    -> {rec['status']} ({rec['duration_sec']:.1f}s)")
            records.append(rec)

    return records


def run_yeast(image_filter: str | None,
              skip_3d: bool, skip_matched: bool, skip_published: bool,
              subset: str = "full",
              ) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    if subset == "one-per-temp":
        zstacks = setup_paths.list_yeast_zstacks_one_per_temp()
        mip_stems = {p.stem for p in zstacks}
        mips = [p for p in setup_paths.list_yeast_mips() if p.stem in mip_stems]
    else:
        # Validation set: the 30 reps that have both a z-stack and a published MIP.
        zstacks = setup_paths.list_yeast_zstacks_subset()
        mips = setup_paths.list_yeast_mips()
    if image_filter:
        zstacks = [p for p in zstacks if image_filter in p.stem]
        mips = [p for p in mips if image_filter in p.stem]
        if not zstacks and not mips:
            raise SystemExit(
                f"--image {image_filter!r} matched none of the yeast "
                f"inputs (0 z-stacks, 0 MIPs); refusing to run 0 images.")

    if not skip_3d:
        print(f"\n=== Yeast 3D pass ({len(zstacks)} images) ===")
        for i, p in enumerate(zstacks, 1):
            print(f"  [{i}/{len(zstacks)}] 3D: {p.name}")
            rec = run_3d_yeast(p)
            print(f"    -> {rec['status']} ({rec['duration_sec']:.1f}s)")
            records.append(rec)

    if not skip_matched:
        print(f"\n=== Yeast matched-2D pass ({len(zstacks)} images) ===")
        for i, p in enumerate(zstacks, 1):
            print(f"  [{i}/{len(zstacks)}] matched-2D: {p.name}")
            rec = matched_2d_yeast(p)
            print(f"    -> {rec['status']} ({rec['duration_sec']:.1f}s)")
            records.append(rec)

    if not skip_published:
        print(f"\n=== Yeast published-2D pass ({len(mips)} MIPs) ===")
        for i, p in enumerate(mips, 1):
            print(f"  [{i}/{len(mips)}] published-2D: {p.name}")
            rec = run_published_2d_yeast(p)
            print(f"    -> {rec['status']} ({rec['duration_sec']:.1f}s)")
            records.append(rec)

    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["mammalian", "yeast", "both"],
                    default="both")
    ap.add_argument("--image", default=None,
                    help="Filter to images containing this substring")
    ap.add_argument("--skip-3d", action="store_true")
    ap.add_argument("--skip-matched", action="store_true")
    ap.add_argument("--skip-published", action="store_true")
    ap.add_argument("--yeast-subset", choices=["full", "one-per-temp"],
                    default="full",
                    help="full=30 manuscript reps; one-per-temp=5 (rep1/temp)")
    ap.add_argument("--resume", action="store_true",
                    help="Skip any pass whose out_dir already has cells.csv")
    args = ap.parse_args()

    global RESUME
    RESUME = args.resume

    # Truncate failures log at the start of each invocation
    if setup_paths.OUT_FAILURES.exists():
        setup_paths.OUT_FAILURES.unlink()

    all_records: list[dict[str, Any]] = []
    if args.dataset in ("mammalian", "both"):
        all_records.extend(run_mammalian(
            args.image, args.skip_3d, args.skip_matched, args.skip_published))
    if args.dataset in ("yeast", "both"):
        all_records.extend(run_yeast(
            args.image, args.skip_3d, args.skip_matched, args.skip_published,
            subset=args.yeast_subset))

    # A run that produced no records processed zero images. That is never a
    # success — it means every enabled pass had an empty work list (missing
    # data dir, over-narrow --image, or everything skipped). Fail loudly rather
    # than writing a summary that reports 0/0 ok and exiting 0.
    if not all_records:
        raise SystemExit(
            "No pipeline steps ran: every enabled pass had an empty work list. "
            "Refusing to report success on 0 images. Check the input "
            "directories in 00_setup_paths.py and any --image / --skip flags.")

    # Save combined timing/status report
    setup_paths.OUT_COMBINED.mkdir(parents=True, exist_ok=True)
    summary_path = setup_paths.OUT_COMBINED / "run_summary.csv"
    df = pd.DataFrame(all_records)
    df.to_csv(summary_path, index=False)
    print(f"\nWrote run summary to {summary_path}")
    n_ok = (df["status"] == "ok").sum()
    print(f"  {n_ok}/{len(df)} steps ok")
    if (df["status"] != "ok").any():
        print(f"  failures logged to {setup_paths.OUT_FAILURES}")


if __name__ == "__main__":
    main()
