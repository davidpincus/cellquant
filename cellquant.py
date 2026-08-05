#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
#   "numpy>=1.24,<2.0",
#   "scipy",
#   "scikit-image>=0.24",
#   "pandas",
#   "matplotlib",
#   "pyyaml",
#   "tifffile",
#   "cellpose>=4.0",
#   "opencv-python-headless<4.10",
#   "torch",
# ]
# ///
"""
cellquant.py — unified 2D + 3D fluorescence image analysis pipeline.

Single-script multi-channel cell quantification: segmentation, puncta detection,
colocalization, spatial proximity, morphometrics, and replicate-level statistics.
Auto-detects whether inputs are 2D MIPs or 3D Z-stacks; same CLI, same outputs.

2D mode (mammalian/yeast/bacteria MIPs) — paper-validated:
  python cellquant.py images/ \
      "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
      --out results/ --cell-type mammalian \
      --filename-pattern "{condition}_rep{replicate}"

3D mode (Z-stacks; auto-detected from input shape):
  python cellquant.py stacks/ \
      "1:Sis1:quantify" "2:Hsp70:quantify" \
      --out results3d/ --cell-type yeast \
      --voxel-size 0.1 0.1 \
      --filename-pattern "{condition}_rep{replicate}"

Mode dispatch:
  --mode auto   default; 2D for ndim==3 inputs (CYX), 3D for ndim==4 (ZCYX/CZYX)
  --mode 2d     force; fails on 3D input
  --mode 3d     force; requires --voxel-size

Per-cell metrics common to both modes (see cells.csv):
  cell_id, n_nuclei, condition, replicate
  {ch}_cell_mean, {ch}_nucleus_mean, {ch}_cytosol_mean (per non-skip channel)
  {ch}_condensate_index_cell, {ch}_condensate_index_cytosol (p95/mean ratio)
  {ch}_puncta_n, {ch}_puncta_mean_intensity, {ch}_diffuse_mean_intensity
  {ch}_frac_intensity_in_puncta, {ch}_puncta_integrated_intensity
  {ch}_fragmentation_index_simple        connected components above LoG-Otsu
  {ch}_fragmentation_index_persistence   threshold-free swept integral [0,1]
  pearson_r_{a}_vs_{b}, manders_m1, manders_m2 (with --colocalization)
  {ch}_mean_distance, {ch}_fraction_proximal (with --nucleolar-proximity)

2D-only: cell_area, *_puncta_area, nucleolar_area/circularity/eccentricity/solidity
3D-only: cell_volume_vox, cell_volume_um3, *_puncta_volume_um3,
         nucleolar_volume_um3, nucleolar_eq_diameter_um
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Any

__version__ = "1.1.0"

# Suppress duplicate-libomp crash on macOS (conda/pip can install conflicting
# copies via MKL, PyTorch, and system llvm-openmp).  Must be set before any
# C-extension that links OpenMP is imported.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

# ---------------------------------------------------------------------------
# Numpy version gate — must run BEFORE importing cellpose, scipy, skimage, or
# matplotlib, because those pull in C extensions (cv2, torch) that crash
# immediately under numpy 2.x instead of producing a useful error message.
# ---------------------------------------------------------------------------
import numpy as np

_np_ver = tuple(int(x) for x in np.__version__.split(".")[:2])
if _np_ver >= (2, 0):
    print(f"ERROR: numpy {np.__version__} detected — cellquant requires numpy <2.0")
    print('Fix: pip install "numpy>=1.24,<2.0" "opencv-python-headless<4.10"')
    print("     If you're on a Mac and using zsh, the quotes are required.")
    sys.exit(1)

from cellpose import models

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pandas as pd
import tifffile as tiff
import yaml

from scipy.ndimage import distance_transform_edt, gaussian_laplace
from scipy.stats import mannwhitneyu
from skimage import filters, measure, morphology, segmentation
from skimage.transform import resize

try:
    import torch
except Exception:
    torch = None

# Multiformat reader — optional, used for ND2/CZI/LIF inputs. We prefer
# bioio (the actively-maintained modular successor); aicsimageio is kept as
# a fallback for legacy environments. Falls back to tifffile-only when
# neither is installed; we surface a clear error if a non-tiff input is
# requested without a backend.
try:
    from bioio import BioImage as _MultiformatReader  # type: ignore
    _MULTIFORMAT_BACKEND = "bioio"
    HAS_MULTIFORMAT = True
except Exception:
    try:
        from aicsimageio import AICSImage as _MultiformatReader  # type: ignore
        _MULTIFORMAT_BACKEND = "aicsimageio"
        HAS_MULTIFORMAT = True
    except Exception:
        _MultiformatReader = None  # type: ignore
        _MULTIFORMAT_BACKEND = None
        HAS_MULTIFORMAT = False

# Backwards-compat alias — external code (tests, sub-tools) may check this.
HAS_AICSIMAGEIO = HAS_MULTIFORMAT

# Non-tiff formats handled exclusively through the multiformat backend.
_AICS_EXTS = {".nd2", ".czi", ".lif", ".lsm"}
_TIFF_EXTS = {".tif", ".tiff", ".ome.tif", ".ome.tiff"}

# ---------------------------------------------------------------------------
# Startup dependency check — scikit-image version
# ---------------------------------------------------------------------------
_ski_ver = tuple(int(x) for x in __import__("skimage").__version__.split(".")[:2])
if _ski_ver < (0, 24):
    print(f"ERROR: scikit-image {__import__('skimage').__version__} detected — cellquant requires >=0.24")
    print('Fix: pip install "scikit-image>=0.24"')
    sys.exit(1)


# ---------------------------------------------------------------------------
# Cell-type presets
# ---------------------------------------------------------------------------
CELL_TYPE_PRESETS: dict[str, dict[str, Any]] = {
    "mammalian": {
        "pretrained_model": "cpsam",
        "nuclei_diameter": None,
        "cell_diameter": None,
        "flow_threshold": 0.4,
        "cellprob_threshold": 0.0,
        "nucleus_dilate_px": 3,
        "log_sigma": 2.0,
        "puncta_min_area_px": 6,
        "puncta_max_area_px": 10_000,
        # Physical-unit erosion of the puncta compartment before LoG. Removes
        # the membrane gradient (high-contrast cytoplasm/background edge)
        # that Otsu otherwise picks up as outline-shaped "puncta." 0.5 µm
        # is the membrane-gradient scale for typical mammalian cells; tune
        # via --puncta-compartment-erode-um.
        "puncta_compartment_erode_um": 0.5,
        "keep_min_nuclei": 1,
        "keep_max_nuclei": 4,
        # 3D-specific overrides — applied only when mode == "3d"
        "_3d": {
            "log_sigma": 2.0,
            "puncta_min_volume_vox": 8,
            "puncta_max_volume_vox": 200_000,
            # Minimum cell volume to keep after 3D segmentation. 30_000 vox
            # ≈ 25 µm³ at typical mammalian voxel size (0.094 µm XY,
            # 0.25 µm Z) — well below the p5 of real-cell volumes (~53k
            # vox) observed across HCT116 z-stacks, but well above the
            # debris/fragment band (median ~2k vox) that the per-Z
            # Cellpose stitch path produces in tight cell clusters. See
            # validation_3d/diagnostics/seg_experiments/ for the
            # histogram-derived justification.
            #
            # Note on tight clusters: for fields with HCT116-style
            # confluent clusters (cells touching with no gaps), the
            # composite cell-seg input occasionally over-fragments. The
            # workaround is `--cell-seg-channel <cytoplasmic_marker>`
            # (e.g. PABPC1); this is flagged for a future preset review
            # rather than changed here because it requires the user to
            # have a named cytoplasmic channel and may not generalize
            # across all mammalian acquisitions.
            "min_cell_volume_vox": 30_000,
            "max_cell_volume_vox": 0,
            "seg_3d_method": "stitch",
            "stitch_threshold": 0.4,
            "nucleolar_opening_radius": 1,
            "nucleolar_min_volume_vox": 10,
        },
    },
    "yeast": {
        "pretrained_model": "cpsam",
        "seg_downsample": 1,
        "cell_diameter": 40.0,
        "flow_threshold": 0.4,
        "cellprob_threshold": -1.0,
        "nucleus_dilate_px": 0,
        "min_cell_area": 200,
        "max_cell_area": 5000,
        "log_sigma": 1.5,
        "puncta_min_area_px": 3,
        "puncta_max_area_px": 300,
        "puncta_compartment": "whole-cell",
        "keep_min_nuclei": 0,
        "keep_max_nuclei": 0,
        "_3d": {
            "log_sigma": 1.5,
            # Yeast cell volume: ~30-100 µm³. At 0.1 µm voxels that's 30k-100k vox.
            "min_cell_volume_vox": 1500,
            "max_cell_volume_vox": 200_000,
            # 3D puncta: cap matches the paper's 2D 300-px area scaled to volume
            # (sphere of radius ~0.9 µm → ~3000 vox at 0.1 µm voxels).
            "puncta_min_volume_vox": 4,
            "puncta_max_volume_vox": 3000,
            # Default to slice-stitching with a strict IoU threshold (0.65).
            # Calibration on 25C_series1_rep1 (native voxel 0.10571/0.23 µm)
            # showed full do_3D does NOT round cells (the residual axial
            # elongation is optical PSF, not a segmentation artifact) yet costs
            # ~2.3x more wall-clock; a stricter stitch threshold links less
            # aggressively across Z and recovers more cells at stitch's cost.
            "seg_3d_method": "stitch",
            "stitch_threshold": 0.65,
            "nucleolar_opening_radius": 1,
            "nucleolar_min_volume_vox": 10,
        },
    },
    "bacteria": {
        "pretrained_model": "cpsam",
        # No nuclei_diameter — bacteria have no nucleus channel;
        # nuclear segmentation is skipped when no channel has role "nucleus".
        "cell_diameter": 15.0,
        "flow_threshold": 0.4,
        "cellprob_threshold": -2.0,
        "nucleus_dilate_px": 0,
        "log_sigma": 1.0,
        "puncta_min_area_px": 2,
        "puncta_max_area_px": 200,
        "keep_min_nuclei": 0,
        "keep_max_nuclei": 0,
        "_3d": {
            "log_sigma": 1.0,
            "puncta_min_volume_vox": 3,
            "puncta_max_volume_vox": 1_000,
            "min_cell_volume_vox": 100,
            "max_cell_volume_vox": 10_000,
            "seg_3d_method": "stitch",
            "stitch_threshold": 0.4,
            "nucleolar_opening_radius": 1,
            "nucleolar_min_volume_vox": 5,
        },
    },
}


# ---------------------------------------------------------------------------
# Hardcoded defaults
# ---------------------------------------------------------------------------
DEFAULTS: dict[str, Any] = {
    # I/O — .tif/.tiff handled by tifffile (paper-validated); .nd2/.czi/.lif
    # handled by aicsimageio when installed.
    "exts": [".tif", ".tiff", ".nd2", ".czi", ".lif"],
    "file_pattern": None,
    "save_masks": True,
    "reuse_masks": False,

    # Segmentation downsampling (smooths subcellular texture before Cellpose)
    "seg_downsample": 3,

    # Threading
    "cpu_threads": 4,

    # Display normalization for QC
    "p_low": 1.0,
    "p_high": 99.8,

    # Cellpose
    "use_gpu": True,
    "pretrained_model": "cpsam",

    # Cellpose thresholds (diameter None avoids the 30/diameter codepath)
    "nuclei_diameter": None,
    "cell_diameter": None,
    "flow_threshold": 0.4,
    "cellprob_threshold": 0.0,

    # Compartmenting
    "nucleus_dilate_px": 3,

    # Puncta detection
    "log_sigma": 2.0,
    "puncta_min_area_px": 6,
    "puncta_max_area_px": 10_000,
    "puncta_threshold_method": "otsu",  # "otsu" | "triangle" | "fixed"
    "puncta_threshold_fixed": 0.20,
    "puncta_filter_round": True,
    "puncta_min_circularity": 0.40,
    "puncta_min_solidity": 0.70,
    "puncta_compartment": "cytosol",  # a built-in or a --compartment NAME
    # User-defined regions from --compartment "NAME = TERM [op TERM]...".
    # Empty by default, so a run with no --compartment behaves exactly as before.
    "compartments": {},
    # Physical-unit inward erosion of the puncta compartment before LoG.
    # 0.0 disables; the mammalian preset overrides to 0.5 µm to suppress the
    # cell-membrane gradient that the LoG-Otsu detector otherwise picks up
    # as outline-shaped "puncta." Anisotropic-aware in 3D (uses XY and Z
    # voxel sizes). See --puncta-compartment-erode-um.
    "puncta_compartment_erode_um": 0.0,

    # QC rendering
    "qc_downsample": 4,
    "qc_dpi": 150,

    # Keep gate
    "keep_min_nuclei": 1,
    "keep_max_nuclei": 4,

    # Metadata parsing
    "parse_from_filename": True,
    "filename_pattern": "MAX_{condition}_rep{replicate}",
    "condition_map": {},
    "condition_order": [],
    "reference_condition": None,

    # Cell area filtering
    "min_cell_area": 0,
    "max_cell_area": 0,

    # Nucleolar segmentation
    "nucleolar_opening_radius": 1,
    "nucleolar_min_area_px": 5,

    # Colocalization
    "colocalization": False,
    "colocalization_compartment": "whole-cell",

    # Nucleolar proximity
    "nucleolar_proximity_channel": None,
    "proximity_threshold_px": 5,

    # Plotting
    "skip_plots": False,
    "trend": False,

    # Specialized metrics (opt-in via flags — off by default so default output
    # stays focused on the general-purpose metrics that apply across cell
    # types and experiments)
    "condensate_index": False,

    # ----- Mode and 3D-specific defaults (only consulted when mode == "3d") -----
    "mode": "auto",                 # "auto" | "2d" | "3d"
    "axes": None,                   # auto-detect; user can override e.g. "ZCYX"
    "voxel_size_xy_um": 1.0,        # microns; only used in 3D
    "voxel_size_z_um": 1.0,         # microns; only used in 3D
    "project_z": None,              # None | "max" | "sum" | "mean"
    "seg_3d_method": "stitch",      # "stitch" | "full"
    "stitch_threshold": 0.4,
    # 3D crop: keep N central slices around the per-image signal peak Z. 0 =
    # disabled (use the whole stack). When >0, cells whose 3D mask touches the
    # top or bottom Z plane of the cropped volume are also excluded.
    "z_crop_center": 0,
    "puncta_min_volume_vox": 8,
    "puncta_max_volume_vox": 100_000,
    "min_cell_volume_vox": 0,
    "max_cell_volume_vox": 0,
    "nucleolar_min_volume_vox": 10,
    "proximity_threshold_um": 0.5,  # 3D proximity threshold (microns, not pixels)

    # Fragmentation index (NEW — works in both modes; threshold-free metric
    # complementing puncta count for diffuse/condensate biology)
    "fragmentation_thresholds": 10,  # 0 disables persistence variant

    # Per-channel puncta detection overrides
    "puncta_rolling_ball": {},          # {channel_name: radius_px_or_vox}
    "puncta_params_per_channel": {},    # {channel_name: {param: value}}
}


# ---------------------------------------------------------------------------
# Channel parsing
# ---------------------------------------------------------------------------
def parse_channels(channel_strs: list[str]) -> list[dict[str, Any]]:
    """Parse "pos:Name:role" strings into channel dicts.

    Roles: nucleus, quantify, cell-boundary, skip.
    """
    channels = []
    for s in channel_strs:
        parts = s.split(":")
        if len(parts) != 3:
            raise ValueError(
                f"Channel spec '{s}' must be 'position:Name:role' "
                f"(e.g. '1:DAPI:nucleus')")
        pos, name, role = int(parts[0]), parts[1], parts[2].lower()
        valid_roles = {"nucleus", "quantify", "cell-boundary", "nucleolus", "skip"}
        if role not in valid_roles:
            raise ValueError(f"Channel role '{role}' not in {valid_roles}")
        channels.append({"position": pos, "name": name, "role": role})
    channels.sort(key=lambda c: c["position"])
    return channels


def get_nucleus_channel(channels: list[dict]) -> dict | None:
    """Return the channel with role 'nucleus', or None if absent."""
    for ch in channels:
        if ch["role"] == "nucleus":
            return ch
    return None


def get_cell_seg_channel(channels: list[dict], cfg: dict) -> dict:
    """Return the channel to use for cell segmentation.

    Priority: explicit --cell-seg-channel > cell-boundary role > first quantify.
    """
    explicit = cfg.get("cell_seg_channel")
    if explicit:
        for ch in channels:
            if ch["name"].lower() == explicit.lower():
                return ch
        raise ValueError(f"--cell-seg-channel '{explicit}' not found in channels")
    for ch in channels:
        if ch["role"] == "cell-boundary":
            return ch
    for ch in channels:
        if ch["role"] == "quantify":
            return ch
    raise ValueError("No channel available for cell segmentation")


def get_quantify_channels(channels: list[dict]) -> list[dict]:
    return [ch for ch in channels if ch["role"] == "quantify"]


def get_puncta_channels(
    channels: list[dict], puncta_names: list[str] | None
) -> list[dict]:
    if not puncta_names:
        return []
    name_map = {ch["name"].lower(): ch for ch in channels}
    result = []
    for name in puncta_names:
        if name.lower() not in name_map:
            raise ValueError(f"--puncta-channels '{name}' not found in channels")
        result.append(name_map[name.lower()])
    return result


def get_nucleolus_channels(channels: list[dict]) -> list[dict]:
    return [ch for ch in channels if ch["role"] == "nucleolus"]


# ---------------------------------------------------------------------------
# Config resolution (4-layer: defaults < cell-type preset < YAML < CLI)
# ---------------------------------------------------------------------------
def _parse_condition_map(pairs: list[str] | None) -> dict[str, str] | None:
    """Parse ["ctrl=control", "as=arsenite"] into a dict."""
    if not pairs:
        return None
    result: dict[str, str] = {}
    for pair in pairs:
        if "=" not in pair:
            raise ValueError(
                f"--condition-map entry '{pair}' must be key=value "
                f"(e.g. 'ctrl=control')")
        k, v = pair.split("=", 1)
        result[k.strip().lower()] = v.strip()
    return result


# Knobs that may legitimately be tuned per channel for puncta detection.
_PUNCTA_PARAM_NUMERIC_KEYS = {
    "log_sigma",
    "puncta_threshold_fixed",
    "puncta_min_circularity",
    "puncta_min_solidity",
}
_PUNCTA_PARAM_INT_KEYS = {
    "puncta_min_area_px",
    "puncta_max_area_px",
    "puncta_min_volume_vox",
    "puncta_max_volume_vox",
}
_PUNCTA_PARAM_STR_KEYS = {
    "puncta_threshold_method",   # otsu | triangle | fixed
    "puncta_compartment",        # cytosol | nucleus | whole-cell
}
_PUNCTA_PARAM_VALID = (
    _PUNCTA_PARAM_NUMERIC_KEYS | _PUNCTA_PARAM_INT_KEYS | _PUNCTA_PARAM_STR_KEYS
)


def _parse_puncta_rolling_ball(specs: list[str] | None) -> dict[str, int]:
    """Parse ["G3BP1:25", "PABPC1:0"] into {channel: radius}."""
    if not specs:
        return {}
    out: dict[str, int] = {}
    for spec in specs:
        if ":" not in spec:
            raise ValueError(
                f"--puncta-rolling-ball entry '{spec}' must be CH:RADIUS "
                f"(e.g. 'G3BP1:25')")
        ch, radius_s = spec.split(":", 1)
        ch = ch.strip()
        if not ch:
            raise ValueError(
                f"--puncta-rolling-ball entry '{spec}' has empty channel name")
        try:
            radius = int(radius_s)
        except ValueError:
            raise ValueError(
                f"--puncta-rolling-ball entry '{spec}' has non-integer radius "
                f"'{radius_s}'")
        if radius < 0:
            raise ValueError(
                f"--puncta-rolling-ball radius for {ch} must be >= 0 "
                f"(0 disables)")
        out[ch] = radius
    return out


def _parse_puncta_params_per_channel(
    specs: list[str] | None,
) -> dict[str, dict[str, Any]]:
    """Parse ["G3BP1:log_sigma=2.5,puncta_threshold_method=triangle", ...].

    Returns {channel_name: {key: typed_value, ...}}.
    """
    if not specs:
        return {}
    out: dict[str, dict[str, Any]] = {}
    for spec in specs:
        if ":" not in spec:
            raise ValueError(
                f"--puncta-params-per-channel entry '{spec}' must be CH:K=V,K=V"
                f" (e.g. 'G3BP1:log_sigma=2.5,puncta_threshold_method=triangle')")
        ch, kv_block = spec.split(":", 1)
        ch = ch.strip()
        if not ch:
            raise ValueError(
                f"--puncta-params-per-channel entry '{spec}' has empty "
                f"channel name")
        params: dict[str, Any] = {}
        for kv in kv_block.split(","):
            kv = kv.strip()
            if not kv:
                continue
            if "=" not in kv:
                raise ValueError(
                    f"--puncta-params-per-channel '{spec}': '{kv}' must be "
                    f"key=value")
            key, val = kv.split("=", 1)
            key = key.strip()
            val = val.strip()
            if key not in _PUNCTA_PARAM_VALID:
                raise ValueError(
                    f"--puncta-params-per-channel '{spec}': unknown key "
                    f"'{key}'. Valid keys: "
                    f"{sorted(_PUNCTA_PARAM_VALID)}")
            if key in _PUNCTA_PARAM_INT_KEYS:
                params[key] = int(val)
            elif key in _PUNCTA_PARAM_NUMERIC_KEYS:
                params[key] = float(val)
            else:
                params[key] = val
        out[ch] = params
    return out


def resolve_per_channel_cfg(cfg: dict, channel_name: str) -> dict:
    """Return a cfg dict with per-channel overrides applied for one channel.

    Layers (later wins):
      1. cfg as-is (global cellquant config)
      2. cfg["puncta_params_per_channel"][channel_name] overrides
      3. cfg["puncta_rolling_ball"][channel_name] radius surfaced as
         "puncta_rolling_ball_radius" (consumed by detect_puncta)
    """
    eff = dict(cfg)
    per_ch = cfg.get("puncta_params_per_channel", {}) or {}
    eff.update(per_ch.get(channel_name, {}))
    rb = cfg.get("puncta_rolling_ball", {}) or {}
    eff["puncta_rolling_ball_radius"] = int(rb.get(channel_name, 0))
    return eff


def build_config(args: argparse.Namespace) -> dict[str, Any]:
    cfg = dict(DEFAULTS)

    # Layer 2: cell-type preset
    preset_3d_overrides: dict[str, Any] = {}
    if args.cell_type:
        preset_name = args.cell_type.lower()
        if preset_name not in CELL_TYPE_PRESETS:
            raise ValueError(
                f"Unknown cell type '{preset_name}'. "
                f"Available: {list(CELL_TYPE_PRESETS.keys())}")
        preset = dict(CELL_TYPE_PRESETS[preset_name])
        # Stash 3D-specific overrides for later application; don't pollute the
        # 2D config with 3D keys.
        preset_3d_overrides = dict(preset.pop("_3d", {}))
        cfg.update(preset)

    # Layer 3: YAML overrides
    if args.config:
        with open(args.config) as fh:
            overrides = yaml.safe_load(fh)
        if isinstance(overrides, dict):
            for key, val in overrides.items():
                cfg[key] = val

    # Layer 4: CLI args (highest priority — always wins)
    cli_map: dict[str, Any] = {
        "seg_downsample": args.seg_downsample,
        "nuclei_diameter": args.nuclei_diameter,
        "cell_diameter": args.cell_diameter,
        "flow_threshold": args.flow_threshold,
        "cellprob_threshold": args.cellprob_threshold,
        "nucleus_dilate_px": args.nucleus_dilate_px,
        "log_sigma": args.log_sigma,
        "puncta_min_area_px": args.puncta_min_area_px,
        "puncta_max_area_px": args.puncta_max_area_px,
        "puncta_threshold_method": args.puncta_threshold_method,
        "puncta_threshold_fixed": args.puncta_threshold_fixed,
        "puncta_min_circularity": args.puncta_min_circularity,
        "puncta_min_solidity": args.puncta_min_solidity,
        "puncta_compartment": args.puncta_compartment,
        "puncta_compartment_erode_um": args.puncta_compartment_erode_um,
        "qc_downsample": args.qc_downsample,
        "qc_dpi": args.qc_dpi,
        "keep_min_nuclei": args.keep_min_nuclei,
        "keep_max_nuclei": args.keep_max_nuclei,
        "pretrained_model": args.pretrained_model,
        "cpu_threads": args.cpu_threads,
        "filename_pattern": args.filename_pattern,
        "condition_order": args.condition_order,
        "reference_condition": args.reference_condition,
        "cell_seg_channel": args.cell_seg_channel,
        "file_pattern": args.file_pattern,
        "min_cell_area": args.min_cell_area,
        "max_cell_area": args.max_cell_area,
        "colocalization_compartment": args.colocalization_compartment,
        "nucleolar_proximity_channel": args.nucleolar_proximity,
        "proximity_threshold_px": args.proximity_threshold,
        # 3D-specific
        "mode": args.mode,
        "axes": args.axes,
        "seg_3d_method": args.seg_3d_method,
        "stitch_threshold": args.stitch_threshold,
        "puncta_min_volume_vox": args.puncta_min_volume_vox,
        "puncta_max_volume_vox": args.puncta_max_volume_vox,
        "min_cell_volume_vox": args.min_cell_volume_vox,
        "max_cell_volume_vox": args.max_cell_volume_vox,
        "fragmentation_thresholds": args.fragmentation_thresholds,
        "project_z": args.project_z,
        "z_crop_center": args.z_crop_center,
    }
    cli_provided: set[str] = set()
    for key, val in cli_map.items():
        if val is not None:
            cfg[key] = val
            cli_provided.add(key)

    # Voxel size from --voxel-size XY [Z]; only meaningful in 3D
    if args.voxel_size is not None:
        if len(args.voxel_size) == 1:
            cfg["voxel_size_xy_um"] = float(args.voxel_size[0])
        elif len(args.voxel_size) == 2:
            cfg["voxel_size_xy_um"] = float(args.voxel_size[0])
            cfg["voxel_size_z_um"] = float(args.voxel_size[1])
        else:
            raise ValueError("--voxel-size takes 1 (2D) or 2 (3D) values")
        # Explicit marker that the user supplied a COMPLETE voxel size.
        # Downstream code used to infer this by comparing against (1.0, 1.0),
        # which silently misreads a legitimate 1 µm isotropic dataset as "not
        # supplied". Only the two-value form asserts a Z; the one-value form
        # leaves Z at its default, so it must NOT be treated as a full 3D
        # assertion or the ladder would reject metadata that supplies the Z.
        if len(args.voxel_size) == 2:
            cfg["_voxel_size_user_set"] = True

    # 3D proximity threshold (microns), if user provided one
    if args.proximity_threshold_um is not None:
        cfg["proximity_threshold_um"] = float(args.proximity_threshold_um)

    # Boolean flags (store_true — only apply when set)
    if args.no_gpu:
        cfg["use_gpu"] = False
    if args.no_save_masks:
        cfg["save_masks"] = False
    if args.reuse_masks:
        cfg["reuse_masks"] = True
    if args.skip_plots:
        cfg["skip_plots"] = True
    if args.colocalization:
        cfg["colocalization"] = True
    if args.trend:
        cfg["trend"] = True
    if args.condensate_index:
        cfg["condensate_index"] = True
    cfg["assume_isotropic"] = bool(args.assume_isotropic)
    cfg["override_metadata"] = bool(args.override_metadata)
    cfg["allow_2d_colocalization"] = bool(args.allow_2d_colocalization)
    cfg["cell_type"] = args.cell_type.lower() if getattr(args, "cell_type", None) else None

    # Condition map from CLI (parsed from key=value pairs)
    cmap = _parse_condition_map(args.condition_map)
    if cmap is not None:
        cfg["condition_map"] = cmap

    # Per-channel puncta detection overrides
    rb_map = _parse_puncta_rolling_ball(args.puncta_rolling_ball)
    if rb_map:
        cfg["puncta_rolling_ball"] = rb_map
    pcm = _parse_puncta_params_per_channel(args.puncta_params_per_channel)
    if pcm:
        cfg["puncta_params_per_channel"] = pcm

    # Compartment definitions. CLI wins over the YAML layer; either may supply a
    # list of "NAME = EXPR" strings. The YAML layer assigns cfg keys without
    # validation (see the --config merge above), so normalise here rather than
    # letting an unparsed list reach the segmentation loop.
    comp_specs = args.compartment
    if comp_specs:
        cli_provided.add("compartments")
    elif isinstance(cfg.get("compartments"), list):
        comp_specs = cfg["compartments"]
    cfg["compartments"] = (parse_compartment_specs(comp_specs)
                           if comp_specs else {})

    # Stash the preset's 3D overrides on the cfg; main() will apply them once
    # mode is finalized (after auto-detection from the first input file).
    cfg["_preset_3d_overrides"] = preset_3d_overrides
    # Stash the set of keys that came from the CLI so the 3D preset overrides
    # don't stomp on them even if their value happens to equal DEFAULTS.
    cfg["_cli_provided_keys"] = cli_provided

    return cfg


def apply_3d_preset_overrides(cfg: dict[str, Any]) -> None:
    """Apply the cell-type preset's 3D-specific keys after mode is finalized.
    Called in main() once we know the input is 3D. CLI-provided keys are
    honored — they win over the preset even when their value equals DEFAULTS."""
    overrides = cfg.pop("_preset_3d_overrides", {})
    cli_provided: set[str] = cfg.pop("_cli_provided_keys", set())
    for key, val in overrides.items():
        if key in cli_provided:
            continue
        if cfg.get(key) == DEFAULTS.get(key):
            cfg[key] = val


# ---------------------------------------------------------------------------
# Filename metadata parsing
# ---------------------------------------------------------------------------
def _build_filename_regex(pattern: str) -> re.Pattern | None:
    """Convert e.g. 'MAX_{condition}_rep{replicate}' to a named-group regex."""
    if not pattern:
        return None
    # Replace placeholders with unique tokens, escape the rest, swap back
    COND_TOK = "__COND__"
    REP_TOK = "__REP__"
    temp = pattern.replace("{condition}", COND_TOK).replace("{replicate}", REP_TOK)
    temp = re.escape(temp)
    temp = temp.replace(COND_TOK, r"(?P<condition>\w+?)")
    temp = temp.replace(REP_TOK, r"(?P<replicate>\d+)")
    return re.compile(temp, re.IGNORECASE)


def _build_filename_regex_no_rep(pattern: str) -> re.Pattern | None:
    """Fallback regex when {replicate} is at the end: capture condition only."""
    if not pattern or not pattern.endswith("{replicate}"):
        return None
    stripped = pattern[: -len("{replicate}")]
    COND_TOK = "__COND__"
    temp = stripped.replace("{condition}", COND_TOK)
    temp = re.escape(temp)
    temp = temp.replace(COND_TOK, r"(?P<condition>\w+)")
    return re.compile(temp + "$", re.IGNORECASE)


def parse_filename_metadata(stem: str, cfg: dict) -> dict[str, str]:
    info: dict[str, str] = {"condition": "", "replicate": ""}
    if not cfg.get("parse_from_filename", False):
        return info
    pattern = cfg.get("filename_pattern", "")
    regex = _build_filename_regex(pattern)
    if regex is None:
        return info
    m = regex.search(stem)
    if m:
        gd = m.groupdict()
        raw_cond = gd.get("condition", "")
        cmap = cfg.get("condition_map", {})
        info["condition"] = cmap.get(raw_cond.lower(), raw_cond) if raw_cond else ""
        info["replicate"] = gd.get("replicate", "")
    else:
        # If {replicate} is trailing and digits are absent, default to "0"
        fallback = _build_filename_regex_no_rep(pattern)
        if fallback:
            m2 = fallback.search(stem)
            if m2:
                raw_cond = m2.group("condition")
                cmap = cfg.get("condition_map", {})
                info["condition"] = cmap.get(raw_cond.lower(), raw_cond) if raw_cond else ""
                info["replicate"] = "0"
    return info


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Configurable multi-channel cell quantification pipeline")

    # Positional args (shorthand for --images and --channels)
    ap.add_argument("positional", nargs="*", default=[],
                    help="Optional: images_dir followed by channel specs")

    # Required (unless supplied via positional)
    ap.add_argument("--images", default=None,
                    help="Directory containing MIP TIFFs")
    ap.add_argument("--out", required=True,
                    help="Output directory")
    ap.add_argument("--channels", nargs="+", default=None,
                    help='Channel definitions: "pos:Name:role" '
                         '(e.g. "1:DAPI:nucleus" "2:G3BP1:quantify")')

    # Optional channel selection
    ap.add_argument("--puncta-channels", nargs="+", default=None,
                    help="Channel names for puncta detection (e.g. G3BP1 PABPC1)")
    ap.add_argument("--cell-seg-channel", default=None,
                    help="Channel name for cell segmentation "
                         "(default: first cell-boundary or first quantify)")

    # Config
    ap.add_argument("--config", default=None,
                    help="YAML config file (overrides presets; CLI args win over this)")
    ap.add_argument("--cell-type", default=None,
                    choices=list(CELL_TYPE_PRESETS.keys()),
                    help="Cell-type preset for default parameters")

    # Segmentation
    ap.add_argument("--seg-downsample", type=int, default=None)
    ap.add_argument("--no-gpu", action="store_true",
                    help="Disable GPU (CPU fallback)")
    ap.add_argument("--pretrained-model", default=None)
    ap.add_argument("--nuclei-diameter", type=float, default=None)
    ap.add_argument("--cell-diameter", type=float, default=None)
    ap.add_argument("--flow-threshold", type=float, default=None)
    ap.add_argument("--cellprob-threshold", type=float, default=None)
    ap.add_argument("--nucleus-dilate-px", type=int, default=None)

    # Puncta detection
    ap.add_argument("--no-puncta", action="store_true",
                    help="Suppress puncta detection (only compute intensity metrics)")
    ap.add_argument("--log-sigma", type=float, default=None)
    ap.add_argument("--puncta-min-area-px", type=int, default=None)
    ap.add_argument("--puncta-max-area-px", type=int, default=None)
    ap.add_argument("--puncta-threshold-method", default=None,
                    choices=["otsu", "triangle", "fixed"])
    ap.add_argument("--puncta-threshold-fixed", type=float, default=None)
    ap.add_argument("--puncta-min-circularity", type=float, default=None)
    ap.add_argument("--puncta-min-solidity", type=float, default=None)
    ap.add_argument("--puncta-compartment", default=None, metavar="REGION",
                    help="Region puncta are detected in: a built-in "
                         "(cytosol, nucleus, whole-cell, cell, nucleolus) or a "
                         "name defined with --compartment. Validated at startup.")
    ap.add_argument("--compartment", action="append", default=None,
                    metavar='"NAME = EXPR"',
                    help='Define a region by set algebra; repeatable. '
                         'EXPR is TERM [op TERM]... with op one of - & + '
                         '(minus, intersect, union), applied strictly left to '
                         'right. Operators MUST have spaces around them. A TERM '
                         'is a built-in (whole-cell, cell, nucleus, cytosol, '
                         'nucleolus), a NAME defined earlier, or '
                         'exclusion(CH,within=PARENT); suffix ~UM to grow it or '
                         '~-UM to shrink it by that many microns. Example: '
                         '--compartment "cytosol_free = cell - nucleus~0.3"')
    ap.add_argument("--puncta-compartment-erode-um", type=float, default=None,
                    metavar="UM",
                    help="Erode the puncta compartment inward by this many "
                         "microns before LoG detection. Removes the cell-"
                         "membrane gradient that the LoG-Otsu detector "
                         "otherwise classifies as outline-shaped puncta. "
                         "Anisotropic-aware in 3D. Default: 0.5 µm for "
                         "mammalian preset, 0 otherwise. Set to 0 to disable.")
    ap.add_argument("--puncta-rolling-ball", nargs="+", default=None,
                    metavar="CH:RADIUS",
                    help="Per-channel rolling-ball (white tophat) background "
                         "subtraction radius applied to the puncta channel "
                         "before LoG detection. e.g. 'G3BP1:25 PABPC1:0'. "
                         "Radius is in pixels (2D) or voxels (3D); 0 disables.")
    ap.add_argument("--puncta-params-per-channel", nargs="+", default=None,
                    metavar="CH:K=V,K=V",
                    help="Per-channel puncta detection overrides. e.g. "
                         "'G3BP1:log_sigma=2.5,puncta_threshold_method=triangle "
                         "PABPC1:log_sigma=1.8'. Supported keys: log_sigma, "
                         "puncta_threshold_method, puncta_threshold_fixed, "
                         "puncta_min_area_px, puncta_max_area_px, "
                         "puncta_min_volume_vox, puncta_max_volume_vox, "
                         "puncta_min_circularity, puncta_min_solidity, "
                         "puncta_compartment.")

    # QC
    ap.add_argument("--qc-downsample", type=int, default=None)
    ap.add_argument("--qc-dpi", type=int, default=None)

    # Keep gate
    ap.add_argument("--keep-min-nuclei", type=int, default=None)
    ap.add_argument("--keep-max-nuclei", type=int, default=None)

    # Metadata
    ap.add_argument("--filename-pattern", default=None,
                    help='Pattern with {condition} and {replicate} placeholders '
                         '(e.g. "MAX_{condition}_rep{replicate}")')
    ap.add_argument("--condition-order", nargs="+", default=None,
                    help="X-axis order for superplots (e.g. control arsenite)")
    ap.add_argument("--reference-condition", default=None,
                    help="Reference condition for pairwise Wilcoxon p-values "
                         "(e.g. NT_NT)")
    ap.add_argument("--condition-map", nargs="+", default=None,
                    help='Map raw condition names: key=value pairs '
                         '(e.g. "ctrl=control" "as=arsenite")')

    # I/O options
    ap.add_argument("--file-pattern", default=None,
                    help='Glob pattern for image files (e.g. "*.tif")')
    ap.add_argument("--no-save-masks", action="store_true",
                    help="Disable writing segmentation/puncta masks to disk")
    ap.add_argument("--reuse-masks", action="store_true",
                    help="Load existing masks from --out masks/ dir instead "
                         "of running Cellpose")
    ap.add_argument("--skip-plots", action="store_true",
                    help="Skip superplot and Prism CSV generation")

    # Cell area filtering
    ap.add_argument("--min-cell-area", type=int, default=None,
                    help="Minimum cell area in pixels")
    ap.add_argument("--max-cell-area", type=int, default=None,
                    help="Maximum cell area in pixels")

    # Colocalization
    ap.add_argument("--colocalization", action="store_true",
                    help="Compute pairwise Pearson + Manders")
    ap.add_argument("--colocalization-compartment", default=None,
                    choices=["whole-cell", "cytosol", "nucleus"],
                    help="Compartment for colocalization (default whole-cell)")

    # Nucleolar proximity
    ap.add_argument("--nucleolar-proximity", default=None, metavar="CHAN",
                    help="Channel name defining nucleolar mask")
    ap.add_argument("--proximity-threshold", type=float, default=None,
                    help="Distance threshold for 'proximal' (default 5 px)")

    # Plotting
    ap.add_argument("--trend", action="store_true",
                    help="Trend line on multi-condition plots")

    # Specialized metrics (opt-in)
    ap.add_argument("--condensate-index", action="store_true",
                    help="Compute the Condensate Index (p95/mean intensity ratio) "
                         "per channel per cell. Off by default — relevant for "
                         "condensate/punctate biology, irrelevant otherwise.")
    ap.add_argument("--allow-2d-colocalization", action="store_true",
                    help="Permit colocalization on 2D / MIP input. Projected "
                         "Pearson/Manders are not defensible on a projection; "
                         "without this flag cellquant refuses colocalization in "
                         "2D mode. When set, colocalization.csv is stamped with "
                         "projection_derived=True.")

    # Threading
    ap.add_argument("--cpu-threads", type=int, default=None)

    # ----- 3D-specific flags (ignored in 2D mode) -----
    ap.add_argument("--mode", default=None, choices=["auto", "2d", "3d"],
                    help="Force pipeline mode (default: auto-detect from input shape)")
    ap.add_argument("--voxel-size", nargs="+", type=float, default=None,
                    metavar=("XY_UM", "Z_UM"),
                    help="Voxel size in microns. 1 value (XY) for 2D, 2 values (XY Z) for 3D. "
                         "Auto-read from OME/ImageJ metadata when present.")
    ap.add_argument("--assume-isotropic", action="store_true",
                    help="In 3D, proceed with 1.0 µm isotropic voxels when no "
                         "voxel size is available from --voxel-size or file "
                         "metadata. Without this flag cellquant aborts rather "
                         "than silently assuming isotropic voxels.")
    ap.add_argument("--override-metadata", action="store_true",
                    help="Allow --voxel-size to override the file's OME/ImageJ "
                         "voxel metadata when the two disagree by more than one "
                         "percent. Without this flag such a disagreement aborts "
                         "the run. The override is recorded in provenance.json.")
    ap.add_argument("--axes", default=None,
                    help="Override axis order for 3D inputs, e.g. ZCYX or CZYX. "
                         "Default: trust the file's declared axes; otherwise heuristic from shape.")
    ap.add_argument("--seg-3d-method", default=None, choices=["stitch", "full"],
                    help="3D Cellpose mode: 'stitch' (fast 2D-per-Z + IoU stitch) "
                         "or 'full' (true do_3D=True).")
    ap.add_argument("--stitch-threshold", type=float, default=None,
                    help="IoU threshold for Z-stitching (default 0.4).")
    ap.add_argument("--puncta-min-volume-vox", type=int, default=None,
                    help="3D puncta volume floor in voxels.")
    ap.add_argument("--puncta-max-volume-vox", type=int, default=None,
                    help="3D puncta volume ceiling in voxels. For yeast at "
                         "0.1 µm voxels, ~3000 (sphere of radius ~0.9 µm).")
    ap.add_argument("--min-cell-volume-vox", type=int, default=None,
                    help="3D cell volume filter (voxels).")
    ap.add_argument("--max-cell-volume-vox", type=int, default=None)
    ap.add_argument("--proximity-threshold-um", type=float, default=None,
                    help="3D proximity threshold in microns (replaces "
                         "--proximity-threshold for 3D inputs).")
    ap.add_argument("--project-z", default=None, choices=["max", "sum", "mean"],
                    help="Project z-stacks to 2D before analysis using the given "
                         "method. When set on 3D input, cellquant reduces the Z "
                         "axis (via max-intensity / sum / mean) and runs the 2D "
                         "pipeline. Use this to generate MIPs inside cellquant "
                         "instead of preprocessing with Fiji/napari. Has no "
                         "effect on inputs that are already 2D.")
    ap.add_argument("--z-crop-center", type=int, default=None,
                    help="Crop 3D inputs to N central slices around the per-image "
                         "peak integrated-signal Z, then segment/measure on the "
                         "cropped volume. Cells whose 3D mask touches the top or "
                         "bottom Z plane of the cropped volume are excluded (so "
                         "truncated half-cells don't corrupt volume). "
                         "Useful on CPU when one stack contains more axial range "
                         "than needed; the dropped tail cells are honest sampling, "
                         "not a measurement bias. 0 (default) = disabled.")

    # Fragmentation index — works in both 2D and 3D
    ap.add_argument("--fragmentation-thresholds", type=int, default=None,
                    help="Number of threshold levels for the persistence-based "
                         "fragmentation index. 0 disables it. Default 10.")

    args = ap.parse_args()

    # Resolve positional args as shorthand for --images and --channels
    if args.positional:
        if args.images is None and args.channels is None and len(args.positional) >= 2:
            args.images = args.positional[0]
            args.channels = args.positional[1:]
        elif args.images is None:
            args.images = args.positional[0]

    if args.images is None:
        ap.error("--images is required (or supply images_dir as first positional arg)")
    if args.channels is None:
        ap.error("--channels is required (or supply channel specs as positional args)")

    return args


# ---------------------------------------------------------------------------
# Utilities (reused verbatim from reference)
# ---------------------------------------------------------------------------
def robust_rescale(img: np.ndarray, p_low: float, p_high: float) -> np.ndarray:
    lo, hi = np.percentile(img, [p_low, p_high])
    if hi <= lo:
        return np.zeros_like(img, dtype=np.float32)
    out = (img - lo) / (hi - lo)
    return np.clip(out, 0.0, 1.0).astype(np.float32, copy=False)


def downsample_nn(img: np.ndarray, factor: int) -> np.ndarray:
    if factor <= 1:
        return img
    return img[::factor, ::factor]


def upsample_labels_nn(
    lab_small: np.ndarray, out_shape: tuple[int, int]
) -> np.ndarray:
    return resize(
        lab_small, out_shape, order=0, preserve_range=True, anti_aliasing=False,
    ).astype(np.int32, copy=False)


def _diam_or_none(x) -> float | None:
    if x is None:
        return None
    try:
        v = float(x)
    except Exception:
        return None
    return v if v > 0 else None


def _set_threads(cfg: dict) -> None:
    n = int(cfg.get("cpu_threads", 0) or 0)
    if n <= 0:
        return
    os.environ.setdefault("OMP_NUM_THREADS", str(n))
    os.environ.setdefault("MKL_NUM_THREADS", str(n))
    if torch is not None:
        try:
            torch.set_num_threads(n)
        except Exception:
            pass


def _safe_p95(im: np.ndarray, mask: np.ndarray) -> float:
    v = im[mask]
    if v.size == 0:
        return np.nan
    return float(np.percentile(v.astype(np.float64), 95))


def _condensate_index(p95: float, mean: float) -> float:
    if np.isnan(p95) or np.isnan(mean) or mean == 0:
        return np.nan
    return float(p95 / mean)


def safe_mean(im: np.ndarray, mask: np.ndarray) -> float:
    v = im[mask]
    return float(v.mean()) if v.size else np.nan


# ---------------------------------------------------------------------------
# Skimage compatibility wrappers (reused verbatim)
# ---------------------------------------------------------------------------
def remove_small_objects_compat(bw: np.ndarray, min_size: int) -> np.ndarray:
    try:
        return morphology.remove_small_objects(bw, max_size=min_size - 1)
    except TypeError:
        return morphology.remove_small_objects(bw, min_size=min_size)


def remove_small_holes_compat(bw: np.ndarray, area_threshold: int) -> np.ndarray:
    try:
        return morphology.remove_small_holes(bw, max_size=area_threshold - 1)
    except TypeError:
        return morphology.remove_small_holes(bw, area_threshold=area_threshold)


# ---------------------------------------------------------------------------
# Image I/O
# ---------------------------------------------------------------------------
def find_images(img_dir: Path, cfg: dict) -> list[Path]:
    file_pattern = cfg.get("file_pattern")
    if file_pattern:
        paths = sorted(img_dir.glob(file_pattern))
    else:
        exts = cfg.get("exts", [".tif", ".tiff"])
        paths = []
        for ext in exts:
            paths.extend(sorted(img_dir.glob(f"*{ext}")))
    return [p for p in paths if p.is_file()]


def _is_aics_format(path: Path) -> bool:
    """True iff this file should be opened via aicsimageio rather than tifffile."""
    return path.suffix.lower() in _AICS_EXTS


def _user_supplied_voxel(cfg: dict) -> bool:
    """True iff the voxel size in cfg came from the user, not from DEFAULTS.

    Single source of truth for "should file metadata be allowed to win?". Every
    loader and the 3D voxel ladder must agree on this: if the TIFF path and the
    multiformat path answer differently, the run reports one voxel size and
    computes with another.

    ``_voxel_size_user_set`` is set by build_config only for the complete
    two-value ``--voxel-size XY Z`` form. The ``!= (1.0, 1.0)`` fallback covers
    values arriving by other routes (a ``--config`` YAML layer, or the
    single-value form, which leaves Z at its default).
    """
    return bool(cfg.get("_voxel_size_user_set")) or (
        (cfg.get("voxel_size_xy_um"), cfg.get("voxel_size_z_um")) != (1.0, 1.0))


def _require_aicsimageio(path: Path) -> None:
    """Backwards-compat name; check that a multiformat backend is installed."""
    if not HAS_MULTIFORMAT:
        raise SystemExit(
            f"{path.name}: reading '{path.suffix}' files requires a "
            f"multiformat backend. Install bioio (recommended):\n"
            f"  pip install bioio bioio-nd2 bioio-czi bioio-lif bioio-ome-tiff\n"
            f"Or the legacy aicsimageio:\n"
            f"  pip install aicsimageio\n"
            f"Alternatively, export your data to TIFF or OME-TIFF."
        )


def _aics_load_2d(path: Path, channels: list[dict]) -> dict[str, np.ndarray]:
    """Load a 2D plane from any aicsimageio-supported format."""
    _require_aicsimageio(path)
    img = _MultiformatReader(str(path))
    arr = img.get_image_data("CYX", T=0, Z=0)  # (C, Y, X) float
    if arr.shape[0] < len(channels):
        raise ValueError(
            f"{path.name}: has {arr.shape[0]} channels but {len(channels)} declared")
    result: dict[str, np.ndarray] = {}
    for ch in channels:
        idx = ch["position"] - 1
        result[ch["name"]] = arr[idx].astype(np.float32, copy=False)
    return result


def _aics_load_3d(
    path: Path, channels: list[dict], cfg: dict
) -> tuple[dict[str, np.ndarray], tuple[float, float]]:
    """Load a 3D z-stack from any aicsimageio-supported format."""
    _require_aicsimageio(path)
    img = _MultiformatReader(str(path))
    arr = img.get_image_data("CZYX", T=0)  # (C, Z, Y, X)

    # Pull voxel size from the file unless the user explicitly supplied one.
    voxel = (cfg["voxel_size_xy_um"], cfg["voxel_size_z_um"])
    if not _user_supplied_voxel(cfg):
        try:
            pps = img.physical_pixel_sizes
            # AICSImage returns PhysicalPixelSizes(Z, Y, X) in microns.
            if pps.Z and pps.X:
                voxel = (float(pps.X), float(pps.Z))
        except Exception:
            pass

    if arr.shape[0] < len(channels):
        raise ValueError(
            f"{path.name}: has {arr.shape[0]} channels but {len(channels)} declared")

    result: dict[str, np.ndarray] = {}
    for ch in channels:
        idx = ch["position"] - 1
        result[ch["name"]] = arr[idx].astype(np.float32, copy=False)
    return result, voxel


def _aics_dimensions(path: Path) -> tuple[int, int, int, int]:
    """Return (T, C, Z, Y, X) → just (C, Z, Y, X) since we ignore T.

    Returns the dimensions reported by aicsimageio; used by detect_input_mode
    to decide 2D vs 3D for non-tiff formats without loading pixel data.
    """
    _require_aicsimageio(path)
    img = _MultiformatReader(str(path))
    dims = img.dims
    return (int(dims.C), int(dims.Z), int(dims.Y), int(dims.X))


def load_tiff(path: Path, channels: list[dict]) -> dict[str, np.ndarray]:
    """Load a multi-channel 2D image (MIP) and return {channel_name: 2D array}.

    Despite the name, dispatches to aicsimageio for ND2/CZI/LIF; tifffile
    handles .tif/.tiff (the paper-validated path).
    """
    if _is_aics_format(path):
        return _aics_load_2d(path, channels)
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 2:
        arr = arr[np.newaxis, :, :]  # single-channel 2D → (1, Y, X)
    if arr.ndim != 3:
        raise ValueError(f"{path.name}: expected 2D or 3D array, got shape {arr.shape}")

    n_ch = len(channels)

    # Determine CYX vs YXC layout
    if arr.shape[0] <= arr.shape[-1] and arr.shape[0] <= 10:
        cyx = arr  # channels-first
    elif arr.shape[-1] <= arr.shape[0] and arr.shape[-1] <= 10:
        cyx = np.moveaxis(arr, -1, 0)  # channels-last → CYX
    else:
        raise ValueError(
            f"{path.name}: cannot determine channel axis for shape {arr.shape}")

    if cyx.shape[0] < n_ch:
        raise ValueError(
            f"{path.name}: has {cyx.shape[0]} channels but {n_ch} defined")

    result: dict[str, np.ndarray] = {}
    for ch in channels:
        idx = ch["position"] - 1  # 1-based → 0-based
        if idx >= cyx.shape[0]:
            raise ValueError(
                f"{path.name}: channel position {ch['position']} out of range "
                f"(image has {cyx.shape[0]} channels)")
        result[ch["name"]] = cyx[idx].astype(np.float32, copy=False)
    return result


# ---------------------------------------------------------------------------
# 3D-specific I/O helpers
# ---------------------------------------------------------------------------
def _read_ome_voxel_size(path: Path) -> tuple[float, float] | None:
    """Try to read XY/Z voxel size in µm from file metadata.

    For .tif/.tiff: parses OME-TIFF Pixels element or ImageJ metadata.
    For .nd2/.czi/.lif: uses aicsimageio's physical_pixel_sizes.
    Returns (xy_um, z_um), or None if metadata is absent or unparseable.
    """
    if _is_aics_format(path):
        if not HAS_AICSIMAGEIO:
            return None
        try:
            img = _MultiformatReader(str(path))
            pps = img.physical_pixel_sizes
            if pps.X and pps.Z:
                return (float(pps.X), float(pps.Z))
        except Exception:
            return None
        return None
    try:
        with tiff.TiffFile(str(path)) as tf:
            if tf.is_ome and tf.ome_metadata:
                from xml.etree import ElementTree as ET
                root = ET.fromstring(tf.ome_metadata)
                ns = (
                    {"o": root.tag.split("}")[0].strip("{")}
                    if "}" in root.tag else {}
                )
                pix = (root.find(".//o:Pixels", ns) if ns
                       else root.find(".//Pixels"))
                if pix is not None:
                    xy = pix.attrib.get("PhysicalSizeX")
                    z = pix.attrib.get("PhysicalSizeZ")
                    if xy and z:
                        return (float(xy), float(z))
            if tf.is_imagej and tf.imagej_metadata:
                meta = tf.imagej_metadata
                spacing = meta.get("spacing")
                page = tf.pages[0]
                xres = page.tags.get("XResolution")
                xy = None
                if xres is not None and xres.value[0] != 0:
                    xy = float(xres.value[1] / xres.value[0])
                if xy is not None and spacing is not None:
                    return (float(xy), float(spacing))
    except Exception:
        return None
    return None


def _read_pixel_size_xy_um(path: Path) -> float | None:
    """Lateral pixel size in µm from file metadata, independent of any Z spacing.

    `_read_ome_voxel_size` deliberately requires BOTH a lateral and an axial
    value, because a 3D analysis needs both. A two-dimensional projection has no
    Z spacing, so that reader returns None for it — and every micron-denominated
    parameter in 2D then fell back to an EDT sampling of 1.0, silently
    reinterpreting microns as pixels. Any value below 1.0 became a no-op: the
    mammalian preset's 0.5 µm puncta-compartment erosion never removed a single
    pixel in a default 2D run.
    """
    if _is_aics_format(path):
        if not HAS_AICSIMAGEIO:
            return None
        try:
            pps = _MultiformatReader(str(path)).physical_pixel_sizes
            return float(pps.X) if pps.X else None
        except Exception:
            return None
    try:
        with tiff.TiffFile(str(path)) as tf:
            if tf.is_ome and tf.ome_metadata:
                from xml.etree import ElementTree as ET
                root = ET.fromstring(tf.ome_metadata)
                ns = ({"o": root.tag.split("}")[0].strip("{")}
                      if "}" in root.tag else {})
                pix = (root.find(".//o:Pixels", ns) if ns
                       else root.find(".//Pixels"))
                if pix is not None and pix.attrib.get("PhysicalSizeX"):
                    return float(pix.attrib["PhysicalSizeX"])
            xres = tf.pages[0].tags.get("XResolution")
            if xres is not None and xres.value[0] != 0:
                unit = tf.pages[0].tags.get("ResolutionUnit")
                ij_unit = (tf.imagej_metadata or {}).get("unit", "") if tf.is_imagej else ""
                # Only trust the resolution tag when it is expressed in microns.
                if ij_unit in ("micron", "um", "µm") or (
                        unit is not None and int(unit.value) == 1 and tf.is_imagej):
                    return float(xres.value[1] / xres.value[0])
    except Exception:
        return None
    return None


def _um_denominated_2d_params(cfg: dict) -> list[tuple[str, str, bool]]:
    """Active parameters expressed in microns, which without a pixel size would
    be silently reinterpreted as pixels in 2D.

    Returns (label, cfg_key, user_asserted). The distinction matters: a value the
    user typed is an assertion about physical units that cannot be honoured, so
    the run should refuse. A value that merely came from a cell-type preset
    should not abort a first run on data that carries no resolution metadata —
    including this tool's own example images — but it must not be applied
    silently either.
    """
    cli_keys = cfg.get("_cli_provided_keys") or set()
    active: list[tuple[str, str, bool]] = []
    if float(cfg.get("puncta_compartment_erode_um", 0.0) or 0.0) > 0:
        active.append(("--puncta-compartment-erode-um",
                       "puncta_compartment_erode_um",
                       "puncta_compartment_erode_um" in cli_keys))
    for _name, d in (cfg.get("compartments") or {}).items():
        if any(pad for _, _, pad in d["terms"]):
            # A '~' pad is always something the user wrote by hand.
            active.append((f"--compartment \"{d['raw']}\" (uses a '~' pad)",
                           "compartments", True))
    return active


def _read_ome_channel_names(path: Path) -> list[str]:
    """Try to read channel names from file metadata.

    For .tif/.tiff: parses OME-TIFF Channel elements.
    For .nd2/.czi/.lif: uses aicsimageio's channel_names.
    Returns a list of channel names in source order, or [] if
    absent/unparseable. Used for cross-checking against user-supplied
    --channels specs and for propagating into the provenance sidecar.
    """
    if _is_aics_format(path):
        if not HAS_AICSIMAGEIO:
            return []
        try:
            img = _MultiformatReader(str(path))
            names = list(img.channel_names or [])
            return [n for n in names if n]
        except Exception:
            return []
    try:
        with tiff.TiffFile(str(path)) as tf:
            if tf.is_ome and tf.ome_metadata:
                from xml.etree import ElementTree as ET
                root = ET.fromstring(tf.ome_metadata)
                ns = (
                    {"o": root.tag.split("}")[0].strip("{")}
                    if "}" in root.tag else {}
                )
                if ns:
                    chans = root.findall(".//o:Channel", ns)
                else:
                    chans = root.findall(".//Channel")
                names = [
                    c.attrib.get("Name") or c.attrib.get("Fluor") or ""
                    for c in chans
                ]
                return [n for n in names if n]
    except Exception:
        pass
    return []


def _voxel_size_is_implausible(xy_um: float, z_um: float) -> str | None:
    """Return a human-readable reason if a voxel size looks wrong, else None.

    Plausible microscopy ranges (very generous):
      XY: 0.01 µm (high-end SIM/STED) to 10 µm (low-mag whole-tissue)
      Z:  0.05 µm (dense optical sectioning) to 50 µm (sparse stack)
    Values outside these ranges almost always indicate a metadata unit
    mismatch (e.g. a stack recorded in nm being read as µm) rather than a
    real acquisition.
    """
    if not (0.01 <= xy_um <= 10.0):
        return f"XY voxel {xy_um} µm is outside plausible 0.01–10 µm range"
    if not (0.05 <= z_um <= 50.0):
        return f"Z voxel {z_um} µm is outside plausible 0.05–50 µm range"
    if z_um < xy_um:
        return (f"Z voxel {z_um} µm < XY voxel {xy_um} µm — Z is rarely "
                f"smaller than XY; check metadata units")
    return None


def _detect_axes(
    arr: np.ndarray, n_channels: int, declared: str | None = None
) -> str:
    """Detect axis order for a 3D/4D TIFF.

    Trusts the file's declared axes string (e.g. 'ZCYX' from ImageJ) when
    present; falls back to a shape-based heuristic. Returns one of:
        'ZYX' (single-channel stack), 'ZCYX', 'CZYX', 'ZYXC', etc.
    """
    if declared:
        d = declared.upper().replace("S", "C")
        kept = "".join(c for c in d if c in "ZCYX")
        if len(kept) == arr.ndim and set(kept) >= {"Y", "X"}:
            return kept

    if arr.ndim == 3:
        return "ZYX"
    if arr.ndim != 4:
        raise ValueError(
            f"Expected 3D or 4D array, got {arr.ndim}D shape {arr.shape}")

    shape = arr.shape
    # Channel axis: dim equal to n_channels and ≤ 10
    cand_c = [i for i, s in enumerate(shape) if s == n_channels and s <= 10]
    if not cand_c:
        small = [(i, s) for i, s in enumerate(shape) if s <= 10]
        if not small:
            raise ValueError(
                f"Cannot identify channel axis in shape {shape} "
                f"for {n_channels} channels")
        cand_c = [min(small, key=lambda x: x[1])[0]]
    c_axis = cand_c[0]

    spatial = sorted(
        [(i, s) for i, s in enumerate(shape) if i != c_axis],
        key=lambda x: x[1], reverse=True,
    )
    yx_axes = {spatial[0][0], spatial[1][0]}
    z_axis = next(i for i in range(4) if i != c_axis and i not in yx_axes)

    layout = [None, None, None, None]
    layout[z_axis] = "Z"
    layout[c_axis] = "C"
    layout[spatial[0][0]] = "Y"
    layout[spatial[1][0]] = "X"
    return "".join(layout)  # type: ignore[arg-type]


def load_zstack(
    path: Path, channels: list[dict], cfg: dict
) -> tuple[dict[str, np.ndarray], tuple[float, float]]:
    """Load a multi-channel 3D Z-stack and reshape to (Z, Y, X) per channel.

    Returns:
        images: dict {channel_name: 3D array shape (Z, Y, X)}
        voxel_size: (xy_um, z_um) — auto-read from metadata unless overridden

    Dispatches to aicsimageio for ND2/CZI/LIF; tifffile handles .tif/.tiff.
    """
    if _is_aics_format(path):
        return _aics_load_3d(path, channels, cfg)
    arr = np.asarray(tiff.imread(str(path)))

    voxel_size = (cfg["voxel_size_xy_um"], cfg["voxel_size_z_um"])
    auto = _read_ome_voxel_size(path)
    # Same rule as the multiformat loader and the 3D ladder: an explicit
    # user-supplied voxel size wins over file metadata. TIFF is the primary
    # input format, so a sentinel that disagreed here would make the run
    # compute with one voxel size while the banner and provenance.json
    # reported another.
    if auto is not None and not _user_supplied_voxel(cfg):
        voxel_size = auto

    declared_axes: str | None = None
    try:
        with tiff.TiffFile(str(path)) as tf:
            if tf.series and tf.series[0].axes:
                declared_axes = tf.series[0].axes
    except Exception:
        pass

    n_ch = len(channels)
    axes = cfg.get("axes")
    if axes is None:
        axes = _detect_axes(arr, n_ch, declared=declared_axes)

    if arr.ndim == 3:
        # (Z, Y, X) — single-channel
        if n_ch != 1:
            raise ValueError(
                f"{path.name}: 3D array but {n_ch} channels declared")
        czyx = arr[np.newaxis, :, :, :]
    else:
        target = "CZYX"
        perm = [axes.index(a) for a in target]
        czyx = np.transpose(arr, perm)

    if czyx.shape[0] < n_ch:
        raise ValueError(
            f"{path.name}: has {czyx.shape[0]} channels, {n_ch} declared")

    result: dict[str, np.ndarray] = {}
    for ch in channels:
        idx = ch["position"] - 1
        result[ch["name"]] = czyx[idx].astype(np.float32, copy=False)
    return result, voxel_size


def detect_input_mode(path: Path, n_channels: int) -> str:
    """Look at a file's shape to decide whether it's 2D (MIP) or 3D (Z-stack).

    Returns "2d" or "3d". Used when --mode auto.
    Dispatches to aicsimageio for ND2/CZI/LIF; tifffile for .tif/.tiff.
    """
    if _is_aics_format(path):
        try:
            C, Z, Y, X = _aics_dimensions(path)
        except SystemExit:
            raise
        except Exception:
            return "2d"  # safe default if probe fails
        return "3d" if Z > 1 else "2d"

    try:
        with tiff.TiffFile(str(path)) as tf:
            if tf.series:
                shape = tf.series[0].shape
            else:
                shape = tf.pages[0].shape
    except Exception:
        # If we can't peek, fall back to imread
        arr = np.asarray(tiff.imread(str(path)))
        shape = arr.shape

    # 2D if ndim <= 3 (channel axis allowed); 3D if ndim == 4 (Z + C + Y + X)
    if len(shape) <= 3:
        # Could be (Y,X), (C,Y,X), or (Z,Y,X) for single-channel stacks.
        # Single-channel stacks are rare and ambiguous — if the smallest dim
        # is much larger than typical channel counts (>10) AND the file looks
        # like a stack from metadata, treat as 3D. Otherwise 2D.
        if len(shape) == 3 and min(shape) > 10:
            return "3d"  # likely (Z, Y, X) single-channel stack
        return "2d"
    return "3d"


def _project_z(stack_zyx: np.ndarray, method: str) -> np.ndarray:
    """Reduce a (Z, Y, X) stack to a (Y, X) projection.

    method ∈ {"max", "sum", "mean"} — corresponds to Fiji's standard MIP /
    sum-projection / mean-projection operations.
    """
    if method == "max":
        return stack_zyx.max(axis=0)
    if method == "sum":
        return stack_zyx.sum(axis=0).astype(np.float32, copy=False)
    if method == "mean":
        return stack_zyx.mean(axis=0).astype(np.float32, copy=False)
    raise ValueError(f"Unknown --project-z method: {method}")


def load_images(
    path: Path, channels: list[dict], cfg: dict
) -> tuple[dict[str, np.ndarray], str, tuple[float, float] | None]:
    """Unified loader. Returns (images, mode, voxel_size).

    voxel_size is None in 2D mode, (xy_um, z_um) in 3D mode.

    Dispatch logic:
      - file is 2D-on-disk → load 2D, return mode="2d"
      - file is 3D-on-disk and --project-z set → load 3D, project, return mode="2d"
      - file is 3D-on-disk and --project-z not set → load 3D, return mode="3d"
      - file is 3D-on-disk and cfg["mode"]=="2d" forced without --project-z →
        fail with a clear error explaining the resolution

    The per-file dimensionality is always trusted; cfg["mode"] acts as a hint
    for pipeline-level behavior and as a force override.
    """
    forced_mode = cfg.get("mode", "auto")
    project_method = cfg.get("project_z")
    file_dim = detect_input_mode(path, len(channels))

    if file_dim == "2d":
        return load_tiff(path, channels), "2d", None

    # file_dim == "3d"
    images, voxel = load_zstack(path, channels, cfg)

    if project_method:
        projected = {name: _project_z(arr, project_method)
                     for name, arr in images.items()}
        return projected, "2d", voxel

    if forced_mode == "2d":
        raise SystemExit(
            f"{path.name}: file is 3D on disk but --mode 2d was forced "
            f"without --project-z. Either drop --mode 2d (cellquant runs 3D "
            f"natively) or add --project-z max/sum/mean to project the stack "
            f"to 2D before analysis."
        )

    return images, "3d", voxel


# ---------------------------------------------------------------------------
# Cellpose init + eval (reused verbatim, parameterized on cfg)
# ---------------------------------------------------------------------------
def _is_mps_system() -> bool:
    """Check if running on Apple Silicon (MPS backend)."""
    if torch is None:
        return False
    return hasattr(torch.backends, "mps") and torch.backends.mps.is_available()


def _describe_device(cfg: dict) -> str:
    """Return a human-readable description of the device cellquant will use."""
    if not cfg.get("use_gpu", True):
        return "CPU (--no-gpu set)"
    if torch is None:
        return "CPU (torch not importable)"
    if torch.cuda.is_available():
        try:
            name = torch.cuda.get_device_name(0)
        except Exception:
            name = "CUDA device"
        return f"CUDA — {name}"
    if _is_mps_system():
        # cellquant forces CPU on MPS because cpsam Transformer ops are
        # not yet supported — see init_model.
        return "CPU (Apple Silicon — MPS not supported by cpsam Transformer)"
    return "CPU (no GPU detected)"


def init_model(cfg: dict) -> models.CellposeModel:
    # Announce the device choice up front so it's visible in run logs and
    # provenance — useful for reproducibility and triage.
    print(f"[device] {_describe_device(cfg)}")

    # MPS (Apple Silicon) doesn't fully support cpsam Transformer ops;
    # force CPU and disable BFloat16 automatically.
    if cfg["use_gpu"] and _is_mps_system():
        print("[warn] MPS GPU not supported by cpsam Transformer; using CPU")
        cfg["use_gpu"] = False
    try:
        return models.CellposeModel(
            gpu=cfg["use_gpu"], pretrained_model=cfg["pretrained_model"],
            use_bfloat16=cfg["use_gpu"])  # BFloat16 only on CUDA
    except TypeError as exc:
        if cfg["use_gpu"] and "BFloat16" in str(exc):
            print(f"[warn] GPU init failed ({exc}); retrying on CPU")
            cfg["use_gpu"] = False
            return models.CellposeModel(
                gpu=False, pretrained_model=cfg["pretrained_model"],
                use_bfloat16=False)
        raise


def eval_masks(
    model: models.CellposeModel,
    img2d: np.ndarray,
    diameter: float | None,
    cfg: dict,
) -> np.ndarray:
    out = model.eval(
        img2d,
        diameter=diameter,
        flow_threshold=cfg["flow_threshold"],
        cellprob_threshold=cfg["cellprob_threshold"],
    )
    masks = out[0] if isinstance(out, (tuple, list)) else out
    return masks.astype(np.int32, copy=False)


def eval_masks_3d(
    model: models.CellposeModel,
    stack_zyx: np.ndarray,
    diameter: float | None,
    cfg: dict,
) -> np.ndarray:
    """Run Cellpose on a 3D stack (Z, Y, X). Method per cfg['seg_3d_method']."""
    method = cfg.get("seg_3d_method", "stitch")
    aniso = cfg["voxel_size_z_um"] / cfg["voxel_size_xy_um"]
    common = dict(
        diameter=diameter,
        flow_threshold=cfg["flow_threshold"],
        cellprob_threshold=cfg["cellprob_threshold"],
        z_axis=0,  # our stacks are (Z, Y, X) so Z is axis 0
    )
    if method == "full":
        out = model.eval(stack_zyx, do_3D=True, anisotropy=aniso, **common)
    elif method == "stitch":
        out = model.eval(
            stack_zyx, do_3D=False,
            stitch_threshold=float(cfg["stitch_threshold"]), **common,
        )
    else:
        raise ValueError(f"Unknown seg_3d_method: {method}")
    masks = out[0] if isinstance(out, (tuple, list)) else out
    return masks.astype(np.int32, copy=False)


def segment_nuclei(
    model: models.CellposeModel, nuc_img: np.ndarray, cfg: dict
) -> np.ndarray:
    if nuc_img.ndim == 3:
        return eval_masks_3d(model, nuc_img,
                             _diam_or_none(cfg.get("nuclei_diameter")), cfg)
    return eval_masks(model, nuc_img, _diam_or_none(cfg["nuclei_diameter"]), cfg)


def segment_cells(
    model: models.CellposeModel, cell_img: np.ndarray, cfg: dict
) -> np.ndarray:
    if cell_img.ndim == 3:
        return eval_masks_3d(model, cell_img,
                             _diam_or_none(cfg.get("cell_diameter")), cfg)
    return eval_masks(model, cell_img, _diam_or_none(cfg["cell_diameter"]), cfg)


def map_nuclei_to_cells(
    nuc_mask: np.ndarray, cell_mask: np.ndarray
) -> dict[int, list[int]]:
    mapping: dict[int, list[int]] = {}
    nuc_ids = np.unique(nuc_mask)
    nuc_ids = nuc_ids[nuc_ids != 0]
    for nid in nuc_ids:
        overlap_cells = cell_mask[nuc_mask == nid]
        overlap_cells = overlap_cells[overlap_cells != 0]
        if overlap_cells.size == 0:
            continue
        cid = int(np.bincount(overlap_cells).argmax())
        mapping.setdefault(cid, []).append(int(nid))
    return mapping


# ---------------------------------------------------------------------------
# Compartment mask creation
# ---------------------------------------------------------------------------
BUILTIN_COMPARTMENTS = ("whole-cell", "cell", "nucleus", "cytosol", "nucleolus")
_COMPARTMENT_OPS = {"-", "&", "+"}


def parse_compartment_specs(specs: list[str] | None) -> dict[str, dict]:
    """Parse repeated ``--compartment "NAME = TERM [op TERM]..."`` definitions.

    Grammar (whitespace-delimited -- operators MUST be surrounded by spaces)::

        NAME = TERM (OP TERM)*
        OP   := '-' (minus) | '&' (intersect) | '+' (union), strictly left to right
        TERM := ATOM ['~' PAD]        PAD in microns; +grows, -shrinks
        ATOM := whole-cell | cell | nucleus | cytosol | nucleolus
              | <a NAME defined earlier>
              | exclusion(CHANNEL,within=PARENT)

    Tokenising on whitespace is deliberate: it is what makes ``whole-cell`` and
    ``cell-boundary`` unambiguous. ``cell - nucleus`` is three tokens and parses;
    ``cell-nucleus`` is one token and is rejected as an unknown region, with a
    message telling the user to put spaces around the operator.

    Returns an ordered dict NAME -> {"terms": [(op, atom, pad_um)], "raw": str}.
    The first term carries op ``None``. Forward references are refused, so
    evaluating definitions in insertion order is always safe (no cycles are
    representable).
    """
    out: dict[str, dict] = {}
    for spec in (specs or []):
        if "=" not in spec:
            raise ValueError(
                f"--compartment {spec!r}: expected 'NAME = EXPRESSION' "
                f"(e.g. --compartment \"cytosol = cell - nucleus\")")
        name, _, expr = spec.partition("=")
        name, expr = name.strip(), expr.strip()
        if not name or not expr:
            raise ValueError(f"--compartment {spec!r}: empty name or expression")
        if name in BUILTIN_COMPARTMENTS:
            raise ValueError(
                f"--compartment {spec!r}: '{name}' is a built-in region and cannot "
                f"be redefined. Built-ins: {', '.join(BUILTIN_COMPARTMENTS)}.")
        if name in out:
            raise ValueError(f"--compartment {spec!r}: '{name}' defined twice")
        if any(c in name for c in " :,=()~"):
            raise ValueError(
                f"--compartment {spec!r}: compartment name {name!r} may not contain "
                f"spaces or any of : , = ( ) ~")

        toks = expr.split()
        if len(toks) % 2 == 0:
            raise ValueError(
                f"--compartment {spec!r}: expression must alternate "
                f"TERM OP TERM ...; got {len(toks)} whitespace-separated tokens. "
                f"Operators need spaces around them ('cell - nucleus', "
                f"not 'cell-nucleus').")
        known = set(BUILTIN_COMPARTMENTS) | set(out)
        terms: list[tuple[str | None, str, float]] = []
        for i, tok in enumerate(toks):
            if i % 2 == 1:
                if tok not in _COMPARTMENT_OPS:
                    raise ValueError(
                        f"--compartment {spec!r}: {tok!r} is not an operator "
                        f"({' '.join(sorted(_COMPARTMENT_OPS))})")
                continue
            atom, pad = tok, 0.0
            if "~" in tok:
                atom, _, pad_s = tok.partition("~")
                try:
                    pad = float(pad_s)
                except ValueError:
                    raise ValueError(
                        f"--compartment {spec!r}: {pad_s!r} after '~' in {tok!r} "
                        f"is not a number of microns") from None
            if atom.startswith("exclusion(") and atom.endswith(")"):
                pass  # validated at resolution time; see _resolve_compartment_expr
            elif atom not in known:
                raise ValueError(
                    f"--compartment {spec!r}: unknown region {atom!r}. "
                    f"Available here: {', '.join(sorted(known))}. "
                    f"(Regions must be defined before they are used, and "
                    f"operators need spaces: 'cell - nucleus', not 'cell-nucleus'.)")
            terms.append((None if i == 0 else toks[i - 1], atom, pad))
        out[name] = {"terms": terms, "raw": spec}
    return out


def _pad_mask_um(mask_bool: np.ndarray, pad_um: float, cfg: dict) -> np.ndarray:
    """Grow (pad_um > 0) or shrink (pad_um < 0) a mask by a physical distance.

    Uses the anisotropic EDT so the pad is uniform in microns rather than in
    voxels -- the same reasoning as `_erode_compartment_to_um`, which handles
    only the shrink direction. In 2D without a voxel size the sampling falls
    back to 1.0, i.e. the pad is in pixels; `validate_compartment_config`
    refuses that case rather than letting it pass silently.
    """
    if not pad_um:
        return mask_bool
    xy = float(cfg.get("voxel_size_xy_um", 1.0) or 1.0)
    z = float(cfg.get("voxel_size_z_um", 1.0) or 1.0)
    sampling = (z, xy, xy) if mask_bool.ndim == 3 else (xy, xy)
    if pad_um > 0:
        return distance_transform_edt(~mask_bool, sampling=sampling) <= pad_um
    return distance_transform_edt(mask_bool, sampling=sampling) >= -pad_um


def _builtin_compartment(
    atom: str, cell_mask: np.ndarray, nuc_mask: np.ndarray, cfg: dict,
    nucleolar_mask: np.ndarray | None,
) -> np.ndarray:
    """The built-in named regions. `cytosol` here is bit-for-bit the historical
    definition: the cell mask minus the GLOBAL union of nuclei dilated by
    cfg['nucleus_dilate_px'] voxels. Not per-cell -- a cell whose mask overlaps
    a neighbour's nucleus has that overlap excluded, and changing this would
    silently move every previously reported cytosol number."""
    if atom in ("cell", "whole-cell"):
        return cell_mask > 0
    if atom == "nucleus":
        return nuc_mask > 0
    if atom == "cytosol":
        radius = int(cfg["nucleus_dilate_px"])
        if radius > 0:
            footprint = (morphology.ball(radius) if cell_mask.ndim == 3
                         else morphology.disk(radius))
            nuc_dil = morphology.dilation((nuc_mask > 0), footprint=footprint)
        else:
            nuc_dil = (nuc_mask > 0)
        return (cell_mask > 0) & (~nuc_dil)
    if atom == "nucleolus":
        if nucleolar_mask is None:
            raise ValueError(
                "compartment expression references 'nucleolus' but no nucleolar "
                "mask exists for this run. Give a channel the 'nucleolus' role, "
                "e.g. \"2:Nsr1:nucleolus\".")
        return nucleolar_mask > 0
    raise ValueError(f"Unknown built-in region: {atom!r}")


def _resolve_compartment_expr(
    name: str, cell_mask: np.ndarray, nuc_mask: np.ndarray, cfg: dict,
    nucleolar_mask: np.ndarray | None, _cache: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    """Evaluate a user-defined compartment, strictly left to right."""
    cache = _cache if _cache is not None else {}
    if name in cache:
        return cache[name]
    defs = cfg.get("compartments") or {}
    terms = defs[name]["terms"]
    result: np.ndarray | None = None
    for op, atom, pad in terms:
        if atom.startswith("exclusion("):
            raise ValueError(
                f"--compartment {defs[name]['raw']!r}: negative-space detection "
                f"(exclusion(...)) is parsed but not yet enabled. The detector "
                f"has no validated defaults: on the yeast test data, dark "
                f"regions below 0.65x the per-cell median are absent from "
                f"maximum-intensity projections entirely, and a single in-focus "
                f"plane yields only ~0.09 candidate objects per cell. Shipping "
                f"defaults calibrated on that would produce plausible masks from "
                f"noise. Calibration against a dedicated vacuole marker "
                f"(FM4-64, Vph1-GFP, CMAC) is required first.")
        if atom in BUILTIN_COMPARTMENTS:
            m = _builtin_compartment(atom, cell_mask, nuc_mask, cfg, nucleolar_mask)
        else:
            m = _resolve_compartment_expr(
                atom, cell_mask, nuc_mask, cfg, nucleolar_mask, cache)
        m = _pad_mask_um(m, pad, cfg)
        if result is None:
            result = m
        elif op == "-":
            result = result & (~m)
        elif op == "&":
            result = result & m
        elif op == "+":
            result = result | m
    cache[name] = result
    return result


def validate_compartment_config(cfg: dict, channels: list[dict]) -> None:
    """Check every compartment REFERENCE resolves, before any segmentation runs.

    Without this a typo in --puncta-compartment raises deep in the per-image
    loop, after Cellpose has already spent minutes. Same reasoning as the
    voxel-size ladder: fail at startup, naming the fix.
    """
    defs = cfg.get("compartments") or {}
    legal = set(BUILTIN_COMPARTMENTS) | set(defs)
    has_nucleolus = any(ch["role"] == "nucleolus" for ch in channels)

    refs = [("--puncta-compartment", cfg.get("puncta_compartment"))]
    if cfg.get("colocalization_compartment"):
        refs.append(("--colocalization-compartment",
                     cfg.get("colocalization_compartment")))
    for ch_name, over in (cfg.get("puncta_params_per_channel") or {}).items():
        if isinstance(over, dict) and over.get("puncta_compartment"):
            refs.append((f"--puncta-params-per-channel {ch_name}",
                         over["puncta_compartment"]))
    for flag, val in refs:
        if val and val not in legal:
            raise SystemExit(
                f"[error] {flag}={val!r} is not a known region.\n"
                f"        Known: {', '.join(sorted(legal))}.\n"
                f"        Define it first with "
                f"--compartment \"{val} = cell - nucleus\".")

    # Micron-denominated pads without a pixel size are caught earlier, in main(),
    # by _um_denominated_2d_params — after the 2D pixel size has been resolved
    # from file metadata, so a file that carries its resolution is not refused.
    for nm, d in defs.items():
        for _, atom, _ in d["terms"]:
            if atom.startswith("exclusion("):
                raise SystemExit(
                    f"[error] --compartment {d['raw']!r}: negative-space detection "
                    f"(exclusion(...)) is parsed but not yet enabled.\n"
                    f"        The detector has no validated defaults. On the yeast "
                    f"test data, dark regions below 0.65x the per-cell median are\n"
                    f"        absent from maximum-intensity projections entirely, "
                    f"and a single in-focus plane yields only ~0.09 candidate\n"
                    f"        objects per cell. Shipping defaults calibrated on "
                    f"that would produce plausible masks out of noise.\n"
                    f"        Calibration against a dedicated vacuole marker "
                    f"(FM4-64, Vph1-GFP, CMAC) is required first.")
            if atom == "nucleolus" and not has_nucleolus:
                raise SystemExit(
                    f"[error] --compartment {d['raw']!r} references 'nucleolus' but "
                    f"no channel has the 'nucleolus' role.\n"
                    f"        Add one, e.g. \"2:Nsr1:nucleolus\".")


def make_compartment_mask(
    cell_mask: np.ndarray,
    nuc_mask: np.ndarray,
    compartment: str,
    cfg: dict,
    *,
    nucleolar_mask: np.ndarray | None = None,
) -> np.ndarray:
    defs = cfg.get("compartments") or {}
    if compartment in defs:
        return _resolve_compartment_expr(
            compartment, cell_mask, nuc_mask, cfg, nucleolar_mask)
    if compartment == "cytosol":
        radius = int(cfg["nucleus_dilate_px"])
        if radius > 0:
            footprint = (morphology.ball(radius) if cell_mask.ndim == 3
                         else morphology.disk(radius))
            nuc_dil = morphology.dilation((nuc_mask > 0), footprint=footprint)
        else:
            nuc_dil = (nuc_mask > 0)
        return (cell_mask > 0) & (~nuc_dil)
    elif compartment == "nucleus":
        return (nuc_mask > 0)
    elif compartment == "whole-cell":
        return (cell_mask > 0)
    elif compartment in ("cell", "nucleolus"):
        return _builtin_compartment(
            compartment, cell_mask, nuc_mask, cfg, nucleolar_mask)
    else:
        legal = sorted(set(BUILTIN_COMPARTMENTS) | set(defs))
        raise ValueError(
            f"Unknown compartment: {compartment!r}. Known: {', '.join(legal)}. "
            f"Define new ones with --compartment \"NAME = cell - nucleus\".")


# ---------------------------------------------------------------------------
# Composite cell segmentation (no-nucleus mode)
# ---------------------------------------------------------------------------
def build_composite_seg_image(
    images: dict[str, np.ndarray],
    channels: list[dict],
    cfg: dict,
) -> np.ndarray:
    """Rescale each non-skip channel to [0,1], sum, renormalize.

    Used when has_nuclei is False and no explicit cell-seg-channel or
    cell-boundary role exists.
    """
    non_skip = [ch for ch in channels if ch["role"] != "skip"]
    if not non_skip:
        raise ValueError("No non-skip channels for composite segmentation image")
    acc = np.zeros_like(next(iter(images.values())), dtype=np.float32)
    for ch in non_skip:
        acc += robust_rescale(images[ch["name"]], cfg["p_low"], cfg["p_high"])
    lo, hi = acc.min(), acc.max()
    if hi > lo:
        acc = (acc - lo) / (hi - lo)
    return acc.astype(np.float32, copy=False)


def filter_cells_by_area(
    cell_mask: np.ndarray, min_area: int, max_area: int
) -> np.ndarray:
    """Zero out cells outside [min_area, max_area]. 0 means disabled."""
    if min_area <= 0 and max_area <= 0:
        return cell_mask
    out = cell_mask.copy()
    for cid in np.unique(out):
        if cid == 0:
            continue
        area = int((out == cid).sum())
        if (min_area > 0 and area < min_area) or (max_area > 0 and area > max_area):
            out[out == cid] = 0
    return out


def find_peak_signal_z(images: dict[str, np.ndarray],
                       channels: list[dict]) -> int:
    """Return the Z slice with the highest integrated intensity summed across
    all non-skip channels. For 3D inputs only — caller must check ndim."""
    summed: np.ndarray | None = None
    for ch in channels:
        if ch["role"] == "skip":
            continue
        arr = images[ch["name"]]
        if arr.ndim != 3:
            continue
        per_z = arr.sum(axis=(1, 2)).astype(np.float64)
        summed = per_z if summed is None else summed + per_z
    if summed is None:
        first = next(iter(images.values()))
        return first.shape[0] // 2 if first.ndim == 3 else 0
    return int(summed.argmax())


def apply_z_crop(images: dict[str, np.ndarray], window_size: int,
                 channels: list[dict]) -> tuple[dict[str, np.ndarray], int, int]:
    """Crop 3D channels to `window_size` slices centered on the per-image
    peak-signal Z. Returns (cropped_images, z_lo, z_hi) — inclusive bounds.

    If window_size >= stack depth, returns the input unchanged."""
    first = next(iter(images.values()))
    if first.ndim != 3:
        return images, 0, 0
    Z = first.shape[0]
    if window_size <= 0 or window_size >= Z:
        return images, 0, Z - 1
    peak_z = find_peak_signal_z(images, channels)
    half = window_size // 2
    z_lo = max(0, peak_z - half)
    z_hi = min(Z - 1, z_lo + window_size - 1)
    # If clipped against the top, slide window down to keep the target size
    if z_hi - z_lo + 1 < window_size:
        z_lo = max(0, z_hi - window_size + 1)
    cropped = {name: arr[z_lo:z_hi + 1] for name, arr in images.items()}
    return cropped, z_lo, z_hi


def exclude_cells_at_z_boundary(cell_mask: np.ndarray) -> tuple[np.ndarray, int]:
    """Zero out cell labels that touch the top (Z=0) or bottom (Z=-1) plane
    of a 3D mask. Use AFTER a z-crop so that cells straddling the crop boundary
    don't contribute truncated half-volumes to the per-cell measurements.

    Returns (filtered_mask, n_excluded_cells)."""
    if cell_mask.ndim != 3 or cell_mask.size == 0:
        return cell_mask, 0
    top_ids = set(int(v) for v in np.unique(cell_mask[0]))
    bot_ids = set(int(v) for v in np.unique(cell_mask[-1]))
    border_ids = (top_ids | bot_ids) - {0}
    if not border_ids:
        return cell_mask, 0
    out = cell_mask.copy()
    border_arr = np.array(sorted(border_ids), dtype=cell_mask.dtype)
    drop = np.isin(out, border_arr)
    out[drop] = 0
    return out, len(border_ids)


# ---------------------------------------------------------------------------
# Nucleolar segmentation
# ---------------------------------------------------------------------------
CELL_SHAPE_EXPECT: dict[Any, float] = {
    "yeast": 1.8,       # near-spherical
    "mammalian": 3.5,   # adherent; can spread / elongate in-plane
    "bacteria": 8.0,    # rods
    None: 3.0,
}


def warn_cell_shape(cell_mask: np.ndarray, cfg: dict) -> None:
    """Warn when segmented cells are implausibly shaped. In 3D this flags two
    specific failure modes — cells flattened along Z with their long axis in the
    imaging plane (the fingerprint of a wrong --voxel-size / anisotropy), and
    gross elongation in any orientation — while deliberately staying silent on
    the mild prolate-along-Z elongation caused by optical axial PSF, which is
    expected and not a segmentation error. In 2D it flags an axis ratio beyond
    the cell-type norm. Read-only heuristic; never raises."""
    try:
        is_3d = cell_mask.ndim == 3
        ct = cfg.get("cell_type")
        max_ratio = CELL_SHAPE_EXPECT.get(ct, CELL_SHAPE_EXPECT[None])
        if is_3d:
            sp = (float(cfg.get("voxel_size_z_um", 1.0) or 1.0),
                  float(cfg.get("voxel_size_xy_um", 1.0) or 1.0),
                  float(cfg.get("voxel_size_xy_um", 1.0) or 1.0))
        ratios: list[float] = []
        zcos: list[float] = []
        for pr in measure.regionprops(cell_mask)[:400]:
            if pr.area < 50:
                continue
            if is_3d:
                T = np.asarray(pr.inertia_tensor)
                C = ((np.trace(T) / 2.0) * np.eye(3) - T) * np.outer(sp, sp)
                w, V = np.linalg.eigh(C)
                w = np.clip(w, 1e-9, None)
                order = np.argsort(w)[::-1]
                semi = np.sqrt(5.0 * w[order])
                ratios.append(float(semi[0] / semi[-1]))
                zcos.append(float(abs(V[0, order[0]])))
            else:
                mn = float(getattr(pr, "axis_minor_length", 0.0) or 0.0)
                mx = float(getattr(pr, "axis_major_length", 0.0) or 0.0)
                if mn > 0:
                    ratios.append(mx / mn)
        if len(ratios) < 5:
            return
        med = float(np.median(ratios))
        if is_3d:
            # Two 3D signatures, measured on the yeast calibration image
            # (25C_series1_rep1, 391 cells, stitch@0.65):
            #   correct voxel 0.10571/0.23 (aniso 2.18): ratio 1.66, |cosZ| 0.98
            #       -> mild PROLATE-along-Z. This is optical axial PSF, is
            #          expected, and must NOT warn.
            #   wrong voxel   0.094/0.10   (aniso 1.06): ratio 1.49, |cosZ| 0.15
            #       -> OBLATE, long axis in the imaging plane: the fingerprint of
            #          an under-scaled Z / wrong --voxel-size anisotropy.
            # So flag oblate cells (low |cosZ|) at a ratio gate (1.3) below the
            # real bug's 1.49, and separately flag gross elongation in any
            # orientation. The orientation flip (|cosZ|), not the ratio
            # magnitude, is what distinguishes the bug from benign PSF.
            medcos = float(np.median(zcos)) if zcos else 1.0
            if med > 1.3 and medcos < 0.3:
                print(f"[warn] cell shape: cells look flattened along Z (median axis "
                      f"ratio {med:.2f}, long-axis |cos Z|={medcos:.2f} — long axis "
                      f"lies in the imaging plane). This is the signature of an "
                      f"under-scaled Z / wrong voxel-size anisotropy (Z spacing too "
                      f"small), not biology. Double-check --voxel-size XY Z against "
                      f"the OME/ImageJ metadata.")
            elif med > 2.5:
                print(f"[warn] cell shape: cells are grossly elongated (median axis "
                      f"ratio {med:.2f}, long-axis |cos Z|={medcos:.2f}). Real cells "
                      f"seldom exceed ~2.5x; check segmentation, --stitch-threshold, "
                      f"and the cell-type preset. (Mild ~1.5-1.7x prolate-along-Z "
                      f"elongation is expected optical PSF and is not flagged.)")
        elif med > max_ratio:
            print(f"[warn] cell shape: median cell axis ratio {med:.2f} exceeds the "
                  f"~{max_ratio:.1f} expected for cell-type '{ct}'. Check "
                  f"segmentation and the cell-type preset.")
    except Exception:
        return


def segment_nucleoli(
    image: np.ndarray, cell_mask: np.ndarray, cfg: dict
) -> np.ndarray:
    """Per-cell Otsu thresholding of nucleolus channel, morphological opening,
    size filtering. Returns combined binary mask (0/1) for the whole image.
    Works for both 2D and 3D images."""
    is_3d = image.ndim == 3
    opening_radius = int(cfg.get("nucleolar_opening_radius", 1))
    if is_3d:
        min_size = int(cfg.get("nucleolar_min_volume_vox", 10))
        footprint = morphology.ball(opening_radius) if opening_radius > 0 else None
    else:
        min_size = int(cfg.get("nucleolar_min_area_px", 5))
        footprint = morphology.disk(opening_radius) if opening_radius > 0 else None
    combined = np.zeros(image.shape, dtype=np.uint8)

    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    for cid in cell_ids:
        roi = (cell_mask == cid)
        vals = image[roi]
        if vals.size < 20:
            continue
        try:
            thr = filters.threshold_otsu(vals)
        except Exception:
            continue
        bw = (image > thr) & roi
        if footprint is not None:
            bw = morphology.opening(bw, footprint=footprint)
        bw = bw & roi  # clip back to cell
        # Size filter (area in 2D, volume in 3D)
        lab = measure.label(bw)
        for prop in measure.regionprops(lab):
            if prop.area >= min_size:
                combined[lab == prop.label] = 1

    return combined


def compute_nucleolar_morphology(
    nucleolar_mask: np.ndarray,
    cell_mask: np.ndarray,
    metadata: dict[str, str],
    image_name: str,
    cfg: dict | None = None,
) -> pd.DataFrame:
    """Per-cell nucleolar morphometrics on the largest nucleolus.

    2D: area + circularity + solidity + eccentricity (paper-validated).
    3D: volume (vox + µm³) + equivalent diameter (µm) with anisotropic spacing.
    """
    is_3d = nucleolar_mask.ndim == 3
    rows: list[dict[str, Any]] = []
    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]
    base = {
        "image": image_name,
        "condition": metadata.get("condition", ""),
        "replicate": metadata.get("replicate", ""),
    }

    if is_3d:
        cfg = cfg or {}
        spacing = (
            cfg.get("voxel_size_z_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
        )
        voxel_um3 = spacing[0] * spacing[1] * spacing[2]

    for cid in cell_ids:
        cell_bin = (cell_mask == cid)
        nol_in_cell = nucleolar_mask * cell_bin
        lab = measure.label(nol_in_cell > 0)
        props = measure.regionprops(lab)
        n_nucleoli = len(props)
        row = dict(base, cell_id=int(cid), n_nucleoli=n_nucleoli)

        if n_nucleoli == 0:
            if is_3d:
                row.update(
                    nucleolar_volume_vox=0,
                    nucleolar_volume_um3=0.0,
                    nucleolar_eq_diameter_um=float("nan"),
                )
            else:
                row.update(
                    nucleolar_area=0,
                    nucleolar_solidity=float("nan"),
                    nucleolar_circularity=float("nan"),
                    nucleolar_eccentricity=float("nan"),
                )
            rows.append(row)
            continue

        largest = max(props, key=lambda p: p.area)

        if is_3d:
            vol_vox = int(largest.area)
            vol_um3 = vol_vox * voxel_um3
            eq_diam_um = (6.0 * vol_um3 / np.pi) ** (1.0 / 3.0)
            row.update(
                nucleolar_volume_vox=vol_vox,
                nucleolar_volume_um3=float(vol_um3),
                nucleolar_eq_diameter_um=float(eq_diam_um),
            )
        else:
            per = measure.perimeter_crofton(largest.image, directions=4)
            circ = (4.0 * np.pi * largest.area) / (per * per) if per > 0 else float("nan")
            sol = (largest.area / largest.area_convex
                   if largest.area_convex > 0 else float("nan"))
            row.update(
                nucleolar_area=int(largest.area),
                nucleolar_solidity=float(sol),
                nucleolar_circularity=float(circ),
                nucleolar_eccentricity=float(largest.eccentricity),
            )
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Puncta detection (parameterized version of reference)
# ---------------------------------------------------------------------------
def _erode_compartment_to_um(
    mask_bool: np.ndarray,
    erode_um: float,
    voxel_xy_um: float,
    voxel_z_um: float,
) -> np.ndarray:
    """Erode a boolean mask inward by `erode_um` microns using anisotropic EDT.

    Inside the input mask, each voxel's distance to the nearest non-mask
    voxel is computed using the supplied voxel sizes as EDT sampling. Voxels
    closer to the boundary than `erode_um` are removed. This produces an
    erosion that is physically uniform (in microns) rather than uniform in
    voxel counts — important when XY and Z voxel sizes differ.

    Used by `detect_puncta` to remove the cell-membrane gradient from the
    LoG detection ROI. See `--puncta-compartment-erode-um`.
    """
    if erode_um <= 0:
        return mask_bool
    is_3d = mask_bool.ndim == 3
    if is_3d:
        sampling = (voxel_z_um, voxel_xy_um, voxel_xy_um)
    else:
        sampling = (voxel_xy_um, voxel_xy_um)
    dist_inside = distance_transform_edt(mask_bool, sampling=sampling)
    return dist_inside >= erode_um


def detect_puncta(
    image: np.ndarray, compartment_mask: np.ndarray, cfg: dict
) -> np.ndarray:
    """LoG-based puncta detection. Auto-dispatches to 2D or 3D.

    2D: paper-validated path — Gaussian → Laplacian, scalar sigma, area filter,
        optional roundness/solidity filter.
    3D: anisotropic gaussian_laplace with sigma_z = sigma_xy · (xy_um / z_um)
        so the kernel is isotropic in physical space; volume filter; no
        roundness filter (sphericity-based filtering would belong here but
        adds cost — the volume cap is the dominant control).
    """
    is_3d = image.ndim == 3

    # Physical-unit compartment erosion before LoG. Removes the cell-membrane
    # gradient that the LoG-Otsu detector otherwise classifies as outline-
    # shaped "puncta." Anisotropic in 3D. Guard against over-erosion that
    # would leave an effectively empty compartment.
    erode_um = float(cfg.get("puncta_compartment_erode_um", 0.0) or 0.0)
    if erode_um > 0:
        voxel_xy = float(cfg.get("voxel_size_xy_um", 1.0) or 1.0)
        voxel_z = float(cfg.get("voxel_size_z_um", 1.0) or 1.0)
        comp_bool = compartment_mask > 0
        eroded = _erode_compartment_to_um(
            comp_bool, erode_um, voxel_xy, voxel_z)
        # Over-erosion guard: if erosion drops the compartment by >95% or
        # leaves fewer than 100 voxels, skip erosion and warn rather than
        # producing an empty compartment.
        before = int(comp_bool.sum())
        after = int(eroded.sum())
        if after < 100 or (before > 0 and after / before < 0.05):
            print(f"[warn] puncta_compartment_erode_um={erode_um} would shrink "
                  f"compartment from {before} to {after} voxels; "
                  f"skipping erosion for this image")
        else:
            compartment_mask = eroded.astype(compartment_mask.dtype)

    # Optional per-channel rolling-ball background subtraction (white tophat
    # with a disk/ball structuring element — fast and standard for punctate
    # fluorescence). Applied before LoG so the LoG response reflects local
    # contrast over a sliding window rather than absolute intensity. Skipped
    # when radius is 0 (the default).
    rb_radius = int(cfg.get("puncta_rolling_ball_radius", 0))
    if rb_radius > 0:
        footprint = (morphology.ball(rb_radius) if is_3d
                     else morphology.disk(rb_radius))
        image = morphology.white_tophat(image, footprint=footprint)

    g = robust_rescale(image, cfg["p_low"], cfg["p_high"])
    g = g * (compartment_mask > 0)

    if is_3d:
        sigma_xy = float(cfg["log_sigma"])
        z_um = cfg.get("voxel_size_z_um", 1.0)
        xy_um = cfg.get("voxel_size_xy_um", 1.0)
        # Z sigma in voxels = sigma_xy_voxels · (xy_um / z_um) — so the kernel
        # is isotropic in microns.
        sigma_z = max(sigma_xy * (xy_um / z_um), 0.5)
        log_im = -gaussian_laplace(g, sigma=(sigma_z, sigma_xy, sigma_xy))
    else:
        log_im = -filters.laplace(filters.gaussian(g, sigma=cfg["log_sigma"]))

    log_im = robust_rescale(log_im, 1.0, 99.9)

    vals = log_im[compartment_mask > 0]
    if vals.size < 50:
        return np.zeros_like(compartment_mask, dtype=np.int32)

    method = str(cfg["puncta_threshold_method"]).lower()
    if method == "otsu":
        thr = filters.threshold_otsu(vals)
    elif method == "triangle":
        thr = filters.threshold_triangle(vals)
    elif method == "fixed":
        thr = float(cfg["puncta_threshold_fixed"])
    else:
        raise ValueError(
            f"Unknown puncta_threshold_method: {cfg['puncta_threshold_method']}")

    bw = (log_im >= thr) & (compartment_mask > 0)
    if is_3d:
        min_size = int(cfg.get("puncta_min_volume_vox", 4))
        max_size = int(cfg.get("puncta_max_volume_vox", 100_000))
    else:
        min_size = int(cfg["puncta_min_area_px"])
        max_size = int(cfg["puncta_max_area_px"])
        bw = remove_small_holes_compat(bw, area_threshold=8)
    bw = remove_small_objects_compat(bw, min_size=min_size)

    lab = measure.label(bw)
    if lab.max() == 0:
        return lab.astype(np.int32, copy=False)

    keep = np.zeros(lab.max() + 1, dtype=bool)
    do_round = (not is_3d) and bool(cfg.get("puncta_filter_round", False))
    min_circ = float(cfg.get("puncta_min_circularity", 0.0))
    min_sol = float(cfg.get("puncta_min_solidity", 0.0))

    for prop in measure.regionprops(lab):
        if not (min_size <= prop.area <= max_size):
            continue
        if do_round:
            per = measure.perimeter_crofton(prop.image, directions=4)
            if per <= 0:
                continue
            circ = (4.0 * np.pi * prop.area) / (per * per)
            if circ < min_circ:
                continue
            if hasattr(prop, "solidity") and prop.solidity < min_sol:
                continue
        keep[prop.label] = True

    lab2 = keep[lab] * lab
    lab2 = measure.label(lab2 > 0)
    return lab2.astype(np.int32, copy=False)


# ---------------------------------------------------------------------------
# Fragmentation index (NEW; works in 2D and 3D natively)
# ---------------------------------------------------------------------------
def compute_fragmentation_indices(
    image: np.ndarray,
    cell_mask: np.ndarray,
    cfg: dict,
) -> dict[int, tuple[int, float]]:
    """For each cell, compute two fragmentation metrics on the LoG response:

      simple = number of connected components above a single LoG-Otsu threshold
               within the cell. Threshold-dependent but interpretable; matches
               the spirit of the paper's puncta count.

      persistence = threshold-free integral over a swept threshold range.
                    Sweep K thresholds spanning the LoG response distribution
                    inside the cell; record component count at each level;
                    return the area under the count-vs-threshold curve,
                    normalized to [0, 1] by the swept range.
                    A uniform cell decays from 1→0 quickly and gives small
                    values; a fragmented cell sustains many components across
                    levels and gives large values.

    Returns: {cell_id: (simple_count, persistence_index)}
    """
    is_3d = image.ndim == 3
    n_levels = max(int(cfg.get("fragmentation_thresholds", 10)), 0)

    # Compute the LoG response over the whole image once
    g = robust_rescale(image, cfg["p_low"], cfg["p_high"])
    if is_3d:
        sigma_xy = float(cfg["log_sigma"])
        z_um = cfg.get("voxel_size_z_um", 1.0)
        xy_um = cfg.get("voxel_size_xy_um", 1.0)
        sigma_z = max(sigma_xy * (xy_um / z_um), 0.5)
        log_im = -gaussian_laplace(g, sigma=(sigma_z, sigma_xy, sigma_xy))
    else:
        log_im = -filters.laplace(filters.gaussian(g, sigma=cfg["log_sigma"]))
    log_im = robust_rescale(log_im, 1.0, 99.9)

    out: dict[int, tuple[int, float]] = {}
    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    for cid in cell_ids:
        cell_bin = (cell_mask == cid)
        vals = log_im[cell_bin]
        if vals.size < 50:
            out[int(cid)] = (0, 0.0)
            continue

        # Simple variant
        try:
            thr = float(filters.threshold_otsu(vals))
        except Exception:
            out[int(cid)] = (0, 0.0)
            continue
        bw = (log_im >= thr) & cell_bin
        lab = measure.label(bw)
        n_simple = int(lab.max())

        # Persistence variant
        if n_levels < 2:
            out[int(cid)] = (n_simple, float("nan"))
            continue
        lo = float(np.percentile(vals, 50))
        hi = float(np.percentile(vals, 99))
        if hi <= lo:
            out[int(cid)] = (n_simple, 0.0)
            continue
        levels = np.linspace(lo, hi, n_levels)
        counts = []
        for t in levels:
            bw_t = (log_im >= t) & cell_bin
            counts.append(int(measure.label(bw_t).max()))
        counts_arr = np.asarray(counts, dtype=np.float64)
        # Normalize: divide integral by (max possible count) × (level range)
        # to get a unitless [0, 1]-ish index. We use the trapezoidal integral
        # of count vs. normalized threshold.
        denom = max(counts_arr.max(), 1.0)
        norm_t = (levels - lo) / (hi - lo)
        persistence = float(np.trapz(counts_arr / denom, norm_t))
        out[int(cid)] = (n_simple, persistence)

    return out


# ---------------------------------------------------------------------------
# Quantification — dynamic per-channel metrics
# ---------------------------------------------------------------------------
def per_cell_metrics(
    image_name: str,
    images: dict[str, np.ndarray],
    cell_mask: np.ndarray,
    nuc_mask: np.ndarray,
    puncta_masks: dict[str, np.ndarray],
    cell_to_nucs: dict[int, list[int]],
    metadata: dict[str, str],
    channels: list[dict],
    cfg: dict,
    fragmentation_per_channel: dict[str, dict[int, tuple[int, float]]] | None = None,
    nucleolar_mask: np.ndarray | None = None,
) -> pd.DataFrame:
    is_3d = cell_mask.ndim == 3
    rows: list[dict[str, Any]] = []
    nuc_union = (nuc_mask > 0)
    radius = int(cfg["nucleus_dilate_px"])
    if radius > 0:
        footprint = morphology.ball(radius) if is_3d else morphology.disk(radius)
        nuc_dil = morphology.dilation(nuc_union, footprint=footprint)
    else:
        nuc_dil = nuc_union

    voxel_um3 = (
        cfg.get("voxel_size_z_um", 1.0)
        * cfg.get("voxel_size_xy_um", 1.0)
        * cfg.get("voxel_size_xy_um", 1.0)
    ) if is_3d else 1.0

    non_skip = [ch for ch in channels if ch["role"] != "skip"]
    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]
    compartment = cfg.get("puncta_compartment", "cytosol")
    fragmentation_per_channel = fragmentation_per_channel or {}

    # User-defined regions are global masks; resolve each once, then intersect
    # with the cell inside the loop. Sorted for a stable cells.csv column order.
    user_comps: dict[str, np.ndarray] = {}
    for _cname in sorted(cfg.get("compartments") or {}):
        user_comps[_cname] = make_compartment_mask(
            cell_mask, nuc_mask, _cname, cfg, nucleolar_mask=nucleolar_mask)

    for cid in cell_ids:
        cell_bin = (cell_mask == cid)
        if int(cell_bin.sum()) < 50:
            continue

        nuc_ids = cell_to_nucs.get(int(cid), [])
        n_nuc = len(nuc_ids)

        nuc_bin = np.zeros_like(cell_bin, dtype=bool)
        for nid in nuc_ids:
            nuc_bin |= (nuc_mask == nid)

        cyt_bin = cell_bin & (~nuc_dil)

        row: dict[str, Any] = {
            "image": image_name,
            "condition": metadata.get("condition", ""),
            "replicate": metadata.get("replicate", ""),
            "cell_id": int(cid),
            "n_nuclei": n_nuc,
        }
        if is_3d:
            cell_vol = int(cell_bin.sum())
            nuc_vol = int(nuc_bin.sum())
            cyt_vol = int(cyt_bin.sum())
            row.update(
                cell_volume_vox=cell_vol,
                cell_volume_um3=cell_vol * voxel_um3,
                nucleus_volume_vox=nuc_vol,
                nucleus_volume_um3=nuc_vol * voxel_um3,
                cytosol_volume_vox=cyt_vol,
                cytosol_volume_um3=cyt_vol * voxel_um3,
            )
        else:
            row.update(
                cell_area_px=int(cell_bin.sum()),
                nucleus_area_px=int(nuc_bin.sum()),
                cytosol_area_px=int(cyt_bin.sum()),
            )
        for _cname, _cmask in user_comps.items():
            _sz = int((_cmask & cell_bin).sum())
            if is_3d:
                row[f"{_cname}_volume_vox"] = _sz
                row[f"{_cname}_volume_um3"] = _sz * voxel_um3
            else:
                row[f"{_cname}_area_px"] = _sz

        # Per non-skip channel: mean intensities + (opt-in) condensate index
        want_ci = bool(cfg.get("condensate_index", False))
        for ch in non_skip:
            name = ch["name"]
            img = images[name]
            row[f"{name}_cell_mean"] = safe_mean(img, cell_bin)
            row[f"{name}_nucleus_mean"] = safe_mean(img, nuc_bin)
            row[f"{name}_cytosol_mean"] = safe_mean(img, cyt_bin)
            for _cname, _cmask in user_comps.items():
                row[f"{name}_{_cname}_mean"] = safe_mean(img, _cmask & cell_bin)

            if want_ci:
                m_cell = row[f"{name}_cell_mean"]
                p95_cell = _safe_p95(img, cell_bin)
                row[f"{name}_condensate_index_cell"] = _condensate_index(p95_cell, m_cell)

                m_cyt = row[f"{name}_cytosol_mean"]
                p95_cyt = _safe_p95(img, cyt_bin)
                row[f"{name}_condensate_index_cytosol"] = _condensate_index(p95_cyt, m_cyt)

        # Per puncta channel: puncta-specific metrics + fragmentation
        for ch_name, p_mask in puncta_masks.items():
            ch_img = images[ch_name]
            puncta_in_cell = p_mask * cell_bin
            puncta_bin = (puncta_in_cell > 0)
            n_puncta = len(np.unique(puncta_in_cell)) - 1  # subtract bg 0

            ch_puncta_vals = ch_img[puncta_bin]
            puncta_int = float(ch_puncta_vals.sum()) if ch_puncta_vals.size else 0.0

            if compartment == "cytosol":
                comp_bin = cyt_bin
            elif compartment == "nucleus":
                comp_bin = nuc_bin
            else:
                comp_bin = cell_bin
            comp_total = float(ch_img[comp_bin].sum()) if np.any(comp_bin) else 0.0
            diffuse_bin = comp_bin & (~puncta_bin)

            row[f"{ch_name}_puncta_n"] = n_puncta
            puncta_size = int(puncta_bin.sum())
            if is_3d:
                row[f"{ch_name}_puncta_volume_vox"] = puncta_size
                row[f"{ch_name}_puncta_volume_um3"] = puncta_size * voxel_um3
            else:
                row[f"{ch_name}_puncta_area_px"] = puncta_size
            row[f"{ch_name}_puncta_mean_intensity"] = (
                float(ch_puncta_vals.mean()) if ch_puncta_vals.size else np.nan)
            row[f"{ch_name}_diffuse_mean_intensity"] = safe_mean(ch_img, diffuse_bin)
            row[f"{ch_name}_frac_intensity_in_puncta"] = (
                puncta_int / comp_total if comp_total > 0 else np.nan)
            row[f"{ch_name}_puncta_integrated_intensity"] = puncta_int

            # Fragmentation indices (NEW; computed in compute_fragmentation_indices)
            frag = fragmentation_per_channel.get(ch_name, {}).get(int(cid))
            if frag is not None:
                simple, persistence = frag
                row[f"{ch_name}_fragmentation_index_simple"] = int(simple)
                row[f"{ch_name}_fragmentation_index_persistence"] = float(persistence)
            else:
                row[f"{ch_name}_fragmentation_index_simple"] = 0
                row[f"{ch_name}_fragmentation_index_persistence"] = float("nan")

        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        df["keep"] = df["n_nuclei"].between(
            int(cfg["keep_min_nuclei"]), int(cfg["keep_max_nuclei"]))
    return df


def per_image_summary(
    cells_df: pd.DataFrame, puncta_chs: list[dict]
) -> pd.DataFrame:
    if cells_df.empty:
        summary: dict[str, Any] = {"n_cells": 0, "n_keep": 0}
        for ch in puncta_chs:
            summary[f"median_{ch['name']}_puncta_n"] = np.nan
            summary[f"median_{ch['name']}_frac_intensity_in_puncta"] = np.nan
        return pd.DataFrame([summary])

    keep = cells_df[cells_df["keep"]].copy()
    summary = {
        "n_cells": int(len(cells_df)),
        "n_keep": int(len(keep)),
    }
    for ch in puncta_chs:
        col_n = f"{ch['name']}_puncta_n"
        col_frac = f"{ch['name']}_frac_intensity_in_puncta"
        summary[f"median_{ch['name']}_puncta_n"] = (
            float(keep[col_n].median())
            if len(keep) and col_n in keep.columns else np.nan)
        summary[f"median_{ch['name']}_frac_intensity_in_puncta"] = (
            float(keep[col_frac].median())
            if len(keep) and col_frac in keep.columns else np.nan)
    return pd.DataFrame([summary])


# ---------------------------------------------------------------------------
# Colocalization module
# ---------------------------------------------------------------------------
def costes_threshold(ch_a: np.ndarray, ch_b: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """Costes auto-threshold: linear regression between channels, step threshold
    down along regression line until sub-threshold Pearson R <= 0.
    Falls back to Otsu on non-convergence. Returns (thr_a, thr_b)."""
    a_vals = ch_a[mask].astype(np.float64)
    b_vals = ch_b[mask].astype(np.float64)
    if a_vals.size < 10:
        return (0.0, 0.0)

    # Linear regression b = slope*a + intercept
    a_mean, b_mean = a_vals.mean(), b_vals.mean()
    cov = ((a_vals - a_mean) * (b_vals - b_mean)).mean()
    var_a = ((a_vals - a_mean) ** 2).mean()
    if var_a == 0:
        return (0.0, 0.0)
    slope = cov / var_a
    intercept = b_mean - slope * a_mean

    # Step threshold from max(a) down
    a_max = float(a_vals.max())
    n_steps = 256
    for step in range(n_steps):
        thr_a = a_max * (1.0 - step / n_steps)
        thr_b = slope * thr_a + intercept
        below = (a_vals < thr_a) | (b_vals < thr_b)
        a_sub = a_vals[below]
        b_sub = b_vals[below]
        if a_sub.size < 5:
            continue
        a_sub_m = a_sub.mean()
        b_sub_m = b_sub.mean()
        a_std = a_sub.std()
        b_std = b_sub.std()
        if a_std == 0 or b_std == 0:
            continue
        r = ((a_sub - a_sub_m) * (b_sub - b_sub_m)).mean() / (a_std * b_std)
        if r <= 0:
            return (thr_a, max(thr_b, 0.0))

    # Fallback to Otsu
    try:
        thr_a = float(filters.threshold_otsu(a_vals))
        thr_b = float(filters.threshold_otsu(b_vals))
    except Exception:
        thr_a, thr_b = 0.0, 0.0
    return (thr_a, thr_b)


def compute_colocalization(
    images: dict[str, np.ndarray],
    coloc_channels: list[dict],
    cell_mask: np.ndarray,
    nuc_mask: np.ndarray,
    compartment: str,
    cfg: dict,
    metadata: dict[str, str],
    image_name: str,
) -> pd.DataFrame:
    """Per cell, per pair of coloc channels: Pearson R, Manders M1/M2."""
    rows: list[dict[str, Any]] = []
    pairs = list(combinations(coloc_channels, 2))
    if not pairs:
        return pd.DataFrame()

    comp_mask = make_compartment_mask(cell_mask, nuc_mask, compartment, cfg)

    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    for ch_a, ch_b in pairs:
        img_a = images[ch_a["name"]].astype(np.float64)
        img_b = images[ch_b["name"]].astype(np.float64)
        pair_name = f"{ch_a['name']}_vs_{ch_b['name']}"

        # Compute Costes thresholds on whole-image masked region
        whole_mask = comp_mask > 0
        thr_a, thr_b = costes_threshold(img_a, img_b, whole_mask)

        for cid in cell_ids:
            cell_comp = (cell_mask == cid) & whole_mask
            if cell_comp.sum() < 10:
                continue

            a_vals = img_a[cell_comp]
            b_vals = img_b[cell_comp]

            # Pearson R
            a_m, b_m = a_vals.mean(), b_vals.mean()
            a_std, b_std = a_vals.std(), b_vals.std()
            if a_std > 0 and b_std > 0:
                pearson_r = float(
                    ((a_vals - a_m) * (b_vals - b_m)).mean() / (a_std * b_std))
            else:
                pearson_r = np.nan

            # Manders M1 and M2 using Costes thresholds
            above_a = a_vals > thr_a
            above_b = b_vals > thr_b
            coloc_ab = above_a & above_b

            total_a = a_vals[above_a].sum()
            total_b = b_vals[above_b].sum()
            m1 = float(a_vals[coloc_ab].sum() / total_a) if total_a > 0 else np.nan
            m2 = float(b_vals[coloc_ab].sum() / total_b) if total_b > 0 else np.nan

            rows.append({
                "image": image_name,
                "condition": metadata.get("condition", ""),
                "replicate": metadata.get("replicate", ""),
                "cell_id": int(cid),
                "pair": pair_name,
                "pearson_r": pearson_r,
                "manders_m1": m1,
                "manders_m2": m2,
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Nucleolar proximity analysis
# ---------------------------------------------------------------------------
def compute_nucleolar_proximity(
    puncta_mask: np.ndarray,
    nucleolar_mask: np.ndarray,
    cell_mask: np.ndarray,
    threshold: float,
    metadata: dict[str, str],
    image_name: str,
    channel_name: str,
    cfg: dict | None = None,
    extra_thresholds: list[float] | None = None,
) -> pd.DataFrame:
    """Per-cell puncta distance to nearest nucleolar boundary.

    2D: Euclidean distance in pixels; threshold compared in pixels.
    3D: Anisotropic EDT with sampling=(z_um, y_um, x_um); distance is in
        microns; threshold (when invoked with cfg) is the µm threshold.

    extra_thresholds (µm) add sensitivity columns fraction_proximal_<t>um
    alongside the primary fraction_proximal, so robustness to the exact
    threshold can be checked from a single run.
    """
    is_3d = puncta_mask.ndim == 3
    rows: list[dict[str, Any]] = []
    cfg = cfg or {}
    extra_cols = [(float(t), f"fraction_proximal_{float(t):g}um")
                  for t in (extra_thresholds or [])]

    if is_3d:
        spacing = (
            cfg.get("voxel_size_z_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
        )
        dist_field = distance_transform_edt(~(nucleolar_mask > 0), sampling=spacing)
    else:
        dist_field = distance_transform_edt(~(nucleolar_mask > 0))

    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    # --- Vectorized punctum-centroid distances -----------------------------
    # One pass over the nonzero puncta voxels replaces the old double loop
    # (re-scanning the whole volume once per cell and again per punctum).
    #
    # Voxels are grouped by (cell_id, punctum_id) *pairs*, not by punctum alone:
    # a punctum straddling two cells therefore contributes a separate centroid
    # to each cell — exactly what the old `puncta_mask * (cell_mask == cid)`
    # masking computed. Collapsing to punctum_id would change the result.
    #
    # This is bit-identical to the old int(round(coords.mean())) per punctum.
    # The reason is exact arithmetic, NOT any traversal-order guarantee:
    # voxel coordinates are integer-valued in float64, and a per-punctum
    # coordinate sum is ~1e7 at most (a few thousand voxels times indices of a
    # few thousand), far below 2**53. Below that bound every partial sum is
    # exactly representable, so np.mean's pairwise summation and np.bincount's
    # sequential accumulation yield the identical float64 sum regardless of the
    # order voxels are visited. The divide by the (identical) voxel count is a
    # single IEEE op on identical operands, and Python's round() and np.rint()
    # both round half to even — so the rounded integer centroid is the same.
    ndim = puncta_mask.ndim
    vox = np.nonzero(puncta_mask)
    pid_all = puncta_mask[vox].astype(np.int64)
    cid_all = cell_mask[vox].astype(np.int64)
    keep = cid_all != 0
    coords = np.stack([axis[keep] for axis in vox], axis=1).astype(np.int64)
    keys = np.stack([cid_all[keep], pid_all[keep]], axis=1)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = inv.ravel()
    counts = np.bincount(inv, minlength=len(uniq))
    centroid = np.empty((len(uniq), ndim), dtype=np.int64)
    for a in range(ndim):
        axis_sum = np.bincount(inv, weights=coords[:, a], minlength=len(uniq))
        centroid[:, a] = np.rint(axis_sum / counts).astype(np.int64)
    if len(uniq):
        group_dist = dist_field[tuple(centroid[:, a] for a in range(ndim))]
    else:
        group_dist = np.empty(0, dtype=dist_field.dtype)

    # uniq is lexicographically sorted by (cell_id, punctum_id), so iterating
    # groups in order and appending preserves the ascending-punctum-id order
    # the old np.unique(punct_ids) inner loop produced for each cell.
    dists_by_cell: dict[int, list[float]] = {}
    for g in range(len(uniq)):
        dists_by_cell.setdefault(int(uniq[g, 0]), []).append(float(group_dist[g]))

    for cid in cell_ids:
        distances = dists_by_cell.get(int(cid))
        if not distances:
            row = {
                "image": image_name,
                "condition": metadata.get("condition", ""),
                "replicate": metadata.get("replicate", ""),
                "cell_id": int(cid),
                "channel": channel_name,
                "n_puncta": 0,
                "mean_distance": np.nan,
                "min_distance": np.nan,
                "fraction_proximal": np.nan,
            }
            for _, col in extra_cols:
                row[col] = np.nan
            rows.append(row)
            continue

        distances_arr = np.array(distances)
        row = {
            "image": image_name,
            "condition": metadata.get("condition", ""),
            "replicate": metadata.get("replicate", ""),
            "cell_id": int(cid),
            "channel": channel_name,
            "n_puncta": len(distances),
            "mean_distance": float(distances_arr.mean()),
            "min_distance": float(distances_arr.min()),
            "fraction_proximal": float((distances_arr <= threshold).mean()),
        }
        for t, col in extra_cols:
            row[col] = float((distances_arr <= t).mean())
        rows.append(row)

    return pd.DataFrame(rows)


def compute_puncta_distance_colors(
    puncta_mask: np.ndarray,
    nucleolar_mask: np.ndarray,
    threshold: float,
    cfg: dict | None = None,
) -> dict[int, str]:
    """Map each punctum label to 'red' (proximal) or 'blue' (distal).

    Distance is in pixels (2D) or microns (3D, with cfg providing voxel size).
    """
    is_3d = puncta_mask.ndim == 3
    cfg = cfg or {}
    if is_3d:
        spacing = (
            cfg.get("voxel_size_z_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
        )
        dist_field = distance_transform_edt(~(nucleolar_mask > 0), sampling=spacing)
    else:
        dist_field = distance_transform_edt(~(nucleolar_mask > 0))
    colors: dict[int, str] = {}
    punct_ids = np.unique(puncta_mask)
    punct_ids = punct_ids[punct_ids != 0]
    for pid in punct_ids:
        coords = np.where(puncta_mask == pid)
        centroid = tuple(int(round(c.mean())) for c in coords)
        d = float(dist_field[centroid])
        colors[int(pid)] = "red" if d <= threshold else "blue"
    return colors


# ---------------------------------------------------------------------------
# QC rendering — dynamic channels
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# QC rendering — dynamic channels (2D path: paper-validated; 3D path: MIP + middle Z)
# ---------------------------------------------------------------------------
def _save_qc_2d(
    out_png, images, channels, cell_mask, nuc_mask, puncta_masks, cfg,
    nucleolar_mask=None, puncta_distance_colors=None,
):
    """Original 2D QC. Single-panel RGB composite with contour overlays."""
    quantify_chs = get_quantify_channels(channels)
    nuc_ch = get_nucleus_channel(channels)
    ref_img = next(iter(images.values()))
    zeros = np.zeros_like(ref_img, dtype=np.float32)

    r = (robust_rescale(images[quantify_chs[0]["name"]], cfg["p_low"], cfg["p_high"])
         if len(quantify_chs) >= 1 else zeros)
    g = (robust_rescale(images[quantify_chs[1]["name"]], cfg["p_low"], cfg["p_high"])
         if len(quantify_chs) >= 2 else zeros)
    b = (robust_rescale(images[nuc_ch["name"]], cfg["p_low"], cfg["p_high"])
         if nuc_ch is not None else zeros)

    f = int(cfg["qc_downsample"])
    rgb = np.stack(
        [downsample_nn(r, f), downsample_nn(g, f), downsample_nn(b, f)], axis=-1)
    rgb = np.clip(rgb, 0.0, 1.0)

    cell_ds = downsample_nn(cell_mask, f)
    nuc_ds = downsample_nn(nuc_mask, f)
    h, w = rgb.shape[:2]
    dpi = int(cfg["qc_dpi"])
    fig, ax = plt.subplots(1, 1, figsize=(w / dpi, h / dpi), dpi=dpi)
    ax.imshow(rgb); ax.set_axis_off()

    for cid in np.unique(cell_ds):
        if cid == 0: continue
        ax.contour(cell_ds == cid, levels=[0.5], colors="cyan", linewidths=1.0)
    if nuc_mask.max() > 0:
        for nid in np.unique(nuc_ds):
            if nid == 0: continue
            ax.contour(nuc_ds == nid, levels=[0.5], colors="yellow", linewidths=0.8)
    if nucleolar_mask is not None and nucleolar_mask.max() > 0:
        nol_ds = downsample_nn(nucleolar_mask, f)
        ax.contour(nol_ds > 0, levels=[0.5], colors="white", linewidths=0.6)
    for ch_name, pmask in puncta_masks.items():
        pmask_ds = downsample_nn(pmask, f)
        if pmask_ds.max() == 0:
            continue
        if puncta_distance_colors and ch_name in puncta_distance_colors:
            colors_map = puncta_distance_colors[ch_name]
            for pid in np.unique(pmask_ds):
                if pid == 0: continue
                color = colors_map.get(int(pid), "magenta")
                ax.contour(pmask_ds == pid, levels=[0.5],
                           colors=color, linewidths=0.6)
        else:
            ax.contour(pmask_ds > 0, levels=[0.5], colors="magenta", linewidths=0.6)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0, dpi=dpi)
    plt.close(fig)


def _save_qc_3d(
    out_png, images, channels, cell_mask, nuc_mask, puncta_masks, cfg,
    nucleolar_mask=None, puncta_distance_colors=None,
):
    """3D QC: side-by-side MIP composite + middle-Z slice with mask outlines."""
    quantify_chs = get_quantify_channels(channels)
    nuc_ch = get_nucleus_channel(channels)
    ref = next(iter(images.values()))
    Z = ref.shape[0]
    z_mid = Z // 2

    def make_rgb(getter):
        proj_shape = getter(ref).shape
        zeros = np.zeros(proj_shape, dtype=np.float32)
        r = (robust_rescale(getter(images[quantify_chs[0]["name"]]),
                            cfg["p_low"], cfg["p_high"])
             if len(quantify_chs) >= 1 else zeros)
        g = (robust_rescale(getter(images[quantify_chs[1]["name"]]),
                            cfg["p_low"], cfg["p_high"])
             if len(quantify_chs) >= 2 else zeros)
        b = (robust_rescale(getter(images[nuc_ch["name"]]),
                            cfg["p_low"], cfg["p_high"])
             if nuc_ch is not None else zeros)
        return np.clip(np.stack([r, g, b], axis=-1), 0, 1)

    mip_rgb = make_rgb(lambda s: s.max(axis=0))
    mid_rgb = make_rgb(lambda s: s[z_mid])
    f = int(cfg["qc_downsample"])

    def proj2d(mask3d, mode):
        return mask3d.max(axis=0) if mode == "mip" else mask3d[z_mid]

    panels = [
        ("MIP", mip_rgb[::f, ::f], proj2d(cell_mask, "mip")[::f, ::f],
         proj2d(nuc_mask, "mip")[::f, ::f], "mip"),
        (f"Z = {z_mid}", mid_rgb[::f, ::f], proj2d(cell_mask, "mid")[::f, ::f],
         proj2d(nuc_mask, "mid")[::f, ::f], "mid"),
    ]
    h, w = panels[0][1].shape[:2]
    dpi = int(cfg["qc_dpi"])
    fig, axes = plt.subplots(1, 2, figsize=(2 * w / dpi, h / dpi + 0.2), dpi=dpi)

    for ax, (title, rgb_p, cell_p, nuc_p, mode) in zip(axes, panels):
        ax.imshow(rgb_p); ax.set_axis_off()
        ax.set_title(title, fontsize=8)
        for cid in np.unique(cell_p):
            if cid == 0: continue
            ax.contour(cell_p == cid, levels=[0.5], colors="cyan", linewidths=0.8)
        if nuc_p.max() > 0:
            for nid in np.unique(nuc_p):
                if nid == 0: continue
                ax.contour(nuc_p == nid, levels=[0.5], colors="yellow", linewidths=0.6)
        if nucleolar_mask is not None and nucleolar_mask.max() > 0:
            nol_p = (nucleolar_mask.max(axis=0) if mode == "mip"
                     else nucleolar_mask[z_mid])[::f, ::f]
            ax.contour(nol_p > 0, levels=[0.5], colors="white", linewidths=0.5)
        for ch_name, pmask in puncta_masks.items():
            pmask_p = (pmask.max(axis=0) if mode == "mip" else pmask[z_mid])[::f, ::f]
            if pmask_p.max() == 0: continue
            if puncta_distance_colors and ch_name in puncta_distance_colors:
                cmap = puncta_distance_colors[ch_name]
                for pid in np.unique(pmask_p):
                    if pid == 0: continue
                    color = cmap.get(int(pid), "magenta")
                    ax.contour(pmask_p == pid, levels=[0.5],
                               colors=color, linewidths=0.5)
            else:
                ax.contour(pmask_p > 0, levels=[0.5], colors="magenta", linewidths=0.5)

    fig.subplots_adjust(left=0, right=1, top=0.95, bottom=0, wspace=0.02)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.05, dpi=dpi)
    plt.close(fig)


def save_qc_png(
    out_png: Path,
    images: dict[str, np.ndarray],
    channels: list[dict],
    cell_mask: np.ndarray,
    nuc_mask: np.ndarray,
    puncta_masks: dict[str, np.ndarray],
    cfg: dict,
    nucleolar_mask: np.ndarray | None = None,
    puncta_distance_colors: dict[str, dict[int, str]] | None = None,
) -> None:
    """Dispatch to 2D or 3D QC based on input dimensionality."""
    if cell_mask.ndim == 3:
        _save_qc_3d(out_png, images, channels, cell_mask, nuc_mask, puncta_masks,
                    cfg, nucleolar_mask, puncta_distance_colors)
    else:
        _save_qc_2d(out_png, images, channels, cell_mask, nuc_mask, puncta_masks,
                    cfg, nucleolar_mask, puncta_distance_colors)


# ---------------------------------------------------------------------------
# Superplots — dynamic metrics + condition order
# ---------------------------------------------------------------------------
def _build_superplot_metrics(
    channels: list[dict],
    puncta_chs: list[dict],
    coloc_pairs: list[tuple[str, str]] | None = None,
    proximity_channels: list[str] | None = None,
    has_nucleolar_morph: bool = False,
    mode: str = "2d",
    cfg: dict | None = None,
) -> list[tuple[str, str, str]]:
    """Return [(column, y_label, title), ...] for superplot generation."""
    metrics: list[tuple[str, str, str]] = []

    for ch in puncta_chs:
        n = ch["name"]
        metrics.append((
            f"{n}_puncta_n",
            f"Puncta per cell ({n})",
            f"Puncta per cell ({n})"))
        metrics.append((
            f"{n}_frac_intensity_in_puncta",
            f"Fraction intensity in puncta ({n})",
            f"Fraction condensed ({n})"))
        # Fragmentation indices (NEW; available in both 2D and 3D)
        metrics.append((
            f"{n}_fragmentation_index_simple",
            f"Fragmentation index simple ({n})",
            f"Fragmentation (simple) — {n}"))
        metrics.append((
            f"{n}_fragmentation_index_persistence",
            f"Fragmentation index persistence ({n})",
            f"Fragmentation (persistence) — {n}"))

    want_ci = bool((cfg or {}).get("condensate_index", False))
    if want_ci:
        for ch in channels:
            if ch["role"] == "skip":
                continue
            n = ch["name"]
            metrics.append((
                f"{n}_condensate_index_cell",
                f"Condensate index (p95/mean) — {n} (cell)",
                f"Condensate index (cell) — {n}"))

    if coloc_pairs:
        for a_name, b_name in coloc_pairs:
            pair = f"{a_name}_vs_{b_name}"
            metrics.append((
                f"pearson_r_{pair}",
                f"Pearson R ({a_name} vs {b_name})",
                f"Pearson R — {a_name} vs {b_name}"))
            metrics.append((
                f"manders_m1_{pair}",
                f"Manders M1 ({a_name} vs {b_name})",
                f"Manders M1 — {a_name} vs {b_name}"))

    if proximity_channels:
        for ch_name in proximity_channels:
            metrics.append((
                f"{ch_name}_mean_distance",
                f"Mean distance to nucleolus ({ch_name})",
                f"Nucleolar distance — {ch_name}"))
            metrics.append((
                f"{ch_name}_fraction_proximal",
                f"Fraction proximal ({ch_name})",
                f"Fraction proximal — {ch_name}"))

    if has_nucleolar_morph:
        if mode == "3d":
            metrics.extend([
                ("nucleolar_volume_um3", "Nucleolar volume (µm³)", "Nucleolar volume"),
                ("nucleolar_eq_diameter_um", "Nucleolar eq. diameter (µm)",
                 "Nucleolar equivalent diameter"),
            ])
        else:
            metrics.extend([
                ("nucleolar_area", "Nucleolar area (px)", "Nucleolar area"),
                ("nucleolar_solidity", "Nucleolar solidity", "Nucleolar solidity"),
                ("nucleolar_circularity", "Nucleolar circularity",
                 "Nucleolar circularity"),
                ("nucleolar_eccentricity", "Nucleolar eccentricity",
                 "Nucleolar eccentricity"),
            ])

    # 3D-only cell volume metric
    if mode == "3d":
        metrics.append(("cell_volume_um3", "Cell volume (µm³)", "Cell volume"))

    return metrics


def compute_pairwise_pvalues(
    df: pd.DataFrame,
    metric: str,
    condition_order: list[str],
    reference_condition: str,
) -> dict[str, dict]:
    """Wilcoxon rank-sum of replicate medians: each condition vs reference.

    Returns {condition: {"pval": float, "pval_corrected": float, "n_ref": int,
    "n_test": int}} for conditions that have ≥3 replicates (ref must too).
    """
    rep_meds = (
        df.groupby(["condition", "image"], as_index=False)[metric]
          .median()
          .rename(columns={metric: "rep_median"})
    )
    ref_meds = rep_meds.loc[
        rep_meds["condition"] == reference_condition, "rep_median"
    ].dropna().values
    if len(ref_meds) < 3:
        return {}

    test_conds = [c for c in condition_order if c != reference_condition]
    # Count testable conditions for Bonferroni
    testable = []
    for cond in test_conds:
        meds = rep_meds.loc[
            rep_meds["condition"] == cond, "rep_median"].dropna().values
        if len(meds) >= 3:
            testable.append((cond, meds))
    n_tests = len(testable)
    if n_tests == 0:
        return {}

    results: dict[str, dict] = {}
    for cond, meds in testable:
        _, pval = mannwhitneyu(meds, ref_meds, alternative="two-sided")
        results[cond] = {
            "pval": pval,
            "pval_corrected": min(pval * n_tests, 1.0),
            "n_test": len(meds),
            "n_ref": len(ref_meds),
            "n_comparisons": n_tests,
        }
    return results


def _format_pval(pval: float) -> str:
    if pval < 0.001:
        return "p < 0.001"
    return f"p = {pval:.2g}"


def _annotate_pvalues(ax, condition_order, pvalues, vals_by_cond,
                      reference_condition=None):
    """Add p-value annotations above each non-reference condition column."""
    if not pvalues:
        return
    all_vals = np.concatenate([v for v in vals_by_cond if v.size > 0])
    if all_vals.size == 0:
        return
    y_max = float(np.nanmax(all_vals))
    y_range = float(np.nanmax(all_vals) - np.nanmin(all_vals))
    if y_range == 0:
        y_range = abs(y_max) * 0.1 or 1.0
    text_y = y_max + 0.06 * y_range
    for i, cond in enumerate(condition_order):
        if cond not in pvalues:
            continue
        info = pvalues[cond]
        p_corr = info["pval_corrected"]
        ax.text(i, text_y, _format_pval(p_corr),
                ha="center", va="bottom", fontsize=7, color="black")
    # Footnote
    if pvalues:
        n_comp = next(iter(pvalues.values()))["n_comparisons"]
        ref_label = (reference_condition or "").capitalize()
        ax.text(0.99, 0.01,
                f"vs {ref_label}, Bonferroni (n={n_comp})",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=6, color="gray")


def _superplot_violin_2(
    df: pd.DataFrame,
    metric: str,
    out_png: Path,
    y_label: str,
    title: str,
    condition_order: list[str],
    jitter_sd: float = 0.06,
    seed: int = 0,
    figsize: tuple[float, float] = (4.4, 6.4),
    pvalues: dict | None = None,
    reference_condition: str | None = None,
) -> None:
    """Violin superplot for <=2 conditions."""
    rng = np.random.default_rng(seed)

    vals_by_cond: list[np.ndarray] = []
    for cond in condition_order:
        vals = (df.loc[df["condition"] == cond, metric]
                .dropna().astype(float).values)
        vals_by_cond.append(vals)

    # Skip if any condition has no data (violinplot crashes on empty arrays)
    if any(v.size == 0 for v in vals_by_cond):
        return

    fig = plt.figure(figsize=figsize)

    plt.violinplot(
        vals_by_cond,
        positions=list(range(len(condition_order))),
        widths=0.72,
        showmeans=False,
        showextrema=False,
        showmedians=False,
    )

    for i, cond in enumerate(condition_order):
        vals = (df.loc[df["condition"] == cond, metric]
                .dropna().astype(float).values)
        x = i + rng.normal(0, jitter_sd, size=len(vals))
        plt.scatter(x, vals, s=10)

    # Replicate medians (per image)
    rep_meds = (
        df.groupby(["condition", "image"], as_index=False)[metric]
          .median()
          .rename(columns={metric: "rep_median"})
    )

    for i, cond in enumerate(condition_order):
        meds = rep_meds.loc[
            rep_meds["condition"] == cond, "rep_median"].values
        if len(meds) == 0:
            continue
        offsets = (np.linspace(-0.12, 0.12, num=len(meds))
                   if len(meds) > 1 else np.array([0.0]))
        plt.scatter(
            i + offsets, meds,
            s=140, marker="D", linewidths=1.2,
            edgecolors="black", facecolors="none", zorder=5,
        )

    # Overall median horizontal line per condition
    for i, cond in enumerate(condition_order):
        vals = (df.loc[df["condition"] == cond, metric]
                .dropna().astype(float).values)
        if len(vals) > 0:
            med = float(np.median(vals))
            plt.hlines(med, i - 0.3, i + 0.3,
                        colors="black", linewidths=1.5, zorder=6)

    ax = plt.gca()
    if pvalues and len(condition_order) == 2 and len(pvalues) == 1:
        # 2-condition bracket style
        info = next(iter(pvalues.values()))
        p_text = _format_pval(info["pval"])
        all_vals = np.concatenate(vals_by_cond)
        y_max = float(np.nanmax(all_vals))
        y_range = float(np.nanmax(all_vals) - np.nanmin(all_vals))
        bar_y = y_max + 0.06 * y_range
        tip_len = 0.02 * y_range
        ax.plot([0, 0, 1, 1],
                [bar_y - tip_len, bar_y, bar_y, bar_y - tip_len],
                color="black", linewidth=1.0)
        ax.text(0.5, bar_y + 0.01 * y_range, p_text,
                ha="center", va="bottom", fontsize=9)
    elif pvalues:
        _annotate_pvalues(ax, condition_order, pvalues, vals_by_cond,
                          reference_condition)

    display_labels = [c.capitalize() for c in condition_order]
    plt.xticks(range(len(condition_order)), display_labels)
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def _superplot_strip_multi(
    df: pd.DataFrame,
    metric: str,
    out_png: Path,
    y_label: str,
    title: str,
    condition_order: list[str],
    show_trend: bool = False,
    jitter_sd: float = 0.08,
    seed: int = 0,
    pvalues: dict | None = None,
    reference_condition: str | None = None,
) -> None:
    """Jittered strip plot for 3+ conditions (or 2 with --trend)."""
    n_cond = len(condition_order)
    figsize = (max(4.4, 1.5 * n_cond), 6.4)
    rng = np.random.default_rng(seed)

    # Check for any data to plot
    has_data = False
    for cond in condition_order:
        if df.loc[df["condition"] == cond, metric].dropna().size > 0:
            has_data = True
            break
    if not has_data:
        return

    fig = plt.figure(figsize=figsize)

    vals_by_cond: list[np.ndarray] = []
    for i, cond in enumerate(condition_order):
        vals = (df.loc[df["condition"] == cond, metric]
                .dropna().astype(float).values)
        vals_by_cond.append(vals)
        x = i + rng.normal(0, jitter_sd, size=len(vals))
        plt.scatter(x, vals, s=8, alpha=0.5)

    # Per-image diamond medians
    rep_meds = (
        df.groupby(["condition", "image"], as_index=False)[metric]
          .median()
          .rename(columns={metric: "rep_median"})
    )

    cond_medians = []
    for i, cond in enumerate(condition_order):
        meds = rep_meds.loc[
            rep_meds["condition"] == cond, "rep_median"].values
        if len(meds) == 0:
            cond_medians.append(np.nan)
            continue
        offsets = (np.linspace(-0.12, 0.12, num=len(meds))
                   if len(meds) > 1 else np.array([0.0]))
        plt.scatter(
            i + offsets, meds,
            s=140, marker="D", linewidths=1.2,
            edgecolors="black", facecolors="none", zorder=5,
        )
        cond_medians.append(float(np.median(meds)))

    # Overall median horizontal line per condition
    for i, cond in enumerate(condition_order):
        vals = (df.loc[df["condition"] == cond, metric]
                .dropna().astype(float).values)
        if len(vals) > 0:
            med = float(np.median(vals))
            plt.hlines(med, i - 0.3, i + 0.3,
                        colors="black", linewidths=1.5, zorder=6)

    # Optional trend line through condition medians
    if show_trend:
        valid = [(i, m) for i, m in enumerate(cond_medians) if not np.isnan(m)]
        if len(valid) >= 2:
            xs, ys = zip(*valid)
            plt.plot(xs, ys, "--", color="gray", linewidth=1.2, zorder=4)

    ax = plt.gca()
    if pvalues:
        _annotate_pvalues(ax, condition_order, pvalues, vals_by_cond,
                          reference_condition)

    display_labels = [c.capitalize() for c in condition_order]
    plt.xticks(range(n_cond), display_labels, rotation=45, ha="right")
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def superplot_violin(
    df: pd.DataFrame,
    metric: str,
    out_png: Path,
    y_label: str,
    title: str,
    condition_order: list[str],
    show_trend: bool = False,
    pvalues: dict | None = None,
    reference_condition: str | None = None,
    **kwargs,
) -> None:
    """Dispatch: violin for <=2 conditions (no trend), strip plot otherwise."""
    if len(condition_order) <= 2 and not show_trend:
        _superplot_violin_2(df, metric, out_png, y_label, title,
                            condition_order, pvalues=pvalues,
                            reference_condition=reference_condition, **kwargs)
    else:
        _superplot_strip_multi(df, metric, out_png, y_label, title,
                               condition_order, show_trend=show_trend,
                               pvalues=pvalues,
                               reference_condition=reference_condition)


# ---------------------------------------------------------------------------
# Prism-ready CSV export — dynamic metrics
# ---------------------------------------------------------------------------
def write_prism_outputs(
    df: pd.DataFrame,
    metrics: list[str],
    prism_dir: Path,
    condition_order: list[str],
) -> None:
    prism_dir.mkdir(parents=True, exist_ok=True)

    base_cols = ["condition", "replicate", "image", "cell_id"]
    keep_cols = [c for c in base_cols if c in df.columns]
    avail = [m for m in metrics if m in df.columns]

    # Long tidy (all cells)
    df[keep_cols + avail].to_csv(prism_dir / "prism_cells_long.csv", index=False)

    # Replicate medians (per image)
    df_rep = df.groupby(
        ["condition", "replicate", "image"], as_index=False)[avail].median()
    df_rep.to_csv(prism_dir / "prism_replicate_medians.csv", index=False)

    # Wide ragged per metric (one column per image replicate)
    for metric in avail:
        cols: dict[str, list[float]] = {}
        for cond in condition_order:
            sub = df[df["condition"] == cond]
            for img, g in sub.groupby("image"):
                colname = f"{cond}__{img}"
                cols[colname] = g[metric].dropna().astype(float).tolist()
        if not cols:
            continue
        maxlen = max(len(v) for v in cols.values())
        wide = {k: v + [np.nan] * (maxlen - len(v)) for k, v in cols.items()}
        pd.DataFrame(wide).to_csv(
            prism_dir / f"wide_{metric}.csv", index=False)


# ---------------------------------------------------------------------------
# Provenance sidecar — version, CLI invocation, input checksums, deps
# ---------------------------------------------------------------------------
def _sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    """Return the hex SHA-256 of a file, streamed in 1 MB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _dep_versions() -> dict[str, str]:
    """Best-effort version pin of key runtime dependencies."""
    out: dict[str, str] = {"python": sys.version.split()[0]}
    for mod_name, attr in [
        ("numpy", "__version__"),
        ("scipy", "__version__"),
        ("pandas", "__version__"),
        ("matplotlib", "__version__"),
        ("skimage", "__version__"),
        ("tifffile", "__version__"),
        ("yaml", "__version__"),
        ("cellpose", "version"),
        ("torch", "__version__"),
        ("cv2", "__version__"),
    ]:
        try:
            mod = __import__(mod_name)
            out[mod_name] = str(getattr(mod, attr, "unknown"))
        except Exception:
            out[mod_name] = "not installed"
    return out


def _resolved_voxel_size(cfg: dict) -> tuple[float | None, float | None, str | None]:
    """Return the (xy, z, source) voxel size a 3D run will ACTUALLY use.

    cfg["voxel_size_*_um"] still holds the (1.0, 1.0) DEFAULTS value when the
    real voxel size is resolved per image from file metadata, so reading cfg
    directly reports a voxel size the run never uses — which is precisely the
    silent 1 µm isotropic fallback this pipeline refuses to perform. The
    startup banner and the provenance sidecar both go through here so they
    cannot drift apart.
    """
    vprov = cfg.get("_voxel_provenance") or {}
    src = vprov.get("resolved_from")
    if src == "metadata":
        mv = vprov.get("metadata_voxel_size_um") or {}
        xy, z = mv.get("xy"), mv.get("z")
        if xy and z:
            return float(xy), float(z), "metadata"
        return None, None, "metadata"
    if src == "assumed_isotropic":
        return 1.0, 1.0, "assumed_isotropic"
    xy, z = cfg.get("voxel_size_xy_um"), cfg.get("voxel_size_z_um")
    return (float(xy) if xy is not None else None,
            float(z) if z is not None else None,
            src)


def _write_provenance(
    out_dir: Path,
    paths: list[Path],
    channels: list[dict],
    cfg: dict,
) -> None:
    """Write a provenance JSON sidecar to out_dir/provenance.json.

    Records cellquant version, the exact CLI invocation, sha256 checksums of
    every input file, the Cellpose model checkpoint name, and pinned versions
    of every key dependency. Together these are sufficient to reproduce a run
    given the same inputs.
    """
    inputs = []
    for p in paths:
        try:
            stat = p.stat()
            entry: dict[str, Any] = {
                "name": p.name,
                "size_bytes": int(stat.st_size),
                "mtime": time.strftime(
                    "%Y-%m-%dT%H:%M:%S", time.localtime(stat.st_mtime)),
                "sha256": _sha256_file(p),
            }
            ome_names = _read_ome_channel_names(p)
            if ome_names:
                entry["ome_channel_names"] = ome_names
            ome_voxel = _read_ome_voxel_size(p)
            if ome_voxel is not None:
                entry["ome_voxel_size_um"] = {"xy": ome_voxel[0], "z": ome_voxel[1]}
            inputs.append(entry)
        except Exception as exc:  # pragma: no cover — best-effort sidecar
            inputs.append({"name": p.name, "error": str(exc)})

    _resolved_voxel = _resolved_voxel_size(cfg)
    record = {
        "cellquant_version": __version__,
        "schema_version": 1,
        "run_at": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime()),
        "cwd": str(Path.cwd()),
        "command_line": " ".join(sys.argv),
        "argv": list(sys.argv),
        "mode": cfg.get("mode", "auto"),
        "channels": [
            {"position": ch["position"], "name": ch["name"], "role": ch["role"]}
            for ch in channels
        ],
        "cell_type": cfg.get("cell_type"),
        "pretrained_model": cfg.get("pretrained_model"),
        "voxel_size_um": {
            "xy": _resolved_voxel[0],
            "z": _resolved_voxel[1],
        } if cfg.get("mode") == "3d" else None,
        "voxel_size_provenance": (cfg.get("_voxel_provenance")
                                  if cfg.get("mode") == "3d" else None),
        "use_gpu": cfg.get("use_gpu"),
        "dependencies": _dep_versions(),
        "inputs": inputs,
    }
    with open(out_dir / "provenance.json", "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=False)


# ---------------------------------------------------------------------------
# Directory setup
# ---------------------------------------------------------------------------
def ensure_dirs(out_root: Path) -> tuple[Path, Path, Path]:
    qc_dir = out_root / "qc"
    mask_dir = out_root / "masks"
    out_root.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    mask_dir.mkdir(parents=True, exist_ok=True)
    return qc_dir, out_root, mask_dir


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    # Parse channel definitions
    channels = parse_channels(args.channels)

    # Build resolved config (4-layer)
    cfg = build_config(args)

    _set_threads(cfg)

    # Resolve channel roles
    nuc_ch = get_nucleus_channel(channels)  # None if no nucleus channel
    has_nuclei = nuc_ch is not None
    quantify_chs = get_quantify_channels(channels)
    nucleolus_chs = get_nucleolus_channels(channels)

    # Auto-populate puncta channels from quantify channels when not specified
    if args.no_puncta:
        puncta_names = []
    elif args.puncta_channels is not None:
        puncta_names = args.puncta_channels
    else:
        # Default: detect puncta in all quantify channels
        puncta_names = [ch["name"] for ch in quantify_chs] if quantify_chs else []
    puncta_chs = get_puncta_channels(channels, puncta_names if puncta_names else None)

    # Colocalization channels = quantify + nucleolus (deduplicated)
    coloc_channels: list[dict] = []
    seen_names: set[str] = set()
    for ch in quantify_chs + nucleolus_chs:
        if ch["name"] not in seen_names:
            coloc_channels.append(ch)
            seen_names.add(ch["name"])

    # --- Validation checks ---
    nuc_prox_ch_name = cfg.get("nucleolar_proximity_channel")
    if nuc_prox_ch_name:
        if not nucleolus_chs:
            raise SystemExit(
                f"--nucleolar-proximity '{nuc_prox_ch_name}' requires at least "
                f"one channel with role 'nucleolus'")
        nuc_prox_match = [ch for ch in nucleolus_chs
                          if ch["name"].lower() == nuc_prox_ch_name.lower()]
        if not nuc_prox_match:
            raise SystemExit(
                f"--nucleolar-proximity '{nuc_prox_ch_name}' must reference a "
                f"channel with role 'nucleolus' (found: "
                f"{[ch['name'] for ch in nucleolus_chs]})")

    do_coloc = bool(cfg.get("colocalization", False))
    if do_coloc and len(coloc_channels) < 2:
        print("[warn] --colocalization requires >=2 eligible channels; skipping")
        do_coloc = False

    compartment = cfg.get("puncta_compartment", "cytosol")
    if not has_nuclei and compartment == "cytosol":
        print("[warn] No nucleus channel — falling back puncta_compartment "
              "from 'cytosol' to 'whole-cell'")
        cfg["puncta_compartment"] = "whole-cell"
        compartment = "whole-cell"


    if not has_nuclei and int(cfg.get("keep_min_nuclei", 0)) > 0:
        print("[warn] No nucleus channel but keep_min_nuclei > 0 — "
              "all cells will fail nuclei gate")

    # Determine cell seg input mode
    # Use composite (sum of all non-skip channels) unless the user specifies
    # a cell-boundary role or explicit --cell-seg-channel.  Composite avoids
    # the problem of a single quantify channel's bright puncta (e.g. stress
    # granules) being mistakenly segmented as tiny "cells" by Cellpose.
    has_cell_boundary = any(ch["role"] == "cell-boundary" for ch in channels)
    explicit_seg_ch = cfg.get("cell_seg_channel")
    use_composite_seg = (not has_cell_boundary and not explicit_seg_ch)

    if not use_composite_seg:
        cell_seg_ch = get_cell_seg_channel(channels, cfg)
    else:
        cell_seg_ch = None  # will use composite

    # I/O paths
    img_dir = Path(args.images)
    out_root = Path(args.out)
    qc_dir, out_dir, mask_dir = ensure_dirs(out_root)

    # Find images
    paths = find_images(img_dir, cfg)
    if not paths:
        # An empty glob is a hard error, never success: name the pattern that
        # matched nothing so a mistyped --filename-pattern or wrong directory
        # is obvious instead of silently producing zero output.
        file_pattern = cfg.get("file_pattern")
        if file_pattern:
            glob_desc = f"glob '{file_pattern}'"
        else:
            exts = cfg.get("exts", [".tif", ".tiff"])
            glob_desc = "extensions " + ", ".join(f"'*{e}'" for e in exts)
        raise SystemExit(
            f"No input images found: {glob_desc} matched no files in "
            f"{img_dir}. Check the images directory and --filename-pattern."
        )

    # ----- Mode detection (auto from first file unless user forced) -----
    requested_mode = cfg.get("mode", "auto")
    if requested_mode == "auto":
        detected = detect_input_mode(paths[0], len(channels))
        cfg["mode"] = detected
        print(f"Mode: {detected} (auto-detected from {paths[0].name})")
    else:
        cfg["mode"] = requested_mode
        print(f"Mode: {requested_mode} (forced)")

    # --project-z forces pipeline-level 2D mode (input is projected per file
    # in load_images before the 2D pipeline runs).
    if cfg.get("project_z") and cfg["mode"] == "3d":
        print(f"--project-z {cfg['project_z']}: projecting z-stacks to 2D "
              f"before analysis")
        cfg["mode"] = "2d"

    # Voxel-size sanity: if we ended up in 3D mode but the user supplied no
    # --voxel-size AND the first file has no OME/ImageJ metadata, downstream
    # 3D metrics (volume in µm³, anisotropic LoG, anisotropic distances) will
    # be wrong. Warn loudly and recommend a fix. This is a one-shot check on
    # the first file — assume the rest of the dataset is homogeneous.
    if cfg["mode"] == "3d":
        # Voxel-size resolution ladder. The CLI value and the file's OME/ImageJ
        # metadata are BOTH read and compared, so a disagreement is surfaced
        # (rather than the CLI silently winning). Both values and any override
        # are recorded in provenance.json.
        user_set_voxel = _user_supplied_voxel(cfg)
        cli_voxel = ((cfg["voxel_size_xy_um"], cfg["voxel_size_z_um"])
                     if user_set_voxel else None)
        ome_voxel = _read_ome_voxel_size(paths[0])
        vprov: dict[str, Any] = {
            "cli_voxel_size_um": ({"xy": cli_voxel[0], "z": cli_voxel[1]}
                                  if cli_voxel else None),
            "metadata_voxel_size_um": ({"xy": ome_voxel[0], "z": ome_voxel[1]}
                                       if ome_voxel else None),
            "override_metadata": bool(cfg.get("override_metadata")),
            "assume_isotropic": bool(cfg.get("assume_isotropic")),
            "cli_metadata_disagreement": None,
            "resolved_from": None,
        }
        bar = "!" * 78
        if cli_voxel is not None and ome_voxel is not None:
            dxy = abs(cli_voxel[0] - ome_voxel[0]) / ome_voxel[0]
            dz = abs(cli_voxel[1] - ome_voxel[1]) / ome_voxel[1]
            vprov["cli_metadata_disagreement"] = {"xy_frac": dxy, "z_frac": dz}
            if max(dxy, dz) > 0.01 and not cfg.get("override_metadata"):
                print(f"\n{bar}")
                print("[error] --voxel-size disagrees with the file's OME/ImageJ voxel")
                print(f"        metadata by >1%: CLI=(XY={cli_voxel[0]:.5f}, "
                      f"Z={cli_voxel[1]:.5f}) µm vs")
                print(f"        metadata=(XY={ome_voxel[0]:.5f}, Z={ome_voxel[1]:.5f}) µm "
                      f"(ΔXY={dxy*100:.1f}%, ΔZ={dz*100:.1f}%).")
                print("        3D volumes, anisotropic LoG and distances all")
                print("        depend on this. Re-run with the correct --voxel-size, or")
                print("        pass --override-metadata to force the CLI value.")
                print(f"{bar}\n")
                sys.exit(1)
            elif max(dxy, dz) > 0.01:
                print(f"\n[warn] --voxel-size overrides file metadata "
                      f"(ΔXY={dxy*100:.1f}%, ΔZ={dz*100:.1f}%); using the CLI value "
                      f"per --override-metadata.\n")
                vprov["resolved_from"] = "cli_override_metadata"
            else:
                vprov["resolved_from"] = "cli"
        elif cli_voxel is not None:
            vprov["resolved_from"] = "cli"
        elif ome_voxel is not None:
            reason = _voxel_size_is_implausible(*ome_voxel)
            if reason:
                print(f"\n[warn] Voxel size read from metadata looks implausible:")
                print(f"         {reason}")
                print(f"       Override with --voxel-size XY_UM Z_UM if this is wrong.\n")
            vprov["resolved_from"] = "metadata"
        elif cfg.get("assume_isotropic"):
            print("\n[warn] 3D mode with no voxel size from --voxel-size or metadata;")
            print("       assuming XY=Z=1.0 µm per --assume-isotropic. 3D physical-unit")
            print("       metrics (volume µm³, distances) will be in voxel")
            print("       units, not microns.\n")
            vprov["resolved_from"] = "assumed_isotropic"
        else:
            print(f"\n{bar}")
            print("[error] 3D mode but no voxel size is available from --voxel-size or")
            print("        from OME/ImageJ metadata. cellquant will not silently assume")
            print("        1.0 µm isotropic voxels — that makes every 3D metric wrong.")
            print("        Pass --voxel-size XY_UM Z_UM, re-export TIFFs with metadata,")
            print("        or pass --assume-isotropic to proceed in voxel units.")
            print(f"{bar}\n")
            sys.exit(1)
        cfg["_voxel_provenance"] = vprov

    # OME channel-name cross-check (one-shot, on the first file): warn the
    # user if the channel ordering they declared on the CLI doesn't match
    # what's in the file's OME metadata. Doesn't override user intent — only
    # surfaces mismatches that often indicate a misconfigured invocation.
    ome_channel_names = _read_ome_channel_names(paths[0])
    if ome_channel_names:
        for ch in channels:
            idx = ch["position"] - 1
            if 0 <= idx < len(ome_channel_names):
                ome_name = ome_channel_names[idx]
                if (ome_name and ch["name"].lower() != ome_name.lower()
                        and ch["role"] != "skip"):
                    print(f"[warn] channel position {ch['position']}: you "
                          f"declared '{ch['name']}' but OME metadata says "
                          f"'{ome_name}' — double-check your --channels order")

    # Apply 3D preset overrides (only if we ended up in 3D mode and user didn't
    # supply explicit CLI args for those keys)
    if cfg["mode"] == "3d":
        apply_3d_preset_overrides(cfg)
    else:
        cfg.pop("_preset_3d_overrides", None)
        # 2D lateral pixel size. The voxel ladder above runs only in 3D, so in
        # 2D cfg["voxel_size_xy_um"] stayed at its 1.0 default and every
        # micron-denominated parameter was silently measured in pixels instead.
        if not _user_supplied_voxel(cfg):
            px = _read_pixel_size_xy_um(paths[0])
            if px:
                cfg["voxel_size_xy_um"] = px
                print(f"Pixel size: {px:.4f} \u00b5m (from file metadata)")
        um_params = _um_denominated_2d_params(cfg)
        if um_params and float(cfg.get("voxel_size_xy_um", 1.0) or 1.0) == 1.0:
            asserted = [u for u in um_params if u[2]]
            bar = "!" * 78
            if asserted:
                print(f"\n{bar}")
                print("[error] you supplied these values in microns, but no pixel size")
                print("        is available, so they would be applied in PIXELS:")
                for label, _, _ in asserted:
                    print(f"          {label}")
                print("        Pass --voxel-size XY_UM, or re-export the images with")
                print("        resolution metadata. (A micron value below 1.0 applied")
                print("        as pixels does nothing at all.)")
                print(f"{bar}\n")
                sys.exit(1)
            # Preset-supplied only: disable rather than apply as pixels, and say so.
            for label, key, _ in um_params:
                print(f"[warn] {label} = {cfg.get(key)} µm comes from the "
                      f"'{cfg.get('cell_type')}' preset, but this 2D input carries "
                      f"no pixel size,")
                print(f"       so the value cannot be interpreted in microns. It is "
                      f"being SKIPPED, not applied as pixels.")
                print(f"       Pass --voxel-size XY_UM to enable it.")
                if key == "puncta_compartment_erode_um":
                    cfg[key] = 0.0

    # Resolve every compartment reference now that mode and the pixel size are
    # final, and before any segmentation runs. A typo used to surface deep in
    # the per-image loop, after Cellpose had already spent minutes.
    validate_compartment_config(cfg, channels)
    if cfg.get("compartments"):
        print("Compartments: " + "; ".join(
            d["raw"] for d in cfg["compartments"].values()))

    # Loud warning when running colocalization on 2D MIP input. Pearson's R
    # and Manders' coefficients are defined on the full 3D voxel distribution;
    # computing them on MIPs is statistically unreliable (overlapping signal
    # from different Z planes inflates apparent colocalization). We do not
    # refuse the run — some users legitimately only have MIPs and want a
    # qualitative comparison — but the warning is conspicuous so the
    # limitation isn't silent.
    if do_coloc and cfg["mode"] == "2d":
        bar = "!" * 78
        if not cfg.get("allow_2d_colocalization"):
            print(f"\n{bar}")
            print("[error] Colocalization was requested on 2D / MIP input.")
            print("        Pearson's R and Manders' M1/M2 are defined on the full 3D")
            print("        voxel distribution; computing them on a projection collapses")
            print("        Z and biases the result. cellquant refuses this by default.")
            print("        Re-run on the source z-stacks (--mode 3d / --voxel-size), or")
            print("        pass --allow-2d-colocalization to force it (outputs will be")
            print("        stamped projection_derived=True).")
            print(f"{bar}\n")
            sys.exit(1)
        cfg["_coloc_projection_derived"] = True
        print(f"\n{bar}")
        print("[warn] Colocalization on a 2D / MIP is projection-derived and")
        print("       statistically unreliable; colocalization.csv is stamped")
        print("       projection_derived=True. Prefer re-running on source z-stacks")
        print("       (cellquant computes 3D colocalization natively).")
        print(f"{bar}\n")

    # Save resolved config (now that mode is locked in)
    save_cfg: dict[str, Any] = dict(cfg)
    save_cfg["_channels"] = [
        f"{ch['position']}:{ch['name']}:{ch['role']}" for ch in channels]
    if puncta_chs:
        save_cfg["_puncta_channels"] = [ch["name"] for ch in puncta_chs]
    save_cfg["_cell_seg_channel"] = (
        cell_seg_ch["name"] if cell_seg_ch else "composite")
    with open(out_dir / "config_used.yml", "w") as fh:
        yaml.dump(save_cfg, fh, default_flow_style=False, sort_keys=False)

    # Write provenance sidecar (FAIR / reproducibility): version, exact CLI
    # invocation, input file checksums, model checkpoint name, key dependency
    # versions. Lets anyone reading the output reconstruct what was run.
    _write_provenance(out_dir, paths, channels, cfg)

    print(f"Found {len(paths)} images in {img_dir}")
    print(f"GPU: {cfg['use_gpu']}")
    if cfg["mode"] == "3d":
        print(f"3D segmentation method: {cfg.get('seg_3d_method', 'stitch')}")
        # Report the voxel size the run will ACTUALLY use, and where it came
        # from. cfg still holds the (1.0, 1.0) default when the value is
        # resolved from file metadata per image, so printing cfg here used to
        # announce "XY=1.0000 µm, Z=1.0000 µm" on runs that went on to read the
        # correct voxel from the file — the exact silent-isotropic-fallback
        # this pipeline refuses to do.
        _vxy, _vz, _vsrc = _resolved_voxel_size(cfg)
        # "first file" is deliberate: the metadata voxel size is latched from
        # the first image and reused for the rest of the run, so this must not
        # claim to be resolved per image.
        _label = {
            "metadata": "from file metadata (first file)",
            "assumed_isotropic": "assumed per --assume-isotropic; 3D metrics "
                                 "are in voxel units, not microns",
            "cli": "user-specified (--voxel-size or --config)",
            "cli_override_metadata": "user-specified, overriding file metadata",
        }.get(_vsrc, "source unrecorded")
        if _vxy and _vz:
            print(f"Voxel size: XY={_vxy:.4f} µm, Z={_vz:.4f} µm "
                  f"(anisotropy={_vz / _vxy:.2f}; {_label})")
        else:
            print(f"Voxel size: {_label}")
    print(f"Channels: "
          + ", ".join(f"{ch['name']}({ch['role']})" for ch in channels))
    print(f"Nuclear segmentation: {'yes' if has_nuclei else 'skipped (no nucleus channel)'}")
    print(f"Cell seg channel: "
          f"{cell_seg_ch['name'] if cell_seg_ch else 'composite (all non-skip)'}")
    if puncta_chs:
        print(f"Puncta channels: "
              + ", ".join(ch["name"] for ch in puncta_chs))
        print(f"Puncta compartment: {cfg['puncta_compartment']}")
    if do_coloc:
        print(f"Colocalization: {[ch['name'] for ch in coloc_channels]}")
    if nuc_prox_ch_name:
        print(f"Nucleolar proximity: {nuc_prox_ch_name}")

    reuse_masks = bool(cfg.get("reuse_masks", False))
    if reuse_masks:
        print("Reusing existing masks from", mask_dir)
    model = None if reuse_masks else init_model(cfg)

    all_cells: list[pd.DataFrame] = []
    all_imgs: list[pd.DataFrame] = []
    all_coloc: list[pd.DataFrame] = []
    all_proximity: list[pd.DataFrame] = []
    all_nol_morph: list[pd.DataFrame] = []
    t_total = time.time()

    is_3d = (cfg["mode"] == "3d")
    if is_3d:
        proximity_threshold = float(cfg.get("proximity_threshold_um", 0.5))
        # Sensitivity band: also emit fraction_proximal at ±0.1 µm around the
        # primary threshold (a robustness check, not a re-derivation of the
        # threshold from pixels). One run yields all three columns.
        proximity_extra_thresholds = [
            round(t, 4)
            for t in (proximity_threshold - 0.1, proximity_threshold + 0.1)
            if t > 0 and abs(t - proximity_threshold) > 1e-9
        ]
    else:
        proximity_threshold = float(cfg.get("proximity_threshold_px", 5))
        proximity_extra_thresholds = []

    for i, p in enumerate(paths, 1):
        print(f"\n=== [{i}/{len(paths)}] {p.name} ===")
        t0 = time.time()

        # Load images via unified loader (mode-aware)
        images, file_mode, voxel_size = load_images(p, channels, cfg)
        if file_mode != cfg["mode"]:
            print(f"  [warn] file looks {file_mode} but pipeline mode is "
                  f"{cfg['mode']}; reading anyway")
        if voxel_size is not None:
            # If metadata gave us a voxel size and the user didn't supply one,
            # adopt it for this file's analysis. Guarding on
            # _user_supplied_voxel rather than a (1.0, 1.0) comparison keeps
            # this consistent with both loaders and the ladder, and stops the
            # message below from claiming "from metadata" for a voxel size the
            # user passed on the command line.
            # The second test matters: the loaders return cfg's own value when
            # a file carries no voxel metadata, so without it this reports
            # "from metadata" for files that have none (e.g. under
            # --assume-isotropic).
            if (not _user_supplied_voxel(cfg)
                    and voxel_size != (cfg["voxel_size_xy_um"],
                                       cfg["voxel_size_z_um"])):
                cfg["voxel_size_xy_um"], cfg["voxel_size_z_um"] = voxel_size
                print(f"  voxel size from metadata: "
                      f"XY={voxel_size[0]:.4f}, Z={voxel_size[1]:.4f} µm")

        # Optional per-image adaptive z-crop. Cells straddling the crop bound-
        # ary are dropped after segmentation (see exclude_cells_at_z_boundary
        # call below). Applied only in 3D mode; 2D inputs are unchanged.
        z_crop_window = int(cfg.get("z_crop_center", 0))
        z_crop_applied = False
        if z_crop_window > 0 and file_mode == "3d":
            orig_z = next(iter(images.values())).shape[0]
            images, z_lo, z_hi = apply_z_crop(images, z_crop_window, channels)
            new_z = next(iter(images.values())).shape[0]
            if new_z < orig_z:
                z_crop_applied = True
                print(f"  z-crop: peak-centered, kept Z=[{z_lo},{z_hi}] "
                      f"({new_z} of {orig_z} slices)")

        first_img = next(iter(images.values()))
        ref_shape = first_img.shape  # (Y, X) in 2D, (Z, Y, X) in 3D

        if reuse_masks:
            cell_mask_path = mask_dir / f"{p.stem}_cellmask.tif"
            if not cell_mask_path.exists():
                print(f"  [warn] mask not found: {cell_mask_path.name}, skipping")
                continue
            cell_mask = np.asarray(
                tiff.imread(str(cell_mask_path))).astype(np.int32)
            nuc_mask_path = mask_dir / f"{p.stem}_nucmask.tif"
            if has_nuclei and nuc_mask_path.exists():
                nuc_mask = np.asarray(
                    tiff.imread(str(nuc_mask_path))).astype(np.int32)
            else:
                nuc_mask = np.zeros(ref_shape, dtype=np.int32)
        else:
            # Determine cell seg input
            if use_composite_seg:
                cell_seg_img = build_composite_seg_image(images, channels, cfg)
            else:
                cell_seg_img = images[cell_seg_ch["name"]]
            ref_shape = cell_seg_img.shape

            if is_3d:
                # 3D path: apply seg_downsample to XY only (preserve Z). The
                # 2D path uses seg_downsample as spatial regularization that
                # lets Cellpose-SAM see whole cells rather than fitting tight
                # around bright nuclear regions; in 3D we get the same effect
                # by skipping every Nth pixel in Y and X while keeping the
                # full Z stack. Yeast/bacteria presets use seg_downsample=1
                # so this is a no-op for those.
                ds = max(1, int(cfg.get("seg_downsample", 1)))
                if has_nuclei:
                    nuc_img = images[nuc_ch["name"]]
                    nuc_seg_in = (nuc_img[:, ::ds, ::ds]
                                  if ds > 1 else nuc_img)
                    nuc_small = segment_nuclei(model, nuc_seg_in, cfg)
                    if ds > 1:
                        nuc_mask = upsample_labels_nn(
                            nuc_small,
                            (nuc_small.shape[0], ref_shape[1], ref_shape[2]))
                    else:
                        nuc_mask = nuc_small
                else:
                    nuc_mask = np.zeros(ref_shape, dtype=np.int32)

                cell_seg_in = (cell_seg_img[:, ::ds, ::ds]
                               if ds > 1 else cell_seg_img)
                cell_small = segment_cells(model, cell_seg_in, cfg)
                if ds > 1:
                    cell_mask = upsample_labels_nn(
                        cell_small,
                        (cell_small.shape[0], ref_shape[1], ref_shape[2]))
                else:
                    cell_mask = cell_small

                # Z-boundary exclusion (only when a z-crop was applied) —
                # drop any cell whose mask touches the top or bottom plane of
                # the cropped volume so truncated half-cells don't enter the
                # per-cell measurements.
                if z_crop_applied:
                    cell_mask, n_zbex = exclude_cells_at_z_boundary(cell_mask)
                    if n_zbex:
                        print(f"  z-boundary exclusion: dropped {n_zbex} "
                              f"straddler cell(s)")

                # Cell volume filtering (3D)
                min_v = int(cfg.get("min_cell_volume_vox", 0))
                max_v = int(cfg.get("max_cell_volume_vox", 0))
                cell_mask = filter_cells_by_area(cell_mask, min_v, max_v)
            else:
                # 2D path: paper-validated, byte-identical to original
                ds = max(1, int(cfg.get("seg_downsample", 1)))
                if has_nuclei:
                    nuc_img = images[nuc_ch["name"]]
                    nuc_seg_in = nuc_img[::ds, ::ds] if ds > 1 else nuc_img
                    nuc_small = segment_nuclei(model, nuc_seg_in, cfg)
                    nuc_mask = (upsample_labels_nn(nuc_small, ref_shape)
                                if ds > 1 else nuc_small)
                else:
                    nuc_mask = np.zeros(ref_shape, dtype=np.int32)

                cell_seg_in = cell_seg_img[::ds, ::ds] if ds > 1 else cell_seg_img
                cell_small = segment_cells(model, cell_seg_in, cfg)
                cell_mask = (upsample_labels_nn(cell_small, ref_shape)
                             if ds > 1 else cell_small)

                # Cell area filtering (2D)
                min_area = int(cfg.get("min_cell_area", 0))
                max_area = int(cfg.get("max_cell_area", 0))
                cell_mask = filter_cells_by_area(cell_mask, min_area, max_area)

        warn_cell_shape(cell_mask, cfg)
        cell_to_nucs = map_nuclei_to_cells(nuc_mask, cell_mask)

        # Parse metadata from filename
        meta = parse_filename_metadata(p.stem, cfg)

        # Nucleolar segmentation
        nucleolar_mask = np.zeros(ref_shape, dtype=np.uint8)
        for nol_ch in nucleolus_chs:
            nucleolar_mask |= segment_nucleoli(
                images[nol_ch["name"]], cell_mask, cfg)

        # Nucleolar morphometrics (cfg now passed for 3D voxel size)
        if nucleolar_mask.max() > 0:
            nol_morph_df = compute_nucleolar_morphology(
                nucleolar_mask, cell_mask, meta, p.name, cfg=cfg)
            all_nol_morph.append(nol_morph_df)

        # Create compartment mask for puncta detection
        compartment_mask = make_compartment_mask(
            cell_mask, nuc_mask, compartment, cfg,
            nucleolar_mask=nucleolar_mask)

        # Detect puncta per channel — each channel sees its own resolved cfg
        # (per-channel overrides applied; rolling-ball radius routed in).
        puncta_masks: dict[str, np.ndarray] = {}
        for pch in puncta_chs:
            ch_name = pch["name"]
            ch_cfg = resolve_per_channel_cfg(cfg, ch_name)
            # Per-channel puncta_compartment can shift which mask is used
            if ch_cfg.get("puncta_compartment") != compartment:
                ch_compartment_mask = make_compartment_mask(
                    cell_mask, nuc_mask,
                    ch_cfg.get("puncta_compartment", compartment), cfg,
                    nucleolar_mask=nucleolar_mask)
            else:
                ch_compartment_mask = compartment_mask
            puncta_masks[ch_name] = detect_puncta(
                images[ch_name], ch_compartment_mask, ch_cfg)

        # Fragmentation indices per puncta channel (NEW)
        fragmentation_per_channel: dict[str, dict[int, tuple[int, float]]] = {}
        for pch in puncta_chs:
            ch_name = pch["name"]
            ch_cfg = resolve_per_channel_cfg(cfg, ch_name)
            fragmentation_per_channel[ch_name] = compute_fragmentation_indices(
                images[ch_name], cell_mask, ch_cfg)

        # Per-cell metrics (now consume fragmentation results)
        cells_df = per_cell_metrics(
            p.name, images, cell_mask, nuc_mask, puncta_masks,
            cell_to_nucs, meta, channels, cfg,
            fragmentation_per_channel=fragmentation_per_channel,
            nucleolar_mask=nucleolar_mask)
        img_df = per_image_summary(cells_df, puncta_chs)
        img_df.insert(0, "image", p.name)
        if meta["condition"]:
            img_df.insert(1, "condition", meta["condition"])
        if meta["replicate"]:
            img_df.insert(2, "replicate", meta["replicate"])

        all_cells.append(cells_df)
        all_imgs.append(img_df)

        # Colocalization
        if do_coloc:
            coloc_compartment = cfg.get(
                "colocalization_compartment", "whole-cell")
            coloc_df = compute_colocalization(
                images, coloc_channels, cell_mask, nuc_mask,
                coloc_compartment, cfg, meta, p.name)
            if cfg.get("_coloc_projection_derived"):
                coloc_df["projection_derived"] = True
            all_coloc.append(coloc_df)

        # Nucleolar proximity (cfg passed for 3D anisotropic EDT)
        puncta_dist_colors: dict[str, dict[int, str]] | None = None
        if nuc_prox_ch_name and nucleolar_mask.max() > 0:
            puncta_dist_colors = {}
            for pch in puncta_chs:
                prox_df = compute_nucleolar_proximity(
                    puncta_masks[pch["name"]], nucleolar_mask, cell_mask,
                    proximity_threshold, meta, p.name, pch["name"], cfg=cfg,
                    extra_thresholds=proximity_extra_thresholds)
                all_proximity.append(prox_df)
                puncta_dist_colors[pch["name"]] = compute_puncta_distance_colors(
                    puncta_masks[pch["name"]], nucleolar_mask,
                    proximity_threshold, cfg=cfg)

        # QC rendering
        save_qc_png(
            qc_dir / f"{p.stem}_qc.png",
            images, channels, cell_mask, nuc_mask, puncta_masks, cfg,
            nucleolar_mask=nucleolar_mask if nucleolar_mask.max() > 0 else None,
            puncta_distance_colors=puncta_dist_colors)

        # Save masks
        if bool(cfg["save_masks"]):
            tiff.imwrite(
                str(mask_dir / f"{p.stem}_cellmask.tif"),
                cell_mask.astype(np.uint16))
            if has_nuclei:
                tiff.imwrite(
                    str(mask_dir / f"{p.stem}_nucmask.tif"),
                    nuc_mask.astype(np.uint16))
            if nucleolar_mask.max() > 0:
                tiff.imwrite(
                    str(mask_dir / f"{p.stem}_nucleolarmask.tif"),
                    nucleolar_mask.astype(np.uint8))
            for ch_name, pmask in puncta_masks.items():
                tiff.imwrite(
                    str(mask_dir / f"{p.stem}_{ch_name}_punctamask.tif"),
                    pmask.astype(np.uint16))

        dt = time.time() - t0
        if not cells_df.empty:
            keep_df = cells_df[cells_df["keep"]]
            parts = [
                f"{dt:.1f}s",
                f"cells: {len(cells_df)}",
                f"keep: {len(keep_df)}",
            ]
            for pch in puncta_chs:
                col = f"{pch['name']}_puncta_n"
                if col in keep_df.columns and len(keep_df):
                    parts.append(
                        f"median {pch['name']} puncta: "
                        f"{float(keep_df[col].median())}")
            print(f"  {' | '.join(parts)}")
        else:
            print(f"  {dt:.1f}s | cells: 0 (segmentation likely failed)")

    # -------------------------------------------------------------------
    # Aggregate and save CSVs
    # -------------------------------------------------------------------
    cells_out = (pd.concat(all_cells, ignore_index=True)
                 if all_cells else pd.DataFrame())
    imgs_out = (pd.concat(all_imgs, ignore_index=True)
                if all_imgs else pd.DataFrame())

    cells_csv = out_dir / "cells.csv"
    imgs_csv = out_dir / "images.csv"
    cells_out.to_csv(cells_csv, index=False)
    imgs_out.to_csv(imgs_csv, index=False)

    # Write colocalization CSV
    if all_coloc:
        coloc_out = pd.concat(all_coloc, ignore_index=True)
        coloc_csv = out_dir / "colocalization.csv"
        coloc_out.to_csv(coloc_csv, index=False)
        print(f"  {coloc_csv}")

        # Pivot colocalization metrics into cells_out
        if not coloc_out.empty and not cells_out.empty:
            for pair_name in coloc_out["pair"].unique():
                pair_df = coloc_out[coloc_out["pair"] == pair_name]
                for metric_col in ["pearson_r", "manders_m1", "manders_m2"]:
                    pivot_col = f"{metric_col}_{pair_name}"
                    merge_df = pair_df[["image", "cell_id", metric_col]].rename(
                        columns={metric_col: pivot_col})
                    cells_out = cells_out.merge(
                        merge_df, on=["image", "cell_id"], how="left")

    # Write nucleolar proximity CSV
    if all_proximity:
        prox_out = pd.concat(all_proximity, ignore_index=True)
        prox_csv = out_dir / "nucleolar_proximity.csv"
        prox_out.to_csv(prox_csv, index=False)
        print(f"  {prox_csv}")

        # Pivot proximity metrics into cells_out
        if not prox_out.empty and not cells_out.empty:
            extra_prox_cols = sorted(
                c for c in prox_out.columns if c.startswith("fraction_proximal_"))
            for ch_name in prox_out["channel"].unique():
                ch_df = prox_out[prox_out["channel"] == ch_name]
                for metric_col in ["mean_distance", "fraction_proximal"] + extra_prox_cols:
                    pivot_col = f"{ch_name}_{metric_col}"
                    merge_df = ch_df[["image", "cell_id", metric_col]].rename(
                        columns={metric_col: pivot_col})
                    cells_out = cells_out.merge(
                        merge_df, on=["image", "cell_id"], how="left")

    # Write nucleolar morphology CSV
    if all_nol_morph:
        nol_morph_out = pd.concat(all_nol_morph, ignore_index=True)
        nol_morph_csv = out_dir / "nucleolar_morphology.csv"
        nol_morph_out.to_csv(nol_morph_csv, index=False)
        print(f"  {nol_morph_csv}")

        # Pivot morphology metrics into cells_out (column set differs by mode)
        if not nol_morph_out.empty and not cells_out.empty:
            if is_3d:
                morph_cols = [
                    "nucleolar_volume_vox", "nucleolar_volume_um3",
                    "nucleolar_eq_diameter_um", "n_nucleoli",
                ]
            else:
                morph_cols = [
                    "nucleolar_area", "nucleolar_solidity",
                    "nucleolar_circularity", "nucleolar_eccentricity",
                    "n_nucleoli",
                ]
            present = [c for c in morph_cols
                       if c in nol_morph_out.columns]
            merge_df = nol_morph_out[["image", "cell_id"] + present]
            cells_out = cells_out.merge(
                merge_df, on=["image", "cell_id"], how="left")

    # Re-save cells.csv with pivoted spatial metrics
    if all_coloc or all_proximity or all_nol_morph:
        cells_out.to_csv(cells_csv, index=False)

    dt_total = time.time() - t_total
    print(f"\n=== Segmentation + quantification done ({dt_total:.1f}s) ===")
    print(f"  {cells_csv}")
    print(f"  {imgs_csv}")

    # -------------------------------------------------------------------
    # Superplots + Prism export (on kept cells with valid conditions)
    # -------------------------------------------------------------------
    if cfg.get("skip_plots", False):
        print("\nSkipping plots (--skip-plots).")
        return

    if cells_out.empty or "condition" not in cells_out.columns:
        print("\nSkipping plots (no cells or no condition column).")
        return

    # Determine condition order
    condition_order: list[str] = cfg.get("condition_order", [])
    if not condition_order:
        condition_order = sorted(
            cells_out["condition"].dropna().unique().tolist())

    plot_df = cells_out[
        cells_out["keep"] & cells_out["condition"].isin(condition_order)
    ].copy()

    if plot_df.empty:
        print(f"\nSkipping plots (no kept cells in {condition_order}).")
        return

    # Build dynamic metric list (including coloc and proximity)
    coloc_pairs: list[tuple[str, str]] | None = None
    if do_coloc and all_coloc:
        coloc_pairs = [(a["name"], b["name"])
                       for a, b in combinations(coloc_channels, 2)]

    proximity_ch_names: list[str] | None = None
    if nuc_prox_ch_name and all_proximity:
        proximity_ch_names = [pch["name"] for pch in puncta_chs]

    metric_defs = _build_superplot_metrics(
        channels, puncta_chs,
        coloc_pairs=coloc_pairs,
        proximity_channels=proximity_ch_names,
        has_nucleolar_morph=bool(all_nol_morph),
        mode=cfg["mode"],
        cfg=cfg)
    avail = [(c, y, t) for c, y, t in metric_defs if c in plot_df.columns]
    missing = [c for c, _, _ in metric_defs if c not in plot_df.columns]
    if missing:
        print(f"\n[warn] Missing columns for plotting: {missing}")

    # Generate superplots
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    ref_cond = cfg.get("reference_condition")
    # For 2-condition datasets, default to using the first condition as
    # reference so p-values are always computed and written to CSV.
    if not ref_cond and len(condition_order) == 2:
        ref_cond = condition_order[0]
    show_trend = bool(cfg.get("trend", False))

    # Compute pairwise p-values per metric
    all_pval_rows: list[dict] = []
    for col, y_label, title in avail:
        pvalues: dict | None = None
        if ref_cond and ref_cond in condition_order:
            pvalues = compute_pairwise_pvalues(
                plot_df, col, condition_order, ref_cond)
            for cond, info in pvalues.items():
                all_pval_rows.append({
                    "metric": col,
                    "condition": cond,
                    "reference": ref_cond,
                    "pval": info["pval"],
                    "pval_corrected": info["pval_corrected"],
                    "n_test_replicates": info["n_test"],
                    "n_ref_replicates": info["n_ref"],
                    "n_comparisons": info["n_comparisons"],
                })
        superplot_violin(
            plot_df, col, plot_dir / f"{col}_superplot.png",
            y_label, title, condition_order, show_trend=show_trend,
            pvalues=pvalues, reference_condition=ref_cond)

    print(f"\n  Wrote {len(avail)} superplots to {plot_dir}/")

    # Write p-value summary CSV
    if all_pval_rows:
        pval_df = pd.DataFrame(all_pval_rows)
        pval_csv = out_dir / "pvalues.csv"
        pval_df.to_csv(pval_csv, index=False)
        print(f"  {pval_csv}")

    # Prism-ready CSVs
    prism_dir = out_dir / "prism"
    write_prism_outputs(
        plot_df, [c for c, _, _ in avail], prism_dir, condition_order)
    print(f"  Wrote Prism CSVs to {prism_dir}/")
    print(f"  {qc_dir}/")


if __name__ == "__main__":
    main()
