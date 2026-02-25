#!/usr/bin/env python3
"""
cellquant.py

Configurable multi-channel cell quantification pipeline: segmentation,
puncta detection, quantification, superplots, and Prism-ready CSVs.

Generalized from condensates_cluster.py to work with arbitrary channels,
cell types, and puncta targets — all driven by CLI args and optional YAML config.

Run:
  python cellquant.py --images MIPs/ --out results/ \
      --channels "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
      --puncta-channels G3BP1

Abani's original 3-channel pipeline (equivalent invocation):
  python cellquant.py --images Abani_SG_MIPs/ --out results/ \
      --channels "1:DAPI:nucleus" "2:PABPC1:quantify" "3:G3BP1:quantify" \
      --puncta-channels G3BP1 PABPC1 \
      --filename-pattern "MAX_{condition}_rep{replicate}" \
      --condition-order control arsenite
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

from scipy.ndimage import distance_transform_edt
from scipy.stats import mannwhitneyu
from skimage import filters, measure, morphology, segmentation
from skimage.transform import resize
from cellpose import models

try:
    import torch
except Exception:
    torch = None


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
        "keep_min_nuclei": 1,
        "keep_max_nuclei": 4,
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
    },
}


# ---------------------------------------------------------------------------
# Hardcoded defaults
# ---------------------------------------------------------------------------
DEFAULTS: dict[str, Any] = {
    # I/O
    "exts": [".tif", ".tiff"],
    "file_pattern": None,
    "save_masks": True,

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
    "puncta_compartment": "cytosol",  # "cytosol" | "nucleus" | "whole-cell"

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


def build_config(args: argparse.Namespace) -> dict[str, Any]:
    cfg = dict(DEFAULTS)

    # Layer 2: cell-type preset
    if args.cell_type:
        preset_name = args.cell_type.lower()
        if preset_name not in CELL_TYPE_PRESETS:
            raise ValueError(
                f"Unknown cell type '{preset_name}'. "
                f"Available: {list(CELL_TYPE_PRESETS.keys())}")
        cfg.update(CELL_TYPE_PRESETS[preset_name])

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
        "qc_downsample": args.qc_downsample,
        "qc_dpi": args.qc_dpi,
        "keep_min_nuclei": args.keep_min_nuclei,
        "keep_max_nuclei": args.keep_max_nuclei,
        "pretrained_model": args.pretrained_model,
        "cpu_threads": args.cpu_threads,
        "filename_pattern": args.filename_pattern,
        "condition_order": args.condition_order,
        "cell_seg_channel": args.cell_seg_channel,
        "file_pattern": args.file_pattern,
        "min_cell_area": args.min_cell_area,
        "max_cell_area": args.max_cell_area,
        "colocalization_compartment": args.colocalization_compartment,
        "nucleolar_proximity_channel": args.nucleolar_proximity,
        "proximity_threshold_px": args.proximity_threshold,
    }
    for key, val in cli_map.items():
        if val is not None:
            cfg[key] = val

    # Boolean flags (store_true — only apply when set)
    if args.no_gpu:
        cfg["use_gpu"] = False
    if args.no_save_masks:
        cfg["save_masks"] = False
    if args.skip_plots:
        cfg["skip_plots"] = True
    if args.colocalization:
        cfg["colocalization"] = True
    if args.trend:
        cfg["trend"] = True

    # Condition map from CLI (parsed from key=value pairs)
    cmap = _parse_condition_map(args.condition_map)
    if cmap is not None:
        cfg["condition_map"] = cmap

    return cfg


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
        # If {replicate} is trailing and digits are absent, default to "1"
        fallback = _build_filename_regex_no_rep(pattern)
        if fallback:
            m2 = fallback.search(stem)
            if m2:
                raw_cond = m2.group("condition")
                cmap = cfg.get("condition_map", {})
                info["condition"] = cmap.get(raw_cond.lower(), raw_cond) if raw_cond else ""
                info["replicate"] = "1"
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
    ap.add_argument("--puncta-compartment", default=None,
                    choices=["cytosol", "nucleus", "whole-cell"])

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
    ap.add_argument("--condition-map", nargs="+", default=None,
                    help='Map raw condition names: key=value pairs '
                         '(e.g. "ctrl=control" "as=arsenite")')

    # I/O options
    ap.add_argument("--file-pattern", default=None,
                    help='Glob pattern for image files (e.g. "*.tif")')
    ap.add_argument("--no-save-masks", action="store_true",
                    help="Disable writing segmentation/puncta masks to disk")
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

    # Threading
    ap.add_argument("--cpu-threads", type=int, default=None)

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
    # max_size removes objects <= threshold; subtract 1 to match old min_size
    # semantics (remove objects strictly < min_size, i.e. keep >= min_size)
    return morphology.remove_small_objects(bw, max_size=min_size - 1)


def remove_small_holes_compat(bw: np.ndarray, area_threshold: int) -> np.ndarray:
    # max_size removes holes <= threshold; old area_threshold removed holes
    # strictly < threshold, so subtract 1 to preserve semantics
    return morphology.remove_small_holes(bw, max_size=area_threshold - 1)


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


def load_tiff(path: Path, channels: list[dict]) -> dict[str, np.ndarray]:
    """Load a multi-channel TIFF and return {channel_name: 2D array}."""
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
# Cellpose init + eval (reused verbatim, parameterized on cfg)
# ---------------------------------------------------------------------------
def _is_mps_system() -> bool:
    """Check if running on Apple Silicon (MPS backend)."""
    if torch is None:
        return False
    return hasattr(torch.backends, "mps") and torch.backends.mps.is_available()


def init_model(cfg: dict) -> models.CellposeModel:
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


def segment_nuclei(
    model: models.CellposeModel, nuc_img: np.ndarray, cfg: dict
) -> np.ndarray:
    return eval_masks(model, nuc_img, _diam_or_none(cfg["nuclei_diameter"]), cfg)


def segment_cells(
    model: models.CellposeModel, cell_img: np.ndarray, cfg: dict
) -> np.ndarray:
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
def make_compartment_mask(
    cell_mask: np.ndarray,
    nuc_mask: np.ndarray,
    compartment: str,
    cfg: dict,
) -> np.ndarray:
    if compartment == "cytosol":
        nuc_dil = morphology.dilation(
            (nuc_mask > 0),
            footprint=morphology.disk(int(cfg["nucleus_dilate_px"])))
        return (cell_mask > 0) & (~nuc_dil)
    elif compartment == "nucleus":
        return (nuc_mask > 0)
    elif compartment == "whole-cell":
        return (cell_mask > 0)
    else:
        raise ValueError(f"Unknown puncta_compartment: {compartment}")


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


# ---------------------------------------------------------------------------
# Nucleolar segmentation
# ---------------------------------------------------------------------------
def segment_nucleoli(
    image: np.ndarray, cell_mask: np.ndarray, cfg: dict
) -> np.ndarray:
    """Per-cell Otsu thresholding of nucleolus channel, morphological opening,
    size filtering. Returns combined binary mask (0/1) for the whole image."""
    opening_radius = int(cfg.get("nucleolar_opening_radius", 1))
    min_area = int(cfg.get("nucleolar_min_area_px", 5))
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
        if opening_radius > 0:
            bw = morphology.opening(
                bw, footprint=morphology.disk(opening_radius))
        bw = bw & roi  # clip back to cell
        # Size filter
        lab = measure.label(bw)
        for prop in measure.regionprops(lab):
            if prop.area >= min_area:
                combined[lab == prop.label] = 1

    return combined


def compute_nucleolar_morphology(
    nucleolar_mask: np.ndarray,
    cell_mask: np.ndarray,
    metadata: dict[str, str],
    image_name: str,
) -> pd.DataFrame:
    """Per-cell nucleolar morphometrics (largest nucleolus by area)."""
    rows: list[dict[str, Any]] = []
    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    for cid in cell_ids:
        cell_bin = (cell_mask == cid)
        nol_in_cell = nucleolar_mask * cell_bin
        lab = measure.label(nol_in_cell > 0)
        props = measure.regionprops(lab)
        n_nucleoli = len(props)

        if n_nucleoli == 0:
            rows.append({
                "image": image_name,
                "condition": metadata.get("condition", ""),
                "replicate": metadata.get("replicate", ""),
                "cell_id": int(cid),
                "nucleolar_area": 0,
                "nucleolar_solidity": np.nan,
                "nucleolar_circularity": np.nan,
                "nucleolar_eccentricity": np.nan,
                "n_nucleoli": 0,
            })
            continue

        # Largest nucleolus by area
        largest = max(props, key=lambda p: p.area)
        per = measure.perimeter_crofton(largest.image, directions=4)
        circ = (4.0 * np.pi * largest.area) / (per * per) if per > 0 else np.nan
        sol = largest.area / largest.area_convex if largest.area_convex > 0 else np.nan

        rows.append({
            "image": image_name,
            "condition": metadata.get("condition", ""),
            "replicate": metadata.get("replicate", ""),
            "cell_id": int(cid),
            "nucleolar_area": int(largest.area),
            "nucleolar_solidity": float(sol),
            "nucleolar_circularity": float(circ),
            "nucleolar_eccentricity": float(largest.eccentricity),
            "n_nucleoli": n_nucleoli,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Puncta detection (parameterized version of reference)
# ---------------------------------------------------------------------------
def detect_puncta(
    image: np.ndarray, compartment_mask: np.ndarray, cfg: dict
) -> np.ndarray:
    g = robust_rescale(image, cfg["p_low"], cfg["p_high"])
    g = g * (compartment_mask > 0)

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
    bw = remove_small_objects_compat(bw, min_size=int(cfg["puncta_min_area_px"]))
    bw = remove_small_holes_compat(bw, area_threshold=8)

    lab = measure.label(bw)
    if lab.max() == 0:
        return lab.astype(np.int32, copy=False)

    keep = np.zeros(lab.max() + 1, dtype=bool)
    minA = int(cfg["puncta_min_area_px"])
    maxA = int(cfg["puncta_max_area_px"])
    do_round = bool(cfg.get("puncta_filter_round", False))
    min_circ = float(cfg.get("puncta_min_circularity", 0.0))
    min_sol = float(cfg.get("puncta_min_solidity", 0.0))

    for prop in measure.regionprops(lab):
        if not (minA <= prop.area <= maxA):
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
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    nuc_union = (nuc_mask > 0)
    nuc_dil = morphology.dilation(
        nuc_union, footprint=morphology.disk(int(cfg["nucleus_dilate_px"])))

    non_skip = [ch for ch in channels if ch["role"] != "skip"]

    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    compartment = cfg.get("puncta_compartment", "cytosol")

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
            "cell_area_px": int(cell_bin.sum()),
            "nucleus_area_px": int(nuc_bin.sum()),
            "cytosol_area_px": int(cyt_bin.sum()),
        }

        # Per non-skip channel: mean intensities + condensate index
        for ch in non_skip:
            name = ch["name"]
            img = images[name]
            row[f"{name}_cell_mean"] = safe_mean(img, cell_bin)
            row[f"{name}_nucleus_mean"] = safe_mean(img, nuc_bin)
            row[f"{name}_cytosol_mean"] = safe_mean(img, cyt_bin)

            m_cell = row[f"{name}_cell_mean"]
            p95_cell = _safe_p95(img, cell_bin)
            row[f"{name}_condensate_index_cell"] = _condensate_index(
                p95_cell, m_cell)

            m_cyt = row[f"{name}_cytosol_mean"]
            p95_cyt = _safe_p95(img, cyt_bin)
            row[f"{name}_condensate_index_cytosol"] = _condensate_index(
                p95_cyt, m_cyt)

        # Per puncta channel: puncta-specific metrics
        for ch_name, p_mask in puncta_masks.items():
            ch_img = images[ch_name]
            puncta_in_cell = p_mask * cell_bin
            puncta_bin = (puncta_in_cell > 0)
            n_puncta = len(np.unique(puncta_in_cell)) - 1  # subtract background 0

            ch_puncta_vals = ch_img[puncta_bin]
            puncta_int = (
                float(ch_puncta_vals.sum()) if ch_puncta_vals.size else 0.0)

            # Compartment for frac_intensity denominator
            if compartment == "cytosol":
                comp_bin = cyt_bin
            elif compartment == "nucleus":
                comp_bin = nuc_bin
            else:
                comp_bin = cell_bin

            comp_total = (
                float(ch_img[comp_bin].sum()) if np.any(comp_bin) else 0.0)

            diffuse_bin = comp_bin & (~puncta_bin)

            row[f"{ch_name}_puncta_n"] = n_puncta
            row[f"{ch_name}_puncta_area_px"] = int(puncta_bin.sum())
            row[f"{ch_name}_puncta_mean_intensity"] = (
                float(ch_puncta_vals.mean()) if ch_puncta_vals.size else np.nan)
            row[f"{ch_name}_diffuse_mean_intensity"] = safe_mean(
                ch_img, diffuse_bin)
            row[f"{ch_name}_frac_intensity_in_puncta"] = (
                puncta_int / comp_total if comp_total > 0 else np.nan)
            row[f"{ch_name}_puncta_integrated_intensity"] = puncta_int

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
) -> pd.DataFrame:
    """Per-cell puncta distance to nearest nucleolar boundary."""
    rows: list[dict[str, Any]] = []
    dist_field = distance_transform_edt(~(nucleolar_mask > 0))

    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    for cid in cell_ids:
        cell_bin = (cell_mask == cid)
        puncta_in_cell = puncta_mask * cell_bin
        punct_ids = np.unique(puncta_in_cell)
        punct_ids = punct_ids[punct_ids != 0]
        n_puncta = len(punct_ids)
        if n_puncta == 0:
            rows.append({
                "image": image_name,
                "condition": metadata.get("condition", ""),
                "replicate": metadata.get("replicate", ""),
                "cell_id": int(cid),
                "channel": channel_name,
                "n_puncta": 0,
                "mean_distance": np.nan,
                "min_distance": np.nan,
                "fraction_proximal": np.nan,
            })
            continue

        distances = []
        for pid in punct_ids:
            ys, xs = np.where(puncta_in_cell == pid)
            cy, cx = ys.mean(), xs.mean()  # centroid
            d = dist_field[int(round(cy)), int(round(cx))]
            distances.append(float(d))

        distances_arr = np.array(distances)
        rows.append({
            "image": image_name,
            "condition": metadata.get("condition", ""),
            "replicate": metadata.get("replicate", ""),
            "cell_id": int(cid),
            "channel": channel_name,
            "n_puncta": n_puncta,
            "mean_distance": float(distances_arr.mean()),
            "min_distance": float(distances_arr.min()),
            "fraction_proximal": float((distances_arr <= threshold).mean()),
        })

    return pd.DataFrame(rows)


def compute_puncta_distance_colors(
    puncta_mask: np.ndarray,
    nucleolar_mask: np.ndarray,
    threshold: float,
) -> dict[int, str]:
    """Map each punctum label to 'red' (proximal) or 'blue' (distal)."""
    dist_field = distance_transform_edt(~(nucleolar_mask > 0))
    colors: dict[int, str] = {}
    punct_ids = np.unique(puncta_mask)
    punct_ids = punct_ids[punct_ids != 0]
    for pid in punct_ids:
        ys, xs = np.where(puncta_mask == pid)
        cy, cx = ys.mean(), xs.mean()
        d = dist_field[int(round(cy)), int(round(cx))]
        colors[int(pid)] = "red" if d <= threshold else "blue"
    return colors


# ---------------------------------------------------------------------------
# QC rendering — dynamic channels
# ---------------------------------------------------------------------------
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
    """Build RGB composite (R=quant1, G=quant2, B=nucleus) + contour overlay.

    Overlays: cyan = cell, yellow = nucleus, white = nucleolus,
    magenta/red/blue = puncta.
    """
    quantify_chs = get_quantify_channels(channels)
    nuc_ch = get_nucleus_channel(channels)
    ref_img = next(iter(images.values()))
    zeros = np.zeros_like(ref_img, dtype=np.float32)

    # R channel: first quantify
    r = (robust_rescale(images[quantify_chs[0]["name"]], cfg["p_low"], cfg["p_high"])
         if len(quantify_chs) >= 1 else zeros)
    # G channel: second quantify
    g = (robust_rescale(images[quantify_chs[1]["name"]], cfg["p_low"], cfg["p_high"])
         if len(quantify_chs) >= 2 else zeros)
    # B channel: nucleus (or zeros if no nucleus channel)
    b = (robust_rescale(images[nuc_ch["name"]], cfg["p_low"], cfg["p_high"])
         if nuc_ch is not None else zeros)

    f = int(cfg["qc_downsample"])
    rgb = np.stack(
        [downsample_nn(r, f), downsample_nn(g, f), downsample_nn(b, f)],
        axis=-1)
    rgb = np.clip(rgb, 0.0, 1.0)

    cell_ds = downsample_nn(cell_mask, f)
    nuc_ds = downsample_nn(nuc_mask, f)

    h, w = rgb.shape[:2]
    dpi = int(cfg["qc_dpi"])
    fig, ax = plt.subplots(1, 1, figsize=(w / dpi, h / dpi), dpi=dpi)
    ax.imshow(rgb)
    ax.set_axis_off()

    for cid in np.unique(cell_ds):
        if cid == 0:
            continue
        ax.contour(cell_ds == cid, levels=[0.5], colors="cyan", linewidths=1.0)

    # Nuclear contours (skip when no nuclei)
    if nuc_mask.max() > 0:
        for nid in np.unique(nuc_ds):
            if nid == 0:
                continue
            ax.contour(nuc_ds == nid, levels=[0.5], colors="yellow", linewidths=0.8)

    # Nucleolar contours (white)
    if nucleolar_mask is not None and nucleolar_mask.max() > 0:
        nol_ds = downsample_nn(nucleolar_mask, f)
        ax.contour(nol_ds > 0, levels=[0.5], colors="white", linewidths=0.6)

    # Puncta contours
    for ch_name, pmask in puncta_masks.items():
        pmask_ds = downsample_nn(pmask, f)
        if pmask_ds.max() == 0:
            continue
        if puncta_distance_colors and ch_name in puncta_distance_colors:
            # Color each punctum individually (red=proximal, blue=distal)
            colors_map = puncta_distance_colors[ch_name]
            for pid in np.unique(pmask_ds):
                if pid == 0:
                    continue
                color = colors_map.get(int(pid), "magenta")
                ax.contour(
                    pmask_ds == pid, levels=[0.5],
                    colors=color, linewidths=0.6)
        else:
            ax.contour(
                pmask_ds > 0, levels=[0.5], colors="magenta", linewidths=0.6)

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0, dpi=dpi)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Superplots — dynamic metrics + condition order
# ---------------------------------------------------------------------------
def _build_superplot_metrics(
    channels: list[dict],
    puncta_chs: list[dict],
    coloc_pairs: list[tuple[str, str]] | None = None,
    proximity_channels: list[str] | None = None,
    has_nucleolar_morph: bool = False,
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

    for ch in channels:
        if ch["role"] == "skip":
            continue
        n = ch["name"]
        metrics.append((
            f"{n}_condensate_index_cell",
            f"Condensate index (p95/mean) — {n} (cell)",
            f"Condensate index (cell) — {n}"))

    # Colocalization metrics
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

    # Proximity metrics
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

    # Nucleolar morphology metrics
    if has_nucleolar_morph:
        metrics.append((
            "nucleolar_area",
            "Nucleolar area (px)",
            "Nucleolar area"))
        metrics.append((
            "nucleolar_solidity",
            "Nucleolar solidity",
            "Nucleolar solidity"))
        metrics.append((
            "nucleolar_circularity",
            "Nucleolar circularity",
            "Nucleolar circularity"))
        metrics.append((
            "nucleolar_eccentricity",
            "Nucleolar eccentricity",
            "Nucleolar eccentricity"))

    return metrics


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

    # Wilcoxon rank-sum on replicate medians (≥3 per condition required)
    if len(condition_order) == 2:
        meds_a = rep_meds.loc[
            rep_meds["condition"] == condition_order[0], "rep_median"
        ].dropna().values
        meds_b = rep_meds.loc[
            rep_meds["condition"] == condition_order[1], "rep_median"
        ].dropna().values
        if len(meds_a) >= 3 and len(meds_b) >= 3:
            _, pval = mannwhitneyu(meds_a, meds_b, alternative="two-sided")
            p_text = f"p < 0.001" if pval < 0.001 else f"p = {pval:.2g}"
            # Bracket annotation
            ax = plt.gca()
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

    for i, cond in enumerate(condition_order):
        vals = (df.loc[df["condition"] == cond, metric]
                .dropna().astype(float).values)
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
    **kwargs,
) -> None:
    """Dispatch: violin for <=2 conditions (no trend), strip plot otherwise."""
    if len(condition_order) <= 2 and not show_trend:
        _superplot_violin_2(df, metric, out_png, y_label, title,
                            condition_order, **kwargs)
    else:
        _superplot_strip_multi(df, metric, out_png, y_label, title,
                               condition_order, show_trend=show_trend)


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

    # Save resolved config
    save_cfg: dict[str, Any] = dict(cfg)
    save_cfg["_channels"] = [
        f"{ch['position']}:{ch['name']}:{ch['role']}" for ch in channels]
    if puncta_chs:
        save_cfg["_puncta_channels"] = [ch["name"] for ch in puncta_chs]
    save_cfg["_cell_seg_channel"] = (
        cell_seg_ch["name"] if cell_seg_ch else "composite")
    with open(out_dir / "config_used.yml", "w") as fh:
        yaml.dump(save_cfg, fh, default_flow_style=False, sort_keys=False)

    # Find images
    paths = find_images(img_dir, cfg)
    if not paths:
        raise SystemExit(f"No TIFFs found in {img_dir}")

    print(f"Found {len(paths)} images in {img_dir}")
    print(f"GPU: {cfg['use_gpu']}")
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

    model = init_model(cfg)

    all_cells: list[pd.DataFrame] = []
    all_imgs: list[pd.DataFrame] = []
    all_coloc: list[pd.DataFrame] = []
    all_proximity: list[pd.DataFrame] = []
    all_nol_morph: list[pd.DataFrame] = []
    t_total = time.time()

    proximity_threshold = float(cfg.get("proximity_threshold_px", 5))

    for i, p in enumerate(paths, 1):
        print(f"\n=== [{i}/{len(paths)}] {p.name} ===")
        t0 = time.time()

        # Load image → dict of channel arrays
        images = load_tiff(p, channels)

        # Determine cell seg input
        if use_composite_seg:
            cell_seg_img = build_composite_seg_image(images, channels, cfg)
        else:
            cell_seg_img = images[cell_seg_ch["name"]]
        ref_shape = cell_seg_img.shape  # (H, W)

        # Downsample before Cellpose, then upsample masks back
        ds = max(1, int(cfg.get("seg_downsample", 1)))

        # Nuclear segmentation (skipped when no nucleus channel)
        if has_nuclei:
            nuc_img = images[nuc_ch["name"]]
            nuc_seg_in = nuc_img[::ds, ::ds] if ds > 1 else nuc_img
            nuc_small = segment_nuclei(model, nuc_seg_in, cfg)
            nuc_mask = (upsample_labels_nn(nuc_small, ref_shape)
                        if ds > 1 else nuc_small)
        else:
            nuc_mask = np.zeros(ref_shape, dtype=np.int32)

        # Cell segmentation
        cell_seg_in = cell_seg_img[::ds, ::ds] if ds > 1 else cell_seg_img
        cell_small = segment_cells(model, cell_seg_in, cfg)
        cell_mask = (upsample_labels_nn(cell_small, ref_shape)
                     if ds > 1 else cell_small)

        # Cell area filtering
        min_area = int(cfg.get("min_cell_area", 0))
        max_area = int(cfg.get("max_cell_area", 0))
        cell_mask = filter_cells_by_area(cell_mask, min_area, max_area)

        cell_to_nucs = map_nuclei_to_cells(nuc_mask, cell_mask)

        # Parse metadata from filename
        meta = parse_filename_metadata(p.stem, cfg)

        # Nucleolar segmentation
        nucleolar_mask = np.zeros(ref_shape, dtype=np.uint8)
        for nol_ch in nucleolus_chs:
            nucleolar_mask |= segment_nucleoli(
                images[nol_ch["name"]], cell_mask, cfg)

        # Nucleolar morphometrics
        if nucleolar_mask.max() > 0:
            nol_morph_df = compute_nucleolar_morphology(
                nucleolar_mask, cell_mask, meta, p.name)
            all_nol_morph.append(nol_morph_df)

        # Create compartment mask for puncta detection
        compartment_mask = make_compartment_mask(
            cell_mask, nuc_mask, compartment, cfg)

        # Detect puncta per channel
        puncta_masks: dict[str, np.ndarray] = {}
        for pch in puncta_chs:
            puncta_masks[pch["name"]] = detect_puncta(
                images[pch["name"]], compartment_mask, cfg)

        # Compute per-cell metrics
        cells_df = per_cell_metrics(
            p.name, images, cell_mask, nuc_mask, puncta_masks,
            cell_to_nucs, meta, channels, cfg)
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
            all_coloc.append(coloc_df)

        # Nucleolar proximity
        puncta_dist_colors: dict[str, dict[int, str]] | None = None
        if nuc_prox_ch_name and nucleolar_mask.max() > 0:
            puncta_dist_colors = {}
            for pch in puncta_chs:
                prox_df = compute_nucleolar_proximity(
                    puncta_masks[pch["name"]], nucleolar_mask, cell_mask,
                    proximity_threshold, meta, p.name, pch["name"])
                all_proximity.append(prox_df)
                # Distance colors for QC
                puncta_dist_colors[pch["name"]] = compute_puncta_distance_colors(
                    puncta_masks[pch["name"]], nucleolar_mask,
                    proximity_threshold)

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
            for ch_name in prox_out["channel"].unique():
                ch_df = prox_out[prox_out["channel"] == ch_name]
                for metric_col in ["mean_distance", "fraction_proximal"]:
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

        # Pivot morphology metrics into cells_out
        if not nol_morph_out.empty and not cells_out.empty:
            morph_cols = ["nucleolar_area", "nucleolar_solidity",
                          "nucleolar_circularity", "nucleolar_eccentricity",
                          "n_nucleoli"]
            merge_df = nol_morph_out[["image", "cell_id"] + morph_cols]
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
        has_nucleolar_morph=bool(all_nol_morph))
    avail = [(c, y, t) for c, y, t in metric_defs if c in plot_df.columns]
    missing = [c for c, _, _ in metric_defs if c not in plot_df.columns]
    if missing:
        print(f"\n[warn] Missing columns for plotting: {missing}")

    # Generate superplots
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    show_trend = bool(cfg.get("trend", False))
    for col, y_label, title in avail:
        superplot_violin(
            plot_df, col, plot_dir / f"{col}_superplot.png",
            y_label, title, condition_order, show_trend=show_trend)

    print(f"\n  Wrote {len(avail)} superplots to {plot_dir}/")

    # Prism-ready CSVs
    prism_dir = out_dir / "prism"
    write_prism_outputs(
        plot_df, [c for c, _, _ in avail], prism_dir, condition_order)
    print(f"  Wrote Prism CSVs to {prism_dir}/")
    print(f"  {qc_dir}/")


if __name__ == "__main__":
    main()
