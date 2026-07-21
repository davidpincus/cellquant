"""validation_3d/00_setup_paths.py

Single source of truth for paths, image lists, and runtime parameters used
across the 3D validation scripts. Import from here; do not re-derive paths
elsewhere.

Run this file directly to print the resolved config and sanity-check that
every expected input is on disk.
"""
from __future__ import annotations

import json
import os
import platform
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Filesystem roots (portable — no hardcoded machine path)
# ---------------------------------------------------------------------------
# The project root is the directory containing cellquant.py, validation_3d/,
# and the data dirs (SG_zstacks/, Tif6_Nsr1_Sis1_6hr/). It is resolved at
# import time so the SAME harness file runs unmodified on both the Mac and
# Midway3 — no divergent hardcoded copies (the old Mac iCloud path and the
# cluster's /scratch/midway3/... path are both gone).
#
# Resolution order:
#   1. $CELLQUANT_VALIDATION_ROOT             (explicit override)
#   2. validation_3d/paths.local.json         ({"project_root": "/abs/path"})
#   3. repo-relative default                  (parent of this validation_3d/)
# The repo-relative default already works on both machines when data sits
# alongside the harness; the env var / json are for data kept elsewhere.
_ENV_VAR = "CELLQUANT_VALIDATION_ROOT"
_ROOT_ANCHORS = ("cellquant.py", "validation_3d")


def _resolve_project_root() -> Path:
    """Resolve the validation project root portably, failing loudly.

    Raises FileNotFoundError naming the resolved path, the env var to set, and
    which expected entry was missing if the resolved directory does not look
    like the project root. (Presence of the large data dirs is validated
    separately, and loudly, when they are actually globbed — see
    ``_glob_required``.)
    """
    here = Path(__file__).resolve().parent  # .../validation_3d
    local_json = here / "paths.local.json"

    env = os.environ.get(_ENV_VAR)
    if env:
        root, source = Path(env).expanduser(), f"${_ENV_VAR}"
    elif local_json.is_file():
        try:
            root = Path(str(json.loads(local_json.read_text())["project_root"]))
        except (ValueError, KeyError, TypeError) as exc:
            raise FileNotFoundError(
                f"{local_json} is malformed ({exc}); expected JSON like "
                '{"project_root": "/abs/path/to/AVCGTIA_Neferkara_Ali_Pincus"}.'
            ) from exc
        root, source = root.expanduser(), str(local_json)
    else:
        root, source = here.parent, "repo-relative default (parent of validation_3d/)"

    root = root.resolve()
    hint = (f"Set {_ENV_VAR}=/abs/path, or create {local_json.name} with "
            '{"project_root": "/abs/path"} (see paths.local.json.example).')
    if not root.is_dir():
        raise FileNotFoundError(
            f"Validation project root does not exist: {root} "
            f"(resolved from {source}). {hint}")
    for anchor in _ROOT_ANCHORS:
        if not (root / anchor).exists():
            raise FileNotFoundError(
                f"Validation project root {root} (resolved from {source}) is "
                f"missing expected entry '{anchor}' — this does not look like "
                f"the AVCGTIA project root. {hint}")
    return root


PROJECT_ROOT = _resolve_project_root()
VALIDATION_ROOT = PROJECT_ROOT / "validation_3d"

CELLQUANT_SCRIPT = PROJECT_ROOT / "cellquant.py"

# ---------------------------------------------------------------------------
# Mammalian inputs (HCT116 stress granules — Abani)
# ---------------------------------------------------------------------------
# 6 z-stacks: control_rep1..3 + arsenite_rep1..3. 5 are 10-slice (10, 3, 1200,
# 1192) ZCYX uint16; control_rep2 is 5-slice (half-depth anomaly — flag in
# report).
MAMMALIAN_ZSTACK_DIR = PROJECT_ROOT / "SG_zstacks"

# 7 MIPs from the published Fig 2 analysis. Naming MAX_{condition}_rep{N}.tif.
# Replicate IDs: arsenite_rep1..4 + control_rep1, control_rep4, control_rep5.
# 4-image overlap with z-stacks: arsenite_rep1..3 + control_rep1.
MAMMALIAN_MIP_DIR = PROJECT_ROOT / "SG_zstacks" / "MIPs"

# Per Methods: 1192x1200 px @ 0.094 µm XY; 10 slices @ 0.25 µm Z (mammalian).
MAMMALIAN_VOXEL_XY_UM = 0.094
MAMMALIAN_VOXEL_Z_UM = 0.25

MAMMALIAN_CHANNELS = ["1:DAPI:nucleus", "2:G3BP1:quantify", "3:PABPC1:quantify"]
MAMMALIAN_CELL_TYPE = "mammalian"
MAMMALIAN_FILENAME_PATTERN_3D = "{condition}_rep{replicate}"
MAMMALIAN_FILENAME_PATTERN_2D = "MAX_{condition}_rep{replicate}"

# ---------------------------------------------------------------------------
# Yeast inputs (temperature series — Asif)
# ---------------------------------------------------------------------------
# 86 per-series z-stacks exported from .sld via Fiji + a patched SlidebookReader.
# Each is (Z=71, C=3, Y=1200, X=1192) uint16, ~580 MB. Named
# {temp}_series1_rep{N}.tif. 30 of them have matching MIPs in MIPs/ — the
# manuscript subset (reps 1-N per temperature; 25C=6, 30C=5, 32C=5, 36C=6,
# 40C=8). Validation runs on this 30-image subset.
YEAST_ZSTACK_DIR = PROJECT_ROOT / "Tif6_Nsr1_Sis1_6hr" / "z-stacks"

# 30 pre-computed published MIP TIFFs (CYX). Shape (3, 1200, 1192) uint16.
YEAST_MIP_DIR: Path | None = PROJECT_ROOT / "Tif6_Nsr1_Sis1_6hr" / "MIPs"

# Voxel size from the files' OME PhysicalSize (authoritative): 1192x1200 px @
# 0.10571 µm XY; 71 slices @ 0.23 µm Z (anisotropy 2.18). The manuscript
# Methods' "0.094 µm XY / 0.1 µm Z" was wrong for yeast, and the original 3D
# run inherited that error via --voxel-size 0.094 0.1. Corrected here to match
# the OME metadata; cellquant's voxel ladder now confirms CLI==OME at run time.
YEAST_VOXEL_XY_UM = 0.10571
YEAST_VOXEL_Z_UM = 0.23

YEAST_CHANNELS = ["1:Tif6:quantify", "2:Nsr1:nucleolus", "3:Sis1:quantify"]
YEAST_CELL_TYPE = "yeast"
YEAST_TEMPERATURES = ["25C", "30C", "32C", "36C", "40C"]
# Both z-stacks and MIPs use the same naming: 25C_series1_rep1.tif
YEAST_FILENAME_PATTERN_3D = "{condition}_series1_rep{replicate}"
YEAST_FILENAME_PATTERN_2D = "{condition}_series1_rep{replicate}"

# ---------------------------------------------------------------------------
# Validation output directories
# ---------------------------------------------------------------------------
OUT_MAMMALIAN = VALIDATION_ROOT / "outputs_mammalian"
OUT_YEAST = VALIDATION_ROOT / "outputs_yeast"
OUT_COMBINED = VALIDATION_ROOT / "outputs_combined"

OUT_MAMM_3D = OUT_MAMMALIAN / "3d"
OUT_MAMM_MATCHED2D = OUT_MAMMALIAN / "matched_2d"
OUT_MAMM_PUBLISHED2D = OUT_MAMMALIAN / "published_2d"
OUT_MAMM_COMPARISONS = OUT_MAMMALIAN / "comparisons"

OUT_YEAST_3D = OUT_YEAST / "3d"
OUT_YEAST_MATCHED2D = OUT_YEAST / "matched_2d"
OUT_YEAST_PUBLISHED2D = OUT_YEAST / "published_2d"
OUT_YEAST_COMPARISONS = OUT_YEAST / "comparisons"

OUT_SYNTH = VALIDATION_ROOT / "synthetic_sphericity"
OUT_REPORT = VALIDATION_ROOT / "REPORT.html"
OUT_FAILURES = VALIDATION_ROOT / "failures.log"

# ---------------------------------------------------------------------------
# Runtime + execution hints
# ---------------------------------------------------------------------------
# Use CPU unless a real CUDA GPU is present. On Apple Silicon there is no
# CUDA (and cpsam Transformer ops aren't supported on MPS), so this resolves
# to True locally and the driver adds --no-gpu. On a CUDA box (e.g. Midway3)
# it resolves to False and cellquant uses the GPU — no manual toggle needed,
# so the identical harness runs on both machines.
def _detect_no_gpu() -> bool:
    try:
        import torch
        return not torch.cuda.is_available()
    except Exception:
        return True


NO_GPU = _detect_no_gpu()

# OMP/libomp on macOS conflicts — cellquant.py also sets this at startup;
# we set it here too in case a subprocess spawn doesn't inherit.
ENV_OVERRIDES = {"KMP_DUPLICATE_LIB_OK": "TRUE"}

# Per-image cellquant timeout (seconds). Yeast 3D on Apple Silicon CPU is the
# binding case: full method took 5h 43min on 25C_series1_rep1; stitch was
# still in Cellpose at 3h 51min when killed. Set headroom to 12h.
PER_IMAGE_TIMEOUT_SEC = 60 * 60 * 12


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _glob_required(directory: Path | None, pattern: str, desc: str) -> list[Path]:
    """Glob ``pattern`` in ``directory``, failing LOUDLY if the directory is
    missing or nothing matches.

    An empty result from these input-listing helpers always means a
    misconfigured data path, never a valid "no work to do" state. Returning
    ``[]`` silently once let a validation run process 0 of 30 images and still
    exit 0 — a whole axis lost, caught only by counting files on disk. Raise
    instead, naming the glob that matched nothing.
    """
    if directory is None:
        raise FileNotFoundError(
            f"{desc}: input directory is not configured (None); "
            f"cannot glob '{pattern}'.")
    if not directory.is_dir():
        raise FileNotFoundError(
            f"{desc}: input directory does not exist: {directory}")
    matches = sorted(directory.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"{desc}: glob '{pattern}' matched no files in {directory}")
    return matches


def list_mammalian_zstacks() -> list[Path]:
    return _glob_required(MAMMALIAN_ZSTACK_DIR, "*.tif", "mammalian z-stacks")


def list_mammalian_mips() -> list[Path]:
    return _glob_required(MAMMALIAN_MIP_DIR, "MAX_*.tif", "mammalian MIPs")


def list_yeast_sld_files() -> list[Path]:
    # Raw acquisition .sld may legitimately be absent after conversion to TIFF,
    # so this one stays tolerant of an empty result (it is a convert-if-present
    # probe, not a validation input list).
    return sorted(YEAST_ZSTACK_DIR.glob("*.sld"))


def list_yeast_zstacks_all() -> list[Path]:
    """All 86 exported yeast z-stacks (full temperature series)."""
    return _glob_required(YEAST_ZSTACK_DIR, "*.tif", "yeast z-stacks")


def list_yeast_mips() -> list[Path]:
    """The 30 published yeast MIP TIFFs (manuscript subset)."""
    return _glob_required(YEAST_MIP_DIR, "*.tif", "yeast MIPs")


def replicate_id_from_zstack(path: Path) -> str:
    """control_rep1.tif -> 'control_rep1'."""
    return path.stem


def replicate_id_from_mip(path: Path) -> str:
    """MAX_control_rep1.tif -> 'control_rep1'."""
    return path.stem.replace("MAX_", "", 1)


def replicate_id_from_yeast_zstack(path: Path) -> str:
    """25C_series1_rep1.tif -> '25C_series1_rep1' (same naming as MIPs)."""
    return path.stem


def replicate_id_from_yeast_mip(path: Path) -> str:
    return path.stem


def overlap_image_ids() -> list[str]:
    """IDs where both a z-stack and a published MIP exist."""
    z_ids = {replicate_id_from_zstack(p) for p in list_mammalian_zstacks()}
    m_ids = {replicate_id_from_mip(p) for p in list_mammalian_mips()}
    return sorted(z_ids & m_ids)


def yeast_overlap_image_ids() -> list[str]:
    """The 30 IDs that have both a z-stack and a published MIP."""
    z_ids = {p.stem for p in list_yeast_zstacks_all()}
    m_ids = {p.stem for p in list_yeast_mips()}
    return sorted(z_ids & m_ids)


def list_yeast_zstacks_subset() -> list[Path]:
    """Z-stacks for the 30 manuscript reps with published-MIP counterparts.

    This is the validation set: 25C 1-6, 30C 1-5, 32C 1-5, 36C 1-6, 40C 1-8.
    """
    overlap = set(yeast_overlap_image_ids())
    return [p for p in list_yeast_zstacks_all() if p.stem in overlap]


def list_yeast_zstacks_one_per_temp() -> list[Path]:
    """Single representative (rep1) per temperature from the manuscript subset.

    Used when full per-image runtime forbids the 30-image set; trades replicate
    power for getting one stack per condition through 3D within a day.
    """
    out: list[Path] = []
    seen: set[str] = set()
    for p in list_yeast_zstacks_subset():
        temp = p.stem.split("_", 1)[0]
        if temp not in seen:
            seen.add(temp)
            out.append(p)
    return out


# Yeast preset defaults to seg_3d_method="full" (3D cpsam) which is ~5 h/image
# on Apple Silicon CPU. Override to "stitch" for the validation runs so the
# batch is feasible. Documented in REPORT.html.
YEAST_SEG_3D_METHOD = "stitch"

# CPU-feasible 3D config for yeast (Part A/B/C/D outcome):
#   z-crop = 41 central slices around per-image peak signal Z (drops top/
#     bottom cell layers; cells straddling the crop boundary are excluded so
#     volumes/sphericity stay valid).
#   XY seg downsample = 2 (Cellpose runs on 596x600 images; masks upsampled
#     back to 1192x1200 for measurement).
YEAST_Z_CROP_CENTER = 41
YEAST_SEG_DOWNSAMPLE = 2


# ---------------------------------------------------------------------------
# Recon / sanity check (run this file directly)
# ---------------------------------------------------------------------------
def _hardware_report() -> None:
    print(f"platform        : {platform.platform()}")
    print(f"python          : {sys.version.split()[0]}")
    print(f"cwd             : {os.getcwd()}")
    try:
        import torch
        cuda = torch.cuda.is_available()
        mps = hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
        print(f"torch           : {torch.__version__}  (cuda={cuda}, mps={mps})")
    except Exception as exc:
        print(f"torch           : NOT importable ({exc})")
    try:
        import psutil  # may not be installed; soft fail
        print(f"cpu_count       : {psutil.cpu_count(logical=True)} "
              f"logical / {psutil.cpu_count(logical=False)} physical")
        vm = psutil.virtual_memory()
        print(f"ram             : {vm.total / 1024**3:.1f} GB total, "
              f"{vm.available / 1024**3:.1f} GB available")
    except Exception:
        pass


def _path_report() -> None:
    print()
    print("Project root    :", PROJECT_ROOT)
    print("Validation root :", VALIDATION_ROOT, "exists" if VALIDATION_ROOT.exists() else "MISSING")
    print("cellquant.py    :", CELLQUANT_SCRIPT, "exists" if CELLQUANT_SCRIPT.exists() else "MISSING")
    print()
    print("=== Mammalian ===")
    print(f"z-stacks dir    : {MAMMALIAN_ZSTACK_DIR}")
    zs = list_mammalian_zstacks()
    print(f"  found {len(zs)} z-stacks:")
    for p in zs:
        print(f"    {p.name}")
    print(f"MIP dir         : {MAMMALIAN_MIP_DIR}")
    ms = list_mammalian_mips()
    print(f"  found {len(ms)} MIPs:")
    for p in ms:
        print(f"    {p.name}")
    overlap = overlap_image_ids()
    print(f"  overlap (both z-stack and MIP available): {overlap}")
    print(f"  voxel size      : XY={MAMMALIAN_VOXEL_XY_UM} um, Z={MAMMALIAN_VOXEL_Z_UM} um")
    print()
    print("=== Yeast ===")
    print(f"z-stacks dir    : {YEAST_ZSTACK_DIR}")
    zs = list_yeast_zstacks_all()
    print(f"  found {len(zs)} z-stack TIFFs (full series)")
    subset = list_yeast_zstacks_subset()
    print(f"  manuscript subset (matches MIPs): {len(subset)} images")
    print(f"MIP dir         : {YEAST_MIP_DIR}")
    mips = list_yeast_mips()
    print(f"  found {len(mips)} MIPs")
    print(f"  voxel size      : XY={YEAST_VOXEL_XY_UM} um, Z={YEAST_VOXEL_Z_UM} um")


if __name__ == "__main__":
    _hardware_report()
    _path_report()
