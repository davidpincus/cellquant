# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cellquant is a single-file Python pipeline (`cellquant.py`, ~4700 lines) for multi-channel
fluorescence microscopy image quantification, in **2D or natively in 3D**. It performs
cell/nuclear/nucleolar segmentation via Cellpose, puncta detection, per-cell intensity
quantification, colocalization, nucleolar morphometry and proximity, and generates QC
overlays, superplots, and Prism-ready CSVs.

The intended user is a biologist who does not write code and drives the tool through an AI
assistant. That shapes the design: one file, no install step, everything configured through
CLI flags, and errors that explain themselves in the user's language.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate cellquant
```

Key dependencies: numpy (>=1.24,<2.0), pandas, scipy, matplotlib, scikit-image (>=0.24),
tifffile, pyyaml, cellpose (>=4.0), opencv-python-headless (<4.10), python 3.11–3.12,
optional torch for GPU, optional bioio readers for ND2/CZI/LIF.

The `numpy<2.0` and `opencv-python-headless<4.10` pins are load-bearing — cellquant checks
numpy at startup (before the cellpose import) and exits with an error if violated.

## Running the Pipeline

```bash
# Positional shorthand (images_dir then channel specs) — the documented form
python cellquant.py example_data/mammalian_SGs/ \
    "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
    --cell-type mammalian --out results/ \
    --filename-pattern "MAX_{condition}_rep{replicate}"

# Equivalent flag form
python cellquant.py --images MIPs/ --out results/ \
    --channels "1:DAPI:nucleus" "2:G3BP1:quantify" --puncta-channels G3BP1

# Native 3D on z-stacks (mode is auto-detected from image shape; --mode forces it)
python cellquant.py example_data/yeast_3d/ \
    "1:Tif6:quantify" "2:Nsr1:nucleolus" "3:Sis1:quantify" \
    --cell-type yeast --out results_3d/ --colocalization

# Reuse existing masks (skip Cellpose, rerun quantification/plotting)
python cellquant.py --images MIPs/ --out results/ --reuse-masks --channels ...
```

## Running Tests

```bash
pytest                                                      # all tests (~50 s)
pytest tests/test_regression.py -q                          # regression suite only
pytest -x                                                   # stop on first failure
```

15 tests. They are **structural** regression tests: they run the full pipeline on the
cropped datasets in `example_data/` and assert the output files exist with the expected
columns and row counts. They deliberately do **not** assert exact numeric values, because
Cellpose is not deterministic across platforms. Session-scoped fixtures mean the pipeline
runs once per dataset.

`pytest.ini` pins `testpaths = tests` so a bare `pytest` cannot wander into `validation_3d/`.

## Architecture

**Single-file design**: everything lives in `cellquant.py`. No packages or modules. This is
intentional — the target user shares and debugs one file with an AI assistant.

**Config resolution is 4-layer** (each overrides the previous):
1. `DEFAULTS` dict (line ~276)
2. Cell-type preset (`CELL_TYPE_PRESETS`, line ~160: mammalian, yeast, bacteria)
3. YAML config file (`--config`)
4. CLI args (highest priority)

In 3D, `apply_3d_preset_overrides` layers mode-specific values on top of the resolved config.

**Channel system**: channels are `"position:Name:role"` strings. Roles: `nucleus`,
`quantify`, `cell-boundary`, `nucleolus`, `skip`. The nucleus channel drives nuclear
segmentation; quantify channels get intensity metrics; puncta detection runs on channels
named in `--puncta-channels`.

**Cell segmentation strategy**: by default a composite image (sum of all non-skip channels)
drives Cellpose cell segmentation. Overridden by `--cell-seg-channel` or a channel with role
`cell-boundary`.

**Compartments**: measurements are made inside a region, not the whole image. Built-ins are
`whole-cell`, `cell`, `nucleus`, `cytosol`, `nucleolus`. Users define their own by set
algebra via `--compartment "perinuc = cell - nucleolus"`; operators apply strictly left to
right and require surrounding spaces. A term can be dilated/eroded in microns
(`nucleolus~0.3`), which requires a known pixel size. See `parse_compartment_specs`,
`_resolve_compartment_expr`, `make_compartment_mask`.

**Guardrails that refuse a run** — these are deliberate, and changing them is a design
decision, not a bug fix:
- 3D with no voxel size in metadata and no `--voxel-size` → abort, rather than assume 1 µm
  cubes. Volumes, micron distances, and the anisotropy correction all depend on it.
- Colocalization on 2D projections → refused by default, because collapsing Z inflates
  apparent overlap. `--allow-2d-colocalization` forces it and stamps `projection_derived=True`.
- Micron-denominated compartment erosion on an image with no pixel size → abort.
- Implausible voxel sizes (`_voxel_size_is_implausible`) and oblate 3D cell shapes are
  warned about, since both are the signature of a wrong voxel size.

**Pipeline flow per image** (`main()`, line ~3872):
1. `find_images` / `load_tiff` → channel dict (TIFF/OME-TIFF natively, ND2/CZI/LIF via bioio)
2. Mode resolution (2D vs 3D) from image shape, plus voxel-size resolution
3. Cellpose segmentation → cell_mask, nuc_mask (or load from `--reuse-masks`)
4. `filter_cells_by_area`, `map_nuclei_to_cells`
5. Optional: `segment_nucleoli`, `compute_nucleolar_morphology`
6. `detect_puncta` per puncta channel (LoG filter + threshold + shape filters)
7. `per_cell_metrics` → DataFrame with intensity, puncta counts, condensate index
8. Optional: `compute_colocalization` (Pearson + Manders with Costes thresholds)
9. Optional: `compute_nucleolar_proximity`
10. `save_qc_png` (or `_save_qc_3d`) → multi-panel QC image
11. Save masks as compressed TIFFs

**Post-loop aggregation**: concatenates per-image DataFrames, pivots colocalization /
proximity / morphology into `cells.csv`, generates superplots with Mann-Whitney U p-values,
writes Prism CSVs and `provenance.json` (versions, resolved config, resolved voxel size).

## Key Functions

- `parse_channels` / `build_config` / `parse_args` — input parsing and config merging
- `segment_nuclei` / `segment_cells` — Cellpose wrappers with downsample/upsample
- `detect_puncta` — LoG filter → threshold (otsu/triangle/fixed) → shape filtering
- `per_cell_metrics` — mean intensity, p95, condensate index, puncta counts per cell
- `costes_threshold` / `compute_colocalization` — Costes auto-threshold + Pearson/Manders
- `compute_nucleolar_proximity` — vectorized single pass over puncta voxels
- `superplot_violin` / `write_prism_outputs` — visualization and export

## Repository Layout

- `cellquant.py` — the pipeline
- `docs/` — user-facing install guide, quickstart, tutorials, CLI reference, concepts,
  troubleshooting, philosophy. **User-facing docs are part of the product**; a behavior
  change that isn't reflected in `docs/CLI_REFERENCE.md` is incomplete.
- `example_data/` — cropped datasets that ship with the repo, plus `expected_output/`
  snapshots the tests compare against
- `tests/` — structural regression tests
- `benchmark/` — cellquant-vs-CellProfiler agreement harness
- `validation_3d/` — the 2D-vs-3D validation harness and the manuscript figure/table
  scripts. Tracked as the reproducibility record; run outputs and large data are gitignored.
  Includes exploratory diagnostics that document dead ends — that is intentional.
- `view_3d_segmentation.py` — standalone 3D mask viewer

Large raw-image directories and the LaTeX manuscript live alongside the repo but are
gitignored; this repository is public, so never commit raw data, working notes, or anything
naming unpublished collaborator experiments.

## Conventions

- Match the existing style: type hints, module-level functions, no classes for pipeline
  stages, docstrings that explain *why* a threshold or default was chosen.
- Error messages are user-facing product surface. Name the flag to set and the value seen.
- When adding a flag: wire it in `parse_args`, give it a `DEFAULTS` entry, document it in
  `docs/CLI_REFERENCE.md`, and consider whether any cell-type preset should override it.
- Keep `README.md`, `CITATION.cff`, `.zenodo.json` and `__version__` in sync on release.
