# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cellquant is a single-file Python pipeline (`cellquant.py`, ~2300 lines) for multi-channel fluorescence microscopy image quantification. It performs cell/nuclear segmentation via Cellpose, puncta detection, per-cell intensity quantification, colocalization analysis, nucleolar morphometry, and generates superplots + Prism-ready CSVs.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate cellquant
```

Key dependencies: numpy (>=1.24,<2.0), pandas, scipy, matplotlib, scikit-image (>=0.24), tifffile, pyyaml, cellpose (>=4.0), opencv-python-headless (<4.10), python 3.11–3.12, optional torch for GPU.

## Running the Pipeline

```bash
# Basic usage
python cellquant.py --images MIPs/ --out results/ \
    --channels "1:DAPI:nucleus" "2:G3BP1:quantify" --puncta-channels G3BP1

# Positional shorthand (images_dir then channel specs)
python cellquant.py MIPs/ "1:DAPI:nucleus" "2:G3BP1:quantify" --out results/

# With cell-type preset (mammalian, yeast, bacteria)
python cellquant.py --images MIPs/ --out results/ --cell-type yeast \
    --channels "1:Med15:quantify" --puncta-channels Med15

# Reuse existing masks (skip Cellpose, rerun quantification/plotting)
python cellquant.py --images MIPs/ --out results/ --reuse-masks --channels ...
```

## Running Tests

```bash
pytest
```

No test files exist yet.

## Architecture

**Single-file design**: Everything lives in `cellquant.py`. No packages or modules.

**Config resolution is 4-layer** (each overrides the previous):
1. `DEFAULTS` dict (hardcoded)
2. Cell-type preset (`CELL_TYPE_PRESETS`: mammalian, yeast, bacteria)
3. YAML config file (`--config`)
4. CLI args (highest priority)

**Channel system**: Channels are specified as `"position:Name:role"` strings. Roles: `nucleus`, `quantify`, `cell-boundary`, `nucleolus`, `skip`. The nucleus channel drives nuclear segmentation; quantify channels get intensity metrics; puncta detection runs on channels named in `--puncta-channels`.

**Cell segmentation strategy**: By default uses a composite image (sum of all non-skip channels) for Cellpose cell segmentation. Overridden by `--cell-seg-channel` or a channel with role `cell-boundary`.

**Pipeline flow per image** (`main()` at line 1834):
1. `load_tiff` → channel dict
2. Cellpose segmentation → cell_mask, nuc_mask (or load from `--reuse-masks`)
3. `filter_cells_by_area`, `map_nuclei_to_cells`
4. Optional: `segment_nucleoli`, `compute_nucleolar_morphology`
5. `detect_puncta` per puncta channel (LoG filter + threshold + shape filters)
6. `per_cell_metrics` → DataFrame with intensity, puncta counts, condensate index
7. Optional: `compute_colocalization` (Pearson + Manders with Costes thresholds)
8. Optional: `compute_nucleolar_proximity`
9. `save_qc_png` → multi-panel QC image
10. Save masks as TIFFs

**Post-loop aggregation**: Concatenates per-image DataFrames, pivots colocalization/proximity/morphology into cells.csv, generates superplots with Mann-Whitney U p-values, writes Prism CSVs.

**Data directories**: `Jenny_FRET/`, `Luke_Hsp70/`, `Luke_Med15/` contain experiment-specific TIFF images and results. These are data, not code.

## Key Functions

- `parse_channels` / `build_config` — input parsing and config merging
- `segment_nuclei` / `segment_cells` — Cellpose wrappers with downsample/upsample
- `detect_puncta` — LoG filter → threshold (otsu/triangle/fixed) → shape filtering
- `per_cell_metrics` — computes mean intensity, p95, condensate index, puncta counts per cell
- `costes_threshold` / `compute_colocalization` — Costes auto-threshold + Pearson/Manders
- `superplot_violin` / `write_prism_outputs` — visualization and export
