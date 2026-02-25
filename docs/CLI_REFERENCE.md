# CLI Reference

Complete reference for all `cellquant.py` command-line arguments.

## Basic usage

```bash
python cellquant.py [images_dir] [channel_specs...] --out [output_dir] [options]
```

Or with named arguments:

```bash
python cellquant.py --images [dir] --channels [specs...] --out [output_dir] [options]
```

Both forms are equivalent.

## Required arguments

### Images directory

The path to a folder containing multi-channel TIFF files.

```bash
python cellquant.py /path/to/images/ ...
# or
python cellquant.py --images /path/to/images/ ...
```

The pipeline scans for `*.tif` and `*.tiff` files. Subdirectories are not searched.

### Channel definitions

Define what each channel is and what to do with it using the format `"position:Name:role"`.

```bash
python cellquant.py /path/ "1:DAPI:nucleus" "2:GFP:quantify" ...
# or
python cellquant.py --channels "1:DAPI:nucleus" "2:GFP:quantify" ...
```

**Roles:**

| Role | Purpose | Used for |
|------|---------|----------|
| `nucleus` | Nuclear segmentation | Cellpose nuclear detection, defining nuclear ROI |
| `quantify` | Puncta detection and intensity measurement | LoG puncta detection, intensity metrics |
| `nucleolus` | Nucleolar mask generation | Otsu thresholding within cells, proximity reference |
| `cell-boundary` | Cell segmentation input | Direct input to Cellpose cell detection |
| `skip` | Ignored | Not used in any analysis |

You can have at most one `nucleus` channel but multiple `quantify` and `nucleolus` channels.

### Output directory

```bash
--out /path/to/output/
```

Created automatically if it doesn't exist.

## Cell type and preset

```bash
--cell-type {mammalian, yeast, bacteria}
```

| Preset | Model | Downsample | Cell diameter | Notes |
|--------|-------|------------|---------------|-------|
| `mammalian` | cpsam | 3× | 120 px | Nuclear segmentation expected |
| `yeast` | cpsam | 1× (none) | 40 px | Area filtering 200–5000 px |
| `bacteria` | cpsam | 1× | 15 px | Specialized for small cells |

All preset values are overridable. For example, `--cell-type yeast --cell-diameter 50` uses the yeast preset but overrides the diameter.

## Filename parsing

```bash
--filename-pattern "MAX_{condition}_rep{replicate}"
# or equivalently:
--file-pattern "MAX_{condition}_rep{replicate}"
```

Curly-brace placeholders are extracted from filenames. Common patterns:

```bash
"MAX_{condition}_rep{replicate}"       # MAX_control_rep1.tif → condition=control, replicate=rep1
"{condition}_{replicate}"              # arsenite_rep3.tif → condition=arsenite, replicate=rep3
"img_{condition}_{field}_{replicate}"  # img_treated_f1_rep2.tif
```

Metadata extracted from filenames appears in all output CSVs and determines how data is grouped for plots and statistics.

### Condition mapping and ordering

```bash
--condition-map ctrl=control ars=arsenite
--condition-order control arsenite
```

`--condition-map` renames extracted condition strings. `--condition-order` sets the display order in plots (default: alphabetical/numeric).

## Segmentation parameters

### Cellpose

```bash
--pretrained-model cpsam        # Cellpose model (default: from preset)
--cell-diameter 40.0            # Expected cell diameter in pixels
--nuclei-diameter 30.0          # Expected nucleus diameter in pixels
--flow-threshold 0.4            # Cellpose flow threshold
--cellprob-threshold -1.0       # Cellpose cell probability threshold
--seg-downsample 3              # Downsample factor before Cellpose (1 = none)
--no-gpu                        # Force CPU mode
--cpu-threads 4                 # Number of CPU threads for Cellpose
```

### Nuclear segmentation

```bash
--nucleus-dilate-px 3           # Dilate nuclear masks by N pixels
```

### Cell filtering

```bash
--min-cell-area 200             # Remove cells smaller than N pixels
--max-cell-area 5000            # Remove cells larger than N pixels
--keep-min-nuclei 1             # Minimum nuclei per cell (0 = no filter)
--keep-max-nuclei 4             # Maximum nuclei per cell (0 = no filter)
```

### Cell segmentation channel

```bash
--cell-seg-channel DAPI         # Which channel to use for cell segmentation
```

By default, the pipeline uses the `nucleus` channel for cell segmentation. If no nucleus channel is defined (e.g., yeast with only `nucleolus`), it creates a composite from all non-skip channels.

## Puncta detection

```bash
--puncta-channels Sis1 Tif6     # Which channels to detect puncta in (by name)
--no-puncta                     # Suppress puncta detection entirely (intensity metrics only)
--log-sigma 1.5                 # Laplacian-of-Gaussian sigma
--puncta-min-area-px 3          # Minimum punctum area (pixels)
--puncta-max-area-px 500        # Maximum punctum area (pixels)
--puncta-threshold-method {otsu, fixed}  # Thresholding method for puncta detection
--puncta-threshold-fixed 500.0  # Fixed threshold value (if method is fixed)
--puncta-min-circularity 0.0    # Minimum punctum circularity (0–1)
--puncta-min-solidity 0.0       # Minimum punctum solidity (0–1)
--puncta-compartment {cytosol, whole-cell, nucleus}
```

**Default behavior:** If `--puncta-channels` is not specified, puncta are automatically detected in all `quantify` channels. Use `--no-puncta` to suppress puncta detection and compute only intensity metrics. Use `--puncta-channels` to override the default and detect puncta in specific channels only.

`--puncta-compartment cytosol` restricts puncta detection to the cytoplasmic region (cell minus nucleus). Requires a `nucleus` channel. Falls back to `whole-cell` if no nucleus is available.

## Colocalization

```bash
--colocalization                                      # Enable pairwise colocalization
--colocalization-compartment {whole-cell, cytosol, nucleus}  # Region for analysis
```

Computes Pearson's R and Manders' M1/M2 (with Costes automatic thresholding) for all pairs of `quantify` + `nucleolus` channels. Requires at least 2 eligible channels.

## Nucleolar proximity

```bash
--nucleolar-proximity Nsr1      # Channel name defining nucleolar mask
--proximity-threshold 5         # Distance threshold for "proximal" (pixels)
```

Measures the distance from each punctum centroid to the nearest nucleolar boundary. Reports per-cell: mean distance, min distance, fraction of puncta within the threshold.

Requires the named channel to have role `nucleolus`. Puncta are measured in all `--puncta-channels`.

## Plotting

```bash
--skip-plots                    # Skip all plot generation
--trend                         # Add trend line to multi-condition plots
```

The pipeline automatically selects plot type based on number of conditions:
- ≤2 conditions: Violin superplots with replicate-level Wilcoxon rank-sum p-value (shown when both conditions have ≥3 replicates)
- 3+ conditions: Jittered strip plots with per-image median diamonds (descriptive only, no statistical test)

## Masks and QC

```bash
--no-save-masks                 # Don't save segmentation masks
--qc-downsample 2               # Downsample factor for QC overlay images
--qc-dpi 150                    # DPI for QC overlay images
```

## Configuration file

```bash
--config config.yml             # Load parameters from a YAML file
```

Config file values override cell-type presets but are overridden by explicit CLI arguments. Priority order: CLI > config file > cell-type preset > defaults.

Example `config.yml`:
```yaml
cell_diameter: 50.0
log_sigma: 2.0
puncta_min_area_px: 5
puncta_compartment: cytosol
```

## Output files

| File | Contents |
|------|----------|
| `cells.csv` | Per-cell measurements (all metrics in one table) |
| `images.csv` | Per-image summaries |
| `colocalization.csv` | Per-cell pairwise colocalization (if `--colocalization`) |
| `nucleolar_proximity.csv` | Per-cell puncta-to-nucleolus distances (if `--nucleolar-proximity`) |
| `nucleolar_morphology.csv` | Per-cell nucleolar shape metrics (if any `nucleolus` channel) |
| `config_used.yml` | Complete parameter record |
| `qc/*.png` | QC overlay images |
| `masks/*.tif` | Segmentation masks (unless `--no-save-masks`) |
| `plots/*.png` | Superplot visualizations (unless `--skip-plots`) |
| `prism/*.csv` | Prism-ready data tables (unless `--skip-plots`) |
