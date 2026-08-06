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

## 3D mode and voxel size

cellquant auto-detects 2D vs 3D from input shape. 2D MIPs run the paper-validated pipeline; 3D z-stacks run a full volumetric pipeline with Cellpose-SAM 3D segmentation, anisotropic LoG puncta detection, and 3D-native Pearson/Manders colocalization.

```bash
--mode {auto, 2d, 3d}                 # Force pipeline mode (default: auto-detect)
--voxel-size XY_UM Z_UM                # Voxel size in microns; read from OME/ImageJ metadata when absent
--assume-isotropic                     # Proceed with 1.0 µm isotropic voxels when no voxel size is available
--override-metadata                    # Let --voxel-size win when it disagrees with the file's metadata
--axes ZCYX                            # Override axis order for unusual TIFF layouts
--seg-3d-method {stitch, full}         # 3D segmentation strategy
--stitch-threshold 0.4                 # IoU threshold for Z-stitching
--puncta-min-volume-vox 8              # 3D puncta volume floor (voxels)
--puncta-max-volume-vox 3000           # 3D puncta volume ceiling (voxels)
--min-cell-volume-vox 0                # 3D cell volume filter (voxels); 0 disables
--max-cell-volume-vox 0
--proximity-threshold-um 0.5           # 3D proximity threshold (microns; replaces --proximity-threshold)
--project-z {max, sum, mean}           # Project 3D inputs to 2D before analysis (inside cellquant)
--z-crop-center 0                      # Analyse only N central Z slices (0 = disabled)
```

**Voxel size resolution order:** `--voxel-size XY Z` > OME-TIFF / ImageJ metadata. If neither supplies a value and you are running in 3D, cellquant **aborts**. All 3D metrics (`cell_volume_um3`, `nucleolar_volume_um3`, anisotropic LoG, anisotropic distances) are built on the voxel size, and a plausible-looking wrong number is worse than no number. Note that cellquant only accepts metadata that carries *both* a lateral and an axial size — an OME-TIFF with `PhysicalSizeX` but no `PhysicalSizeZ` counts as no metadata.

**`--assume-isotropic`** is the escape hatch for that abort: the run proceeds with 1.0 µm isotropic voxels and warns that your 3D sizes and distances are then in voxel units, not microns. Use it to get a run out of unlabelled data, not to publish from it.

**`--override-metadata`.** When you pass `--voxel-size` *and* the file carries voxel metadata, the two are compared axis by axis and a disagreement of more than 1 % aborts the run — that is nearly always a typo or the wrong folder. `--override-metadata` forces your value through, and the override is recorded in `provenance.json`.

**`--seg-3d-method stitch`** runs Cellpose 2D per-Z and stitches by IoU; a higher `--stitch-threshold` links less aggressively across Z (fewer over-extended cells). **`--seg-3d-method full`** runs Cellpose with `do_3D=True` and anisotropy-aware processing — more conservative but ~2.3× slower, and on our yeast calibration it did not measurably improve cell roundness (the residual axial elongation is optical PSF, not a segmentation artifact). All cell-type presets therefore default to `stitch`: yeast uses a strict `--stitch-threshold 0.65`; mammalian and bacteria use 0.4.

**`--project-z`** is for users who have 3D acquisition but want 2D analysis (e.g. legacy pipelines, lower-spec hardware). When set on 3D input, cellquant projects each stack with the chosen reduction (max-intensity, sum, or mean) and runs the 2D pipeline. No external Fiji/napari step needed.

## 2D pixel size

Two dimensions need no Z spacing, but they still need a **lateral** pixel size, because some parameters are expressed in microns. cellquant reads it from OME/ImageJ metadata; `--voxel-size XY_UM` supplies it directly (in 2D the single lateral value is enough).

Exactly two things depend on it:

- `--puncta-compartment-erode-um`
- any `~UM` pad inside a `--compartment` definition

If neither metadata nor `--voxel-size` supplies a pixel size, what happens depends on **who asked for the micron value**:

| Source of the value | Behaviour |
|---|---|
| You typed it on the command line | The run **aborts**. It is an assertion about physical units that cannot be honoured. |
| It came from a cell-type preset | The run **warns loudly and skips it** (treated as 0), so a first run on unlabelled data is not blocked. |

The asymmetry matters because a micron value silently reinterpreted as pixels is not merely imprecise — anything below 1.0 becomes a no-op, since the smallest nonzero distance on a pixel grid is one pixel. The shipped 2D example data carries no resolution metadata, so the mammalian preset's 0.5 µm erosion is skipped there and says so.

## Input formats

```bash
# Always supported (via tifffile, paper-validated):
.tif, .tiff (single- or multi-channel, 2D or 3D)
OME-TIFF (metadata is read for voxel size, channel names)

# Supported when bioio is installed:
.nd2 (Nikon)
.czi (Zeiss)
.lif (Leica)
.lsm (Zeiss legacy)
```

Install the non-tiff support:

```bash
pip install bioio bioio-nd2 bioio-czi bioio-lif bioio-ome-tiff
```

cellquant prefers `bioio` (the actively-maintained modular successor to aicsimageio); legacy `aicsimageio` is also accepted if already installed. cellquant raises a clear error if you pass a non-tiff file without a backend installed.

## Filename parsing

```bash
--filename-pattern "MAX_{condition}_rep{replicate}"
```

`--filename-pattern` and `--file-pattern` are **not** the same flag. `--file-pattern` is a discovery glob that replaces the default extension sweep (`--file-pattern "*.tif"`); it does not understand `{condition}` placeholders, and passing one there matches no files and aborts the run.

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

```bash
--reference-condition control
```

`--condition-map` renames extracted condition strings. `--condition-order` sets the display order in plots (default: alphabetical/numeric). `--reference-condition` names the condition every other condition is tested against; without it, comparisons are anchored on the first condition in the order.

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
--puncta-compartment REGION     # a built-in or a name defined with --compartment
--puncta-compartment-erode-um 0 # shrink that region by N microns before detecting
```

**Default behavior:** If `--puncta-channels` is not specified, puncta are automatically detected in all `quantify` channels. Use `--no-puncta` to suppress puncta detection and compute only intensity metrics. Use `--puncta-channels` to override the default and detect puncta in specific channels only.

`--puncta-compartment cytosol` restricts puncta detection to the cytoplasmic region (cell minus nucleus). Requires a `nucleus` channel. Falls back to `whole-cell` if no nucleus is available.

`--puncta-compartment-erode-um` shrinks that region inward by a given distance before detection, which keeps a bright membrane rim from inflating the threshold and suppressing real puncta. It is **in microns** (default 0; the mammalian preset sets 0.5), so it needs a known pixel size — in 3D the erosion is anisotropy-aware. If you type a value and no pixel size is available the run aborts; if the value came from a preset it warns and is skipped. See [2D pixel size](#2d-pixel-size) and [TROUBLESHOOTING](TROUBLESHOOTING.md).

Note that `{channel}_frac_intensity_in_puncta` is normalised against whichever region puncta were detected in, so changing `--puncta-compartment` changes both the numerator and the denominator.

### Defining your own regions

Built-in regions: `whole-cell`, `cell`, `nucleus`, `cytosol`, `nucleolus`. You can define more by set algebra, and use the name anywhere a region is accepted (`--puncta-compartment`, `--colocalization-compartment`, per-channel overrides).

```bash
--compartment "NAME = TERM [op TERM]..."   # repeatable
```

- `op` is `-` (minus), `&` (intersect) or `+` (union), applied **strictly left to right** — there is no operator precedence and no parentheses.
- **Operators must have spaces around them.** `cell - nucleus` is three tokens and parses; `cell-nucleus` is one token and is rejected, because otherwise it could not be told apart from the region named `whole-cell`.
- A `TERM` is a built-in, a name you defined **earlier** (forward references are refused, so a definition can never be circular), or `exclusion(CHANNEL,within=PARENT)`.
- Suffix a term with `~UM` to grow it or `~-UM` to shrink it by that many microns, measured with the anisotropic distance transform so the pad is physically uniform in 3D. In 2D this requires `--voxel-size`; without it the run aborts rather than silently applying the pad in pixels.

```bash
# cytoplasm excluding a generously-padded nucleus
--compartment "cytosol_wide = cell - nucleus~0.3" --puncta-compartment cytosol_wide

# everything in the cell that is not nucleolus, then trim the rim further
--compartment "peri = cell - nucleolus" \
--compartment "peri_tight = peri - nucleolus~0.3" --voxel-size 0.10571
```

Each defined region adds `{channel}_{NAME}_mean` and `{NAME}_area_px` (2D) or `{NAME}_volume_vox` / `{NAME}_volume_um3` (3D) to `cells.csv`. Every region reference is checked at startup, before segmentation runs, so a typo fails immediately instead of after Cellpose.

`exclusion(CHANNEL,within=PARENT)` is intended for regions defined by *absence* of signal — a vacuole seen as a hole in a cytosolic marker. The syntax parses, but the detector is **not yet enabled**: it has no defaults validated against a dedicated vacuole marker, and on projection data the signal it would key on is largely absent. Using it aborts with an explanation.

### Per-channel puncta tuning

When different channels in the same image need different puncta parameters (e.g. one diffuse-condensation channel with a low threshold and one tight-puncta channel with a high threshold), use these:

```bash
--puncta-rolling-ball G3BP1:25 PABPC1:0
# Per-channel rolling-ball (white tophat) background subtraction radius.
# Pixels in 2D, voxels in 3D. 0 disables for that channel.

--puncta-params-per-channel G3BP1:log_sigma=2.5,puncta_threshold_method=triangle \
                            PABPC1:log_sigma=1.8,puncta_min_area_px=4
# Per-channel overrides. Each spec is CH:KEY=VAL,KEY=VAL,...
# Supported keys:
#   log_sigma, puncta_threshold_method, puncta_threshold_fixed,
#   puncta_min_area_px, puncta_max_area_px,
#   puncta_min_volume_vox, puncta_max_volume_vox,
#   puncta_min_circularity, puncta_min_solidity,
#   puncta_compartment
```

Per-channel overrides supersede the global flag values (`--log-sigma`, `--puncta-threshold-method`, etc.) for the named channel only; other channels still see the global values.

### Fragmentation indices

cellquant also computes two threshold-handling regimes for the same underlying puncta biology — useful when discrete-puncta detection is borderline:

| Column | What it is |
|---|---|
| `{ch}_fragmentation_index_simple` | Connected components above LoG-Otsu within the cell |
| `{ch}_fragmentation_index_persistence` | Threshold-free integral over a swept threshold range |

```bash
--fragmentation-thresholds 10   # Sweep depth for the persistence variant (0 disables)
```

See [CONCEPTS.md](CONCEPTS.md) for when each is preferable to `puncta_n`.

## Colocalization

```bash
--colocalization                                      # Enable pairwise colocalization
--colocalization-compartment {whole-cell, cytosol, nucleus}  # Region for analysis
```

Computes Pearson's R and Manders' M1/M2 (with Costes automatic thresholding) for all pairs of `quantify` + `nucleolus` channels. Requires at least 2 eligible channels.

In 3D mode cellquant computes Pearson and Manders natively over the full voxel distribution. Pearson's R and Manders' M1/M2 are defined on the 3D voxel distribution; collapsing Z first changes apparent overlap, so **in 2D mode cellquant refuses the run and exits** rather than reporting a number it cannot stand behind. Either re-run on the source z-stacks (cellquant runs 3D colocalization natively), or pass `--allow-2d-colocalization` to force it — in which case every row of `colocalization.csv` is stamped `projection_derived=True`.

See [CONCEPTS.md — Colocalization](CONCEPTS.md#colocalization-opt-in-via---colocalization) for a longer discussion.

## Opt-in metrics

Some specialized metrics are off by default because they're irrelevant to most workflows. Enable them per-flag:

```bash
--condensate-index              # p95/mean intensity ratio per channel per cell
```

The Condensate Index quantifies internal contrast within the cell — useful for studying phase-separating proteins. It was always-on in previous versions of cellquant; as of the 2026 revision it's opt-in so the default output stays focused on general-purpose metrics. See [CONCEPTS.md — Condensate index](CONCEPTS.md#condensate-index-opt-in-via---condensate-index) for interpretation.

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
--reuse-masks                   # Load saved masks instead of re-running Cellpose
--qc-downsample 2               # Downsample factor for QC overlay images
--qc-dpi 150                    # DPI for QC overlay images
```

Masks are written zlib-compressed, which keeps a large 3D run from filling a disk.

`--reuse-masks` loads `masks/{stem}_cellmask.tif` and `{stem}_nucmask.tif` from the output directory instead of running Cellpose, which makes re-quantification with different puncta or colocalization settings very fast. Puncta and nucleolar masks are still recomputed, so changing those parameters takes effect. Any image with no saved mask is skipped with a warning — and `--no-save-masks` leaves nothing to reuse.

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

## GPU device

cellquant auto-detects the runtime device and prints what it chose at startup (look for `[device]` in the log). The actual decision tree:

- CUDA is available → use CUDA
- Apple Silicon MPS is available → fall back to CPU (Cellpose's cpsam Transformer ops are not yet supported on MPS)
- Otherwise → CPU

```bash
--no-gpu                        # Force CPU even if a GPU is available
```

For 3D workflows, a CUDA GPU is roughly an order of magnitude faster than CPU. See [INSTALL.md — GPU setup](INSTALL.md#step-6-gpu-setup-optional-recommended-for-3d) for driver and CUDA wheel installation per platform.

## Output files

| File | Contents |
|------|----------|
| `cells.csv` | Per-cell measurements (all metrics in one table) |
| `images.csv` | Per-image summaries |
| `colocalization.csv` | Per-cell pairwise colocalization (if `--colocalization`) |
| `nucleolar_proximity.csv` | Per-cell puncta-to-nucleolus distances (if `--nucleolar-proximity`) |
| `nucleolar_morphology.csv` | Per-cell nucleolar shape metrics (if any `nucleolus` channel) |
| `config_used.yml` | Complete parameter record |
| `provenance.json` | cellquant version, exact CLI invocation, input file SHA-256 checksums, model checkpoint, dependency versions, OME channel names + voxel sizes (FAIR-aligned reproducibility record) |
| `qc/*.png` | QC overlay images (2D: single panel; 3D: MIP + middle-Z side-by-side) |
| `masks/*.tif` | Segmentation masks (unless `--no-save-masks`); 3D masks are multi-page TIFFs |
| `plots/*.png` | Superplot visualizations (unless `--skip-plots`) |
| `prism/*.csv` | Prism-ready data tables (unless `--skip-plots`) |
