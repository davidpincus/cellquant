# Tutorial 2: Yeast Temperature Series

**Time:** ~45 minutes (plus ~25 minutes pipeline runtime)
**Prerequisites:** Complete [Tutorial 1](TUTORIAL_1_mammalian_SGs.md) first.
**What you'll learn:** Nucleolar segmentation, colocalization, spatial proximity, nucleolar morphology, multi-condition plotting.

> **A note before you start.** cellquant is not a substitute for understanding what the analysis is doing. This tutorial uses colocalization (Pearson, Manders) and 3D-derived nucleolar morphology — both of which have failure modes that look like normal output. Before publishing a number from this pipeline, read the relevant sections of [CONCEPTS.md](CONCEPTS.md). In particular: if you run `--colocalization` on a MIP, the resulting Pearson and Manders values are statistically unreliable, and **cellquant refuses the run by default**. The commands below pass `--allow-2d-colocalization` to force it so the tutorial can show you the output — and every row of `colocalization.csv` is then stamped `projection_derived=True`. [Step 8](#step-8-the-same-experiment-in-3d) repeats the analysis on z-stacks, where no such flag is needed.

> **Using the example data subset?** The repository ships 2 cropped images (25°C and 40°C only) for quick testing. This tutorial describes the full dataset (6 temperatures, 1 replicate each), but all commands work identically on the subset — you'll just see fewer conditions in the output and the trend lines won't be as informative with only 2 points. The cropped images also have fewer cells (~20–50 vs ~200–600 per image) and run much faster (~2 minutes total vs ~25 minutes).

## The biology

Budding yeast expressing three fluorescent markers are imaged across a temperature gradient from 25°C (permissive) to 40°C (severe heat stress):

1. **Tif6-mCherry** (channel 1) — a late ribosome biogenesis factor, normally nucleolar
2. **Nsr1-BFP** (channel 2) — a dense fibrillar component (DFC) nucleolar marker, defines nucleolar position
3. **Sis1-GFP** (channel 3) — an Hsp40 chaperone that forms cytoplasmic condensates under heat stress

> **Note:** Both Tif6 and Sis1 change localization at high temperatures, but for different biological reasons. Tif6 is a ribosome biogenesis factor that normally localizes to the nucleolus; under heat stress, ribosome biogenesis shuts down and Tif6 redistributes to the cytoplasm. Sis1 is a cytosolic Hsp40 chaperone that forms condensates (visible puncta) under stress. If you see nuclear-to-cytoplasmic redistribution in channel 1, that is the expected Tif6 behavior, not a channel labeling error.

We expect to see:
- Sis1 condensates appearing at ≥36°C
- Nucleolar shape changing from crescent (growing) to round (arrested)
- Tif6 potentially redistributing away from the nucleolus at extreme temperatures (ribosome biogenesis shutdown)

## What's different from Tutorial 1?

| Feature | Tutorial 1 (Mammalian) | Tutorial 2 (Yeast) |
|---------|----------------------|---------------------|
| Cell size | ~20–50 μm | ~5 μm |
| Nuclear stain | DAPI (dedicated channel) | None available |
| Nucleolus | Not analyzed | Nsr1 marks it |
| Conditions | 2 (± arsenite) | 6 (temperature series) |
| Replicates | 4–5 per condition | 1 per temperature |
| Segmentation | 3x downsampled | Full resolution |
| Extra analyses | None | Colocalization, proximity, morphology |

## Step 1: Understand the channel roles

This dataset introduces a new channel role: **`nucleolus`**.

```
"1:Tif6:quantify"    — measure Tif6 puncta and intensity
"2:Nsr1:nucleolus"   — use Nsr1 to define the nucleolar region
"3:Sis1:quantify"    — measure Sis1 puncta and intensity
```

The `nucleolus` role tells the pipeline:
- Do NOT use this channel for cell segmentation
- Instead, threshold it within each cell to create a nucleolar mask
- Make this mask available for proximity measurements and morphometrics

Since there's no DAPI or nuclear stain, the pipeline will segment cells from a **composite image** (sum of all channels) rather than from a dedicated nuclear channel. This is handled automatically by the yeast preset.

## Step 2: Construct the command

```bash
python cellquant.py example_data/yeast_temperature/ \
  "1:Tif6:quantify" "2:Nsr1:nucleolus" "3:Sis1:quantify" \
  --cell-type yeast \
  --out example_data/yeast_temperature/output/ \
  --colocalization \
  --allow-2d-colocalization \
  --nucleolar-proximity Nsr1 \
  --puncta-channels Sis1 Tif6 \
  --trend \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

New flags compared to Tutorial 1:

| Flag | What it does |
|------|-------------|
| `--cell-type yeast` | Small cells, no downsampling, area filtering |
| `--colocalization` | Compute Pearson's R and Manders' M1/M2 for all channel pairs |
| `--nucleolar-proximity Nsr1` | Measure distance from puncta to nucleolar boundary (defined by Nsr1) |
| `--puncta-channels Sis1 Tif6` | Which channels to detect puncta in (by name) |
| `--trend` | Add trend lines to multi-condition plots |

## Step 3: Run it

```bash
conda activate cellquant
python cellquant.py example_data/yeast_temperature/ \
  "1:Tif6:quantify" "2:Nsr1:nucleolus" "3:Sis1:quantify" \
  --cell-type yeast \
  --out example_data/yeast_temperature/output/ \
  --colocalization \
  --allow-2d-colocalization \
  --nucleolar-proximity Nsr1 \
  --puncta-channels Sis1 Tif6 \
  --trend \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

With the full dataset, this takes ~25 minutes on CPU because yeast images are processed at full resolution (no downsampling), there are more cells per image (~200–600 vs ~10–15 for mammalian), and the additional analyses (colocalization, proximity, morphology) add computation. The 2-image subset finishes in ~2 minutes.

You'll see a warning on Apple Silicon Macs:
```
[warn] MPS GPU not supported by cpsam Transformer; using CPU
```
This is expected and normal.

## Step 4: Check QC overlays

The yeast QC overlays look different from mammalian:

- **Cyan outlines:** Cell boundaries
- **White outlines:** Nucleolar boundaries (from Nsr1 thresholding)
- **Red dots:** Puncta close to the nucleolus (≤5 pixels)
- **Blue dots:** Puncta far from the nucleolus (>5 pixels)
- No yellow nuclear outlines (because there's no nucleus channel)

**What to check:**
- Are individual yeast cells separated, including mother-daughter pairs?
- Are the white nucleolar outlines capturing the bright Nsr1 signal inside cells?
- At 36–40°C, do you see red/blue puncta dots appearing (Sis1 condensates)?

**If cell segmentation is wrong:**
- Cells merging: try `--cell-diameter 30`
- Cells being split: try `--cell-diameter 50`
- Too much debris: try `--min-cell-area 300`
- Merged clusters passing filter: try `--max-cell-area 3000`

## Step 5: Understand the additional output files

Beyond the standard `cells.csv` and `images.csv`, you now have:

### colocalization.csv

One row per cell per channel pair. Three pairs: Tif6 vs Sis1, Tif6 vs Nsr1, Sis1 vs Nsr1.

| Column | Meaning |
|--------|---------|
| `pair` | Which two channels (e.g., "Tif6_vs_Nsr1") |
| `pearson_r` | Pearson's correlation coefficient (-1 to 1) |
| `manders_m1` | Fraction of channel A signal overlapping channel B |
| `manders_m2` | Fraction of channel B signal overlapping channel A |

High Nsr1-Tif6 colocalization at 25°C confirms both proteins are nucleolar. A decrease at 40°C suggests Tif6 is leaving the nucleolus.

### nucleolar_proximity.csv

One row per cell per puncta channel (Sis1 and Tif6).

| Column | Meaning |
|--------|---------|
| `channel` | Which puncta channel |
| `n_puncta` | Number of puncta in this cell |
| `mean_distance` | Average distance from puncta to nearest nucleolar boundary (pixels) |
| `min_distance` | Distance of closest punctum to nucleolus (pixels) |
| `fraction_proximal` | Fraction of puncta within 5 pixels of nucleolus |

### nucleolar_morphology.csv

One row per cell, describing the largest nucleolus in that cell.

| Column | Meaning |
|--------|---------|
| `nucleolar_area` | Area of the largest nucleolus (pixels) |
| `nucleolar_solidity` | Area / convex hull area. Low = crescent, high = round. |
| `nucleolar_circularity` | 4π × area / perimeter². 1.0 = perfect circle. |
| `nucleolar_eccentricity` | Ellipse fit. 0 = circle, →1 = elongated. |
| `n_nucleoli` | Number of nucleolar objects detected in this cell |

**The key metric for nucleolar shape:** Solidity captures the crescent-to-round transition. Actively growing cells have crescent-shaped nucleoli (solidity ~0.6–0.7), while heat-stressed cells have round nucleoli (solidity ~0.9+).

## Step 6: Interpret multi-condition plots

With 3 or more conditions, the pipeline generates strip plots instead of violins:
- Each dot is one cell
- Large diamonds are per-image medians (one per temperature, since n=1 replicate)
- Dashed line shows the trend through condition medians (because you used `--trend`)

**No p-values are shown.** With one replicate per condition, there is no meaningful statistical test. The data are presented descriptively. This is honest — the temperature series shows trends, not statistically confirmed differences. (If you're running the 2-image subset, you'll see the same: strip plots with 2 points and a trend line connecting them.)

**What to look for in the plots:**
- Sis1 puncta count increasing at ≥36°C
- Nucleolar circularity increasing with temperature
- Nucleolar solidity increasing with temperature
- Nsr1-Tif6 Pearson's R decreasing at 40°C
- Cell area increasing at high temperatures (growth arrest → bigger cells)

## Step 7: Use the data in other software

### GraphPad Prism
The `prism/` folder contains pre-formatted CSVs that can be directly imported into Prism for custom plotting.

### R / Python
Load `cells.csv` directly:
```r
# R
library(readr)
cells <- read_csv("output/cells.csv")
```
```python
# Python
import pandas as pd
cells = pd.read_csv("output/cells.csv")
```

All colocalization, proximity, and morphology metrics are pivoted into `cells.csv` as additional columns, so you can work from a single file.

## Step 8: the same experiment in 3D

Everything above was measured on maximum-intensity projections, and two things about that should have bothered you. You had to pass `--allow-2d-colocalization` to get Pearson and Manders at all. And the nucleolar shape metrics describe the outline of a projected shadow, not the shape of the nucleolus.

`example_data/yeast_3d/` is the same yeast temperature experiment as `example_data/yeast_temperature/` — same strain, same three markers, same question — but as z-stacks instead of projections. Two OME-TIFFs, `25C_rep1.tif` and `40C_rep1.tif`, each 41 z-planes × 3 channels × 256 × 256 pixels (~8 MB). At 0.1057 μm laterally and 0.23 μm axially, that is a 27 × 27 μm field, 9.4 μm deep.

### Run it

```bash
python cellquant.py example_data/yeast_3d/ \
  "1:Tif6:quantify" "2:Nsr1:nucleolus" "3:Sis1:quantify" \
  --cell-type yeast \
  --out results_3d/ \
  --seg-downsample 2 \
  --puncta-channels Tif6 Sis1 \
  --colocalization \
  --nucleolar-proximity Nsr1 \
  --condensate-index \
  --filename-pattern "{condition}_rep{replicate}"
```

The interesting part of this command is what it does **not** contain.

**No `--mode 3d`.** cellquant looks at the shape of the first file and decides:
```
Mode: 3d (auto-detected from 25C_rep1.tif)
```

**No `--voxel-size`.** These files carry their physical pixel size in their OME metadata, so cellquant reads it from the file:
```
Voxel size: XY=0.1057 µm, Z=0.2300 µm (anisotropy=2.18; from file metadata (first file))
```
This is the entire argument for exporting OME-TIFF from your acquisition software rather than a plain TIFF. In 3D a voxel size is not optional — every volume, every distance in microns, and the anisotropy correction on the puncta filter depend on it — so if cellquant can find neither metadata nor `--voxel-size XY_UM Z_UM`, it aborts instead of guessing. (If you do type the values yourself, pass **both**, XY first and Z second.)

**No `--allow-2d-colocalization`.** This is the point of the gate you hit back in Step 2. On a projection, everything along the light path is collapsed into one pixel, so proteins that never share a plane still overlap perfectly and Pearson's R is inflated by geometry rather than biology. In 3D there is nothing to collapse: colocalization is computed over the voxel distribution, so it is allowed by default and `colocalization.csv` carries no `projection_derived` stamp.

Two smaller changes: `--filename-pattern` drops the `MAX_` prefix because these files are stacks, not projections; and `--seg-downsample 2` halves the XY resolution for segmentation (Z is never decimated). The yeast preset segments at full resolution, which in 3D means running Cellpose once per z-plane — 41 times per image. Downsampling keeps the demo short.

**Runtime: about 65 seconds per image**, 129 seconds for the pair, on an Apple Silicon laptop with no GPU. You will still see the `[warn] MPS GPU not supported by cpsam Transformer; using CPU` line from Step 3.

What you should *not* see is a `[warn] cell shape` line. cellquant checks whether your segmented yeast cells look physically plausible for the voxel size you gave it; cells that come out flattened along Z are the signature of an under-scaled Z spacing. Silence here means the metadata was right.

### What changes in the output

The run writes the same file set as before — `cells.csv`, `images.csv`, `colocalization.csv`, `nucleolar_proximity.csv`, `nucleolar_morphology.csv`, `config_used.yml`, `provenance.json`, plus `masks/`, `qc/`, `plots/` (24 superplots) and `prism/`. The columns inside them are what differ.

**Areas become volumes.** Every 2D size column has a 3D counterpart reported twice: once as a raw voxel count, once converted to microns.

| 2D column | 3D columns |
|-----------|------------|
| `cell_area_px` | `cell_volume_vox`, `cell_volume_um3` |
| `nucleus_area_px` | `nucleus_volume_vox`, `nucleus_volume_um3` |
| `cytosol_area_px` | `cytosol_volume_vox`, `cytosol_volume_um3` |
| `Sis1_puncta_area_px` | `Sis1_puncta_volume_vox`, `Sis1_puncta_volume_um3` |
| `nucleolar_area` | `nucleolar_volume_vox`, `nucleolar_volume_um3` |

The `_vox` column is the raw count and cannot be wrong for units reasons. The `_um3` column is that count multiplied by the voxel volume, so it inherits whatever voxel size the run resolved — which is why the banner above is worth reading every time.

**Nucleolar shape metrics go away, and are not replaced.** `nucleolar_solidity`, `nucleolar_circularity` and `nucleolar_eccentricity` do not appear in a 3D run. There is no 3D counterpart to them. Solidity is what Step 5 called "the key metric for nucleolar shape", so this is a real loss, and it is deliberate: the 3D shape descriptors you would reach for instead (sphericity, surface area) depend on surface reconstruction, which is not reliable at 0.23 μm axial sampling with a PSF several times that. Reporting them would produce numbers that look precise and are not. What you get instead is honest size: `nucleolar_volume_um3` and `nucleolar_eq_diameter_um`, the diameter of a sphere with the same volume.

**Proximity switches to microns.** The 2D cutoff is `--proximity-threshold 5` (pixels); in 3D it is `--proximity-threshold-um 0.5` and distances are measured with an anisotropy-aware transform, so a voxel step in Z counts for more than a voxel step in XY. `cells.csv` and `nucleolar_proximity.csv` also gain two sensitivity columns at ±0.1 μm around your cutoff — `Sis1_fraction_proximal_0.4um` and `Sis1_fraction_proximal_0.6um`. Check them. If your conclusion flips between 0.4 and 0.6 μm, it was never about proximity.

**Tuning flags change too.** The Step 4 segmentation fixes `--min-cell-area` and `--max-cell-area` are silently ignored in 3D; use `--min-cell-volume-vox` and `--max-cell-volume-vox` (the yeast preset sets 1500 and 200000). Likewise `--puncta-min-area-px` / `--puncta-max-area-px` become `--puncta-min-volume-vox` / `--puncta-max-volume-vox` (yeast: 4 and 3000), and the 2D roundness filters `--puncta-min-circularity` / `--puncta-min-solidity` do nothing at all in 3D.

### The numbers, and what they are worth

The run finds 45 cells at 25°C and 18 at 40°C — a sparser field, not a segmentation failure — for a `cells.csv` of 63 rows and 66 columns.

| Median per cell | 25°C | 40°C |
|-----------------|------|------|
| `Tif6_puncta_n` | 9 | 35 |
| `Sis1_puncta_n` | 11 | 44 |
| `cell_volume_um3` | 27.77 | 86.69 |
| `nucleolar_volume_um3` | 2.58 | 3.88 |

The directions all agree with the 2D analysis: more Sis1 condensates under heat stress, larger cells under growth arrest, a swollen nucleolus. That agreement is the useful result — the same biology survives the change in dimensionality.

**The magnitudes are not.** These crops are small, dense windows chosen so the tutorial runs in two minutes, and a 27 μm field is a handful of cells, not a field of view. Absolute values here differ from what the same command produces on full frames. Treat the table as a demonstration that the columns are populated and point the right way.

**And there are no statistics here at all.** One image per condition means one median per condition, and cellquant needs at least three images per condition on both sides of a comparison before it will run a test — so no `pvalues.csv` is written, exactly as in Step 6. With two conditions and no `--trend`, you also get violin plots rather than the strip plots you saw earlier. This section shows you the workflow. It does not show you a result.

## Adapting this for your own yeast images

**Different markers?** Change the channel definitions:
```bash
"1:Hsp104:quantify" "2:DAPI:nucleus" "3:Sec63:quantify"
```
If you have DAPI, use `nucleus` instead of `nucleolus` — the pipeline will segment nuclei normally.

**Different number of channels?** The pipeline handles any number of channels (2+). Just define each one.

**Want colocalization but not proximity?** Drop the `--nucleolar-proximity` flag.

**Longer temperature series or drug titration?** The multi-condition plotting handles any number of conditions. It sorts them automatically if they're numeric (e.g., "25deg", "30deg" → ordered numerically).

## Summary of what you've learned

In this tutorial, you:
- Used the `nucleolus` channel role to define a subnuclear ROI
- Ran colocalization analysis across all channel pairs
- Measured spatial proximity from puncta to nucleolar boundaries
- Quantified nucleolar shape changes with morphometrics
- Interpreted multi-condition strip plots without statistical testing
- Understood when honest data presentation means *not* computing p-values

You now have the tools to analyze virtually any multi-channel fluorescence microscopy dataset with punctate structures, spatial relationships, or organelle morphology changes.
