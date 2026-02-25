# Tutorial 2: Yeast Temperature Series

**Time:** ~45 minutes (plus ~25 minutes pipeline runtime)
**Prerequisites:** Complete [Tutorial 1](TUTORIAL_1_mammalian_SGs.md) first.
**What you'll learn:** Nucleolar segmentation, colocalization, spatial proximity, nucleolar morphology, multi-condition plotting.

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
