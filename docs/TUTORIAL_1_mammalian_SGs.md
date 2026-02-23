# Tutorial 1: Mammalian Stress Granules

**Time:** ~30 minutes (plus ~5 minutes pipeline runtime)
**What you'll learn:** How to run `cellquant`, interpret QC overlays, and understand the output files.

> **Using the example data subset?** The repository ships 2 cropped images (1 control, 1 arsenite) for quick testing. This tutorial describes the full dataset (9 images, 4–5 replicates per condition), but all commands work identically on the subset — you'll just see fewer cells and conditions in the output, and no p-values will appear (the statistical test requires ≥3 replicates per condition).

## The biology

U2OS cells were treated with 500 μM sodium arsenite to induce stress granules — cytoplasmic condensates of RNA and RNA-binding proteins that form under stress. The images have three channels:

1. **DAPI** — stains DNA (nuclei)
2. **G3BP1** — stress granule marker (green)
3. **PABPC1** — poly(A)-binding protein, also recruited to stress granules (magenta)

We expect arsenite-treated cells to have bright puncta (stress granules) in both G3BP1 and PABPC1 channels that are absent in untreated controls.

## Step 0: Prerequisites

Make sure you've completed the [Installation Guide](INSTALL.md). You should be able to run:

```bash
conda activate cellquant
python cellquant.py --help
```

## Step 1: Organize your images

Your images should be in a single folder with a consistent naming pattern. The full example dataset looks like this:

```
mammalian_SGs/
├── MAX_control_rep1.tif
├── MAX_control_rep2.tif
├── MAX_control_rep3.tif
├── MAX_control_rep4.tif
├── MAX_arsenite_rep1.tif
├── MAX_arsenite_rep2.tif
├── MAX_arsenite_rep3.tif
├── MAX_arsenite_rep4.tif
└── MAX_arsenite_rep5.tif
```

(The repository subset includes only `MAX_control_rep1.tif` and `MAX_arsenite_rep1.tif`.)

Key points:
- All files are multi-channel TIFF maximum intensity projections (MIPs)
- The filename encodes the condition (`control` or `arsenite`) and replicate number
- The naming pattern is consistent: `MAX_{condition}_rep{replicate}.tif`

**Using your own images?** They need to be multi-channel TIFFs. If your images are in other formats (`.nd2`, `.czi`, `.lif`), use Fiji/ImageJ to export them as TIFFs first. If you have z-stacks, create maximum intensity projections in Fiji (`Image > Stacks > Z Project > Max Intensity`).

## Step 2: Know your channels

Before running the pipeline, you need to know:
- How many channels are in your images?
- What is each channel?
- What order are they in (channel 1, 2, 3...)?

For this dataset:
- Channel 1: DAPI (nuclei)
- Channel 2: G3BP1 (stress granule marker — we want to quantify this)
- Channel 3: PABPC1 (another marker — we also want to quantify this)

If you're not sure about your channel order, open one image in Fiji and check `Image > Color > Channels Tool`.

## Step 3: Construct the command

The command has four parts:

```bash
python cellquant.py [images_dir] [channel_definitions] [options] --out [output_dir]
```

For this dataset:

```bash
python cellquant.py example_data/mammalian_SGs/ \
  "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
  --cell-type mammalian \
  --out example_data/mammalian_SGs/output/ \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

Let's break this down:

| Part | Meaning |
|------|---------|
| `example_data/mammalian_SGs/` | Where your images are |
| `"1:DAPI:nucleus"` | Channel 1 is DAPI, use it for nuclear segmentation |
| `"2:G3BP1:quantify"` | Channel 2 is G3BP1, quantify puncta in it |
| `"3:PABPC1:quantify"` | Channel 3 is PABPC1, quantify puncta in it |
| `--cell-type mammalian` | Use mammalian cell defaults (large cells, 3x downsampling) |
| `--out ...` | Where to save results |
| `--filename-pattern ...` | How to extract condition and replicate from filenames |

**Don't understand a flag?** Ask your AI assistant: "What does `--cell-type mammalian` do in cellquant?" It will explain.

## Step 4: Run it

```bash
conda activate cellquant
python cellquant.py example_data/mammalian_SGs/ \
  "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
  --cell-type mammalian \
  --out example_data/mammalian_SGs/output/ \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

You'll see progress output:

```
Found 9 images in example_data/mammalian_SGs
GPU: False (or True if you have a CUDA GPU)
Channels: DAPI(nucleus), G3BP1(quantify), PABPC1(quantify)
Nuclear segmentation: yes
Cell seg channel: DAPI

=== [1/9] MAX_control_rep1.tif ===
  Cellpose: cells 14, nuclei 15 ...
  Kept 12 cells (nuclei filter: 1-4)
  Puncta: G3BP1 → 8, PABPC1 → 3
  ...
```

(If you're using the 2-image subset, you'll see "Found 2 images" instead.)

**How long will it take?** About 1–5 minutes per image on CPU, depending on your computer and image size. On a GPU, ~10–30 seconds per image. The cropped subset images are faster (~1 minute total).

## Step 5: Check the QC overlays

This is the most important step. Open the `output/qc/` folder and look at the QC overlay images.

```bash
open example_data/mammalian_SGs/output/qc/   # Mac
# or navigate there in your file browser
```

Each QC image shows:
- **Cyan outlines:** Cell boundaries (from Cellpose)
- **Yellow outlines:** Nuclear boundaries (from DAPI)
- **Magenta dots:** Detected puncta

**What to check:**
- Are cells correctly segmented? Each cell should have its own cyan outline.
- Are nuclei inside cells? Yellow outlines should be within cyan outlines.
- Are puncta real? Magenta marks should correspond to visible bright spots, not background noise.

**If something looks wrong:**
- Cells merging together → try `--cell-diameter 80` (smaller diameter)
- Cells being split → try `--cell-diameter 150` (larger diameter)
- Too many false puncta → try `--log-sigma 2.0` (stricter detection)
- Missing real puncta → try `--log-sigma 1.0` (more sensitive detection)

Ask your AI assistant: "The cells in my QC overlay are being merged together. Here's what the image looks like: [describe what you see]. How do I fix this?"

## Step 6: Understand the output files

```
output/
├── config_used.yml          # Exact parameters used (for reproducibility)
├── cells.csv                # Per-cell measurements (main data file)
├── images.csv               # Per-image summaries
├── qc/                      # QC overlay images
│   ├── MAX_control_rep1_qc.png
│   └── ...
├── masks/                   # Segmentation masks (for further analysis)
│   ├── MAX_control_rep1_cellmask.tif
│   └── ...
├── plots/                   # Superplot visualizations
│   ├── n_puncta_G3BP1_superplot.png
│   └── ...
└── prism/                   # Prism-ready CSV files
    ├── n_puncta_G3BP1_prism.csv
    └── ...
```

### cells.csv — your main data file

Each row is one cell. Key columns:

| Column | Meaning |
|--------|---------|
| `image` | Source image filename |
| `condition` | Extracted from filename (e.g., "control", "arsenite") |
| `replicate` | Extracted from filename (e.g., "rep1") |
| `cell_id` | Unique cell identifier within each image |
| `cell_area_px` | Cell area in pixels |
| `n_puncta_G3BP1` | Number of G3BP1 puncta in this cell |
| `n_puncta_PABPC1` | Number of PABPC1 puncta in this cell |
| `frac_condensed_G3BP1` | Fraction of G3BP1 signal in puncta |
| `frac_condensed_PABPC1` | Fraction of PABPC1 signal in puncta |

This CSV can be opened in Excel, imported into GraphPad Prism, or loaded in R/Python for further analysis.

## Step 7: Interpret the superplots

Open the `output/plots/` folder. You'll see superplots for each metric.

The superplot shows:
- **Small colored dots:** Individual cells (each dot = one cell)
- **Large black diamonds:** Replicate medians (each diamond = the median of all cells from one image)
- **Horizontal line:** Overall median per condition

**About the p-values:** When both conditions have ≥3 biological replicates, the pipeline computes a Wilcoxon rank-sum test on the replicate medians (the black diamonds), NOT on the individual cells. This is statistically correct — cells from the same image are not independent observations. The p-value is shown as a bracket above the two conditions.

If either condition has fewer than 3 replicates, no p-value is shown. The data are presented descriptively, and that's the honest thing to do.

With 4–5 biological replicates, you may see p-values of 0.05–0.15 even for visually obvious effects. This is expected. It does not mean the effect isn't real — it means you need more biological replicates to reach conventional statistical significance. The per-cell data clearly show the effect; the replicate-level test honestly reports that the sample size is small.

## Step 8: Next steps

**Using your own data:**
1. Organize your images in a folder with consistent naming
2. Figure out your channel order and assign roles
3. Ask your AI assistant to help construct the command
4. Run it, check QC overlays, iterate on parameters

**Ready for more?** Continue to [Tutorial 2: Yeast Temperature Series](TUTORIAL_2_yeast_temp.md) to learn about colocalization, nucleolar proximity, and multi-condition analysis.

## Common questions

**Q: My images have 4 channels but I only care about 3. What do I do?**
Mark the unwanted channel as `skip`: `"4:BF:skip"`

**Q: My images are not z-projections. Can I use them?**
The pipeline expects 2D images (single-plane or MIPs). If you have z-stacks, create MIPs in Fiji first.

**Q: I get a "no images found" error.**
Check that your image directory contains `.tif` files and the path is correct. The pipeline looks for files matching `*.tif` and `*.tiff`.

**Q: The puncta detection is finding spots in the DAPI channel. Why?**
The `nucleus` role channel is not searched for puncta. Only `quantify` channels are. If you're seeing this, double-check your channel role assignments.

**Q: Can I run this on a cluster?**
Yes. `cellquant.py` works on any system with the right Python environment. On an HPC cluster with SLURM, you can submit it as a job. Ask your AI assistant for a SLURM submission script.
