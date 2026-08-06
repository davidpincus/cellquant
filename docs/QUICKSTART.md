# Quick Start

Get results in 5 minutes. For detailed explanations, see the [full tutorials](TUTORIAL_1_mammalian_SGs.md).

## Prerequisites

Pick one (see [INSTALL.md](INSTALL.md) for details):

```bash
# Recommended: install uv once, then no further setup is needed
# (uv reads cellquant.py's PEP 723 header on every run)
curl -LsSf https://astral.sh/uv/install.sh | sh      # Mac / Linux / WSL

# OR: traditional conda environment
conda activate cellquant
```

Throughout this guide, `uv run cellquant.py …` and `python cellquant.py …` are interchangeable. Use whichever matches your setup.

## Try the included example data

The repository ships a small subset of each dataset (2 images per dataset, cropped to 400×400 pixels) for quick testing. Run this to verify everything works:

```bash
uv run cellquant.py example_data/mammalian_SGs/ \
  "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
  --cell-type mammalian \
  --out example_data/mammalian_SGs/test_output/ \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

This should finish in under a minute. Check the results:

```bash
open example_data/mammalian_SGs/test_output/qc/      # Mac
# or navigate to the qc/ folder in your file browser
```

If you see images with cyan cell outlines and yellow nuclear outlines, it's working.

### The warning you will see — it's expected

```
[warn] --puncta-compartment-erode-um = 0.5 µm comes from the 'mammalian' preset,
       but this 2D input carries no pixel size, so the value cannot be
       interpreted in microns. It is being SKIPPED, not applied as pixels.
```

The example TIFFs carry no pixel-size metadata, so cellquant cannot convert 0.5 µm into pixels. Rather than guess — or silently treat "0.5 µm" as "0.5 pixels", which would be a no-op that looks like it worked — it disables the erosion and says so. On your own images, pass `--voxel-size XY_UM` (your lateral pixel size in microns) to enable it.

The distinction matters: a micron value that came from a **preset** warns and is skipped, as here. A micron value that **you typed** is an assertion about physical units cellquant cannot honour, so it aborts instead.

## Use it on your own data

### Option 1: Mammalian cells with nuclear stain

```bash
python cellquant.py /path/to/images/ \
  "1:DAPI:nucleus" "2:MarkerA:quantify" "3:MarkerB:quantify" \
  --cell-type mammalian \
  --out /path/to/output/ \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

### Option 2: Yeast cells

```bash
python cellquant.py /path/to/images/ \
  "1:ChannelA:quantify" "2:ChannelB:quantify" \
  --cell-type yeast \
  --out /path/to/output/ \
  --filename-pattern "{condition}_{replicate}"
```

### Option 3: Yeast with nucleolar analysis

```bash
python cellquant.py /path/to/images/ \
  "1:Protein1:quantify" "2:NucMarker:nucleolus" "3:Protein2:quantify" \
  --cell-type yeast \
  --out /path/to/output/ \
  --colocalization \
  --allow-2d-colocalization \
  --nucleolar-proximity NucMarker \
  --puncta-channels Protein1 Protein2 \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

## After running

1. **Check QC overlays:** `open /path/to/output/qc/` — do segmentation boundaries look right?
2. **Get your data:** `open /path/to/output/cells.csv` — one row per cell, all metrics
3. **See the plots:** `open /path/to/output/plots/` — superplots with statistics

## Something wrong?

- **Cells merging:** add `--cell-diameter 30` (smaller)
- **Cells splitting:** add `--cell-diameter 150` (bigger)
- **Too many false puncta:** add `--log-sigma 1.0` (smaller = more selective)
- **Apple Silicon warning:** normal, ignore it

See [Troubleshooting](TROUBLESHOOTING.md) or ask your AI assistant.
