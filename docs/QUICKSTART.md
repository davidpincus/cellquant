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

The repository ships a small cropped subset of each dataset for quick testing — 2 images each: `mammalian_SGs/` at 800×800 pixels, `yeast_temperature/` at 400×400 pixels, and `yeast_3d/` at 256×256 pixels × 41 z-planes. Run this to verify everything works:

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

### Option 4: 3D z-stacks

You do not have to tell cellquant your images are stacks — it detects that from the first file and switches to the 3D pipeline. This command runs on the shipped example data, so you can paste it right now:

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

Budget about 65 seconds per image on a Mac CPU. Two things worth noticing:

**No `--allow-2d-colocalization`.** Compare with Option 3. In 3D, colocalization is measured over the actual voxel distribution, so there is nothing to force — the guardrail exists only because a projection collapses Z and makes unrelated proteins look like they overlap.

**No `--voxel-size` either — and that is the point.** These files are OME-TIFFs that carry their own voxel size (XY = 0.1057 µm, Z = 0.23 µm), so cellquant reads it. The startup banner tells you where the number came from:

```
Voxel size: XY=0.1057 µm, Z=0.2300 µm (anisotropy=2.18; from file metadata (first file))
```

and `provenance.json` records `resolved_from: metadata`. When your microscope wrote correct metadata, the tool needs no help from you.

**If your own files lack it, you must pass `--voxel-size XY_UM Z_UM` — both values, lateral first, axial second.** Z spacing is almost never equal to XY spacing, and every physical measurement in 3D depends on getting the ratio right: `cell_volume_um3`, nucleolar volume, and every proximity distance in microns. cellquant will not assume 1 µm cubic voxels on your behalf; it aborts instead (see below).

3D output columns are volumes rather than areas — `cell_volume_um3` where a 2D run gives `cell_area_px`.

## After running

1. **Check QC overlays:** `open /path/to/output/qc/` — do segmentation boundaries look right?
2. **Get your data:** `open /path/to/output/cells.csv` — one row per cell, all metrics
3. **See the plots:** `open /path/to/output/plots/` — superplots with statistics

## Something wrong?

- **Cells merging:** add `--cell-diameter 30` (smaller)
- **Cells splitting:** add `--cell-diameter 150` (bigger)
- **Too many false puncta:** add `--log-sigma 1.0` (smaller = more selective)
- **Apple Silicon warning:** normal, ignore it
- **`[error] 3D mode but no voxel size is available`:** your stacks carry no voxel metadata, and cellquant refuses to invent one. Three ways forward: pass `--voxel-size XY_UM Z_UM` (lateral first, axial second — look the numbers up in your acquisition software); or `--project-z max` to flatten each stack and run the 2D pipeline instead; or `--assume-isotropic`, which proceeds with 1 µm cubic voxels and puts every volume and distance in **voxel units, not microns**
- **Running low on disk:** cellquant writes segmentation masks for every image by default — cell, nuclear, nucleolar, and one per puncta channel. They are zlib-compressed now (roughly 85× smaller than uncompressed), but a large 3D dataset still adds up. Add `--no-save-masks` if you do not need them — but then there is nothing for `--reuse-masks` to load, so a re-run has to re-segment from scratch

See [Troubleshooting](TROUBLESHOOTING.md) or ask your AI assistant.
