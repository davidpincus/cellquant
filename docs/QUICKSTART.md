# Quick Start

Get results in 5 minutes. For detailed explanations, see the [full tutorials](TUTORIAL_1_mammalian_SGs.md).

## Prerequisites

```bash
conda activate cellquant    # or set up env first: see INSTALL.md
```

## Option 1: Mammalian cells with nuclear stain

```bash
python cellquant.py /path/to/images/ \
  "1:DAPI:nucleus" "2:MarkerA:quantify" "3:MarkerB:quantify" \
  --cell-type mammalian \
  --out /path/to/output/ \
  --filename-pattern "MAX_{condition}_rep{replicate}"
```

## Option 2: Yeast cells

```bash
python cellquant.py /path/to/images/ \
  "1:ChannelA:quantify" "2:ChannelB:quantify" \
  --cell-type yeast \
  --out /path/to/output/ \
  --filename-pattern "{condition}_{replicate}"
```

## Option 3: Yeast with nucleolar analysis

```bash
python cellquant.py /path/to/images/ \
  "1:Protein1:quantify" "2:NucMarker:nucleolus" "3:Protein2:quantify" \
  --cell-type yeast \
  --out /path/to/output/ \
  --colocalization \
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
- **Too many false puncta:** add `--log-sigma 2.0`
- **Apple Silicon warning:** normal, ignore it

See [Troubleshooting](TROUBLESHOOTING.md) or ask your AI assistant.
