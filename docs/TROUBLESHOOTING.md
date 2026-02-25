# Troubleshooting

Common issues and how to fix them. If your problem isn't here, copy the error message and paste it to your AI assistant — it will know what to do.

## Installation issues

### "conda: command not found"

You either don't have conda installed, or it's not in your PATH.

**Fix:** Close and reopen your terminal. If it still doesn't work, reinstall Miniforge (see [INSTALL.md](INSTALL.md)).

### "ModuleNotFoundError: No module named 'cellpose'"

The cellquant environment isn't activated.

**Fix:**
```bash
conda activate cellquant
```

### Cellpose installation is very slow

PyTorch is a large download (~2 GB). This is normal. Wait for it to finish.

### "ERROR: Could not find a version that satisfies the requirement..."

Your Python version may be too old.

**Fix:** Check your version:
```bash
python --version
```
You need 3.11 or higher. Create a new environment with the right version:
```bash
conda create -n cellquant python=3.11 -y
```

## Runtime errors

### "[warn] MPS GPU not supported by cpsam Transformer; using CPU"

**This is not an error.** This is an informational warning on Apple Silicon Macs. The Cellpose Transformer model doesn't support Apple's GPU backend, so the pipeline automatically uses CPU. Everything works correctly, just slower.

### "pretrained model ... not found, using default model"

Cellpose couldn't find the specified model file. The pipeline will fall back to the default model, which should work fine.

**If you want a specific model:** Check `~/.cellpose/models/` to see what's installed. You can download models from within Python:
```python
from cellpose import models
models.Cellpose(model_type="cyto3")  # Downloads if not present
```

### "No images found in [path]"

The pipeline looks for `.tif` and `.tiff` files in the specified directory.

**Fix:**
- Check the path is correct: `ls /your/path/`
- Make sure files end in `.tif` or `.tiff` (not `.TIF`, `.nd2`, `.czi`, etc.)
- The pipeline does not search subdirectories

### "Channel spec must be 'position:Name:role'"

Your channel definition format is wrong.

**Correct format:** `"1:DAPI:nucleus"` (position:Name:role, quotes required)
**Common mistakes:**
- Missing quotes: `1:DAPI:nucleus` → `"1:DAPI:nucleus"`
- Wrong separator: `"1-DAPI-nucleus"` → `"1:DAPI:nucleus"`
- Invalid role: `"1:DAPI:nuclear"` → `"1:DAPI:nucleus"`

### "--nucleolar-proximity requires at least one channel with role 'nucleolus'"

You asked for nucleolar proximity but no channel is assigned the `nucleolus` role.

**Fix:** Change your channel definition to include a nucleolus channel:
```bash
"2:Nsr1:nucleolus"   # instead of "2:Nsr1:quantify"
```

### "Pipeline seems stuck / very slow"

Cellpose CPU mode is slow. Expected times:
- Mammalian cells (with 3x downsampling): 1–3 minutes per image
- Yeast cells (full resolution): 3–5 minutes per image
- Bacteria (full resolution): 2–4 minutes per image

For faster processing, use a computer with an NVIDIA GPU and don't pass `--no-gpu`.

## Segmentation issues

### Cells are being merged together

The Cellpose diameter is too large — it thinks adjacent cells are one big cell.

**Fix:** Decrease the cell diameter:
```bash
--cell-diameter 30     # yeast (default: 40)
--cell-diameter 80     # mammalian (default: 120)
```

### Cells are being split in half

The Cellpose diameter is too small — it thinks one cell is multiple cells.

**Fix:** Increase the cell diameter:
```bash
--cell-diameter 50     # yeast
--cell-diameter 150    # mammalian
```

### Too many small debris objects detected as cells

Area filtering isn't aggressive enough.

**Fix:**
```bash
--min-cell-area 300    # increase minimum (default varies by preset)
```

### Mother-daughter yeast pairs are being merged

This is a common challenge. Try:
```bash
--cell-diameter 30     # smaller diameter helps separate pairs
--flow-threshold 0.6   # stricter flow threshold
```

### Some cells are not being segmented at all

**Symptom:** QC overlays show cells with no boundary drawn. In puncta channels, the puncta mask may appear to outline the missing cell's boundary — this happens because a neighboring cell's mask expanded into the undetected cell's space.

**Likely cause:** Cellpose's cell probability threshold is too conservative, so low-contrast cells are ignored.

**Fix:** Lower the cell probability threshold (default 0.0):
```bash
--cellprob-threshold -2.0    # more aggressive detection
--cellprob-threshold -4.0    # even more aggressive (may over-segment)
```

If cells are much smaller or larger than the preset expects, also set the diameter explicitly:
```bash
--cell-diameter 40           # adjust to match your cells in pixels
```

### Multi-nucleated cells are being removed

The nuclei filter is too strict.

**Fix:**
```bash
--keep-max-nuclei 6    # allow more nuclei per cell (default: 4 for mammalian)
```

For yeast without a nuclear channel, set both to 0:
```bash
--keep-min-nuclei 0 --keep-max-nuclei 0
```

## Puncta detection issues

### Too many false puncta (noise being detected)

**Fix:** Lower the LoG sigma (smaller sigma responds to smaller/sharper structures, making detection more selective):
```bash
--log-sigma 1.0        # default: 1.5
```

Or increase minimum punctum area:
```bash
--puncta-min-area-px 5  # default: 3
```

### Real puncta are being missed

**Fix:** Raise the LoG sigma (larger sigma responds to broader structures, making detection more permissive):
```bash
--log-sigma 2.0
```

Or decrease minimum punctum area:
```bash
--puncta-min-area-px 2
```

### Puncta detected in the nucleus when they should only be cytoplasmic

**Fix:**
```bash
--puncta-compartment cytosol   # restrict to cytoplasm only
```

Requires a `nucleus` channel.

## Output issues

### "nucleolar_proximity.csv is empty / not generated"

Most likely you forgot `--puncta-channels`:
```bash
--puncta-channels Sis1 Tif6    # specify which channels to measure proximity for
```

### Condition names are wrong in plots

Use `--condition-map` to rename:
```bash
--condition-map 25deg=25°C 30deg=30°C 36deg=36°C
```

### Conditions are in the wrong order

Use `--condition-order`:
```bash
--condition-order 25deg 28deg 30deg 32deg 36deg 40deg
```

### Superplots show only one condition

The filename pattern isn't extracting conditions correctly. Check:
```bash
--filename-pattern "MAX_{condition}_rep{replicate}"
```

Make sure the pattern matches your actual filenames. The `{condition}` placeholder must capture the part of the filename that differs between groups.

## How to ask an AI for help

When asking an AI assistant about a `cellquant` error:

1. **Include the full error message** — copy everything from "Traceback" to the end
2. **Include your command** — the exact command you ran
3. **Describe what you see** — "The cells are being merged" is more helpful than "it's not working"
4. **Mention your system** — Mac/Windows/Linux, Apple Silicon or Intel

Example prompt:

> I'm running cellquant.py on yeast images and getting this error:
> ```
> [paste error here]
> ```
> My command was:
> ```
> [paste command here]
> ```
> I'm on a Mac with Apple Silicon. The QC overlays show cells being merged together. How do I fix this?
