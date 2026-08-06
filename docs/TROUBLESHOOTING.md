# Troubleshooting

Common issues and how to fix them. If your problem isn't here, copy the error message and paste it to your AI assistant — it will know what to do.

## Installation issues

### "OMP: Error #15: Initializing libomp.dylib, but found libomp already initialized"

This is a common macOS issue where multiple copies of the OpenMP runtime are pulled in by different packages (e.g., NumPy via MKL, PyTorch, system llvm-openmp). The pipeline now sets `KMP_DUPLICATE_LIB_OK=TRUE` automatically, so **updating to the latest cellquant.py should fix this**.

If you still see the error (e.g., running an older copy of the script), set the environment variable yourself before running:
```bash
export KMP_DUPLICATE_LIB_OK=TRUE
python cellquant.py --help
```

For a permanent fix, add that `export` line to your `~/.zshrc` (Mac) or `~/.bashrc` (Linux) file.

**If you want to eliminate the duplicate entirely** (optional):
```bash
# In your active conda environment, replace MKL with a non-conflicting backend
conda install nomkl
```

### "numpy.core.multiarray failed to import"

This error appears during `from cellpose import models` and looks like a cellpose bug, but it's actually a numpy version issue. Cellpose imports OpenCV (cv2), which tries to load numpy C extensions that changed in numpy 2.0 — the import chain is cellpose → cv2 → numpy C extensions → crash.

**Fix:**
```bash
pip install "numpy>=1.24,<2.0" "opencv-python-headless<4.10"
```
> **zsh users (default Mac shell):** The quotes around the version specs are required.

The latest version of cellquant detects numpy 2.x at startup and prints this fix automatically. If you're seeing this error, update your copy of `cellquant.py`.

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
You need **3.11 or 3.12** (3.13 is not supported). Rebuild the environment from the spec file — this installs the dependencies too, which a bare `conda create` would not:
```bash
conda env remove -n cellquant
conda env create -f environment.yml
conda activate cellquant
```

### "TypeError: remove_small_objects() got an unexpected keyword argument 'max_size'"

Your scikit-image version is too old. cellquant requires scikit-image >= 0.24.

**Fix:**
```bash
pip install "scikit-image>=0.24"
```

### Connection timeout or hang on first run

On the first real run (not `--help`), Cellpose needs to download its segmentation model (~500 MB). This requires an internet connection and may take several minutes on slow connections.

**Fix:**
- Check your internet connection.
- If you're behind a proxy or firewall, ensure Python can reach the Cellpose model server.
- Be patient — the download will show a progress bar. Once complete, subsequent runs work offline.

## Runtime errors

### "[warn] MPS GPU not supported by cpsam Transformer; using CPU"

**This is not an error.** This is an informational warning on Apple Silicon Macs. The Cellpose Transformer model doesn't support Apple's GPU backend, so the pipeline automatically uses CPU. Everything works correctly, just slower.

### "pretrained model ... not found, using default model"

Cellpose couldn't find the specified model file. The pipeline will fall back to the default model, which should work fine.

**If you want a specific model:** cellpose 4.x ships exactly one model, `cpsam`, and the older `model_type` names (`cyto3`, `nuclei`, …) no longer exist. `--pretrained-model` accepts `cpsam` (the default) or an absolute path to a checkpoint you trained yourself; anything else falls back to cpsam with the warning above. Installed checkpoints live in `~/.cellpose/models/`.

### "[error] Colocalization was requested on 2D / MIP input."

cellquant refuses `--colocalization` on 2D input and exits. This is deliberate, not a bug. Pearson's R and Manders' M1/M2 are defined on the 3D voxel distribution; a maximum-intensity projection collapses Z, so two proteins that never share a plane can overlap perfectly in the projection and inflate the coefficient by geometry rather than biology.

**Fix, in order of preference:**
- Re-run on the source z-stacks. cellquant computes colocalization natively in 3D and does not gate it there.
- If you only have projections, add `--allow-2d-colocalization`. The run proceeds and every row of `colocalization.csv` is stamped `projection_derived=True`. Treat those values as qualitative.

### "[error] 3D mode but no voxel size is available from --voxel-size or from OME/ImageJ metadata."

Your input is a z-stack, but nothing told cellquant how big a voxel is. It will not assume 1.0 µm cubic voxels, because that silently makes every volume, every micron distance, and the anisotropy correction on the puncta filter wrong while still producing believable-looking numbers.

**Fix — pick one:**
- `--voxel-size XY_UM Z_UM` — **both values, lateral first, axial second.** Look them up in your acquisition software.
- Re-export the images as OME-TIFF (or ImageJ TIFF) carrying the pixel size.
- `--project-z max` — flatten each stack and run the 2D pipeline instead.
- `--assume-isotropic` — proceed with 1 µm cubic voxels. Volumes and distances are then **voxel counts wearing micron labels**. Use it to get a run out, not to publish from.

Note that cellquant only accepts metadata carrying *both* a lateral and an axial size. An OME-TIFF with `PhysicalSizeX` but no `PhysicalSizeZ` counts as no metadata and lands you here.

### "[error] --voxel-size disagrees with the file's OME/ImageJ voxel metadata by >1%"

You passed `--voxel-size` and the file also carries voxel metadata, and they differ by more than 1% on some axis. Almost always this is a typo or the wrong input folder — for example typing `0.10` for a file whose metadata says `0.10571` is 5.4% off.

**Fix:** drop `--voxel-size` and let the metadata win, or correct your numbers. If you genuinely know the metadata is wrong, `--override-metadata` forces your values through; the override is recorded in `provenance.json`.

### "[error] you supplied these values in microns, but no pixel size is available"

A 2D run where you typed a micron-denominated value that cellquant cannot convert to pixels, because the image carries no lateral pixel size. It covers exactly two things: `--puncta-compartment-erode-um`, and any `~UM` pad inside a `--compartment` definition.

**Fix:** pass `--voxel-size XY_UM` (in 2D, the lateral size alone is enough), or drop the micron-denominated argument.

A micron value that came from a **cell-type preset** rather than from you does not abort — it warns and is skipped. That is what you see on the shipped mammalian example data, which carries no pixel size:

```
[warn] --puncta-compartment-erode-um = 0.5 µm comes from the 'mammalian' preset,
       but this 2D input carries no pixel size ... It is being SKIPPED
```

See [QUICKSTART](QUICKSTART.md) for the longer explanation. The asymmetry is intentional: a value *you* typed is an assertion about physical units that cannot be honoured, so the run stops; a preset default should not block a first run on unlabelled data.

### "file is 3D on disk but --mode 2d was forced without --project-z"

You forced `--mode 2d` on a z-stack without saying how to flatten it.

**Fix:** drop `--mode 2d` and let auto-detection run the 3D pipeline, or add `--project-z max` (or `sum`/`mean`) to project each stack before analysis.

### "[warn] cell shape: cells look flattened along Z"

Not fatal, but take it seriously. Your segmented cells came out oblate with their long axis lying in the imaging plane. That is the fingerprint of an under-scaled Z spacing — usually `--voxel-size` given in the wrong order, or a Z step that is too small.

**Fix:** check `--voxel-size XY Z` against the file's metadata. **The order is lateral first, axial second.** On the yeast calibration image, the correct voxel size (0.10571 / 0.23 µm) gives a median axis ratio of 1.66 with |cos Z| 0.98 and stays silent — that mild prolate-along-Z elongation is optical PSF and is expected. A wrong voxel size (0.094 / 0.10) gives ratio 1.41 with |cos Z| 0.17, and warns.

This check runs only for cell types it has been calibrated on (currently the `yeast` preset). Adherent mammalian cells are genuinely much wider than they are tall, and bacterial rods genuinely lie in the imaging plane, so for those the same signature appears at the *correct* voxel size and the check would be pure noise. Its silence under `--cell-type mammalian` is therefore not evidence that your voxel size is right — verify it against `provenance.json`, which records the value used and where it came from.

### "No input images found: ... matched no files in [path]"

The pipeline looks for `.tif` and `.tiff` files in the specified directory.

**Fix:**
- Check the path is correct: `ls /your/path/`
- Make sure files end in a supported, lower-case extension: `.tif`/`.tiff` (always), or `.nd2`/`.czi`/`.lif` when the `bioio` readers are installed (see [INSTALL.md](INSTALL.md)). Upper-case extensions like `.TIF` are not matched.
- The pipeline does not search subdirectories
- Check you have not passed a `{condition}` pattern to `--file-pattern`. That flag is a discovery **glob** (`--file-pattern "*.tif"`) and replaces the default extension sweep; the flag that parses conditions and replicates out of filenames is `--filename-pattern`. Passing a brace pattern to `--file-pattern` matches nothing and produces exactly this error.

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
