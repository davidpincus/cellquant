# cellquant vs CellProfiler agreement benchmark (HCT116 stress granules)

**Part A of the Figure 5 replacement / response to reviewer point R2.4.** Runs
`cellquant` and CellProfiler *independently* on the same HCT116 MIPs, matches cells
by centroid, and reports per-cell agreement on three measures: G3BP1 puncta count,
puncta area fraction, and mean G3BP1 intensity. The claim under test is
**equivalence** (cellquant offers accessibility, not an accuracy advantage), so the
benchmark succeeds by showing agreement — and reports disagreement straight when it
finds it.

> **No parameter here is tuned to improve agreement.** Detection/segmentation
> parameters are set once from the image biology and from cellquant's own mammalian
> preset (see "Parameters & provenance"), then whatever agreement results is reported.

Everything is versioned; **only code is committed, never data** (`.gitignore` keeps
`out_cq/`, `out_cp/`, `cp_input/`, `compare_out/`, masks, CSVs and figures out of git).

---

## Files

| File | What it is |
|---|---|
| `build_cppipe.py` | Builds `cellprofiler_hct116.cppipe` programmatically via `cellprofiler_core` (run inside the container). Parameters + provenance live at the top. |
| `cellprofiler_hct116.cppipe` | The generated CellProfiler pipeline (committed so it can be run without rebuilding). |
| `compare_tools.py` | The comparison harness: centroid matching, correlations, OLS, Bland-Altman, figures. Reusable / CLI-driven. |
| `run_benchmark.sh` | One-command orchestration: cellquant → channel split → CellProfiler (Docker) → compare. Fails loud. |
| `.gitignore` | Keeps all outputs/data out of git. |

## Requirements

- **Docker** with the pinned CellProfiler image (primary path — native install on
  Apple Silicon is not viable):
  ```bash
  docker pull cellprofiler/cellprofiler:4.2.6
  ```
  *Fallback (documented, not used here):* conda env with Python **3.9** (CellProfiler
  4.2.x does not support 3.11/3.12), `conda-forge openjdk=11`, then
  `pip install cellprofiler==4.2.6`. On Apple Silicon prefer Docker.
- A **cellquant** Python env (this repo's `environment.yml`) — also used by
  `compare_tools.py` (needs numpy, scipy, pandas, matplotlib, scikit-image, tifffile;
  **no scikit-learn**).

## Inputs

- `HCT116_DIR` = the **full-resolution** 3-channel `MAX_*.tif` (1192×1200), e.g.
  `SG_zstacks/MIPs/` — **not** the 800×800 `example_data/` crops. Both tools see the
  identical pixel data.
- **Channel order** is `1=DAPI, 2=G3BP1, 3=PABPC1`. The TIFFs carry no embedded
  channel *names* (ImageJ composite; `spacing=0.5` µm matches the paper's HCT116), so
  this order follows the paper's validated `cellquant` invocation and prior analysis.
  CellProfiler cannot easily consume a 3-channel stack, so `run_benchmark.sh` splits
  each stack into byte-identical single-channel files (`MAX_..._{DAPI,G3BP1,PABPC1}.tif`).

## Run it

```bash
./run_benchmark.sh
# or override any of:
HCT116_DIR=/path/to/MIPs OUT_ROOT=/tmp/bench ./run_benchmark.sh
```

This runs, in order: cellquant (fresh) → split channels → (re)build the `.cppipe` in
the container → CellProfiler headless → `compare_tools.py`. Outputs land in
`out_cq/`, `out_cp/`, and `compare_out/`.

Run the comparison alone (after both tools have produced output):

```bash
python compare_tools.py \
  --cellquant-cells out_cq/cells.csv \
  --cellquant-masks out_cq/masks \
  --cp-cells out_cp/cp_Cells.csv \
  --cp-puncta out_cp/cp_Puncta.csv \
  --cp-image out_cp/cp_Image.csv \
  --out-dir compare_out
```

`--match-tol-px` overrides the default tolerance; `--cq-*-col` / `--cp-*-col` override
column names. **`--cellquant-masks` is required** because cellquant's `cells.csv` has
no centroid column — centroids are derived from the integer cell masks
(`{image}_cellmask.tif`, label == `cell_id`). `--cp-image` maps CellProfiler's
`ImageNumber` → filename.

## Pipeline (native segmentation — no RunCellpose)

`Images → Metadata → NamesAndTypes → Groups → IdentifyPrimaryObjects (Nuclei ← DAPI,
Otsu) → IdentifySecondaryObjects (Cells ← Nuclei, propagation, PABPC1 cell body) →
EnhanceOrSuppressFeatures (Speckles on G3BP1) → IdentifyPrimaryObjects (Puncta ←
enhanced G3BP1) → RelateObjects (Puncta→Cells) → MeasureObjectSizeShape →
MeasureObjectIntensity → ExportToSpreadsheet`.

CellProfiler segments cells with a *classical* Otsu-nuclei + propagation pipeline —
genuinely **independent** of cellquant's Cellpose, which is the point of the benchmark.

### Parameters & provenance (set once; not tuned to agreement)

| Parameter | Value | Source |
|---|---|---|
| Speckle feature size | **6 px** | cellquant mammalian `log_sigma = 2.0`; a LoG of σ maximally responds to a blob of diameter 2√2·σ = 5.66 → 6. |
| Nuclei diameter | 60–250 px | measured HCT116 nuclei ~110–185 px here |
| Puncta diameter | 3–30 px | lower bound = cellquant `puncta_min_area_px` 6 → equiv-diam ~3 |
| Puncta threshold | Manual **0.003** on the enhanced image | fixed/absolute (adaptive/Otsu detect cytoplasmic texture as "puncta" in control → biologically wrong). Calibrated *once* to the enhanced-image distribution (granule residuals ~0.004–0.006, background ~10× lower), not to cellquant. |
| Cell-body image | **PABPC1** | brighter/uniform cytoplasm marker → robust propagation; also cellquant's own mammalian seg channel (Table S4). (Spec suggested G3BP1; changed after checkpoint because G3BP1 is diffuse in control cells.) |

## Matching & statistics (`compare_tools.py`)

- Cells matched **within each image** by **mutual nearest neighbour** (`scipy cKDTree`):
  a pair is kept only if each is the other's nearest and within tolerance. Default
  tolerance = **0.25 × median cell equivalent diameter** (`d = 2√(area/π)`),
  overridable with `--match-tol-px`. Unmatched cells are counted, never dropped
  silently; images with zero matches warn.
- Per measure, on matched pairs: **Pearson r (+p), Spearman ρ (+p), OLS slope/intercept
  (cellquant ~ cellprofiler), Bland-Altman bias + ±1.96 SD LoA**.
- Puncta area fraction: cellquant `G3BP1_puncta_area_px / cell_area_px`; CellProfiler
  Σ(child-puncta `AreaShape_Area`) / cell `AreaShape_Area` via `Parent_Cells`.
- **Intensity caveat:** CellProfiler rescales intensity to 0–1; cellquant reports raw
  16-bit counts. Pearson r and the OLS slope are scale-invariant and are the primary
  intensity metrics; the intensity Bland-Altman bias reflects the unit difference and
  is labelled as such. Count and area fraction are on shared units (BA is direct).
- Sanity check: `Pearson(intensity) ≥ Pearson(area) ≥ Pearson(count)` roughly — a
  violation is flagged as a likely harness bug (intensity is a mask-average and should
  agree best; count is noisiest).

Outputs: `benchmark_summary.csv` + `.json`, `matches.csv` (per-pair raw evidence), and
`agreement_combined.{png,pdf}` (+ standalone scatter / Bland-Altman), plus a stdout
summary table.

## Result (full 7-image HCT116 set, n = 62 matched cells)

| measure | Pearson r | Spearman ρ | OLS slope | Bland-Altman |
|---|---|---|---|---|
| **mean cell G3BP1 intensity** | **0.95** (p≈9e-32) | 0.95 | 66 233 (≈2¹⁶ rescale) | +809 (unit diff — read r/slope) |
| puncta area fraction | 0.22 | 0.34 | 0.20 | bias +8.4e-4, LoA [−5.9e-3, +7.6e-3] |
| G3BP1 puncta count | 0.12 | 0.20 | 0.08 | bias +1.84/cell, LoA [−7.0, +10.7] |

Segmentation match rate: cellquant 0.78 / CellProfiler 0.72 (62/79, 62/86), median
centroid distance 13 px.

### How to read this (equivalence, honestly)

- **Mean cell intensity agrees tightly** across tools (r = 0.95) — a pixel average in a
  matched mask is essentially the same measurement in both. The intensity result also
  validates the harness end-to-end (matching + column mapping + alignment).
- **Both tools recover the population-level biology**: arsenite ≫ control puncta in
  each (cellquant 5.1 vs 2.35, CellProfiler 3.6 vs 0.21 puncta/cell).
- **Per-cell puncta count and area fraction do *not* correlate** (r ≈ 0.1–0.2). The two
  tools *define a punctum differently* — cellquant detects Laplacian-of-Gaussian blobs
  in the nucleus-excluded, 0.5 µm-eroded cytoplasm; CellProfiler thresholds
  tophat-enhanced speckles over the whole cell. With low per-cell counts (median 3),
  those definitional differences dominate single-cell values while cancelling in the
  population mean. **Neither tool is ground truth**, so this is a difference in
  definition, not an accuracy gap.

**Bottom line for the manuscript:** equivalence holds for what a fixed mask measures
(intensity) and for the population-level condensate response; per-cell puncta *counts*
are method-dependent in both tools and should be compared at the replicate/population
level, not cell by cell. No accuracy claim is made over CellProfiler.

## Caveats

- **Not tuned to agreement.** Parameters are fixed from biology / cellquant's preset;
  the puncta threshold was calibrated once to the enhanced-image distribution before
  any tool comparison.
- **Small n** (7 images, 62 matched cells). Per-cell count/area correlations are
  low-power; the population-level and intensity results are the robust signals.
- **Channel order** rests on the paper's validated invocation (no embedded channel
  names in the TIFFs); G3BP1 vs PABPC1 cannot be distinguished from pixels alone.
- **cellquant Cellpose is non-deterministic across platforms**, so exact cell counts
  can vary run to run; the agreement statistics are stable to that.
- CellProfiler counts whole-cell puncta (`RelateObjects`); cellquant counts
  cytoplasm-compartment puncta — one contributor to the per-cell count difference.
