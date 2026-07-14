# What cellquant Measures (plain language)

This document is for users who want to know what each cellquant metric *means* before they trust the number it spits out. It is deliberately non-technical — no code, no equations unless absolutely necessary. If you're an experienced quantitative image-analysis person, you already know all of this and can skip to [CLI_REFERENCE.md](CLI_REFERENCE.md).

cellquant is a tool, not an oracle. The numbers it produces are only meaningful if (a) the segmentation behind them was sensible, (b) the assumptions baked into each metric match your biology, and (c) you are comparing like with like across conditions. The sections below walk through every metric cellquant reports, what it tries to capture, what would have to be true for the number to be trustworthy, and the situations in which it can mislead.

## How to read this guide

Every metric below has the same three subsections:

- **What it is** — one sentence saying what the number represents.
- **When it's reliable** — the conditions under which the number actually means what it sounds like.
- **When it can mislead** — the situations where the number looks normal but the underlying biology is being misrepresented. **Read this one carefully.**

If you finish a section and you're not sure your data meets the "reliable" conditions, treat the metric as exploratory and validate independently before drawing conclusions.

## Before any of this matters: the segmentation

Every per-cell number cellquant reports starts with a cell mask from Cellpose. If the segmentation is wrong, every downstream number is wrong, no matter how internally consistent the pipeline is. **Always look at the `qc/*.png` overlays first.** If cells are obviously merged together, obviously split in half, or systematically missed, fix the segmentation (tune `--cell-diameter`, `--flow-threshold`, `--cellprob-threshold`, or for 3D `--seg-3d-method`) before you trust any cells.csv values.

Visually inspect at least one QC image per condition. If the segmentation looks great on the negative control and broken on the positive condition, that's a confound, not biology.

---

## Cell-level basics

### `cell_area_px` (2D) / `cell_volume_um3` (3D)

**What it is** — the size of the cell mask in pixels (2D) or in cubic microns (3D, derived from voxel counts × voxel size).

**When it's reliable** — when the segmentation matches what you'd outline by hand. In 3D, also when the voxel size is correct (cellquant reads it from OME/ImageJ metadata or you pass `--voxel-size XY Z`; if neither is set, the value is in raw voxel counts and you should not interpret `_um3`).

**When it can mislead** — when Cellpose over- or under-segments. Merged cells inflate volume; over-split cells deflate it. A bimodal volume distribution across conditions is often a segmentation artifact, not biology.

### `n_nuclei`

**What it is** — the number of nuclei detected inside this cell's mask. Used as a "real cell" filter: by default cellquant keeps cells with `n_nuclei` between `--keep-min-nuclei` and `--keep-max-nuclei`.

**When it's reliable** — when there's a clean nuclear stain (DAPI, Hoechst) and a `nucleus`-role channel is declared.

**When it can mislead** — when nuclei are clustered tightly enough to merge, or when out-of-focus nuclei from neighboring cells leak into a cell's mask. For yeast and bacteria, nuclear segmentation is usually skipped entirely and `n_nuclei = 0` is normal.

---

## Intensity

### `{channel}_cell_mean`, `{channel}_nucleus_mean`, `{channel}_cytosol_mean`

**What it is** — the average pixel intensity of `{channel}` over the cell mask, the nuclear mask, or the cytosolic mask (cell minus dilated nucleus), respectively.

**When it's reliable** — when image acquisition is consistent across conditions: same exposure, same laser power, same gain, same day if possible. Means are sensitive to all of these.

**When it can mislead** — across imaging sessions where any acquisition parameter differs, means become apples-to-oranges. cellquant cannot detect this; you have to ensure it on the microscope side. Bleed-through from another channel is also invisible to the mean — if your green channel picks up red signal, the mean goes up but the protein didn't change.

---

## Puncta detection

cellquant detects puncta with a **Laplacian-of-Gaussian (LoG) filter** — a method that responds to small, locally-bright spots above a smoothed background. The Gaussian smoothing scale is set by `--log-sigma`; the threshold by `--puncta-threshold-method` (Otsu by default).

### `{channel}_puncta_n`

**What it is** — the number of discrete LoG-detected puncta inside this cell.

**When it's reliable** — when puncta are well-separated, roughly circular (or spherical in 3D), and brighter than the local cytosolic background. The mammalian stress granule and yeast nucleolus systems in the cellquant examples satisfy these conditions.

**When it can mislead** — in three common failure modes:
1. **Diffuse condensation regime.** If your protein is starting to phase-separate but hasn't formed discrete puncta yet — e.g. an "early stress" condition — `puncta_n` flickers between 0 and 5 essentially at random. Use the **fragmentation indices** (below) instead.
2. **Threshold sensitivity.** Otsu is data-driven; if a condition shifts the overall intensity distribution, Otsu's threshold shifts with it, and counts move for reasons that have nothing to do with discrete puncta. Try `--puncta-threshold-method triangle` or `--puncta-threshold-method fixed --puncta-threshold-fixed 0.30` and check whether your conclusion survives.
3. **Saturated images.** If puncta saturate the detector, LoG can't resolve them and reports too few. Re-acquire below saturation.

### `{channel}_puncta_area_px` (2D) / `{channel}_puncta_volume_um3` (3D)

**What it is** — total area (or volume) of detected puncta in the cell.

**When it's reliable** — same conditions as `puncta_n`. Use it as a sanity check: if `puncta_n` goes up but `puncta_area` doesn't move much, the new puncta are smaller; if both go up together, the cell has more *and* larger puncta.

**When it can mislead** — same failure modes as `puncta_n`. Additionally, in 2D, MIP-based puncta area underestimates true 3D extent of puncta that span multiple Z planes.

### `{channel}_puncta_mean_intensity`, `{channel}_puncta_integrated_intensity`, `{channel}_diffuse_mean_intensity`, `{channel}_frac_intensity_in_puncta`

**What they are** — for the pixels inside detected puncta: their average intensity (`puncta_mean`), their summed intensity (`puncta_integrated`), and the fraction of total channel signal in the cell that's in puncta (`frac_intensity_in_puncta`). For pixels in the puncta-detection compartment but outside puncta: their mean (`diffuse_mean`).

**When they're reliable** — same as `puncta_n`: discrete, well-separated puncta on a clean background.

**When they can mislead** — `frac_intensity_in_puncta` is sensitive to how aggressive the puncta segmentation is. A more permissive threshold pulls in more pixels and pushes this value up. Compare conditions only when puncta segmentation parameters are identical across them.

### `{channel}_fragmentation_index_simple`

**What it is** — number of connected components above a single LoG-Otsu threshold within the cell. Threshold-dependent. A drop-in replacement for `puncta_n` when the puncta volume cap is the limiting factor — useful when one or two puncta dominate the cell but you care about how *fragmented* the signal is.

**When it's reliable** — same regimes as `puncta_n`.

**When it can mislead** — same as `puncta_n`; threshold sensitivity in particular.

### `{channel}_fragmentation_index_persistence`

**What it is** — a **threshold-free** measure of how fragmented the LoG response is. cellquant sweeps a range of threshold levels, counts connected components at each level, and integrates the count-vs-threshold curve into a single number normalized to roughly [0, 1]. A spatially uniform cell decays from 1 → 0 quickly across thresholds and gives small values; a cell with discrete spots sustains many components across thresholds and gives large values; a cell with diffuse-but-clumpy condensation lands in between.

**When it's reliable** — this is the metric cellquant recommends for **chronic-stress and diffuse-condensation regimes** where `puncta_n` flickers between 0 and 5. It also degrades gracefully when puncta count is well-defined: it's monotonic with `puncta_n` for clean punctate biology, just on a different (continuous) scale.

**When it can mislead** — less than the simpler metrics, but it inherits the LoG filter's spatial scale (`log_sigma`) sensitivity. Holding sigma constant across conditions is essential.

---

## Condensate index (opt-in via `--condensate-index`)

### `{channel}_condensate_index_cell`, `{channel}_condensate_index_cytosol`

**What it is** — the ratio of the 95th-percentile intensity to the mean intensity, computed across the cell mask (cell) or the cytosolic mask (cytosol). Higher values indicate a "clumpier" intensity distribution — a few bright pixels lifting the percentile above the average. Lower values indicate uniform signal.

**When it's reliable** — when you genuinely care about the *internal contrast* of a channel as a proxy for condensation (a phase-separating protein concentrating into condensates raises this; the same protein staying soluble keeps it flat). It's most useful as a complement to `puncta_n`/`frac_intensity_in_puncta`, not a replacement.

**When it can mislead** — three regimes:
1. **Diffuse vs punctate confound.** A cell with many small dim puncta and a cell with one large bright punctum can give the same condensate index even though the biology is different.
2. **Cross-cell-type comparisons.** Yeast and mammalian cells differ by orders of magnitude in cytosolic volume; the percentile-to-mean ratio is not directly comparable.
3. **Sample size dependence.** Cells with very few pixels (over-segmented fragments) give noisy percentile estimates. cellquant skips cells with fewer than 50 pixels; tighter filters may be appropriate for your data.

The Condensate Index is **off by default** as of the 2026 revision — it was previously always computed, which surfaced an irrelevant column for users not studying condensates. Add `--condensate-index` to bring it back.

---

## Colocalization (opt-in via `--colocalization`)

### `pearson_r_{A}_vs_{B}`

**What it is** — Pearson's correlation coefficient between the per-voxel intensities of channels A and B, computed inside the selected compartment (whole-cell by default, override with `--colocalization-compartment`).

**When it's reliable** — when computed on 3D voxel data with a clean, segmented compartment, and with both channels acquired at non-saturating exposure. cellquant computes it natively on z-stacks when given 3D input.

**When it can mislead** — **most importantly: when computed on maximum intensity projections.** A MIP collapses Z and can make two distinct structures from different Z planes appear to overlap. cellquant emits a loud runtime warning when `--colocalization` is run in 2D mode for this reason. Treat MIP-based Pearson values as advisory and re-run on z-stacks for publication-grade numbers.

Pearson is also confounded by low signal-to-noise: noisy channels report low R even when the underlying biology is colocalized. And it's a linear measure — a non-linear relationship (e.g. one channel saturating while the other still rises) drives R toward zero even if the structures coincide.

### `manders_m1_{A}_vs_{B}`, `manders_m2_{A}_vs_{B}`

**What they are** — Manders' M1 is the fraction of channel-A signal (above its Costes threshold) that occurs in pixels where channel B is also above its threshold. M2 is the symmetric statement for channel B. Together they say "how much of A is in regions where B is also there" and vice versa.

**When they're reliable** — same conditions as Pearson: 3D-native computation, clean compartments, non-saturating exposure. Manders is less sensitive than Pearson to brightness mismatches.

**When they can mislead** — the Costes threshold is a heuristic and can land in a bad spot when one channel has very different dynamic range from the other. The same MIP warning applies — running Manders on MIPs is statistically unsound and cellquant warns when you do it. Manders also depends on the segmentation of the compartment: if the compartment mask includes a lot of background, both M1 and M2 trend low for the wrong reason.

---

## Nucleolar morphometrics (only when a `nucleolus`-role channel is declared)

### 2D: `nucleolar_area`, `nucleolar_circularity`, `nucleolar_solidity`, `nucleolar_eccentricity`

**What they are** — geometric descriptors of the largest detected nucleolus in each cell: pixel area, perimeter-based circularity (1 = perfect circle), solidity (area / convex-hull area; 1 = no concavities), and eccentricity (0 = circle, 1 = line).

**When they're reliable** — when the nucleolar marker (e.g. Nsr1, fibrillarin) gives a clean signal contained inside the nucleus.

**When they can mislead** — when the marker bleeds into the nucleoplasm or when there are multiple nucleoli that get merged by the segmentation. Reported "circularity" of a multi-lobed nucleolus is geometric, not biological.

### 3D: `nucleolar_volume_um3`, `nucleolar_sphericity`, `nucleolar_eq_diameter_um`

**What they are** — volume in cubic microns; sphericity (1 = perfect sphere, lower = more irregular surface) computed from a marching-cubes surface mesh; equivalent diameter (diameter of a sphere with the same volume).

**When they're reliable** — when voxel size is correct (set `--voxel-size XY Z` or rely on OME metadata), and when the nucleolar segmentation is clean. Sphericity is well-defined for closed surfaces.

**When they can mislead** — voxel-size mismatch is the most common failure: a value that looks normal but is in units of voxels rather than microns. cellquant warns if XY=Z=1.0 µm (i.e. no metadata and no `--voxel-size`); take that warning seriously.

---

## Nucleolar proximity (`--nucleolar-proximity`)

### `{channel}_mean_distance`, `{channel}_min_distance`, `{channel}_fraction_proximal`

**What they are** — for each detected punctum in `{channel}` in a cell, the distance from the punctum's centroid to the nearest nucleolar boundary. cellquant reports the mean, the minimum, and the fraction of puncta below `--proximity-threshold` (pixels in 2D) or `--proximity-threshold-um` (microns in 3D).

**When they're reliable** — when the nucleolar mask is clean, the puncta detection is reliable, and (for 3D) the voxel size is correct so distances are in real microns. 3D distances are anisotropy-aware.

**When they can mislead** — if puncta detection is shaky, these distances are computed for the wrong objects. A cell with zero detected puncta reports `nan`, which is honest, but if you compare conditions where one has 1 punctum per cell and the other has 50 puncta per cell, the means describe different populations.

---

## Replicate-level statistics

### Pairwise p-values in `pvalues.csv`

**What they are** — Mann-Whitney U test on **replicate medians**, comparing each non-reference condition to the reference. cellquant requires at least 3 replicates per condition to compute a p-value; otherwise the comparison is skipped (the cell is empty in the output table). p-values are Bonferroni-corrected for the number of conditions tested.

**When they're reliable** — when you have at least 3 biologically independent replicates per condition, when the metric you're testing is well-defined for the cell type, and when the segmentation is comparable across replicates.

**When they can mislead** — three regimes:
1. **Under-powered tests.** With 3 replicates per condition, only large effect sizes reach significance. A non-significant p-value with N=3 is **not** evidence of no effect — it's evidence of insufficient data.
2. **Treating the cell as the replicate.** The right unit for replicate statistics is the experimental replicate (the well, the day, the mouse), not the individual cell. cellquant correctly aggregates to replicate medians before testing, but if all your "replicates" are technical replicates from the same well, the test is misleading.
3. **Multiple metrics, no correction.** Bonferroni corrects within-metric across conditions but not across the many metrics cellquant reports. If you scan 20 metrics looking for a significant one, the family-wise error rate is much higher than the per-metric p suggests.

---

## What cellquant does *not* do

To be explicit:

- cellquant does **not** validate that the metric you're computing is the right one for your biology. If you ask for `puncta_n` and your protein is doing something subtler than discrete punctation, the number is meaningful but the question is wrong.
- cellquant does **not** detect bleed-through, photobleaching during acquisition, focal drift, or other acquisition artifacts. These produce coherent-looking outputs that are nevertheless wrong.
- cellquant does **not** decide whether your replicate count is sufficient. It will compute a Mann-Whitney U on 3 replicates and report a p-value; it is your responsibility to know that with N=3 you are only powered to detect large effects.
- cellquant does **not** know which compartment is biologically relevant. The default puncta compartment is `cytosol` (in mammalian) or `whole-cell` (in yeast/bacteria); if your puncta are nuclear, override with `--puncta-compartment nucleus`.
- cellquant does **not** substitute for understanding the analysis. Treat each metric as a hypothesis the tool helped you test, not as a verdict.

If you read a number from `cells.csv` and you cannot explain in one sentence what the number means and what it would have to be true for it to be trustworthy, **stop and re-read the relevant section above**. That is the failure mode this document is trying to prevent.

---

## Where to go next

- [TUTORIAL_1_mammalian_SGs.md](TUTORIAL_1_mammalian_SGs.md) — worked example end-to-end on stress granule data
- [TUTORIAL_2_yeast_temp.md](TUTORIAL_2_yeast_temp.md) — yeast temperature series with colocalization and proximity
- [TUTORIAL_3_four_condition.md](TUTORIAL_3_four_condition.md) — multi-condition experiments with positive/negative controls
- [CLI_REFERENCE.md](CLI_REFERENCE.md) — every flag, with examples
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) — common errors and how to diagnose them
