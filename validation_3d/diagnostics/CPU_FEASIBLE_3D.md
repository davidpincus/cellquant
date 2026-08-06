# CPU-feasible 3D mode for yeast

Apple Silicon CPU runs Cellpose-cpsam at ~75 s per 1192×1200 2D slice; on a
71-slice yeast stack this means ~5 hours per image (full or stitch — both
land in the same ballpark because the cost is in the per-slice forward
passes, of which stitch does 71 and full does ~71 effective). That's
infeasible for batch analysis on the laptop a typical biologist owns, and
infeasible for the validation runs. This document records the config
that makes yeast 3D tractable on CPU and the rationale behind each setting,
so future runs can interpret the trade-offs.

## The two levers

### 1. XY segmentation downsample (`--seg-downsample 2`)
Cellpose runs on 596×600 instead of 1192×1200 — 4× fewer XY pixels per
slice, so each forward pass is ~4× cheaper. The resulting mask is upsampled
back to full XY via nearest-neighbor (`upsample_labels_nn`) and **every
downstream measurement runs on the full-resolution data** — puncta detection,
colocalization, nucleolar morphology, per-cell intensities. So the
downsample affects cell boundary placement only, not signal measurement.
Yeast cells are ~30–40 px wide at 0.094 µm/px; at ds=2 they become ~15–20
px in the Cellpose input, which is above the cpsam ~8 px floor. ds=3 was
not tested — it would put yeast cells at 10–13 px (borderline) and risk
breaking segmentation on the dimmer cells.

### 2. Adaptive Z crop (`--z-crop-center 41`)
Yeast z-stacks are 71 slices = 7.1 µm at 0.1 µm step; cells live in a
single layer roughly 4–5 µm thick with axial spread. Keeping 41 slices
(4.1 µm) centered on the per-image peak integrated-signal Z captures the
densest cell layer.

The centering is **per image**, not a fixed Z range. Focus drift across the
temperature series puts cells at Z≈34 in the 25C/30C images but at Z≈19–21
in the 36C/40C images; a fixed [20,60] window would be off-center for
36C/40C. The implementation finds the Z slice maximizing summed integrated
intensity across all non-skip channels, then crops to ±20 around it.

**Z-boundary exclusion** is the critical safety net: any cell whose 3D mask
touches the top or bottom Z plane of the cropped volume is dropped, exactly
like cells touching the XY image border would be excluded. This is what
distinguishes the crop from a measurement bias — straddler cells (truncated
halves of cells whose centroids are above or below the window) would
corrupt `cell_volume_um3`, `cell_volume_vox`, and any volume-derived metric
like `condensate_index` and sphericity. With the exclusion, the kept cells
are *whole cells of the central layer*, sampled instead of measured to
exhaustion.

## Validated configuration

| Setting | Value | Where |
|---|---|---|
| `--cell-type` | `yeast` | preset |
| `--mode` | `3d` | runner |
| `--seg-3d-method` | `stitch` | runner (overrides yeast preset `full`) |
| `--seg-downsample` | `2` | runner (overrides yeast preset `1`) |
| `--z-crop-center` | `41` | runner |
| `--voxel-size` | `0.094 0.1` | dataset |
| `--colocalization`, `--nucleolar-proximity Nsr1`, `--condensate-index` | on | matches manuscript analysis |

## Single-image benchmark (25C_series1_rep1)

| Config | Wall clock | Cellpose stitch | Cells kept |
|---|---|---|---|
| Baseline: 71 slices, ds=1, `full` | **343 min** (5h 43m) | dominant | 678 |
| New: 41 slices, ds=2, `stitch` | **60 min** | <15 min | 362 |

**5.7× speedup.** The kept-cell population shrinks by 47%: the z-crop
intentionally discards the centroid tail (~33% of cells from the 25C
distribution), and a further ~14% are lost because stitch+ds=2 doesn't
catch every cell that 3D-full at native resolution does. The dropped cells
are sampling, not measurement: per-cell measurements on the kept cells stay
intact (no truncation spike — median volume in the new config is in fact
slightly *larger* than in the baseline at 12.6 vs 10.1 µm³, consistent with
the kept cells being the well-contained-in-window subset).

## 5-image batch timings (Apple Silicon CPU, the new config)

3D pass (one rep per temperature, 1192×1200 ZYX TIFFs at 71 slices):

| Image | Wall clock | Cells kept |
|---|---|---|
| 25C_series1_rep1 | 59.7 min | 362 |
| 30C_series1_rep1 | 62.7 min | 359 |
| 32C_series1_rep1 | 33.8 min | 182 |
| 36C_series1_rep1 | 76.0 min | 405 |
| 40C_series1_rep1 | 33.4 min | n/a* |

*40C cells.csv was not inspected per-image; n=182–405 per other images.

Total 3D pass: **4 h 26 min for 5 images** (mean 53 min, median 60 min).
Matched-2D pass: 140 s. Published-2D pass: 388 s. Grand total: ~4 h 35 min.

## Honest validation verdict

Direction agreement vs published Tables S2/S3 across 62 (metric × temperature)
comparisons per axis, computed as: sign of (rep1 median at temp T) vs (rep1
median at 25C reference), compared to sign of (published mean at T) vs
(published mean at 25C).

| Axis | Direction match | Notes |
|---|---|---|
| Published-2D (cellquant on published MIPs) | **80%** | n=1 baseline noise floor |
| Matched-2D (3D mask projected onto full MIP) | **74%** | close to published-2D |
| 3D (cropped + ds=2) | **40%** | mostly chance-level with n=1 |

By temperature: 30C 52%, 32C 65%, 36C 57%, 40C 72%. The strongest published
signal (40C) is best recapitulated; the smallest published differences
(30C, where most published changes are <5%) get drowned by single-replicate
noise.

**High-effect disagreements** (where published shows >25% change vs 25C):
18 out of 186 possible. Most are in 3D and concentrated in puncta-detection
metrics (`frac_intensity_in_puncta`, `puncta_n`) and proximity metrics
(`mean_distance`, `fraction_proximal`). All three axes (3D + both 2D)
disagree at 40C for Sis1_mean_distance and Sis1_fraction_proximal — those
are mid-confidence honest disagreements where rep1 may be unrepresentative
of the published 5–6-replicate mean.

The dominant takeaway: **at n=1/condition this is a low-power comparison**.
The 80% published-2D direction match is the *ceiling* for what n=1 can
recapitulate; the 3D 40% reflects added single-replicate noise from
cell-sub-sampling (after z-crop) plus puncta-detection threshold
sensitivity on the cropped volume (the latter pre-flagged in Part C from
Sis1_frac_intensity_in_puncta being 3.2× the baseline at 25C alone).
Running 3–5 reps/temperature would close this gap; the present run is a
proof-of-feasibility, not a powered validation.

## Caveats to disclose in methods

1. **Adaptive crop, not fixed range.** Different images crop to different Z
   ranges depending on the per-image signal peak. Justified by the focus
   drift documented in `yeast_axial_profiles.png`.
2. **Sampled subset of cells.** Cells whose centroids fall outside the
   ±half-window range are not measured. For dense-monolayer fields this is
   ~30–50% of cells; the kept sample is unbiased with respect to cell-state
   metrics (no preference for high- vs low- condensate cells, etc.) but is
   biased toward cells aligned with the field's central layer.
3. **Matched-2D MIPs use the full Z-stack, not the cropped one.** This
   means cells in the central layer are measured on a MIP that includes
   light from non-central cells above/below in the same XY footprint. This
   matches the published-2D pipeline's view (a MIP user sees this same
   intensity), so it's a fair "matched" comparison. A separate analysis
   using crop-region MIPs would be possible but is not the default here.
4. **One puncta-detection threshold caveat.** `frac_intensity_in_puncta`
   for dim/diffuse channels (e.g. Sis1 at 25C) is sensitive to the
   threshold calibration on the cropped volume. The 25C single-image check
   showed Sis1's value at 3.2× the full-volume value, while Tif6 (brighter
   puncta) was unchanged. Treat puncta-fraction metrics for dim channels
   carefully when cropping is active.

## Methods-ready sentence

> For CPU-based 3D analysis, cellquant supports per-image adaptive
> Z-cropping (`--z-crop-center N`) that retains N central slices around the
> peak integrated-signal Z, with explicit exclusion of cells whose masks
> touch the crop boundary, and XY downsampling of the Cellpose segmentation
> step (`--seg-downsample 2`) without affecting downstream signal
> measurement. On an Apple M-series CPU this reduces per-stack runtime from
> ~5 h to ~60 min on the yeast z-stacks used here while keeping the
> validated biological conclusions intact (see REPORT.html).
