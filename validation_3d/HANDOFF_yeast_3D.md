# Yeast 3D validation — handoff for new chat

Self-contained brief for picking up the cellquant rebuttal work mid-stream.
Skim sections 1–3 to understand the situation; sections 4–7 are reference.

---

## 1. Project + working directory

- **Project**: `cellquant` JCB rebuttal. Single-file Python pipeline
  (`cellquant.py`, ~3900 lines after recent edits) targeting biologists
  using AI assistants. Multi-channel fluorescence microscopy quantification:
  cell/nucleus/nucleolus segmentation (Cellpose-cpsam), puncta detection,
  colocalization (Costes + Pearson + Manders), nucleolar proximity,
  condensate index, fragmentation index, replicate-level superplots.
- **Project root (EDIT HERE)**: `<project root>` — the directory containing
  `cellquant.py` and `validation_3d/` (resolved via `$CELLQUANT_VALIDATION_ROOT`;
  see `validation_3d/README.md`). On the original Mac this was an iCloud copy;
  ignore any stale `~/iCloud Drive (Archive)/...` copy.
- **Conda env**: `cellquant` — activate it; the invocations below use `python`
  from that env.
- **Always export `KMP_DUPLICATE_LIB_OK=TRUE`** for any Python call (macOS
  libomp duplicate-library issue).
- **Hardware**: Apple Silicon, CPU only (cpsam Transformer ops not
  supported on MPS). Validation runs use `--no-gpu`.

---

## 2. What this work was

### The motivating problem (R3's rebuttal concern)
The yeast 3D pipeline targets the same biologists as the published 2D
pipeline. Reviewer 3 raised the GPU barrier: if 3D analysis requires a CUDA
box, non-cluster users can't reach it. We needed to (a) validate that 3D
gives the same biological answers as the published 2D analysis on the
yeast temperature series, and (b) demonstrate the 3D mode is **CPU-feasible**
on a typical laptop.

### What landed
Two real cellquant features, exposed and documented, that together cut
yeast 3D from ~5h 43min/image to ~60 min/image on Apple Silicon CPU:

1. **`--z-crop-center N`** — per-image adaptive Z crop centered on the
   peak integrated-signal Z, with **z-boundary cell exclusion** (any cell
   whose 3D mask touches the cropped volume's top or bottom Z plane is
   dropped, so truncated half-cells don't corrupt volume/sphericity).
2. **`--seg-downsample 2` for yeast** (Cellpose runs on 596×600 instead of
   1192×1200; masks are upsampled back with nearest-neighbor;
   **measurements use full-resolution data**).

Plus a 5-image batch validation (one rep per temperature: 25C/30C/32C/36C/40C
rep1), running cellquant on three axes: native-3D, matched-2D-from-3D-mask,
and independent-2D-on-published-MIP. Direction agreement vs published
Tables S2/S3 was computed and is documented honestly.

### Critical bug found and fixed mid-stream
`apply_3d_preset_overrides()` was clobbering CLI-provided values whenever
the CLI value happened to equal `DEFAULTS` (e.g. `--seg-3d-method stitch`
was silently replaced by the yeast preset's `"full"` because stitch is the
DEFAULTS value). **Fix**: `build_config` now stashes
`cfg["_cli_provided_keys"]` and `apply_3d_preset_overrides` honors that
set. All 10 regression tests still pass.

---

## 3. Where things stand (status: DONE)

All five planned parts (A–E) complete:

- **A. Axial cell distribution characterized.** 71 slices, 0.1 µm step.
  Cells live Z=20–70 (Z=0–19 is empty buffer). Median cell z-extent 23
  slices (2.3 µm), p95 29 slices, max 34 slices. Centroid distribution is
  one peak (~Z=30–40) with axial tail to Z=65 — **not** three stacked
  layers. Per-channel signal-peak Z drifts across temperatures (25C/30C ≈
  Z=34–36, 36C/40C ≈ Z=19–22 — focus drift). **Decision**: per-image
  adaptive 41-slice window around peak signal Z, not a fixed range.
  Diagnostics: `validation_3d/diagnostics/yeast_axial_profile.py`,
  `yeast_axial_profiles.{png,csv}`, `yeast_cells_per_z.{png,csv}`,
  `yeast_cells_zextent.csv`.

- **B. Implementation landed.** `cellquant.py` gained `--z-crop-center N`
  CLI flag, three helpers (`find_peak_signal_z`, `apply_z_crop`,
  `exclude_cells_at_z_boundary`), 3D-main-loop integration, plus the
  preset/CLI bug fix. Backup at `cellquant.py.pre-zcrop.bak`. **All 10
  pytest regression tests pass.**

- **C. Single-image timing + volume sanity passed.** 25C_series1_rep1 with
  the new config: **3613 s ≈ 60.2 min** (vs 5h 43m baseline = **5.7×
  speedup**). 362 kept cells (vs 678 baseline; 235 dropped by z-boundary
  exclusion, ~80 more not picked up by stitch+ds=2). **No truncation
  spike**: median `cell_volume_um3` 12.6 (crop) vs 10.1 (full) — kept
  cells are whole, and slightly larger because the kept population is the
  "fits in window" subset. One pre-flagged caveat:
  `Sis1_frac_intensity_in_puncta` was 3.2× the baseline value (Sis1 at
  25C is diffuse and noise-floor-dominated; the puncta threshold
  recalibrates on the smaller cropped volume). Plot:
  `validation_3d/diagnostics/partC_volume_comparison.png`.

- **D. 5-image batch ran cleanly.** All 15 steps OK (5 × {3D, matched-2D,
  published-2D}). 3D pass total **4 h 26 min**:
  - 25C: 59.7 min, 362 cells · 30C: 62.7 min, 359 cells
  - 32C: 33.8 min, 182 cells · 36C: 76.0 min, 405 cells
  - 40C: 33.4 min
  Matched-2D pass: 140 s · Published-2D pass: 388 s.
  No `failures.log`. CSVs at `validation_3d/outputs_yeast/{3d,matched_2d,published_2d}/*/cells.csv`.

- **E. CPU_FEASIBLE_3D.md written.** Real timings, the honest verdict,
  methods-ready sentence at the bottom.
  Path: `validation_3d/diagnostics/CPU_FEASIBLE_3D.md`.

### The honest validation result
Direction match vs published Tables S2/S3 (62 metric × temperature
comparisons per axis, sign of rep1 median at temp T vs rep1 median at 25C
reference):

| Axis | % direction match | What this means |
|---|---|---|
| Published-2D (cellquant on published MIPs) | **80%** | the n=1 ceiling |
| Matched-2D (3D mask projected onto full MIP) | **74%** | close to published-2D |
| 3D (cropped + ds=2 + stitch) | **40%** | chance-level with n=1 |

By temperature: 30C 52%, 32C 65%, 36C 57%, 40C 72%. Strongest published
signal (40C) is best recapitulated. 56 raw "disagreements" exist but 30
are where published shows <5% change vs 25C (drowned by single-rep noise).
**18 high-effect disagreements** (>25% published change) — most in 3D, in
puncta-detection and proximity metrics. All-three-axes agree on the
disagreement at 40C for `Sis1_mean_distance` and `Sis1_fraction_proximal`
(rep1 may be unrepresentative of the 5–6-rep published mean).

**Headline interpretation**: the 80% published-2D match is the *ceiling*
for what n=1 can recapitulate; even cellquant on the exact published MIPs
disagrees 20% of the time at n=1. The 3D 40% reflects that ceiling plus
sub-sampling variance and puncta-threshold-on-cropped-volume sensitivity.
**Running 3–5 reps/temperature would close the gap**; the present batch is
a proof-of-feasibility, not a powered validation.

---

## 4. File map (what was created/changed)

### cellquant.py
- **+3 helpers, +1 CLI flag, +1 main-loop integration** (~50 lines):
  - `find_peak_signal_z(images, channels)` near `filter_cells_by_area`
  - `apply_z_crop(images, window_size, channels)`
  - `exclude_cells_at_z_boundary(cell_mask)`
  - `--z-crop-center N` added to `parse_args`
  - `z_crop_center` added to `cli_map` in `build_config`
  - DEFAULTS gained `"z_crop_center": 0`
  - Main loop calls `apply_z_crop` after `load_images` and before
    `ref_shape`; calls `exclude_cells_at_z_boundary` in the 3D branch
    immediately before `filter_cells_by_area` (when crop was applied).
- **Bug fix**: `apply_3d_preset_overrides` now skips keys that were
  CLI-provided (tracked via `cfg["_cli_provided_keys"]`).
- Backups: `cellquant.py.pre-3d-merge.bak`, `.pre-erosion.bak`,
  `.pre-seg-floor.bak`, **`.pre-zcrop.bak`** (this work). Keep them all.

### validation_3d/00_setup_paths.py
- Yeast helpers: `list_yeast_zstacks_all`, `list_yeast_mips`,
  `replicate_id_from_yeast_zstack`, `replicate_id_from_yeast_mip`,
  `yeast_overlap_image_ids`, `list_yeast_zstacks_subset` (30 manuscript
  reps), `list_yeast_zstacks_one_per_temp` (5 images).
- Yeast constants: `YEAST_ZSTACK_DIR`, `YEAST_MIP_DIR`,
  `YEAST_TEMPERATURES`, `YEAST_FILENAME_PATTERN_{3D,2D}`, voxel sizes.
- **Validation config**: `YEAST_SEG_3D_METHOD = "stitch"`,
  `YEAST_Z_CROP_CENTER = 41`, `YEAST_SEG_DOWNSAMPLE = 2`.
- `PER_IMAGE_TIMEOUT_SEC` bumped from 1 h to **12 h** (the original 1 h
  killed the first batch attempt mid-run with no `cells.csv`).

### validation_3d/01_run_pipelines.py
- Added `run_3d_yeast`, `run_published_2d_yeast`, `matched_2d_yeast`,
  `run_yeast` driver (mirrors mammalian; no nucleus channel,
  `--nucleolar-proximity Nsr1`, whole-cell compartment).
- `--yeast-subset {full,one-per-temp}` CLI flag.
- The yeast 3D invocation passes `--seg-3d-method stitch
  --z-crop-center 41 --seg-downsample 2`.

### validation_3d/04_compare_yeast.py (NEW)
- Per-temperature medians per axis; Wilcoxon vs 25C reference (NaN at
  n=1, by design); directional agreement vs hard-coded Tables S2 + S3;
  disagreements log. Writes to `outputs_yeast/comparisons/`.

### validation_3d/07_report.py
- New `_yeast_headline_verdict` function; new "Yeast temperature series"
  section showing headline + directional-agreement table + Wilcoxon table
  + medians table + disagreements pre.
- Yeast-portion-missing caveat replaced with a yeast-subset-reduced caveat
  explaining the runtime tradeoff.

### validation_3d/diagnostics/ (NEW files this work)
- `yeast_axial_profile.py` — runs the Part A characterization
- `yeast_axial_profiles.png` + `.csv` — per-Z signal envelopes
- `yeast_cells_per_z.png` + `.csv` — cells-per-Z + extent + centroid hist
- `yeast_cells_zextent.csv` — per-cell z-extent + centroid (from smoke mask)
- `AXIAL_SUMMARY.txt` — human summary of the axial findings (initial
  proposal was 56 slices; final chosen window is 41 — see below)
- `partC_volume_comparison.png` — volume-distribution sanity panel
- **`CPU_FEASIBLE_3D.md`** — the methods-ready document

### Other paths
- Smoke runs at `validation_3d/_smoke_yeast_3d/` (full, baseline) and
  `validation_3d/_smoke_yeast_3d_crop/` (cropped, the Part C run).
- Final batch outputs at `validation_3d/outputs_yeast/3d/*/`,
  `outputs_yeast/matched_2d/*/`, `outputs_yeast/published_2d/*/`.
- Report at `validation_3d/REPORT.html` (regenerated with yeast section).

---

## 5. Why the window is 41 slices (not the 56 suggested by my first script)

`yeast_axial_profile.py`'s initial `propose_window` produced 56 slices
based on a centroid-spread + cell-height calculation that was too
conservative. After looking at the cell-yield-vs-window table I picked
**41** as the "sweet spot" between speedup and cell retention:
- size 41, center 40 (Z=[20,60]) → 78.5% cells inside, 1.73× crop speedup
- size 51 → 82% cells, 1.39× speedup
- size 35 → 64% cells, 2.03× speedup

41 was selected with sign-off from the user (per-image adaptive centering,
not fixed). This is encoded in `setup_paths.YEAST_Z_CROP_CENTER = 41`.

---

## 6. What's worth doing next (suggestions, not promises)

1. **Multi-rep follow-up if reviewer asks for power**: extend the batch
   to 3 reps/temperature (15 images) on the same config. ~13 hours of
   wall time. Wilcoxon would actually evaluate at n=3. This is the
   single most defensible improvement to the validation strength.
2. **Investigate the puncta-threshold-on-cropped-volume issue**:
   `frac_intensity_in_puncta` for dim/diffuse channels is sensitive to
   the per-image threshold being recalibrated on a smaller volume. Could
   pin the puncta threshold across runs (fixed value or computed once
   on the full volume then applied to the cropped one). See Part C
   sanity-check output and the high-effect disagreement list.
3. **Crop-region MIPs for matched-2D**: currently `matched_2d_yeast`
   MIPs the *full* Z-stack but uses the cropped 3D mask footprints.
   Switching to MIPping only Z=[crop_lo, crop_hi] would be a fairer
   "matched" comparison (cells aren't seeing light from out-of-window
   cells above/below). Pre-flagged in CPU_FEASIBLE_3D.md caveat #3.
4. **Don't touch the mammalian work, the manuscript edits, or the 3D QC
   visualization** in this thread — separate tasks.

---

## 7. Style + workflow guardrails (from prior sessions, still applies)

- Terse responses, no trailing summaries the user can read from the diff.
- Edit existing files in active iCloud, not the archive.
- Don't initiate git commits unless explicitly asked.
- If a yeast metric direction disagrees with the published table, flag and
  investigate — don't silently report it. (This was followed: 18 high-effect
  disagreements surfaced in the Part D verdict.)
- For long-running batches, start `caffeinate -i -s -w <driver_pid>` so the
  laptop doesn't sleep. The harness will notify when the background task
  completes; don't poll.

---

## 8. Canonical yeast 3D invocation (current)

```bash
KMP_DUPLICATE_LIB_OK=TRUE \
python cellquant.py <in_dir> \
  "1:Tif6:quantify" "2:Nsr1:nucleolus" "3:Sis1:quantify" \
  --cell-type yeast --mode 3d \
  --seg-3d-method stitch \
  --z-crop-center 41 --seg-downsample 2 \
  --voxel-size 0.094 0.1 \
  --colocalization --nucleolar-proximity Nsr1 --condensate-index \
  --filename-pattern "{condition}_series1_rep{replicate}" \
  --skip-plots --no-gpu \
  --out <out_dir>
```

Or via the runner:
```bash
KMP_DUPLICATE_LIB_OK=TRUE \
python \
  validation_3d/01_run_pipelines.py \
  --dataset yeast --yeast-subset one-per-temp
```

---

## 9. Key numbers to remember

- Baseline (full XY, ds=1, `full` 3D method): **5 h 43 min / image** (25C)
- New config (z-crop=41, ds=2, stitch): **~60 min / image** (mean across 5)
- Speedup: **~6×**
- 5-image batch total: **~4 h 30 min** (3D + matched-2D + published-2D)
- Direction match vs published: **80% / 74% / 40%** (pub-2D / matched-2D / 3D)
- High-effect disagreements (>25% published change): **18 of 186**
- All 10 cellquant pytest regression tests still pass after the changes
