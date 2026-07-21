# Segmentation experiments — results

**Date:** 2026-06-10 (post-erosion-fix cellquant.py)
**Cellquant change committed:** mammalian `_3d` preset `min_cell_volume_vox` 0 → **30,000**
**Hypotheses tested:** B1 (stitch_threshold), B2 (volume floor), B3 (cell-seg input), C (nuclei_diameter)
**Result:** one default change (B2), one flagged-for-review CLI workaround (B3), two ruled-out (B1, C)

---

## Baseline (existing state on disk)

| Image | n_3d_total | n_3d_keep | n_0_nuc | median_kept_vox | median_rej_vox | n_2d_MIP_keep |
|---|---|---|---|---|---|---|
| arsenite_rep1 | 48 | 14 | 34 | 201,726 | 481 | 11 |
| arsenite_rep2 | 17 | 11 | 6 | 271,956 | 405 | 7 |
| **arsenite_rep3** | **47** | **17** | **30** | **54,528** | 2,425 | 10 |
| control_rep1 | 15 | 13 | 2 | 293,697 | 176,187 | 14 |
| control_rep2 | 21 | 19 | 2 | 138,744 | 31,726 | — |
| control_rep3 | 11 | 6 | 5 | 166,908 | 18,513 | — |

`baseline.csv` (preserved on disk). Real cells across 5 correctly-segmented images: min volume 39,636 vox; p5 52,848; median 213,987. Rejected debris p75 across all 6: 3,087.

---

## B1 — `stitch_threshold` sweep on arsenite_rep3

| stitch_threshold | n_total | n_keep | n_0nuc | kept_vol_med | rej_vol_med |
|---|---|---|---|---|---|
| 0.1 | 45 | 17 | 28 | 59,130 | 2,426 |
| 0.2 | 46 | 17 | 29 | 54,528 | 2,511 |
| 0.3 | 46 | 17 | 29 | 54,528 | 2,511 |
| **0.4 (baseline)** | 47 | 17 | 30 | 54,528 | 2,426 |

**Verdict: no meaningful effect.** Going 0.4 → 0.1 drops 2 objects, kept count flat. The fragmentation is not driven by per-Z IoU linkage — Cellpose is finding different mask outlines in each Z slice, and lower stitch thresholds don't merge those across slices because they're not slight variations of the same object. **Hypothesis #2 ruled out.**

CSV: `stitch_sweep_arsenite_rep3.csv`.

---

## B2 — `min_cell_volume_vox` floor

### Floor selection from histograms

The data shows a clean gap between real cells and fragments/debris (computed across all 6 images):

- **Real cells (5 correctly-segmented images, n=63):** min 39,636; p5 52,848; p10 89,026; p25 140,202; median 213,987
- **Rejected debris (all 6):** p25 ~500; median ~3,000; p75 3,087
- **arsenite_rep3 "kept" fragmented cells:** min 5,859; p10 27,040; median 54,528

A floor at **30,000 vox** sits comfortably below any real cell (p5 ≈ 53k) and well above the debris/fragment band (p75 ≈ 3k).

### Floor=30,000 on all 6 images

| Image | n_keep base | n_keep @30k | Δkeep | n_total base | n_total @30k | n_0nuc base | n_0nuc @30k |
|---|---|---|---|---|---|---|---|
| arsenite_rep1 | 14 | 14 | 0 | 48 | 15 | 34 | 1 |
| arsenite_rep2 | 11 | 11 | 0 | 17 | 11 | 6 | 0 |
| **arsenite_rep3** | 17 | **15** | **−2** | 47 | 16 | 30 | 1 |
| control_rep1 | 13 | 13 | 0 | 15 | 15 | 2 | 2 |
| control_rep2 | 19 | 19 | 0 | 21 | 20 | 2 | 1 |
| control_rep3 | 6 | 6 | 0 | 11 | 6 | 5 | 0 |

**Verdict: locked in as new mammalian-preset default.**

- Zero regression on the 4 good images and control_rep3 (n_keep identical to baseline)
- arsenite_rep3 loses 2 fragments (17 → 15); the lost fragments are the smallest "kept" cells (min 5,859 vox, p10 27k vox), well below the 39k minimum of real cells in the other images
- Debris cleanup is dramatic across all images (n_0nuc 34→1, 6→0, 30→1, 5→0)

CSV: `volume_floor_arsenite_rep3.csv`, `volume_floor_all6.csv`.

---

## B3 — cell-seg channel: composite vs cytoplasm-only

### arsenite_rep3 probe

| label | cell_seg_channel | min_vol_vox | n_total | n_keep | n_0nuc | median_kept_vol |
|---|---|---|---|---|---|---|
| baseline | composite | 0 | 47 | 17 | 30 | 54,528 |
| B2 | composite | 30,000 | 16 | 15 | 1 | 59,130 |
| **B3** | **PABPC1** | 0 | 21 | 17 | 4 | **248,535** |
| B3 | PABPC1 | 30,000 | 20 | 17 | 3 | **248,535** |
| B3 | G3BP1 | 30,000 | 18 | 16 | 2 | 258,096 |

**Striking finding:** using PABPC1 (cytoplasmic marker) as cell-seg input instead of the composite (sum of DAPI + G3BP1 + PABPC1) lifts the kept-cell volume **4.5×** for arsenite_rep3 (54k → 248k). The composite is dominated by DAPI in the tight cluster, blowing out cell boundaries; the cytoplasmic-only channel preserves them.

### B3 all-6 validation (PABPC1 + floor=30k)

| Image | n_keep base | n_keep B3 | Δkeep | kept_vol base | kept_vol B3 |
|---|---|---|---|---|---|
| arsenite_rep1 | 14 | 14 | 0 | 201,726 | 233,734 |
| arsenite_rep2 | 11 | 11 | 0 | 271,956 | 278,484 |
| **arsenite_rep3** | 17 | **17** | **0** | 54,528 | **248,535** |
| control_rep1 | 13 | **14** | +1 | 293,697 | 282,212 |
| **control_rep2** | 19 | **18** | **−1** | 138,744 | 141,423 |
| control_rep3 | 6 | 6 | 0 | 166,908 | 313,524 |

**Verdict: flagged for Asif's review; NOT committed as default.**

- arsenite_rep3 fragmentation fully resolved (kept_vol 55k → 249k, matching other arsenite images)
- 4 of 5 other images: unchanged or improved
- **control_rep2 lost 1 cell** (19 → 18) — a regression that fails the "no degradation of good images" bar for a default change

**Practical guidance for users with tight HCT116-style clusters:** add `--cell-seg-channel PABPC1` (or your cytoplasmic marker of choice) to the cellquant invocation. The CLI flag is already exposed; no code change is needed. The comment in the mammalian `_3d` preset now documents this workaround.

**For Asif:** if you can find a way to default this without the control_rep2 regression — e.g. "prefer the brightest cytoplasmic-style channel" with an automatic heuristic — that would be a clean follow-up. The fragmentation fix is significant enough that the workaround is worth surfacing in documentation regardless.

CSV: `cellseg_input_arsenite_rep3.csv`, `cellseg_input_all6.csv`.

---

## C — nuclei_diameter on control_rep1

| nuclei_diameter | n_total | n_keep | n_0nuc | median_kept_vol | median_rej_vol |
|---|---|---|---|---|---|
| auto (default) | 15 | 13 | 2 | 293,697 | 176,187 |
| 100 | 15 | 13 | 2 | 293,697 | 176,187 |
| 130 | 15 | 13 | 2 | 293,697 | 176,187 |

**Verdict: not a fix.** Setting explicit nuclei_diameter values (100, 130 — both reasonable for HCT116) produces bit-identical output to the auto default. The 2 missing nuclei in control_rep1 are below the Cellpose-SAM detection threshold or in some configuration the model can't recover via diameter alone.

**Document as a known minor limitation** — under low-DAPI conditions, the nuclei gate can drop occasional real cells. Mitigations a user could apply (not changing defaults):

- Lower `--keep-min-nuclei` to 0 (accepts cells without segmented nuclei — but loses the debris filter)
- Manually patch the nuclear mask in those edge cases

CSV: `nuclei_diam_control_rep1.csv`.

---

## Committed change

`cellquant.py` mammalian `_3d` preset:

```python
"_3d": {
    ...
    # min_cell_volume_vox = 30_000 (was 0)
    # Justified by histogram analysis: real cell p5 ≈ 53k; debris p75 ≈ 3k.
    # 30k sits in the gap. Cleans the per-slice Cellpose fragmentation debris
    # in tight HCT116 clusters without removing any real cells from the
    # validated good images.
    "min_cell_volume_vox": 30_000,
    ...
}
```

Pre-change backup: `cellquant.py.pre-seg-floor.bak`. pytest 10/10 passing post-change.

---

## Final all-6 results with new default

Re-running `01_run_pipelines.py --dataset mammalian --skip-published` after the preset change reproduces B2 numbers exactly:

| Image | n_total | n_keep | n_0_nuc | kept_vol_median |
|---|---|---|---|---|
| arsenite_rep1 | 15 | 14 | 1 | 201,726 |
| arsenite_rep2 | 11 | 11 | 0 | 271,956 |
| arsenite_rep3 | 16 | 15 | 1 | 59,130 |
| control_rep1 | 15 | 13 | 2 | 293,697 |
| control_rep2 | 20 | 19 | 1 | 138,744 |
| control_rep3 | 6 | 6 | 0 | 166,908 |

Comparison report regenerated; HTML report at `validation_3d/REPORT.html` (4.7 MB).

---

## Honest summary

- **The mammalian 3D-segmentation correctness improvement is real but partial.** The volume floor cleans up debris and brings cell counts closer to the 2D-MIP truth without harming any image. But it does not fix the underlying per-slice Cellpose fragmentation in tight clusters — arsenite_rep3's kept cells are still half the volume of their counterparts in other images.

- **The deeper fix (cytoplasm-only cell-seg input) was probed and works** (4.5× kept-volume recovery on arsenite_rep3), but it costs one cell in control_rep2, so it does not meet the bar for an automatic default. It IS available to users now via `--cell-seg-channel PABPC1`, documented in the preset comment.

- **The minor nuclei-gating failure on control_rep1 is unfixable by parameter tuning** — those 2 cells stay rejected unless the user explicitly lowers the nuclei gate.

- **No regression** to the 2D path (10/10 regression tests pass).

- **All knobs remain exposed** as CLI flags: `--min-cell-volume-vox`, `--max-cell-volume-vox`, `--stitch-threshold`, `--cell-seg-channel`, `--nuclei-diameter`.

## Files

- `baseline.csv` — pre-change state
- `stitch_sweep_arsenite_rep3.csv` — B1
- `volume_floor_arsenite_rep3.csv`, `volume_floor_all6.csv` — B2
- `cellseg_input_arsenite_rep3.csv`, `cellseg_input_all6.csv` — B3
- `nuclei_diam_control_rep1.csv` — C
- This document (`RESULTS.md`)
- `../../REPORT.html` — full validation report with the new defaults
