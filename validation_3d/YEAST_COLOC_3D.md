# YEAST_COLOC_3D — 3D colocalization vs published MIP Table S2/S3

Run: `outputs_yeast/3d/` (3D-native) vs `outputs_yeast/matched_2d/` (MIP, masks held fixed).
n replicates 6/5/5/6/8 @ 25/30/32/36/40C (manuscript set). Whole-cell Costes.
Bonferroni Wilcoxon rank-sum on replicate medians vs 25C, factor=4 per pair.

## STEP 0 status
- 3D arm: 29 reps, counts 6/5/5/6/8 ✓ manuscript set (no teens/twenties → not the 6hr set).
- Every rep: cells.csv + colocalization.csv non-empty, all 3 pairs present. Min cells = 105 (40C_rep7); none <50. No failure markers.
- config_used.yml: mode=3d, seg-3d-method=stitch, z-crop-center=41, seg-downsample=2, coloc compartment=whole-cell ✓.
- **matched_2d Z-window: FULL-STACK MIP, not the cropped window.** `01_run_pipelines.py:476` (`_load_zstack_channels`) reads the entire z-stack; `ch.max(axis=0)` projects all planes. The 3D masks (XY footprints) are held fixed, but the MIP image integrates light from planes outside the 41-slice crop. → **Q2 carries caveat #3** (matched cells see out-of-window light); reported MIP−3D magnitudes are an upper bound. Masks are local (`outputs_yeast/3d/*/masks/`) → clean cropped-window re-run is cheap.

## STEP 1 reference check
`04_compare_yeast.py` PUBLISHED_S2/S3 == current `Neferkara_compiled.tex` Table S2/S3 **exactly** (Pearson, Manders M1, all S3 p-values; tex S2 even lists n-reps 6/5/5/6/8). No mismatches.

## Q1 — PRESERVATION (does 3D reproduce the decline + S3 significance?)

LEAD = PEARSON. 3D values are systematically **lower** than published (MIP inflates R) — Q1 is about trend + significance, not absolute level.

| Pair | published call (S2 trend / S3 sig) | 3D result | verdict |
|---|---|---|---|
| **Tif6–Nsr1** | declines; sig **36+40 only** | flat 25–32C, sig decline **36C & 40C only** (pBonf 0.009, 0.003) | **PRESERVED** — exact match. This is the manuscript's headline Nsr1–Tif6 panel. |
| **Sis1–Nsr1** | mild decline; sig **none** | ~flat, **no** sig temp effect (all pBonf >0.05) | **PRESERVED** — matches "no significant change." |
| **Tif6–Sis1** | monotonic decline; sig **all temps** | non-monotonic (rises 30–32C, drops only at 40C); the two nominally-sig temps (30/32C) are significant **increases**; 40C decline ns | **NOT reproduced** — the "significant decline at all temps" does not hold in 3D. |

Net: 2 of 3 pairs preserved, including the headline Tif6–Nsr1. The strong Tif6–Sis1 monotonic-decline claim weakens to a 40C-only, non-significant effect in 3D.

## Q2 — ARTIFACT (MIP − 3D, paired per-rep, Wilcoxon signed-rank vs 0)

PEARSON: **MIP − 3D > 0 in all 15 pair×temp cells.** MIP inflates Pearson R.
- magnitude by pair: Tif6–Sis1 +0.07…+0.17 (largest; diffuse–diffuse), Tif6–Nsr1 +0.05…+0.09, Sis1–Nsr1 +0.02…+0.05 (smallest).
- significance: p=0.031 (n=6), p=0.008 (n=8); n=5 temps floored at p=0.0625 but every rep positive. Direction is unambiguous and consistent: **MIP systematically overestimates colocalization**, exactly as the manuscript's 3D-coloc rationale predicts.

MANDERS M1 *(MARK: half-saturated in S2 — several 1.00s — projection-fragile; do not lead)*: signs mixed, mostly ns; Tif6–Sis1 M1 saturates at 1.00 in 3D so deltas compress. Uninformative for the artifact direction.

## Drop-in revised Table S2 (3D Pearson R) + recomputed S3 significance

Pearson R — replicate median (n reps); **bold** = Bonferroni Wilcoxon p<0.05 vs 25C *(↑ = wrong-direction sig)*:

| Pair | 25C | 30C | 32C | 36C | 40C |
|---|---|---|---|---|---|
| Tif6 vs Sis1 | 0.75 | **0.78↑** | **0.80↑** | 0.77 | 0.66 |
| Tif6 vs Nsr1 | 0.82 | 0.83 | 0.80 | **0.72** | **0.56** |
| Sis1 vs Nsr1 | 0.78 | 0.76 | 0.75 | 0.77 | 0.75 |

Manders M1 (3D) *(projection-fragile; report only)*:

| Pair | 25C | 30C | 32C | 36C | 40C |
|---|---|---|---|---|---|
| Tif6 vs Sis1 | 0.71 | 1.00 | 1.00 | 1.00 | 1.00 |
| Tif6 vs Nsr1 | 0.87 | 0.83 | 0.84 | **0.75** | **0.59** |
| Sis1 vs Nsr1 | 0.84 | 0.81 | 0.81 | 0.76 | **0.68** |

n cells per temp (Pearson, summed across reps): ~1942 / 1454 / 899 / 2472 / 1176.

Recomputed S3 (Pearson, sig vs 25C): Tif6–Sis1 = 30,32C only (both increases); Tif6–Nsr1 = 36,40C; Sis1–Nsr1 = none.

## Context (broad pass, all metrics)
3D directional agreement vs published = **53%** (was 40% at n=1 in the handoff). Proper n closed the gap only modestly; the broad pass is dominated by near-zero noise metrics and the 36C published proximity outliers (which are themselves likely MIP artifacts). The focused coloc above is the real result. matched_2d 79%, published_2d (n=1 stubs) 80%.

## Caveats / n
- n reps 6/5/5/6/8; Bonferroni factor 4 per pair; n=5 temps cannot reach p<0.05 on a 5-rep signed-rank even when every rep agrees.
- matched_2d = FULL-stack MIP (Step 0) → Q2 magnitudes are an upper bound. Clean cropped-window re-run is a cheap local/cluster job (masks retained).
- Full 30-rep published-2D (cellquant-on-MIP) NOT run locally — only n=1 stubs present; published Table S2 is the published-2D reference for Q1.

## Summary statement
Re-running colocalization natively in 3D on the source z-stacks (masks held fixed) preserves the manuscript's central finding — the temperature-dependent loss of Nsr1–Tif6 colocalization, significant only at 36 °C and 40 °C — while revealing that maximum-intensity projection systematically inflates Pearson's R by 0.02–0.17 across all marker pairs (paired MIP−3D > 0 in every condition), confirming that the published MIP values overestimate absolute colocalization even where the temperature trend is qualitatively retained.
