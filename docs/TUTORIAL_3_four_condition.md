# Tutorial 3: Four-condition arsenite dose-response

This tutorial covers a multi-condition experiment with explicit positive and negative controls. The structure here generalizes to most dose-response or panel comparisons — drug concentrations, time points, mutant alleles, temperatures — that you'd display as a four-panel superplot with statistics anchored to a chosen reference.

> **This tutorial uses your own data.** Unlike [Tutorial 1](TUTORIAL_1_mammalian_SGs.md) and [Tutorial 2](TUTORIAL_2_yeast_temp.md), no `arsenite_doseresponse/` dataset ships with cellquant. The commands below are a template to adapt to a four-condition experiment of your own, and the numbers shown are illustrative rather than reproducible.

> **Reminder.** cellquant is not a substitute for understanding what the analysis is doing. The numbers it produces are only meaningful if the segmentation is good, the metric matches your biology, and the conditions are comparable. If you haven't read [CONCEPTS.md](CONCEPTS.md), skim the puncta-detection section before trusting the puncta counts you'll see in this tutorial.

## The experiment

We're imaging U2OS cells expressing GFP-G3BP1 (stress granule marker) and stained with PABPC1 (cytoplasmic mRNA marker) and DAPI (nucleus). Cells are treated with sodium arsenite at four concentrations:

| Condition | Arsenite (µM) | Role | Expected G3BP1 phenotype |
|-----------|---------------|------|--------------------------|
| `mock` | 0 | Negative control | Diffuse cytosolic GFP, no granules |
| `low` | 50 | Sub-threshold dose | Few granules, mixed response |
| `medium` | 100 | Standard SG-inducing dose | Robust granules in most cells |
| `high` | 500 | Positive control | Maximal granulation |

Three biological replicates per condition, three fields per replicate. So 4 conditions × 3 replicates × 3 fields = 36 images total.

## File organization

Put all 36 TIFFs (max-intensity projections of z-stacks) in one folder, named with condition and replicate so cellquant can parse them:

```
arsenite_doseresponse/
├── MAX_mock_rep1_field1.tif
├── MAX_mock_rep1_field2.tif
├── MAX_mock_rep1_field3.tif
├── MAX_mock_rep2_field1.tif
├── ...
├── MAX_low_rep1_field1.tif
├── ...
├── MAX_high_rep3_field3.tif
```

Each file is a 3-channel TIFF: channel 1 = DAPI, channel 2 = G3BP1, channel 3 = PABPC1.

cellquant doesn't care about `_field{N}` — that suffix just keeps your filenames unique. The pattern `MAX_{condition}_rep{replicate}` matches the part it cares about.

## The command

```bash
uv run cellquant.py arsenite_doseresponse/ \
    "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
    --out arsenite_doseresponse/output/ \
    --cell-type mammalian \
    --filename-pattern "MAX_{condition}_rep{replicate}" \
    --condition-order mock low medium high \
    --reference-condition mock \
    --colocalization \
    --allow-2d-colocalization \
    --puncta-channels G3BP1 PABPC1
```

Two flags do most of the multi-condition work:

- **`--condition-order mock low medium high`** sets the left-to-right order of bars in the superplots. Without this, conditions appear in alphabetical order, which here would put `high` before `low` — confusing.
- **`--reference-condition mock`** sets the comparison baseline for the Mann-Whitney U p-values. Each other condition (`low`, `medium`, `high`) is tested against `mock`. Bonferroni correction is applied automatically for the three comparisons.

## What you'll see

After the run completes:

```
arsenite_doseresponse/output/
├── cells.csv                # one row per cell, all metrics
├── images.csv               # one row per image (n_cells, n_keep, medians)
├── colocalization.csv       # one row per cell per channel pair
├── pvalues.csv              # all pairwise Mann-Whitney U tests
├── provenance.json          # version, CLI, input checksums, deps
├── config_used.yml
├── qc/                      # 36 PNG overlays — check these first
├── masks/                   # cell/nuc/puncta TIFFs per image
├── plots/                   # superplots per metric
└── prism/                   # Prism-ready CSVs
```

### Step 1: validate segmentation

Open `qc/` and look at one image from each condition. You want cyan cell outlines around real cells, yellow outlines around nuclei, and magenta outlines around granules in the `high` condition. If `mock` cells have spurious magenta blobs, your puncta threshold is too permissive — tune `--log-sigma` or `--puncta-threshold-method` and re-run.

If `high`-condition cells look segmented as huge merged blobs (granules can confuse Cellpose into treating the cell boundary as the granule boundary), you may need to switch to a composite cell segmentation. With DAPI declared as `nucleus`, cellquant already builds a composite from all non-skip channels for cell segmentation — this is usually sufficient.

### Step 2: read the superplots

In `plots/`, look at `G3BP1_puncta_n_superplot.png`. With four conditions you'll get a strip-plot layout: each condition has a jittered scatter of cells (small dots), per-image median diamonds, and a Mann-Whitney p-value above each non-reference condition. The `mock` condition is the reference; p-values appear above `low`, `medium`, and `high`.

A clean dose-response looks like:

| Condition | Median G3BP1_puncta_n | p vs mock |
|-----------|-----------------------|-----------|
| mock | 0–1 | (reference) |
| low | 1–3 | p ≈ 0.1 — not significant |
| medium | 5–10 | p < 0.01 |
| high | 10–30 | p < 0.001 |

The non-significant `low` vs `mock` comparison is biologically meaningful (sub-threshold dose) — don't read that as "no effect," read it as "the dose curve has a knee somewhere between 50 and 100 µM." If you wanted to nail that down you'd add a `25` µM and `75` µM condition.

### Step 3: read `pvalues.csv`

```
metric,condition,reference,pval,pval_corrected,n_test_replicates,n_ref_replicates,n_comparisons
G3BP1_puncta_n,low,mock,0.0863,0.2589,3,3,3
G3BP1_puncta_n,medium,mock,0.0080,0.0240,3,3,3
G3BP1_puncta_n,high,mock,0.0049,0.0147,3,3,3
G3BP1_frac_intensity_in_puncta,low,mock,...
...
```

`pval_corrected` is the Bonferroni-corrected value (multiplied by `n_comparisons = 3`). Use `pval_corrected` when reporting — that's the family-wise error-rate-controlled number.

**Important.** With three replicates per condition, only large effect sizes reach significance even after correction. A non-significant `pval_corrected` in this dataset is not evidence of no effect — it's a power problem. See the "Replicate-level statistics" section of [CONCEPTS.md](CONCEPTS.md) for the underlying tradeoff.

### Step 4: cross-check with the colocalization data

G3BP1 and PABPC1 are both stress-granule markers; they should co-localize at high arsenite doses (both proteins inside the same granules) and be distributed independently in `mock` cells. In `cells.csv`:

| Condition | Median pearson_r_G3BP1_vs_PABPC1 |
|-----------|----------------------------------|
| mock | ~0.2 (background correlation) |
| low | ~0.3 |
| medium | ~0.5 |
| high | ~0.7 (strong colocalization) |

If `mock` cells show R ≈ 0.7, that's suspicious — likely bleed-through or a co-staining artifact, not biology. Acquisition-quality issues like this are invisible to cellquant; you have to check them.

The same caveat about MIP-based colocalization from [CONCEPTS.md](CONCEPTS.md#colocalization-opt-in-via---colocalization) applies. If you have the z-stacks (not just the MIPs), re-run on the source 3D data — cellquant will pick up the z-stacks automatically and run 3D Pearson/Manders natively, producing more reliable numbers.

## Variations

### Different reference condition

If your positive control is the natural baseline (e.g. comparing different rescue mutants against the wild-type-treated condition):

```bash
--reference-condition WT_arsenite \
--condition-order WT_arsenite mut1_arsenite mut2_arsenite mut3_arsenite
```

### More conditions, denser dose-response

cellquant doesn't cap the number of conditions. A 7-condition dose-response renders the same way; the figure just gets wider. Use `--trend` to add a dashed line through condition medians — useful for monotonic dose-responses:

```bash
--condition-order ctrl 10uM 25uM 50uM 100uM 250uM 500uM \
--reference-condition ctrl \
--trend
```

### Two-factor experiment (e.g. drug × time)

cellquant treats each unique `condition` string as a separate group, so encode two-factor designs in the filename:

```
MAX_DMSO_15min_rep1.tif
MAX_DMSO_30min_rep1.tif
MAX_DMSO_60min_rep1.tif
MAX_drug_15min_rep1.tif
MAX_drug_30min_rep1.tif
MAX_drug_60min_rep1.tif
```

Then `--filename-pattern "MAX_{condition}_rep{replicate}"` and `--condition-order DMSO_15min drug_15min DMSO_30min drug_30min DMSO_60min drug_60min` lays them out drug-vs-DMSO at each time point.

For more involved factorial designs (e.g. 3 drugs × 4 doses × 2 cell lines), pull the resulting `cells.csv` into Prism or R for a proper two-way analysis. cellquant's superplots are best for 1-D condition strings.

## Common mistakes

**Reference condition with fewer than 3 replicates.** The reference must have at least 3 replicates for p-values to be computed at all; cellquant silently skips conditions that fail this. Check `pvalues.csv` — if it's empty, you don't have enough replicates somewhere.

**Treating fields as replicates.** If your "3 replicates" are actually 3 fields from the same well, the statistical test is misleading — you have N=1 biological replicate, three technical fields. Acquire from different wells / days / mice for real replication.

**Forgetting `--reference-condition`** when you have 3+ conditions. With ≤2 conditions cellquant defaults to using the first as reference; with 3+ it doesn't guess. No reference → no p-values.

**Reading p < 0.05 as "the effect is real."** With three replicates and N family-wise corrected comparisons, p < 0.05 still represents a moderate amount of evidence. Report effect sizes (medians, fold changes) alongside the p-value and let readers calibrate.

## Next steps

- Validate replication by checking that the three replicates within each condition look similar — open the per-image diamond medians on the superplot. If one replicate sits far from the other two, look at that day's QC images for an acquisition issue.
- If your `mock` condition has too many false-positive puncta, tighten `--puncta-min-area-px` or switch to `--puncta-threshold-method triangle`. Check that the change doesn't suppress real granules in `high`.
- For very subtle dose-responses where `puncta_n` flickers near zero, switch to `{ch}_fragmentation_index_persistence` — it's a threshold-free continuous measure of clumpiness that doesn't depend on detecting individual puncta.

For background on what each metric means, see [CONCEPTS.md](CONCEPTS.md).
