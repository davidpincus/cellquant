# PCA_3D_YEAST — yeast multivariate temperature phenotype, rebuilt on 3D metrics

## Published recipe (pinned)
- Source: `AVCGTIA_manuscript/workups/pca_umap_temperature.py` on `Asif_Tif6_Nsr1_Sis1_cellquant.csv` (MIP run).
- Samples: **per-cell**, keep==True (9355 cells = published S2 totals 2913+1603+879+2836+1124).
- 39 features = all numeric cols except {image,condition,replicate,cell_id,n_nuclei,keep} with >50% non-NaN.
- Preprocessing: fillna(0) → StandardScaler (z-score, unit variance). **No log, no other imputation.**
- PCA: `sklearn.decomposition.PCA()`, all components (SVD). Plot **PC1 vs PC2**, color = temperature, cmap=coolwarm, order 25/30/32/36/40C.
- Existing `pca_*.csv` are the **PUBLISHED MIP** PCA (not a 3D attempt).
- (Methods text says ">10% missing excluded" + UMAP n_neighbors=15/min_dist=0.1; the committed script uses >50% non-NaN + n_neighbors=30/0.3. Followed the **script** that produced the published 39 features.)

## Feature mapping (39 published → 37 3D)
30 map 1:1 (all coloc pearson/manders, all channel intensities, condensate indices, n_nucleoli).
7 are 2D→3D analogs: nucleolar_circularity→**sphericity**; cell/cytosol/nucleolar/Sis1puncta/Tif6puncta `_area_px`→`_volume_vox`; nucleus_area_px→nucleus_volume_vox (degenerate, no nucleus channel, ~0 in both).
2 dropped (no 3D analog; pipeline computes neither): **nucleolar_solidity, nucleolar_eccentricity**.
→ 37 comparable features, 7943 cells (keep==True).

## STEP 4a — does it hold up
- **Same features dominate PC1** (3D): Sis1 & Tif6 condensate indices (cell+cytosol) + channel intensities + pearson_Sis1_Nsr1 — matches published PC1 ("driven primarily by Sis1/Tif6 condensate indices and cytoplasmic variants"). ✓
- **Variance**: 3D PC1=36.8% / PC2=12.6% vs published 33.1% / 15.1%; cum-5 70% vs 67%. ✓ close.
- **Temperature trajectory preserved.** In BOTH runs the temperature gradient maps onto **PC2**, not PC1: corr(temp-index, PC2)=**+0.91 (3D)** vs **+0.83 (published)**; corr with PC1 is weak in both (+0.43 / −0.31, sign arbitrary). Per-temp PC2 centroid trajectory matches (25/30C low, 32C mid, 36C rising, 40C extreme outlier), with the largest 36→40C step in both. Centroid-separation ratio (between-centroid / within-spread) = **0.72 (3D)** vs **0.68 (published)** — 3D separation is marginally *stronger*. (The manuscript's "progressively along PC1" wording is loose: the published per-cell scores themselves are non-monotonic on PC1 and carry temperature on PC2.)

## STEP 4b — robustness (drop fragile puncta/condensate features)
Dropped 14 fragile features (Sis1/Tif6/Nsr1 condensate_index cell+cytosol, Sis1/Tif6 puncta_n, frac_intensity_in_puncta, puncta_integrated_intensity, puncta_volume_vox) → 23 robust features (intensities, coloc, nucleolar morphology, n_nucleoli).
- Robust PC1=40.3%, PC2=16.2%.
- corr(temp, PC2)=**+0.92** (≈ full +0.91); PC2 trajectory near-identical.
- Separation ratio **rises to 0.82** (from 0.72 full).
→ **Temperature separation is driven by the robust features, not the fragile puncta ones.** Removing puncta-derived features *improves* separation, so the multivariate temperature phenotype is not an artifact of fragile puncta detection.

## UMAP (published panel B)
Added on the same standardized 37-feature matrix, published params (`umap-learn`, n_neighbors=30, min_dist=0.3, random_state=42). Merge + per-temperature small multiples. 3D UMAP recovers the published topology: 25/30 °C cluster together, progressive shift through 32/36 °C, 40 °C separated — each temperature a distinct contiguous region, adjacent temps overlapping, extremes largely non-overlapping (matches manuscript Fig 4B narrative).

## Outputs
- `validation_3d/pca_3d_yeast.png` (A: PC1×PC2 scatter + centroid trajectory; B: UMAP merge + per-temp; C: top PC1 loadings; D: PC1–5 loadings heatmap — mirrors published Fig 4A,B,C,D)
- `validation_3d/pca_3d_loadings.csv`, `pca_3d_variance_explained.csv`, `pca_3d_scores.csv` (now incl. UMAP1/UMAP2)
- script: `validation_3d/09_pca_3d_yeast.py` (needs sklearn+umap env, e.g. `miniforge3/envs/bm_seq`; umap-learn installed via conda-forge)

## Rebuttal sentence
Re-deriving the 39-feature single-cell PCA on the native 3D metrics (37 comparable features) recovers the published multivariate temperature phenotype — the same condensate-index/intensity loadings on PC1 and a continuous temperature trajectory along PC2 culminating in a separated 40 °C state — with marginally stronger temperature separation than the MIP analysis, and the separation survives (indeed strengthens) when all fragile puncta-derived features are removed.
