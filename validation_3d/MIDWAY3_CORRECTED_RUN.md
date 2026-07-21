# Corrected yeast 3D re-run — Midway3 checklist

The local (Mac/CPU) run was stopped in favor of Midway3 (GPU). Everything below
reproduces the **corrected** yeast validation run on the cluster. The core fix
is the voxel size: the original run used `--voxel-size 0.094 0.1`, but the
files' OME PhysicalSize (authoritative) is **0.10571 µm XY / 0.23 µm Z**
(anisotropy 2.18).

Scope: 30 manuscript reps (25C=6, 30C=5, 32C=5, 36C=6, 40C=8), three axes
(3D, matched-2D-from-3D, published-2D).

---

## 1. Update `cellquant.py` (via git)

Branch `native-3d-pipeline` (tip contains: yeast stitch@0.65 default, voxel/
shape guardrails, 2D-coloc gate, nucleolar sphericity/surface removed,
proximity ±0.1 µm sensitivity band).

**On the Mac** (push the branch — tell me and I can do it, or run):
```
git push -u origin native-3d-pipeline
```

**On Midway3**, in the repo at `/scratch/midway3/pincus/cellquant_run/`:
```
git fetch origin && git checkout native-3d-pipeline
```
(If that path is not a git clone, just `scp` the Mac `cellquant.py` over it.)

---

## 2. Apply the harness corrections (gitignored — not carried by git)

These live in `validation_3d/` and must be transferred by hand. Apply the
edits; **keep Midway3's own `PROJECT_ROOT`** in `00_setup_paths.py`.

### `00_setup_paths.py` — voxel + GPU
```
YEAST_VOXEL_XY_UM = 0.10571      # was 0.094
YEAST_VOXEL_Z_UM  = 0.23         # was 0.1
```
And enable the GPU. Either set `NO_GPU = False`, or (recommended, so the same
file works on both machines) replace `NO_GPU = True` with auto-detect:
```python
def _detect_no_gpu() -> bool:
    try:
        import torch
        return not torch.cuda.is_available()
    except Exception:
        return True

NO_GPU = _detect_no_gpu()
```

### `01_run_pipelines.py` — un-break the published-2D pass (ESSENTIAL)
cellquant now aborts 2D colocalization unless you opt in. The published-2D
(MIP) pass is intentional — it generates the MIP-vs-3D artifact comparison — so
add the flag after `"--colocalization",` in **`run_published_2d_yeast`**:
```python
        "--colocalization",
        "--allow-2d-colocalization",   # <-- add
```
(Optional: `run_published_2d_mammalian` too, if you run mammalian; and the
`--resume` skip-if-cells.csv-exists support is nice for restarts but not
required on GPU.)

### `04_compare_yeast.py` — drop the removed metric
Nucleolar sphericity/surface were removed from the pipeline. Change the 3D
morphology list:
```python
    "3d": ["nucleolar_volume_um3", "nucleolar_eq_diameter_um"],   # was [..., "nucleolar_sphericity"]
```

---

## 3. Archive the old (wrong-voxel) outputs, then run

**On Midway3**, from `validation_3d/`:
```
mv outputs_yeast outputs_yeast_wrongvoxel_0.094_0.1_archive   # preserve, don't delete
python 01_run_pipelines.py --dataset yeast                    # add --resume for restart safety
python 04_compare_yeast.py                                    # aggregate vs published S2/S3
python 04_compare_yeast.py --coloc                            # focused 3D-vs-MIP coloc (Q1/Q2)
```

Sanity checks while it runs:
- Each 3D log should show `Voxel size: XY=0.1057 µm, Z=0.2300 µm (anisotropy=2.18)`.
- The corrected voxel now **agrees with the files' OME metadata**, so the
  cellquant voxel ladder proceeds (it would abort on the old 0.094/0.1).
- A `[warn] cell shape` line firing as *oblate* (low |cosZ|) would signal a
  voxel problem; the benign ~1.5–1.7× prolate-along-Z elongation is expected
  optical PSF and no longer warns.

---

## 4. Later (not needed for the run/compare above)

Downstream scripts still assume sphericity exists — reconcile before using them:
- `09_pca_3d_yeast.py` (line ~40 maps `nucleolar_circularity`→`nucleolar_sphericity`)
- `07_report.py` (the "Synthetic sphericity ground truth" section)
- `04_synthetic_sphericity.py` (now obsolete — tested the removed calculation)

Separate, out of scope for the yeast re-derivation: `MAMMALIAN_VOXEL_Z_UM = 0.25`
may have the same class of metadata error — worth checking the mammalian OME
before trusting mammalian 3D physical-unit metrics.
