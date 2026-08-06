# 3D validation harness

Scripts that reproduce the 2D-vs-3D validation analysis for cellquant. Run
outputs, logs, and large image data are gitignored; only the
`.py`/`.md` code and this README are tracked as the reproducibility record.

## Configuring the data root

Paths are resolved portably in `00_setup_paths.py` — there is no hardcoded
machine path, so the identical harness runs on the Mac and on Midway3. The
project root (the directory containing `cellquant.py`, `validation_3d/`, and the
data dirs `SG_zstacks/` and `Tif6_Nsr1_Sis1_6hr/`) is resolved in this order:

1. `$CELLQUANT_VALIDATION_ROOT` — explicit override, e.g.
   `export CELLQUANT_VALIDATION_ROOT=/scratch/midway3/pincus/cellquant_run`
2. `validation_3d/paths.local.json` — copy `paths.local.json.example` and set
   `project_root` (this file is gitignored).
3. Repo-relative default — the parent of `validation_3d/`. This already works on
   both machines when the data sits alongside the harness, so most setups need
   no configuration at all.

Resolution fails loudly (`FileNotFoundError`) if the resolved root does not
exist or is missing `cellquant.py` / `validation_3d/`, naming the resolved path
and the env var to set. Missing image data is reported separately, and loudly,
when a dataset is actually listed (`_glob_required`).
