"""B2 validation: test both floor=15k and floor=30k on all 6 images."""
from __future__ import annotations

import importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
runner = importlib.import_module("seg_experiments_runner")
setup_paths = importlib.import_module("00_setup_paths")

import pandas as pd

IMAGES = ['arsenite_rep1','arsenite_rep2','arsenite_rep3',
          'control_rep1','control_rep2','control_rep3']
FLOORS = [15_000, 30_000]

records = []
print("=== B2 all-6 validation ===")
for floor in FLOORS:
    label = f"B2_floor{floor}"
    print(f"\n=== floor = {floor} ===")
    for stem in IMAGES:
        # arsenite_rep3 already done in seg_exp_B2.py; skip
        out = setup_paths.VALIDATION_ROOT / "diagnostics" / "seg_experiments" / label / stem
        if (out / "cells.csv").exists():
            print(f"  {stem}: already on disk, summarizing")
            s = runner.summarize_cells(label, stem)
        else:
            print(f"  {stem}: running")
            r = runner.run(label, stem, min_cell_volume_vox=floor)
            print(f"    rc={r['returncode']} dt={r['duration_sec']:.0f}s")
            s = runner.summarize_cells(label, stem)
        records.append({**s, "min_cell_volume_vox": floor})

# Baseline from existing runs
print("\n=== baseline (floor=0, existing on-disk) ===")
for stem in IMAGES:
    cf = f"validation_3d/outputs_mammalian/3d/{stem}/cells.csv"
    df = pd.read_csv(cf)
    kept = df[df["keep"]]; rej = df[~df["keep"]]
    records.append({
        "label": "baseline_floor0", "stem": stem,
        "n_total": int(len(df)), "n_keep": int(len(kept)),
        "n_0nuc": int((df["n_nuclei"]==0).sum()),
        "median_kept_vol": float(kept["cell_volume_vox"].median()) if len(kept) else float("nan"),
        "median_rej_vol": float(rej["cell_volume_vox"].median()) if len(rej) else float("nan"),
        "min_cell_volume_vox": 0,
    })

df = pd.DataFrame(records)
df.to_csv("validation_3d/diagnostics/seg_experiments/volume_floor_all6.csv", index=False)

# Pivot for comparison
pd.set_option('display.width', 200)
print()
print("=== Per-image n_keep by floor ===")
print(df.pivot(index="stem", columns="min_cell_volume_vox", values="n_keep").to_string())
print()
print("=== Per-image n_total by floor ===")
print(df.pivot(index="stem", columns="min_cell_volume_vox", values="n_total").to_string())
print()
print("=== Per-image n_0nuc by floor ===")
print(df.pivot(index="stem", columns="min_cell_volume_vox", values="n_0nuc").to_string())
print()
print("=== Per-image kept-cell median vol by floor ===")
print(df.pivot(index="stem", columns="min_cell_volume_vox", values="median_kept_vol").to_string())
