"""B2: min_cell_volume_vox floor sweep on arsenite_rep3."""
from __future__ import annotations

import importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
runner = importlib.import_module("seg_experiments_runner")

import pandas as pd

FLOORS = [15_000, 30_000]

records = []
print("=== B2: min_cell_volume_vox floor sweep on arsenite_rep3 ===")
for floor in FLOORS:
    label = f"B2_floor{floor}"
    print(f"\n--> {label} on arsenite_rep3")
    r = runner.run(label, "arsenite_rep3", min_cell_volume_vox=floor)
    print(f"    rc={r['returncode']} dt={r['duration_sec']:.0f}s")
    s = runner.summarize_cells(label, "arsenite_rep3")
    print(f"    total={s['n_total']} keep={s['n_keep']} 0nuc={s['n_0nuc']} "
          f"kept_vol={s['median_kept_vol']:.0f} rej_vol={s['median_rej_vol']:.0f}")
    records.append({**s, "min_cell_volume_vox": floor,
                    "duration_sec": r["duration_sec"]})

# baseline (floor=0) from existing run
base = pd.read_csv("validation_3d/outputs_mammalian/3d/arsenite_rep3/cells.csv")
kept = base[base["keep"]]; rej = base[~base["keep"]]
records.append({"label": "baseline_floor0", "stem": "arsenite_rep3",
                "n_total": len(base), "n_keep": len(kept),
                "n_0nuc": int((base["n_nuclei"]==0).sum()),
                "median_kept_vol": float(kept["cell_volume_vox"].median()),
                "median_rej_vol": float(rej["cell_volume_vox"].median()),
                "min_cell_volume_vox": 0, "duration_sec": float("nan")})

df = pd.DataFrame(records).sort_values("min_cell_volume_vox").reset_index(drop=True)
df.to_csv("validation_3d/diagnostics/seg_experiments/volume_floor_arsenite_rep3.csv", index=False)
pd.set_option('display.width', 160)
print()
print(df[["min_cell_volume_vox", "n_total", "n_keep", "n_0nuc",
          "median_kept_vol", "median_rej_vol", "duration_sec"]].to_string(index=False))
