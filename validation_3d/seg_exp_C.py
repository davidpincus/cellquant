"""C: nuclei_diameter on control_rep1 — does an explicit value
recover the 2 missed nuclei without harming the other images?

For HCT116 at 0.094 µm/px, typical nucleus is ~12-15 µm = 130-160 px.
Test with diameter 100 (small estimate) and 130 (typical).
"""
from __future__ import annotations

import importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
runner = importlib.import_module("seg_experiments_runner")

import pandas as pd

DIAMS = [100, 130]

records = []
print("=== C: nuclei_diameter probe on control_rep1 ===")
for d in DIAMS:
    label = f"C_nucdiam{d}"
    print(f"\n--> {label} on control_rep1")
    r = runner.run(label, "control_rep1", nuclei_diameter=d)
    print(f"    rc={r['returncode']} dt={r['duration_sec']:.0f}s")
    s = runner.summarize_cells(label, "control_rep1")
    print(f"    total={s['n_total']} keep={s['n_keep']} 0nuc={s['n_0nuc']} "
          f"kept_vol={s['median_kept_vol']:.0f} rej_vol={s['median_rej_vol']:.0f}")
    records.append({**s, "nuclei_diameter": d})

# Baseline
base = pd.read_csv("validation_3d/outputs_mammalian/3d/control_rep1/cells.csv")
kept = base[base["keep"]]; rej = base[~base["keep"]]
records.append({
    "label": "baseline_nucdiam_auto", "stem": "control_rep1",
    "n_total": len(base), "n_keep": int(len(kept)),
    "n_0nuc": int((base["n_nuclei"]==0).sum()),
    "median_kept_vol": float(kept["cell_volume_vox"].median()),
    "median_rej_vol": float(rej["cell_volume_vox"].median()),
    "nuclei_diameter": "auto",
})

df = pd.DataFrame(records)
df.to_csv("validation_3d/diagnostics/seg_experiments/nuclei_diam_control_rep1.csv", index=False)
pd.set_option('display.width', 200)
print()
cols = ["label", "nuclei_diameter", "n_total", "n_keep", "n_0nuc",
        "median_kept_vol", "median_rej_vol"]
print(df[cols].to_string(index=False))
print()
print("Target: recover the 2 rejected real cells (vol 137-215k) by detecting their nuclei.")
print("If n_keep rises from 13 -> 15 with reasonable diameter, that's the fix.")
