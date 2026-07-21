"""B3 probe: cytoplasm-only cell-seg input on arsenite_rep3.
Tests whether using PABPC1 (cytoplasmic) instead of composite reduces
the per-slice fragmentation in the tight HCT116 cluster.

Tests with floor=0 and floor=30k.
"""
from __future__ import annotations

import importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
runner = importlib.import_module("seg_experiments_runner")

import pandas as pd

CONFIGS = [
    ("B3_PABPC1_floor0",     {"cell_seg_channel": "PABPC1"}),
    ("B3_PABPC1_floor30k",   {"cell_seg_channel": "PABPC1", "min_cell_volume_vox": 30_000}),
    ("B3_G3BP1_floor30k",    {"cell_seg_channel": "G3BP1",  "min_cell_volume_vox": 30_000}),
]

records = []
print("=== B3 probe: cytoplasm-channel cell-seg on arsenite_rep3 ===")
for label, kwargs in CONFIGS:
    print(f"\n--> {label}")
    r = runner.run(label, "arsenite_rep3", **kwargs)
    print(f"    rc={r['returncode']} dt={r['duration_sec']:.0f}s")
    if r["returncode"] != 0:
        print(f"    stderr: {r['stderr_tail'][:200]}")
    s = runner.summarize_cells(label, "arsenite_rep3")
    print(f"    total={s['n_total']} keep={s['n_keep']} 0nuc={s['n_0nuc']} "
          f"kept_vol={s['median_kept_vol']:.0f} rej_vol={s['median_rej_vol']:.0f}")
    records.append({**s, **kwargs})

# Add baseline (composite, no floor) and best B2 (composite, floor=30k)
base = pd.read_csv("validation_3d/outputs_mammalian/3d/arsenite_rep3/cells.csv")
kept = base[base["keep"]]; rej = base[~base["keep"]]
records.append({
    "label": "baseline_composite_floor0", "stem": "arsenite_rep3",
    "n_total": len(base), "n_keep": int(len(kept)),
    "n_0nuc": int((base["n_nuclei"]==0).sum()),
    "median_kept_vol": float(kept["cell_volume_vox"].median()),
    "median_rej_vol": float(rej["cell_volume_vox"].median()),
    "cell_seg_channel": "composite", "min_cell_volume_vox": 0,
})
b2 = pd.read_csv("validation_3d/diagnostics/seg_experiments/B2_floor30000/arsenite_rep3/cells.csv")
kept2 = b2[b2["keep"]]; rej2 = b2[~b2["keep"]]
records.append({
    "label": "B2_composite_floor30k", "stem": "arsenite_rep3",
    "n_total": len(b2), "n_keep": int(len(kept2)),
    "n_0nuc": int((b2["n_nuclei"]==0).sum()),
    "median_kept_vol": float(kept2["cell_volume_vox"].median()),
    "median_rej_vol": float(rej2["cell_volume_vox"].median()) if len(rej2) else float("nan"),
    "cell_seg_channel": "composite", "min_cell_volume_vox": 30_000,
})

df = pd.DataFrame(records)
df.to_csv("validation_3d/diagnostics/seg_experiments/cellseg_input_arsenite_rep3.csv", index=False)
pd.set_option('display.width', 200)
print()
cols = ["label", "cell_seg_channel", "min_cell_volume_vox",
        "n_total", "n_keep", "n_0nuc", "median_kept_vol", "median_rej_vol"]
print(df[cols].to_string(index=False))
