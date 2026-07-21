"""B1: stitch_threshold sweep on arsenite_rep3."""
from __future__ import annotations

import importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
runner = importlib.import_module("seg_experiments_runner")
setup_paths = importlib.import_module("00_setup_paths")

import pandas as pd

THRESHOLDS = [0.1, 0.2, 0.3]  # 0.4 baseline already on disk

records = []
print("=== B1: stitch_threshold sweep on arsenite_rep3 ===")
for thr in THRESHOLDS:
    label = f"B1_stitch{thr}"
    print(f"\n--> {label} on arsenite_rep3")
    r = runner.run(label, "arsenite_rep3", stitch_threshold=thr)
    print(f"    rc={r['returncode']} dt={r['duration_sec']:.0f}s")
    s = runner.summarize_cells(label, "arsenite_rep3")
    print(f"    total={s['n_total']} keep={s['n_keep']} 0nuc={s['n_0nuc']} "
          f"kept_vol={s['median_kept_vol']:.0f} rej_vol={s['median_rej_vol']:.0f}")
    records.append({**s, "stitch_threshold": thr,
                    "duration_sec": r["duration_sec"]})

# Add the baseline (0.4) from the existing run
import pandas as pd
base_df = pd.read_csv("validation_3d/outputs_mammalian/3d/arsenite_rep3/cells.csv")
kept = base_df[base_df["keep"]]; rej = base_df[~base_df["keep"]]
records.append({"label": "baseline_stitch0.4", "stem": "arsenite_rep3",
                "n_total": len(base_df), "n_keep": int(kept.shape[0]),
                "n_0nuc": int((base_df["n_nuclei"]==0).sum()),
                "median_kept_vol": float(kept["cell_volume_vox"].median()),
                "median_rej_vol": float(rej["cell_volume_vox"].median()),
                "stitch_threshold": 0.4, "duration_sec": float("nan")})

df = pd.DataFrame(records).sort_values("stitch_threshold").reset_index(drop=True)
df.to_csv("validation_3d/diagnostics/seg_experiments/stitch_sweep_arsenite_rep3.csv", index=False)
pd.set_option('display.width', 160)
print()
print(df[["stitch_threshold", "n_total", "n_keep", "n_0nuc",
          "median_kept_vol", "median_rej_vol", "duration_sec"]].to_string(index=False))
