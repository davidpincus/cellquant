"""B3 validation: PABPC1 + floor=30k on all 6 images."""
from __future__ import annotations

import importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
runner = importlib.import_module("seg_experiments_runner")
setup_paths = importlib.import_module("00_setup_paths")

import pandas as pd

IMAGES = ['arsenite_rep1','arsenite_rep2','arsenite_rep3',
          'control_rep1','control_rep2','control_rep3']

records = []
label = "B3_PABPC1_floor30k"
print(f"=== {label} on all 6 images ===")
for stem in IMAGES:
    out = setup_paths.VALIDATION_ROOT / "diagnostics" / "seg_experiments" / label / stem
    if (out / "cells.csv").exists():
        s = runner.summarize_cells(label, stem)
        print(f"  {stem}: cached → total={s['n_total']} keep={s['n_keep']} 0nuc={s['n_0nuc']}")
    else:
        print(f"  {stem}: running…")
        r = runner.run(label, stem, cell_seg_channel="PABPC1",
                       min_cell_volume_vox=30_000)
        print(f"    rc={r['returncode']} dt={r['duration_sec']:.0f}s")
        s = runner.summarize_cells(label, stem)
        print(f"    total={s['n_total']} keep={s['n_keep']} 0nuc={s['n_0nuc']} "
              f"kept_vol={s['median_kept_vol']:.0f}")
    records.append(s)

# Also pull baseline + B2_floor30k for comparison
for stem in IMAGES:
    base = pd.read_csv(f"validation_3d/outputs_mammalian/3d/{stem}/cells.csv")
    kept_b = base[base["keep"]]
    rej_b = base[~base["keep"]]
    records.append({
        "label": "baseline_composite_floor0", "stem": stem,
        "n_total": len(base), "n_keep": int(len(kept_b)),
        "n_0nuc": int((base["n_nuclei"]==0).sum()),
        "median_kept_vol": float(kept_b["cell_volume_vox"].median()) if len(kept_b) else float("nan"),
        "median_rej_vol": float(rej_b["cell_volume_vox"].median()) if len(rej_b) else float("nan"),
    })
    b2cf = f"validation_3d/diagnostics/seg_experiments/B2_floor30000/{stem}/cells.csv"
    if Path(b2cf).exists():
        b2df = pd.read_csv(b2cf)
        keptb2 = b2df[b2df["keep"]]; rejb2 = b2df[~b2df["keep"]]
        records.append({
            "label": "B2_composite_floor30k", "stem": stem,
            "n_total": len(b2df), "n_keep": int(len(keptb2)),
            "n_0nuc": int((b2df["n_nuclei"]==0).sum()),
            "median_kept_vol": float(keptb2["cell_volume_vox"].median()) if len(keptb2) else float("nan"),
            "median_rej_vol": float(rejb2["cell_volume_vox"].median()) if len(rejb2) else float("nan"),
        })

df = pd.DataFrame(records)
df.to_csv("validation_3d/diagnostics/seg_experiments/cellseg_input_all6.csv", index=False)
pd.set_option('display.width', 200)
print()
print("=== Per-image n_keep ===")
piv = df.pivot_table(index="stem", columns="label", values="n_keep", aggfunc="first")
print(piv.to_string())
print()
print("=== Per-image kept_vol median ===")
piv = df.pivot_table(index="stem", columns="label", values="median_kept_vol", aggfunc="first")
print(piv.to_string())
print()
print("=== Per-image n_total ===")
piv = df.pivot_table(index="stem", columns="label", values="n_total", aggfunc="first")
print(piv.to_string())
