"""validation_3d/15_figureS1_hct116_supp.py

Figure S1B -- the HCT116 metrics cellquant emits that main Figure 2B does not
show. Figure 2B carries G3BP1 puncta count, G3BP1 puncta area, and fraction
condensed for both channels; this fills in PABPC1's puncta measurements, both
condensate indices, and both fragmentation indices.

Everything -- panel style, point sizes, symlog rules, effect-size annotation --
is imported from 13_figure2_mammalian_2d.py rather than reimplemented, so the
supplement cannot drift from the main figure.

Metric definitions, for the legend:
  condensate index          p95 / mean of the per-cell intensity distribution
                            (cellquant.py:1116). > 1 means non-uniform signal.
  fragmentation persistence threshold-free swept integral on [0, 1]; higher
                            means the signal stays broken into more pieces
                            across the whole intensity sweep, so it does not
                            depend on picking one threshold.

Usage:
    python 15_figureS1_hct116_supp.py [--outputs-dir DIR] [--outdir DIR]
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "fig2", HERE / "13_figure2_mammalian_2d.py")
fig2 = importlib.util.module_from_spec(spec); sys.modules["fig2"] = fig2
spec.loader.exec_module(fig2)

# (column, panel title, y label)
PANELS = [
    ("PABPC1_puncta_n",                       "PABPC1 puncta per cell",
     "puncta / cell"),
    ("PABPC1_puncta_area_px",                 "PABPC1 puncta area",
     "area (px)"),
    ("G3BP1_condensate_index_cell",           "G3BP1 condensate index",
     "p95 / mean"),
    ("PABPC1_condensate_index_cell",          "PABPC1 condensate index",
     "p95 / mean"),
    ("G3BP1_fragmentation_index_persistence", "G3BP1 fragmentation",
     "persistence index"),
    ("PABPC1_fragmentation_index_persistence","PABPC1 fragmentation",
     "persistence index"),
]

# PABPC1 puncta area spans the same ~100x range as its G3BP1 twin in Fig 2B, so
# it gets the same symlog treatment. Register it on the imported module's table
# so draw() picks it up. The other five are bounded quantities and stay linear.
fig2.YSCALE["PABPC1_puncta_area_px"] = ("symlog", 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path,
                    default=(HERE.parent / "jcb_revision_midway3_cellquant_analysis"
                             / "validation_3d" / "outputs_mammalian" / "published_2d_v1.1"))
    ap.add_argument("--outdir", type=Path, default=HERE / "figure2_2d")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = fig2.load_cells(args.outputs_dir)
    cells["condition"] = cells["condition"].astype(str)
    reps = {c: sorted(cells.loc[cells["condition"] == c, "rep_id"].unique())
            for c in fig2.ORDER}
    n1, n2 = len(reps["control"]), len(reps["arsenite"])
    print("replicates:", {c: len(v) for c, v in reps.items()})

    fig = plt.figure(figsize=(11.5, 7.4))
    gs = GridSpec(2, 3, figure=fig, wspace=0.42, hspace=0.52)

    rows = []
    for i, (col, title, ylab) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        if col not in cells.columns:
            ax.text(0.5, 0.5, f"{col}\nabsent", ha="center", va="center",
                    fontsize=fig2.FONT_SIZE, transform=ax.transAxes)
            ax.set_axis_off(); print(f"  !! {col} absent"); continue
        rm = fig2.rep_medians(cells, col)
        p, mc, ma = fig2.test(rm)
        es = fig2.effect_size(rm)
        fig2.draw(ax, cells, rm, col, title, ylab, reps, es["fold"])
        print(f"  {title:28s} control={mc:.4g}  arsenite={ma:.4g}  "
              f"fold={es['fold']:.2f}  delta={es['cliffs_delta']:+.2f}  (p={p:.4f})")
        rows.append({"metric": col, "panel": title, "control_median": mc,
                     "arsenite_median": ma, "fold_change": es["fold"],
                     "cliffs_delta": es["cliffs_delta"], "p_exact_two_sided": p,
                     "n_control": n1, "n_arsenite": n2})

    stem = args.outdir / "FigureS1B_hct116_supp"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix("." + ext), bbox_inches="tight", dpi=300)
    print(f"\n  wrote {stem.name}.svg / .pdf / .png")
    pd.DataFrame(rows).to_csv(args.outdir / "figureS1b_stats.csv", index=False)
    print("  wrote figureS1b_stats.csv")


if __name__ == "__main__":
    np.random.seed(0)
    main()
