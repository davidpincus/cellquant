"""validation_3d/16_figureS2_yeast_supp.py

Figure S2 -- the yeast temperature-series metrics cellquant emits that main
Figure 3B does not show. One page, 15 panels.

This REPLACES the published 20-panel S2 rather than updating it. Four of those
panels were nucleolar area, circularity, eccentricity and solidity, which are
properties of a projection and have no 3D counterpart, and several others were
duplicates that only looked distinct because the 2D run had a nucleus channel.

WHAT WAS LEFT OUT, AND WHY. Of 47 numeric columns not already in Figure 3B, 17
are redundant by construction and were never candidates:
  - *_condensate_index_cytosol == *_condensate_index_cell  (verified identical)
  - *_cytosol_mean == *_cell_mean, cytosol_volume == cell_volume  (identical:
    cytosol is cell-minus-nucleus, and this series has no nucleus channel)
  - nucleolar_eq_diameter_um is monotone in nucleolar_volume_um3 (Spearman 1.000)
  - *_volume_vox are the pre-calibration twins of *_volume_um3
  - fraction_proximal_0.4um / _0.6um are a sensitivity band around the 0.5 default
Of the 30 that remain, the six Manders coefficients are omitted because main
Figure 3B already shows colocalization for all three channel pairs; Manders M1
and M2 for every pair are in the deposited cells.csv. Tif6 puncta volume,
puncta/diffuse mean intensity and integrated intensity are likewise in the data
but are dominated by the counts and fractions shown here.

Panel style and statistics are imported from 11_figure3_yeast_3d.py so the
supplement cannot drift from the main figure.

Usage:
    python 16_figureS2_yeast_supp.py [--outputs-dir DIR] [--outdir DIR]
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
    "fig3", HERE / "11_figure3_yeast_3d.py")
fig3 = importlib.util.module_from_spec(spec); sys.modules["fig3"] = fig3
spec.loader.exec_module(fig3)

# 15 panels at ~2.3 x 2.0 in each: a 12 pt "p = 0.009" is taller than the default
# 0.075 stacking pitch, so four brackets in one panel overlap. Drop the bracket
# labels to 9 pt and open the pitch. Everything a reader actually reads -- titles,
# axis labels, ticks -- stays at 12 pt.
fig3.BRACKET_FONT_SIZE = 9
fig3.BRACKET_PITCH = 0.105

# 5 rows x 3, grouped so a reader can scan by theme:
#   size -> condensation -> Tif6 puncta -> Sis1 puncta -> nucleolar proximity
PANELS = [
    ("cell_volume_um3",            "Cell volume",             "volume (µm$^3$)"),
    ("nucleolar_volume_um3",       "Nucleolar volume",        "volume (µm$^3$)"),
    ("nucleolar_volume_fraction",  "Nucleolar volume fraction", "nucleolus / cell"),

    ("Tif6_condensate_index_cell", "Tif6 condensate index",   "p95 / mean"),
    ("Nsr1_condensate_index_cell", "Nsr1 condensate index",   "p95 / mean"),
    ("Sis1_condensate_index_cell", "Sis1 condensate index",   "p95 / mean"),

    ("Tif6_puncta_n",                       "Tif6 puncta per cell",  "puncta / cell"),
    ("Tif6_fragmentation_index_persistence","Tif6 fragmentation",    "persistence index"),
    ("Tif6_frac_intensity_in_puncta",       "Tif6 fraction condensed", "fraction of signal"),

    ("Sis1_puncta_n",                       "Sis1 puncta per cell",  "puncta / cell"),
    ("Sis1_fragmentation_index_persistence","Sis1 fragmentation",    "persistence index"),
    ("Sis1_frac_intensity_in_puncta",       "Sis1 fraction condensed", "fraction of signal"),

    ("Sis1_mean_distance",         "Sis1 puncta to nucleolus", "distance (µm)"),
    ("Sis1_fraction_proximal",     "Sis1 puncta proximal",    "fraction $\\leq$0.5 µm"),
    ("n_nucleoli",                 "Nucleoli per cell",       "count"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path,
                    default=(HERE.parent / "jcb_revision_midway3_cellquant_analysis"
                             / "validation_3d" / "outputs_yeast" / "3d"))
    ap.add_argument("--outdir", type=Path, default=HERE / "figure3_3d")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = fig3.load_cells(args.outputs_dir)
    reps = {t: sorted(cells.loc[cells["condition"] == t, "rep_id"].unique())
            for t in fig3.TEMP_ORDER}
    print("replicates per condition:", {t: len(v) for t, v in reps.items()})

    # Sized to the printable area of a portrait page (~6.5 x 9 in) so that at
    # \textwidth it scales by ~0.8 and the 12 pt type lands near 10 pt in print.
    fig = plt.figure(figsize=(8.0, 11.5))
    gs = GridSpec(5, 3, figure=fig, hspace=0.62, wspace=0.46)

    rows = []
    for i, (col, title, ylab) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        if col not in cells.columns:
            ax.text(0.5, 0.5, f"{col}\nnot in 3D output", ha="center", va="center",
                    fontsize=fig3.FONT_SIZE, color="0.45", transform=ax.transAxes)
            ax.set_title(title, fontsize=fig3.FONT_SIZE, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            print(f"  !! {col} absent"); continue
        rm = fig3.replicate_medians(cells, col)
        st = fig3.stats_vs_ref(rm)
        fig3.draw_panel(ax, cells, rm, st, col, title, ylab, reps)
        med = {t: float(np.nanmedian(rm.loc[rm["condition"] == t, "value"].values))
               for t in fig3.TEMP_ORDER}
        sig = [t for t in fig3.TEMP_ORDER if t != fig3.REF and st[t]["sig"]]
        print(f"  {title:28s} " + " ".join(f"{t}={med[t]:.3g}" for t in fig3.TEMP_ORDER)
              + (f"   sig vs 25C: {','.join(sig)}" if sig else "   sig: none"))
        for t in fig3.TEMP_ORDER:
            rows.append({"panel": title, "metric": col, "temperature": t,
                         "replicate_median_of_medians": med[t],
                         "n_replicates": len(reps[t]),
                         "p_raw_vs_25C": st[t]["p_raw"] if t != fig3.REF else np.nan,
                         "p_bonferroni_vs_25C": st[t]["p_adj"] if t != fig3.REF else np.nan})

    stem = args.outdir / "FigureS2_yeast_metrics_3D"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix("." + ext), bbox_inches="tight", dpi=300)
    print(f"\n  wrote {stem.name}.svg / .pdf / .png")
    pd.DataFrame(rows).to_csv(args.outdir / "figureS2_3d_stats.csv", index=False)
    print("  wrote figureS2_3d_stats.csv")


if __name__ == "__main__":
    np.random.seed(0)
    main()
