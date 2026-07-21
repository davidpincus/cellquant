"""validation_3d/08_yeast_coloc_figure.py

Pearson R, 3 pairs x 5 temps, 3D vs matched-2D(MIP) overlaid. Superplot style:
small per-cell points colored by replicate, large replicate-median diamonds,
published Table S2 value marked per pair/temp.

Output: validation_3d/yeast_coloc_3D_vs_MIP.png
"""
from __future__ import annotations
import importlib, os, sys
from pathlib import Path
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")

PAIRS = ["Tif6_vs_Sis1", "Tif6_vs_Nsr1", "Sis1_vs_Nsr1"]
TEMPS = ["25C", "30C", "32C", "36C", "40C"]
PUB = {
    "Tif6_vs_Sis1": [0.87, 0.85, 0.84, 0.79, 0.71],
    "Tif6_vs_Nsr1": [0.84, 0.83, 0.82, 0.72, 0.54],
    "Sis1_vs_Nsr1": [0.77, 0.74, 0.73, 0.76, 0.59],
}
ARM_COL = {"3d": "#2166ac", "matched_2d": "#b2182b"}
ARM_LBL = {"3d": "3D (native)", "matched_2d": "MIP (matched-2D)"}


def load_coloc(axis):
    root = {"3d": setup_paths.OUT_YEAST_3D,
            "matched_2d": setup_paths.OUT_YEAST_MATCHED2D}[axis]
    frames = []
    for sub in sorted(root.iterdir()):
        cf = sub / "colocalization.csv"
        if not cf.exists() or cf.stat().st_size == 0:
            continue
        df = pd.read_csv(cf)
        stem = sub.name
        cond, rep = stem.split("_series1_rep", 1)
        df["condition"] = cond
        df["rep_id"] = stem
        df["rep_num"] = int(rep)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def main():
    data = {a: load_coloc(a) for a in ARM_COL}
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)
    cmap = plt.get_cmap("tab10")

    for pi, pair in enumerate(PAIRS):
        ax = axes[pi]
        for ti, temp in enumerate(TEMPS):
            # two sub-positions per temp: 3D left, MIP right
            for ai, axis in enumerate(["3d", "matched_2d"]):
                xc = ti + (-0.18 if ai == 0 else 0.18)
                df = data[axis]
                sub = df[(df["pair"] == pair) & (df["condition"] == temp)].copy()
                sub["pearson_r"] = pd.to_numeric(sub["pearson_r"], errors="coerce")
                # per-cell points colored by rep, jittered
                reps = sorted(sub["rep_num"].unique())
                rng = np.linspace(-0.13, 0.13, max(len(reps), 1))
                rep_meds = []
                for ri, rep in enumerate(reps):
                    rs = sub[sub["rep_num"] == rep]["pearson_r"].dropna().values
                    if len(rs) == 0:
                        continue
                    jit = xc + rng[ri] + np.random.uniform(-0.015, 0.015, len(rs))
                    ax.scatter(jit, rs, s=3, color=cmap(ri % 10), alpha=0.18,
                               linewidths=0, zorder=1)
                    rep_meds.append(np.median(rs))
                # replicate-median diamonds
                for ri, rm in enumerate(rep_meds):
                    ax.scatter(xc, rm, marker="D", s=55,
                               facecolor=cmap(ri % 10), edgecolor=ARM_COL[axis],
                               linewidth=1.6, zorder=4)
                # arm grand median line
                if rep_meds:
                    gm = np.median(rep_meds)
                    ax.plot([xc - 0.16, xc + 0.16], [gm, gm], color=ARM_COL[axis],
                            lw=2.5, zorder=5,
                            label=ARM_LBL[axis] if (pi == 0 and ti == 0) else None)
            # published marker (star) centered on temp
            ax.scatter(ti, PUB[pair][ti], marker="*", s=240, color="gold",
                       edgecolor="black", linewidth=0.9, zorder=6,
                       label="Published (MIP, Table S2)" if (pi == 0 and ti == 0) else None)

        ax.set_title(pair.replace("_vs_", " vs "), fontsize=13, fontweight="bold")
        ax.set_xticks(range(len(TEMPS)))
        ax.set_xticklabels(TEMPS)
        ax.set_xlabel("Temperature")
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(axis="y", alpha=0.25)
        if pi == 0:
            ax.set_ylabel("Pearson's R (whole-cell, Costes)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=11,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Yeast colocalization: 3D-native vs MIP (matched masks)  "
                 "— diamonds=replicate medians, bars=arm median, star=published",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    out = Path(__file__).resolve().parent / "yeast_coloc_3D_vs_MIP.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    np.random.seed(0)
    main()
