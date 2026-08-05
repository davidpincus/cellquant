"""validation_3d/12_tables_s2_s3_3d.py

Regenerate supplementary Tables S2 (values) and S3 (p-values) for the yeast
temperature series from the native 3D metrics, under ONE statistical convention.

THE CONVENTION (used here, in Figure 2B, and in Figure 3B):
  - unit of analysis: the replicate (image) median
  - reported value:   the median of those replicate medians
  - test:             two-sided Mann-Whitney U (scipy default: exact when n is small and untied), each temperature vs 25 C
  - correction:       Bonferroni across the four temperature comparisons
At n = 5-8 replicates the exact test and the normal approximation
(scipy.stats.ranksums) are not interchangeable -- e.g. Tif6-Nsr1 at 36 C is
0.0022 exact vs 0.0039 approximate, which after x4 is 0.009 vs 0.016. Mixing
the two within one table is what produced the inconsistency this replaces.

Nucleolar shape descriptors (circularity, solidity, eccentricity) are properties
of a projection and have no 3D counterpart; nucleolar volume and equivalent
diameter take their place. Distances are in microns, not pixels.

Usage:
    python 12_tables_s2_s3_3d.py [--outputs-dir DIR] [--outdir DIR]
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

TEMPS = ["25C", "30C", "32C", "36C", "40C"]
REF = "25C"
NCOMP = 4

BLOCKS = [
    ("Condensate metrics", [
        ("Sis1_condensate_index_cell",      "Sis1 condensate index",   "{:.2f}"),
        ("Sis1_frac_intensity_in_puncta",   "Sis1 frac.\\ condensed",  "{:.3f}"),
        ("Sis1_puncta_n",                   "Sis1 puncta/cell",        "{:.2f}"),
        ("Tif6_condensate_index_cell",      "Tif6 condensate index",   "{:.2f}"),
        ("Tif6_frac_intensity_in_puncta",   "Tif6 frac.\\ condensed",  "{:.3f}"),
        ("Tif6_puncta_n",                   "Tif6 puncta/cell",        "{:.2f}"),
        ("Nsr1_condensate_index_cell",      "Nsr1 condensate index",   "{:.2f}"),
    ]),
    ("Colocalization", [
        ("pearson_r_Tif6_vs_Sis1",  "Pearson $R$ (Tif6 vs Sis1)", "{:.2f}"),
        ("pearson_r_Tif6_vs_Nsr1",  "Pearson $R$ (Tif6 vs Nsr1)", "{:.2f}"),
        ("pearson_r_Sis1_vs_Nsr1",  "Pearson $R$ (Sis1 vs Nsr1)", "{:.2f}"),
        ("manders_m1_Tif6_vs_Sis1", "Manders M1 (Tif6 vs Sis1)",  "{:.2f}"),
        ("manders_m1_Tif6_vs_Nsr1", "Manders M1 (Tif6 vs Nsr1)",  "{:.2f}"),
        ("manders_m1_Sis1_vs_Nsr1", "Manders M1 (Sis1 vs Nsr1)",  "{:.2f}"),
    ]),
    ("Nucleolar proximity", [
        ("Sis1_fraction_proximal", "Sis1 fraction proximal",        "{:.2f}"),
        ("Sis1_mean_distance",     "Sis1 mean distance (\\textmu m)", "{:.2f}"),
        ("Tif6_fraction_proximal", "Tif6 fraction proximal",        "{:.2f}"),
        ("Tif6_mean_distance",     "Tif6 mean distance (\\textmu m)", "{:.2f}"),
    ]),
    ("Cell and nucleolar size", [
        ("cell_volume_um3",          "Cell volume (\\textmu m$^3$)",      "{:.1f}"),
        ("nucleolar_volume_um3",     "Nucleolar volume (\\textmu m$^3$)", "{:.2f}"),
        ("nucleolar_eq_diameter_um", "Nucleolar eq.\\ diameter (\\textmu m)", "{:.2f}"),
    ]),
]


def load_cells(d: Path) -> pd.DataFrame:
    files = sorted(d.glob("*/cells.csv"))
    if not files:
        raise FileNotFoundError(f"no */cells.csv under {d}")
    fr = []
    for f in files:
        df = pd.read_csv(f)
        rid = f.parent.name
        if "condition" not in df.columns or df["condition"].isna().all():
            cond, rep = rid.split("_series1_rep", 1)
            df["condition"], df["replicate"] = cond, rep
        df["rep_id"] = rid
        fr.append(df)
    c = pd.concat(fr, ignore_index=True)
    if "keep" in c.columns:
        c = c[c["keep"] == True].copy()   # noqa: E712
    return c.reset_index(drop=True)


def stats_for(cells: pd.DataFrame, col: str):
    rm = (cells.groupby(["condition", "rep_id"], as_index=False)[col]
               .median().rename(columns={col: "v"}))
    med = {t: float(np.nanmedian(rm.loc[rm["condition"] == t, "v"].values))
           if (rm["condition"] == t).any() else np.nan for t in TEMPS}
    ref = rm.loc[rm["condition"] == REF, "v"].dropna().values
    padj = {}
    for t in TEMPS:
        if t == REF:
            continue
        v = rm.loc[rm["condition"] == t, "v"].dropna().values
        if len(v) < 2 or len(ref) < 2:
            padj[t] = np.nan
        else:
            p = mannwhitneyu(v, ref, alternative="two-sided").pvalue
            padj[t] = min(p * NCOMP, 1.0)
    return med, padj


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument("--outputs-dir", type=Path,
                    default=(here.parent / "jcb_revision_midway3_cellquant_analysis"
                             / "validation_3d" / "outputs_yeast" / "3d"))
    ap.add_argument("--outdir", type=Path, default=here / "tables_3d")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(args.outputs_dir)
    ncells = {t: int((cells["condition"] == t).sum()) for t in TEMPS}
    nreps = {t: int(cells.loc[cells["condition"] == t, "rep_id"].nunique()) for t in TEMPS}
    print("cells:", ncells, "\nreps :", nreps)

    s2 = [r"\small", r"\begin{tabular}{lccccc}", r"\toprule",
          r"Metric & " + " & ".join(f"{t[:-1]}\\textdegree C" for t in TEMPS) + r" \\",
          r"\midrule",
          r"$n$ cells & " + " & ".join(str(ncells[t]) for t in TEMPS) + r" \\",
          r"$n$ replicates & " + " & ".join(str(nreps[t]) for t in TEMPS) + r" \\"]
    s3 = [r"\small", r"\begin{tabular}{lcccc}", r"\toprule",
          r"Metric & " + " & ".join(f"{t[:-1]}\\textdegree C" for t in TEMPS[1:]) + r" \\",
          r"\midrule"]
    rows = []

    for block, metrics in BLOCKS:
        s2 += [r"\addlinespace", r"\multicolumn{6}{l}{\textit{%s}} \\" % block]
        s3 += [r"\addlinespace", r"\multicolumn{5}{l}{\textit{%s}} \\" % block]
        for col, label, fmt in metrics:
            if col not in cells.columns:
                print(f"  !! {col} absent — skipped"); continue
            med, padj = stats_for(cells, col)
            s2.append(label + " & " + " & ".join(
                (fmt.format(med[t]) if np.isfinite(med[t]) else "--") for t in TEMPS) + r" \\")
            cellsS3 = []
            for t in TEMPS[1:]:
                p = padj[t]
                if not np.isfinite(p):
                    cellsS3.append("--")
                elif p < 0.05:
                    cellsS3.append(r"\textbf{%s}" % (f"{p:.3f}" if p < 0.01 else f"{p:.2f}"))
                else:
                    cellsS3.append(f"{p:.3f}" if p < 0.01 else f"{p:.2f}")
            s3.append(label + " & " + " & ".join(cellsS3) + r" \\")
            for t in TEMPS:
                rows.append({"block": block, "metric": col, "label": label,
                             "temperature": t, "median_of_replicate_medians": med[t],
                             "n_replicates": nreps[t],
                             "p_bonferroni_vs_25C": padj.get(t, np.nan)})
            sig = [t for t in TEMPS[1:] if np.isfinite(padj[t]) and padj[t] < 0.05]
            print(f"  {label:38s} " + "  ".join(
                f"{t}={med[t]:.3g}" for t in TEMPS) + f"   sig: {','.join(sig) or 'none'}")

    s2 += [r"\bottomrule", r"\end{tabular}"]
    s3 += [r"\bottomrule", r"\end{tabular}"]
    (args.outdir / "table_s2_body.tex").write_text("\n".join(s2) + "\n")
    (args.outdir / "table_s3_body.tex").write_text("\n".join(s3) + "\n")
    pd.DataFrame(rows).to_csv(args.outdir / "tables_s2_s3_3d.csv", index=False)
    print("\n  wrote table_s2_body.tex, table_s3_body.tex, tables_s2_s3_3d.csv")


if __name__ == "__main__":
    main()
