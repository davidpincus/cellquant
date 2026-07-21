"""validation_3d/07_report.py

Assembles outputs from steps 1-6 into a single self-contained HTML report.
PNGs base64-embedded; CSVs rendered as tables inline.
"""
from __future__ import annotations

import base64
import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")


def _png_to_b64_img(path: Path, width_pct: int = 90) -> str:
    if not path.exists():
        return f"<p><em>missing: {path}</em></p>"
    data = path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return (f'<img src="data:image/png;base64,{b64}" '
            f'style="max-width:{width_pct}%; display:block; margin:1em auto;" />')


def _csv_to_html(path: Path, max_rows: int = 200) -> str:
    if not path.exists():
        return f"<p><em>missing: {path}</em></p>"
    df = pd.read_csv(path)
    if len(df) > max_rows:
        df = df.head(max_rows)
    return df.to_html(index=False, classes="data-table",
                      float_format=lambda v: f"{v:.4g}" if pd.notna(v) else "")


def _headline_verdict() -> str:
    """Read aggregate_recapitulation.csv + per_cell_summary.csv to draft a verdict."""
    agg_path = setup_paths.OUT_MAMM_COMPARISONS / "aggregate_recapitulation.csv"
    perc_path = setup_paths.OUT_MAMM_COMPARISONS / "per_cell_summary.csv"
    lines = []
    if agg_path.exists():
        agg = pd.read_csv(agg_path)
        n_metrics = len(agg)
        if "all_axes_agree_direction" in agg.columns:
            n_agree = int(agg["all_axes_agree_direction"].sum())
            lines.append(
                f"Across {n_metrics} mammalian metrics, the arsenite effect "
                f"direction agreed across all three axes (3D / matched-2D / "
                f"published-2D) for <strong>{n_agree}/{n_metrics}</strong>.")
        # significant in 3D, in 2D-published?
        if "3d_p" in agg.columns:
            sig_3d = int((agg["3d_p"] < 0.05).sum())
            sig_pub = int((agg["published_2d_p"] < 0.05).sum()) if "published_2d_p" in agg.columns else None
            lines.append(
                f"At p &lt; 0.05 (Wilcoxon rank-sum on replicate medians), "
                f"<strong>{sig_3d}/{n_metrics}</strong> metrics are significant in 3D"
                + (f"; <strong>{sig_pub}/{n_metrics}</strong> in independent 2D-on-MIP."
                   if sig_pub is not None else "."))
    if perc_path.exists():
        perc = pd.read_csv(perc_path)
        if "agreement" in perc.columns:
            n_high = (perc["agreement"] == "high agreement").sum()
            n_mod = (perc["agreement"] == "moderate agreement").sum()
            n_dis = (perc["agreement"] == "discrepant").sum()
            lines.append(
                f"Per-cell paired comparison (3D-native vs 2D-from-MIP using "
                f"the same cells): <strong>{n_high}</strong> metrics in high agreement, "
                f"<strong>{n_mod}</strong> in moderate, <strong>{n_dis}</strong> discrepant.")
    if not lines:
        return "<p><em>Validation outputs not yet available.</em></p>"
    return "<p>" + "</p><p>".join(lines) + "</p>"


def _yeast_headline_verdict() -> str:
    """Read 04_compare_yeast outputs to draft the yeast verdict."""
    out_dir = setup_paths.OUT_YEAST_COMPARISONS
    agree_path = out_dir / "directional_agreement.csv"
    disagree_path = out_dir / "disagreements.txt"
    if not agree_path.exists():
        return "<p><em>Yeast comparison outputs not yet generated.</em></p>"
    df = pd.read_csv(agree_path)
    lines = []
    n_rows = len(df)
    if n_rows == 0:
        return "<p><em>Yeast comparison ran but produced no rows.</em></p>"
    # Direction agreement vs published S2 by axis (only count rows where
    # published has a non-empty sign).
    has_pub = df["published_sign"].isin(["+", "-"])
    base = df[has_pub]
    n_compare = len(base)
    for axis in ("3d", "matched_2d", "published_2d"):
        sign_col = f"{axis}_sign"
        if sign_col not in base.columns:
            continue
        agreeing = ((base[sign_col].isin(["+", "-"]))
                    & (base[sign_col] == base["published_sign"])).sum()
        lines.append(
            f"<strong>{axis.replace('_', '-')}</strong>: direction matches "
            f"published Table S2 for "
            f"<strong>{agreeing}/{n_compare}</strong> "
            f"(metric, temperature) cells.")
    # Disagreement count
    if disagree_path.exists():
        txt = disagree_path.read_text()
        n_dis = sum(1 for ln in txt.splitlines()
                    if ":" in ln and "published" in ln)
        if n_dis:
            lines.append(
                f"<strong>{n_dis} direction disagreement(s)</strong> flagged — "
                f"see the disagreements section below.")
        else:
            lines.append("No direction disagreements vs Table S2 detected.")
    return "<p>" + "</p><p>".join(lines) + "</p>"


def main() -> None:
    out_dir = setup_paths.VALIDATION_ROOT
    parts = []
    parts.append("<!DOCTYPE html>")
    parts.append("<html><head><meta charset='utf-8'>")
    parts.append("<title>cellquant 3D validation report</title>")
    parts.append("""<style>
body { font-family: -apple-system, Helvetica, Arial, sans-serif;
       max-width: 1100px; margin: 2em auto; padding: 0 2em; line-height: 1.5;
       color: #222; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.3em; }
h2 { border-bottom: 1px solid #ccc; margin-top: 2em; padding-bottom: 0.2em; }
h3 { color: #555; }
.headline { background: #f7f7e8; padding: 1em 1.5em; border-left: 4px solid #c80;
            margin: 1em 0 2em 0; }
.data-table { border-collapse: collapse; margin: 1em 0; font-size: 0.9em; }
.data-table th, .data-table td { border: 1px solid #ccc; padding: 4px 8px; }
.data-table th { background: #eee; }
.caveat { background: #fdf5f5; padding: 0.6em 1em; border-left: 3px solid #c33;
          margin: 1em 0; font-size: 0.95em; }
.muted { color: #888; font-size: 0.9em; }
</style></head><body>""")
    parts.append("<h1>cellquant 3D validation report</h1>")
    parts.append(f"<p class='muted'>Generated by validation_3d/07_report.py</p>")

    # Headline verdict
    parts.append("<h2>Headline verdict</h2>")
    parts.append("<div class='headline'>")
    parts.append(_headline_verdict())
    parts.append("</div>")

    # Caveats up front
    parts.append("<h2>Caveats</h2>")
    parts.append("<div class='caveat'>")
    parts.append("<strong>Yeast subset is reduced.</strong> Per-image 3D runtime on "
                 "Apple Silicon CPU is ~5 hours with <code>seg_3d_method=full</code> "
                 "and ~4–5 hours with <code>stitch</code> (Cellpose-cpsam is the "
                 "dominant cost on 71-slice 1192×1200 stacks). Running the full "
                 "30-image manuscript subset would take ~7 days. The yeast validation "
                 "below uses a <strong>one-rep-per-temperature subset (5 images)</strong> "
                 "with <code>--seg-3d-method stitch</code>. Replicate-level Wilcoxon "
                 "is undefined at n=1/condition for this subset; the comparison "
                 "instead leans on directional agreement vs published Tables S2/S3.")
    parts.append("</div>")
    parts.append("<div class='caveat'>")
    parts.append("<strong>Mammalian z-stacks and published MIPs are a partial overlap.</strong> "
                 "Six z-stacks are available (control_rep1–3 + arsenite_rep1–3); seven "
                 "published MIPs are available (control_rep1, control_rep4, control_rep5 + "
                 "arsenite_rep1–4). Only four images have both axes available for "
                 "per-cell paired comparison. Aggregate comparisons use all available "
                 "images per axis but are not strictly matched.")
    parts.append("</div>")

    # Per-cell paired comparison
    parts.append("<h2>Per-cell paired comparison (3D-native vs 2D-from-3D-mask)</h2>")
    parts.append("<p>For the 4-image overlap, each 3D cell mask is projected to its XY "
                 "footprint and used to compute matched-2D metrics from the MIP. "
                 "This gives a paired dataset where the same cell appears in both axes "
                 "and any systematic 2D-vs-3D bias is visible per cell.</p>")
    parts.append(_csv_to_html(setup_paths.OUT_MAMM_COMPARISONS / "per_cell_summary.csv"))

    # Per-cell PNGs
    parts.append("<h3>Per-metric paired scatters + Bland-Altman plots</h3>")
    metric_dir = setup_paths.OUT_MAMM_COMPARISONS
    if metric_dir.exists():
        png_files = sorted(metric_dir.glob("per_cell_*.png"))
        for png in png_files:
            parts.append(f"<h4>{png.stem.replace('per_cell_', '')}</h4>")
            parts.append(_png_to_b64_img(png, width_pct=85))

    # Aggregate
    parts.append("<h2>Aggregate recapitulation (does 3D recover the published conclusion?)</h2>")
    parts.append("<p>Replicate-level Wilcoxon rank-sum on each metric for arsenite vs "
                 "control, computed independently on the 3D pipeline output, the matched-2D "
                 "derivation, and the independent cellquant-on-MIP run. "
                 "Signs (+/-) indicate effect direction; <code>p</code> is the two-sided "
                 "p-value.</p>")
    parts.append(_csv_to_html(setup_paths.OUT_MAMM_COMPARISONS / "aggregate_recapitulation.csv"))
    parts.append(_png_to_b64_img(setup_paths.OUT_MAMM_COMPARISONS / "aggregate_bars.png"))

    # Synthetic
    parts.append("<h2>Synthetic sphericity ground truth</h2>")
    parts.append("<p>cellquant's marching-cubes sphericity calculation, tested against "
                 "analytic values for spheres, prolate/oblate spheroids, plus qualitative "
                 "low-sphericity shapes (dumbbell, crescent).</p>")
    parts.append(_csv_to_html(setup_paths.OUT_SYNTH / "results.csv"))
    parts.append(_png_to_b64_img(setup_paths.OUT_SYNTH / "results.png"))

    # Algorithmic sanity
    parts.append("<h2>Algorithmic sanity panels</h2>")
    sanity_dir = setup_paths.OUT_MAMM_COMPARISONS / "sanity"
    if sanity_dir.exists():
        for png in sorted(sanity_dir.glob("*.png")):
            parts.append(f"<h3>{png.stem}</h3>")
            parts.append(_png_to_b64_img(png, width_pct=92))
        # ROI coloc three-way summary
        roi_csv = sanity_dir / "sanity_e_summary.csv"
        if roi_csv.exists():
            parts.append("<h4>Three-way colocalization summary (one cell)</h4>")
            parts.append(_csv_to_html(roi_csv))

    # Yeast temperature series
    parts.append("<h2>Yeast temperature series (3D + matched-2D + published-2D)</h2>")
    parts.append("<p>One-rep-per-temperature subset (25C/30C/32C/36C/40C, rep1 each) "
                 "with stitch 3D segmentation. For each metric and temperature, "
                 "median is taken across replicates (n=1 for this subset; reported "
                 "for completeness only). Direction (+/−) and significance (p &lt; 0.05) "
                 "are compared against published Tables S2 and S3.</p>")
    parts.append("<div class='headline'>")
    parts.append(_yeast_headline_verdict())
    parts.append("</div>")

    yeast_cmp = setup_paths.OUT_YEAST_COMPARISONS
    parts.append("<h3>Directional agreement vs published Table S2</h3>")
    parts.append("<p>Sign (+/−) of each axis's median vs 25C reference, alongside the "
                 "published S2 means' direction. <code>any_direction_disagreement</code> "
                 "= true flags an honest disagreement worth investigating.</p>")
    parts.append(_csv_to_html(yeast_cmp / "directional_agreement.csv"))

    parts.append("<h3>Wilcoxon rank-sum vs 25C reference</h3>")
    parts.append("<p>Replicate-level Wilcoxon test on per-image medians, by axis and "
                 "temperature. P-values are NaN where n &lt; 2 per side — expected for the "
                 "one-rep-per-temperature subset.</p>")
    parts.append(_csv_to_html(yeast_cmp / "wilcoxon_vs_25C.csv"))

    parts.append("<h3>Per-temperature medians</h3>")
    parts.append(_csv_to_html(yeast_cmp / "per_temperature_medians.csv"))

    parts.append("<h3>Direction disagreements (if any)</h3>")
    disagree_path = yeast_cmp / "disagreements.txt"
    if disagree_path.exists():
        parts.append("<pre>" + disagree_path.read_text() + "</pre>")
    else:
        parts.append("<p><em>Comparison not yet run.</em></p>")

    # Timing
    parts.append("<h2>Timing baseline</h2>")
    timing_md = setup_paths.OUT_COMBINED / "timing_summary.md"
    if timing_md.exists():
        # render markdown table-of-contents inline (basic)
        for line in timing_md.read_text().splitlines():
            if line.startswith("# "):
                parts.append(f"<h3>{line[2:]}</h3>")
            elif line.startswith("## "):
                parts.append(f"<h4>{line[3:]}</h4>")
            elif line.startswith("- "):
                parts.append(f"<p>{line[2:]}</p>")
            elif line.startswith("|"):
                parts.append(f"<pre>{line}</pre>")
            elif line.strip() == "":
                parts.append("")
            else:
                parts.append(f"<p>{line}</p>")

    # Failures
    parts.append("<h2>Failures and skipped steps</h2>")
    if setup_paths.OUT_FAILURES.exists():
        parts.append("<pre>" + setup_paths.OUT_FAILURES.read_text() + "</pre>")
    else:
        parts.append("<p>No failures logged.</p>")

    parts.append("</body></html>")
    out_path = setup_paths.OUT_REPORT
    out_path.write_text("".join(parts))
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
