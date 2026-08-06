"""validation_3d/06_timing_report.py

Aggregates per-image timing from outputs_combined/run_summary.csv into a
markdown report suitable for the manuscript's methods note.
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
setup_paths = importlib.import_module("00_setup_paths")


def _stat_row(label: str, vals: list[float]) -> str:
    if not vals:
        return f"| {label} | — | — | — | — |"
    arr = np.array(vals)
    return (f"| {label} | {len(arr)} | {arr.mean():.1f} | "
            f"{np.median(arr):.1f} | {np.percentile(arr, 90):.1f} |")


def main() -> None:
    summary_csv = setup_paths.OUT_COMBINED / "run_summary.csv"
    if not summary_csv.exists():
        print(f"no run summary at {summary_csv}; nothing to report")
        return
    df = pd.read_csv(summary_csv)
    ok = df[df["status"] == "ok"].copy()
    if ok.empty:
        print("no successful runs; nothing to report")
        return
    ok["dataset"] = ok["label"].str.split("/").str[0]
    ok["axis"] = ok["label"].str.split("/").str[1]
    ok["image"] = ok["label"].str.split("/").str[2]

    md_lines = []
    md_lines.append("# 3D validation timing report")
    md_lines.append("")
    md_lines.append("All times are wall-clock seconds for one invocation on this hardware:")
    md_lines.append("")
    import platform
    md_lines.append(f"- Platform: {platform.platform()}")
    try:
        import psutil
        md_lines.append(
            f"- CPU: {psutil.cpu_count(logical=False)} physical / "
            f"{psutil.cpu_count(logical=True)} logical")
        md_lines.append(
            f"- RAM: {psutil.virtual_memory().total / 1024**3:.0f} GB")
    except Exception:
        pass
    try:
        import torch
        md_lines.append(
            f"- torch device chosen: cuda={torch.cuda.is_available()}, "
            f"mps={hasattr(torch.backends, 'mps') and torch.backends.mps.is_available()}")
    except Exception:
        pass
    md_lines.append(f"- --no-gpu enforced: {setup_paths.NO_GPU}")
    md_lines.append("")

    md_lines.append("## Per-axis wall-clock seconds (per image)")
    md_lines.append("")
    md_lines.append("| Dataset/Axis | n | mean | median | p90 |")
    md_lines.append("|---|---|---|---|---|")
    for dataset in ok["dataset"].unique():
        for axis in sorted(ok.loc[ok["dataset"] == dataset, "axis"].unique()):
            mask = (ok["dataset"] == dataset) & (ok["axis"] == axis)
            md_lines.append(_stat_row(
                f"{dataset} / {axis}",
                ok.loc[mask, "duration_sec"].tolist()))
    md_lines.append("")

    md_lines.append("## Total wall-clock seconds")
    md_lines.append("")
    md_lines.append("| Dataset/Axis | total (s) | total (min) |")
    md_lines.append("|---|---|---|")
    for dataset in ok["dataset"].unique():
        for axis in sorted(ok.loc[ok["dataset"] == dataset, "axis"].unique()):
            mask = (ok["dataset"] == dataset) & (ok["axis"] == axis)
            total = ok.loc[mask, "duration_sec"].sum()
            md_lines.append(
                f"| {dataset} / {axis} | {total:.0f} | {total / 60:.1f} |")
    md_lines.append("")

    md_lines.append("## Headline extrapolation")
    md_lines.append("")
    for dataset in ok["dataset"].unique():
        m3 = ok[(ok["dataset"] == dataset) & (ok["axis"] == "3d")]
        if not m3.empty:
            m = m3["duration_sec"].mean()
            md_lines.append(
                f"- A {dataset} z-stack runs through cellquant 3D in "
                f"**{m / 60:.1f} min** mean (CPU only).")
    md_lines.append("")

    out_md = setup_paths.OUT_COMBINED / "timing_summary.md"
    out_md.write_text("\n".join(md_lines))
    print(f"wrote {out_md}")
    print()
    print("\n".join(md_lines))


if __name__ == "__main__":
    main()
