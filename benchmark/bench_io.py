"""Shared IO for the cellquant-vs-CellProfiler benchmark.

Loaders extracted verbatim from compare_tools.py so replicate_compare.py and
puncta_diagnostic.py reuse the exact same parsing (column mapping, cellquant
centroid-from-mask derivation, CellProfiler per-cell puncta-area aggregation).
compare_tools.py imports these back, so its behaviour is unchanged.

Plus puncta-geometry helpers (Part 2): cellquant exposes a LABELED puncta mask
({image}_G3BP1_punctamask.tif, one label per punctum), so per-punctum areas and
per-cell assignment are available directly.

numpy + pandas + scikit-image(regionprops) + tifffile only.
"""
import os

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops


# ---- Verified column mapping (fill/override from the tools' actual outputs) ----
COLMAP = {
    "cellquant": {
        "image": "image",              # e.g. MAX_control_rep1.tif
        "cell_id": "cell_id",          # == integer label in {image}_cellmask.tif
        "keep": "keep",                # analysis-set filter (1-4 nuclei); optional
        "count": "G3BP1_puncta_n",
        "puncta_area_px": "G3BP1_puncta_area_px",
        "cell_area_px": "cell_area_px",
        "intensity": "G3BP1_cell_mean",
    },
    "cp_cells": {
        "imagenumber": "ImageNumber",
        "object": "ObjectNumber",
        "cx": "Location_Center_X",
        "cy": "Location_Center_Y",
        "cell_area": "AreaShape_Area",
        "count": "Children_Puncta_Count",
        "intensity": "Intensity_MeanIntensity_G3BP1",
    },
    "cp_puncta": {
        "imagenumber": "ImageNumber",
        "parent_cell": "Parent_Cells",
        "area": "AreaShape_Area",
        "cx": "Location_Center_X",
        "cy": "Location_Center_Y",
    },
    "cp_image": {
        "imagenumber": "ImageNumber",
        "filename_g3bp1": "FileName_G3BP1",   # -> base image key
    },
}


def _require_cols(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: {what} is missing expected column(s) {missing}.\n"
            f"       Present columns: {list(df.columns)}\n"
            f"       Fix the mapping (COLMAP / --*-col overrides)."
        )


def _nonempty(df, path):
    if df is None or len(df) == 0:
        raise SystemExit(f"ERROR: {path} is empty (0 rows) — refusing to proceed (fail loud).")
    return df


def base_of(image_name):
    """MAX_control_rep1.tif -> MAX_control_rep1 ; strip a trailing _<CHANNEL>."""
    b = str(image_name)
    for ext in (".tif", ".tiff"):
        if b.endswith(ext):
            b = b[: -len(ext)]
    for tok in ("_G3BP1", "_DAPI", "_PABPC1"):
        if b.endswith(tok):
            b = b[: -len(tok)]
    return b


def load_cellquant(cells_csv, masks_dir, cm):
    if not os.path.exists(cells_csv):
        raise SystemExit(f"ERROR: cellquant cells.csv not found: {cells_csv}")
    df = _nonempty(pd.read_csv(cells_csv), cells_csv)
    _require_cols(df, [cm["image"], cm["cell_id"], cm["count"],
                       cm["puncta_area_px"], cm["cell_area_px"], cm["intensity"]], "cellquant cells.csv")
    if cm["keep"] in df.columns:
        n0 = len(df)
        df = df[df[cm["keep"]].astype(bool)].copy()
        print(f"[cellquant] kept {len(df)}/{n0} cells passing the '{cm['keep']}' filter (analysis set).")
    df["base"] = df[cm["image"]].map(base_of)
    df["measure_count"] = df[cm["count"]].astype(float)
    df["measure_area_fraction"] = df[cm["puncta_area_px"]].astype(float) / df[cm["cell_area_px"]].astype(float)
    df["measure_intensity"] = df[cm["intensity"]].astype(float)

    # centroids from the integer cell masks (label == cell_id)
    cx, cy = [], []
    cache = {}
    n_missing_mask = 0
    for _, r in df.iterrows():
        b = r["base"]
        if b not in cache:
            mpath = os.path.join(masks_dir, f"{b}_cellmask.tif")
            if not os.path.exists(mpath):
                raise SystemExit(
                    f"ERROR: cellquant cell mask not found: {mpath}\n"
                    f"       centroids are derived from the masks (cells.csv has no centroid);\n"
                    f"       point --cellquant-masks at the run's masks/ dir."
                )
            m = tifffile.imread(mpath)
            cache[b] = {p.label: (float(p.centroid[1]), float(p.centroid[0])) for p in regionprops(m)}
        c = cache[b].get(int(r[cm["cell_id"]]))
        if c is None:
            n_missing_mask += 1
            cx.append(np.nan); cy.append(np.nan)
        else:
            cx.append(c[0]); cy.append(c[1])
    df["cx"], df["cy"] = cx, cy
    if n_missing_mask:
        print(f"[cellquant] WARNING: {n_missing_mask} cells in cells.csv had no matching label in the mask.")
    df["equiv_d"] = 2.0 * np.sqrt(df[cm["cell_area_px"]].astype(float) / np.pi)
    df = df.dropna(subset=["cx", "cy"]).reset_index(drop=True)
    df["uid"] = df["base"] + "#cq" + df[cm["cell_id"]].astype(int).astype(str)
    return df


def load_cellprofiler(cp_cells_csv, cp_puncta_csv, cp_image_csv, cm_c, cm_p, cm_i):
    if not os.path.exists(cp_cells_csv):
        raise SystemExit(f"ERROR: CellProfiler cp_Cells.csv not found: {cp_cells_csv}")
    cells = _nonempty(pd.read_csv(cp_cells_csv), cp_cells_csv)
    punc = pd.read_csv(cp_puncta_csv) if os.path.exists(cp_puncta_csv) else pd.DataFrame()
    _require_cols(cells, [cm_c["imagenumber"], cm_c["object"], cm_c["cx"], cm_c["cy"],
                          cm_c["cell_area"], cm_c["count"], cm_c["intensity"]], "cp_Cells.csv")

    # ImageNumber -> base, via cp_Image.csv
    if not os.path.exists(cp_image_csv):
        raise SystemExit(
            f"ERROR: cp_Image.csv not found: {cp_image_csv}\n"
            f"       needed to map ImageNumber -> image name (cp_Cells has no filename).\n"
            f"       pass --cp-image or keep it beside --cp-cells."
        )
    imgdf = _nonempty(pd.read_csv(cp_image_csv), cp_image_csv)
    _require_cols(imgdf, [cm_i["imagenumber"], cm_i["filename_g3bp1"]], "cp_Image.csv")
    imgdf["base"] = imgdf[cm_i["filename_g3bp1"]].map(base_of)
    num2base = dict(zip(imgdf[cm_i["imagenumber"]], imgdf["base"]))
    cells["base"] = cells[cm_c["imagenumber"]].map(num2base)

    # per-cell puncta area from parent links (sum child puncta area / cell area)
    punc_area_by_cell = {}
    if len(punc):
        _require_cols(punc, [cm_p["imagenumber"], cm_p["parent_cell"], cm_p["area"]], "cp_Puncta.csv")
        g = punc[punc[cm_p["parent_cell"]] > 0].groupby(
            [cm_p["imagenumber"], cm_p["parent_cell"]])[cm_p["area"]].sum()
        punc_area_by_cell = g.to_dict()

    def cell_punc_area(row):
        return float(punc_area_by_cell.get((row[cm_c["imagenumber"]], row[cm_c["object"]]), 0.0))

    cells["punc_area"] = cells.apply(cell_punc_area, axis=1)
    cells["measure_count"] = cells[cm_c["count"]].astype(float)
    cells["measure_area_fraction"] = cells["punc_area"] / cells[cm_c["cell_area"]].astype(float)
    cells["measure_intensity"] = cells[cm_c["intensity"]].astype(float)
    cells["cx"] = cells[cm_c["cx"]].astype(float)
    cells["cy"] = cells[cm_c["cy"]].astype(float)
    cells["equiv_d"] = 2.0 * np.sqrt(cells[cm_c["cell_area"]].astype(float) / np.pi)
    cells["uid"] = cells["base"] + "#cp" + cells[cm_c["object"]].astype(int).astype(str)
    return cells


# ---------------------------------------------------------------------------
# Puncta-geometry helpers (Part 2 diagnostic)
# ---------------------------------------------------------------------------
def cellquant_punctamask(masks_dir, base, puncta_channel="G3BP1"):
    """Return (punctamask_array, cellmask_array_or_None); (None, None) if no punctamask.

    cellquant writes a LABELED puncta mask per image ({base}_{ch}_punctamask.tif),
    so per-punctum geometry is available directly (no fabrication needed).
    """
    ppath = os.path.join(masks_dir, f"{base}_{puncta_channel}_punctamask.tif")
    if not os.path.exists(ppath):
        return None, None
    pm = tifffile.imread(ppath)
    cpath = os.path.join(masks_dir, f"{base}_cellmask.tif")
    cmk = tifffile.imread(cpath) if os.path.exists(cpath) else None
    return pm, cmk


def cellquant_puncta_table(masks_dir, base, puncta_channel="G3BP1"):
    """Per-punctum table from the labeled punctamask: label, cx, cy, area, cell_id, bbox.

    cell_id = cell mask value at the punctum centroid (reproduces G3BP1_puncta_n
    counted per cell). Returns None if cellquant exposes no punctamask.
    """
    pm, cmk = cellquant_punctamask(masks_dir, base, puncta_channel)
    if pm is None:
        return None
    rows = []
    for p in regionprops(pm):
        cyy, cxx = p.centroid
        cid = int(cmk[int(round(cyy)), int(round(cxx))]) if cmk is not None else -1
        rows.append({"label": int(p.label), "cx": float(cxx), "cy": float(cyy),
                     "area": float(p.area), "cell_id": cid, "bbox": tuple(int(v) for v in p.bbox)})
    return pd.DataFrame(rows, columns=["label", "cx", "cy", "area", "cell_id", "bbox"])


def load_cp_puncta(cp_puncta_csv, cp_image_csv, cm_p=None, cm_i=None):
    """CellProfiler puncta table with a 'base' image column (via cp_Image.csv)."""
    cm_p = cm_p or COLMAP["cp_puncta"]
    cm_i = cm_i or COLMAP["cp_image"]
    if not os.path.exists(cp_puncta_csv):
        raise SystemExit(f"ERROR: cp_Puncta.csv not found: {cp_puncta_csv}")
    punc = _nonempty(pd.read_csv(cp_puncta_csv), cp_puncta_csv)
    _require_cols(punc, [cm_p["imagenumber"], cm_p["parent_cell"], cm_p["area"],
                         cm_p["cx"], cm_p["cy"]], "cp_Puncta.csv")
    imgdf = _nonempty(pd.read_csv(cp_image_csv), cp_image_csv)
    imgdf["base"] = imgdf[cm_i["filename_g3bp1"]].map(base_of)
    num2base = dict(zip(imgdf[cm_i["imagenumber"]], imgdf["base"]))
    punc = punc.copy()
    punc["base"] = punc[cm_p["imagenumber"]].map(num2base)
    return punc


def read_channel(images_dir, base, channel_index):
    """Read one channel of a multichannel source MIP as a 2-D array. Fails loud."""
    for ext in (".tif", ".tiff"):
        p = os.path.join(images_dir, base + ext)
        if os.path.exists(p):
            a = tifffile.imread(p)
            if a.ndim != 3 or channel_index >= a.shape[0]:
                raise SystemExit(f"ERROR: {p} shape {a.shape}; cannot take channel index {channel_index}.")
            return a[channel_index]
    raise SystemExit(f"ERROR: source image not found for base '{base}' in {images_dir}")
