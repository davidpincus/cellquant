"""Equivalence regression test for the vectorized compute_nucleolar_proximity.

The per-cell/per-punctum double loop was replaced by a single pass over
nonzero puncta voxels grouped by (cell_id, punctum_id) pairs. This test pins
that the vectorized implementation is *bit-identical* to the original loop by
running a verbatim copy of the old code (`_compute_nucleolar_proximity_old`)
alongside the current one and asserting equality with check_exact=True.

Coverage:
  - masks derived deterministically from the real example TIFFs (thousands of
    puncta voxels — validates the exact-summation claim at realistic scale, and
    naturally produces puncta that straddle cell borders / fall outside cells);
  - a hand-built 2D case with a punctum straddling two cells (must contribute a
    separate centroid to each cell);
  - a hand-built 3D case (anisotropic EDT + straddle across z).

Exactness holds because voxel coordinates are integer-valued in float64 and
per-punctum coordinate sums stay far below 2**53, so pairwise (np.mean) and
sequential (np.bincount) summation agree, and round()/np.rint() both round half
to even. It is NOT a consequence of traversal order.
"""

import importlib.util
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi
from pandas.testing import assert_frame_equal

ROOT = Path(__file__).resolve().parent.parent

# Import cellquant.py by path (it is a single script, not an installed package).
_spec = importlib.util.spec_from_file_location("cellquant", ROOT / "cellquant.py")
cellquant = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cellquant)

distance_transform_edt = cellquant.distance_transform_edt


# --------------------------------------------------------------------------
# Verbatim copy of the ORIGINAL loop-based implementation (pre-vectorization).
# Kept here as the reference oracle; do not "simplify" it.
# --------------------------------------------------------------------------
def _compute_nucleolar_proximity_old(
    puncta_mask: np.ndarray,
    nucleolar_mask: np.ndarray,
    cell_mask: np.ndarray,
    threshold: float,
    metadata: dict,
    image_name: str,
    channel_name: str,
    cfg: dict | None = None,
    extra_thresholds: list | None = None,
) -> pd.DataFrame:
    is_3d = puncta_mask.ndim == 3
    rows: list[dict[str, Any]] = []
    cfg = cfg or {}
    extra_cols = [(float(t), f"fraction_proximal_{float(t):g}um")
                  for t in (extra_thresholds or [])]

    if is_3d:
        spacing = (
            cfg.get("voxel_size_z_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
            cfg.get("voxel_size_xy_um", 1.0),
        )
        dist_field = distance_transform_edt(~(nucleolar_mask > 0), sampling=spacing)
    else:
        dist_field = distance_transform_edt(~(nucleolar_mask > 0))

    cell_ids = np.unique(cell_mask)
    cell_ids = cell_ids[cell_ids != 0]

    for cid in cell_ids:
        cell_bin = (cell_mask == cid)
        puncta_in_cell = puncta_mask * cell_bin
        punct_ids = np.unique(puncta_in_cell)
        punct_ids = punct_ids[punct_ids != 0]
        n_puncta = len(punct_ids)
        if n_puncta == 0:
            row = {
                "image": image_name,
                "condition": metadata.get("condition", ""),
                "replicate": metadata.get("replicate", ""),
                "cell_id": int(cid),
                "channel": channel_name,
                "n_puncta": 0,
                "mean_distance": np.nan,
                "min_distance": np.nan,
                "fraction_proximal": np.nan,
            }
            for _, col in extra_cols:
                row[col] = np.nan
            rows.append(row)
            continue

        distances = []
        for pid in punct_ids:
            coords = np.where(puncta_in_cell == pid)
            centroid = tuple(int(round(c.mean())) for c in coords)
            d = float(dist_field[centroid])
            distances.append(d)

        distances_arr = np.array(distances)
        row = {
            "image": image_name,
            "condition": metadata.get("condition", ""),
            "replicate": metadata.get("replicate", ""),
            "cell_id": int(cid),
            "channel": channel_name,
            "n_puncta": n_puncta,
            "mean_distance": float(distances_arr.mean()),
            "min_distance": float(distances_arr.min()),
            "fraction_proximal": float((distances_arr <= threshold).mean()),
        }
        for t, col in extra_cols:
            row[col] = float((distances_arr <= t).mean())
        rows.append(row)

    return pd.DataFrame(rows)


META = {"condition": "test", "replicate": "1"}


def _assert_equivalent(puncta_mask, nucleolar_mask, cell_mask, threshold,
                       cfg=None, extra_thresholds=None):
    kwargs = dict(metadata=META, image_name="img", channel_name="chan",
                  cfg=cfg, extra_thresholds=extra_thresholds)
    old = _compute_nucleolar_proximity_old(
        puncta_mask, nucleolar_mask, cell_mask, threshold, **kwargs)
    new = cellquant.compute_nucleolar_proximity(
        puncta_mask, nucleolar_mask, cell_mask, threshold, **kwargs)
    assert_frame_equal(old, new, check_exact=True)
    return new


# --------------------------------------------------------------------------
# Deterministic masks derived from the real example TIFFs (no Cellpose).
# --------------------------------------------------------------------------
def _masks_from_tiff(path: Path):
    """Threshold + connected-component labels → cell/nucleolar/puncta masks.

    Fully deterministic: fixed percentile thresholds, no randomness. Produces
    realistically shaped, overlapping regions from real image geometry, so many
    puncta land across cell borders and outside cells (exercising the straddle
    and out-of-cell branches at thousands-of-voxels scale).
    """
    import tifffile
    stack = tifffile.imread(str(path)).astype(np.float64)  # (C, Y, X)
    c0, c1, c2 = stack[0], stack[1], stack[2]

    cell_bin = c0 > np.percentile(c0, 60)
    cell_mask = ndi.label(cell_bin)[0].astype(np.int32)

    nucleolar_mask = (c1 > np.percentile(c1, 90)).astype(np.int32)

    puncta_bin = c2 > np.percentile(c2, 97)
    puncta_mask = ndi.label(puncta_bin)[0].astype(np.int32)

    return puncta_mask, nucleolar_mask, cell_mask


@pytest.mark.parametrize("tiff", [
    "example_data/yeast_temperature/MAX_25deg_rep1.tif",
    "example_data/mammalian_SGs/MAX_control_rep1.tif",
])
def test_equivalence_on_example_data(tiff):
    path = ROOT / tiff
    if not path.exists():
        pytest.skip(f"missing example data: {tiff}")
    puncta_mask, nucleolar_mask, cell_mask = _masks_from_tiff(path)
    # Sanity: the derived masks actually exercise the interesting branches.
    assert cell_mask.max() > 1 and puncta_mask.max() > 1
    df = _assert_equivalent(
        puncta_mask, nucleolar_mask, cell_mask, threshold=5.0,
        extra_thresholds=[2.0, 10.0])
    assert (df["n_puncta"] > 0).any()


def test_equivalence_2d_straddling_punctum():
    """Punctum 7 spans the border of cells 1 and 2 → a centroid in each cell."""
    cell_mask = np.zeros((20, 40), dtype=np.int32)
    cell_mask[:, :20] = 1
    cell_mask[:, 20:] = 2

    nucleolar_mask = np.zeros((20, 40), dtype=np.int32)
    nucleolar_mask[8:12, 0:3] = 1  # nucleolus near left edge

    puncta_mask = np.zeros((20, 40), dtype=np.int32)
    puncta_mask[5:8, 17:23] = 7    # straddles the x=20 boundary between cells
    puncta_mask[2:4, 30:33] = 8    # wholly inside cell 2
    puncta_mask[15:17, 2:5] = 9    # wholly inside cell 1

    df = _assert_equivalent(
        puncta_mask, nucleolar_mask, cell_mask, threshold=6.0,
        extra_thresholds=[3.0])
    # The straddling punctum must be counted once in each cell.
    assert int(df.loc[df.cell_id == 1, "n_puncta"].iloc[0]) == 2
    assert int(df.loc[df.cell_id == 2, "n_puncta"].iloc[0]) == 2


def test_equivalence_3d_anisotropic():
    """3D anisotropic EDT path with a punctum straddling two cells across z."""
    cell_mask = np.zeros((6, 16, 16), dtype=np.int32)
    cell_mask[:3] = 1
    cell_mask[3:] = 2

    nucleolar_mask = np.zeros((6, 16, 16), dtype=np.int32)
    nucleolar_mask[0:2, 0:3, 0:3] = 1

    puncta_mask = np.zeros((6, 16, 16), dtype=np.int32)
    puncta_mask[2:4, 5:9, 5:9] = 3   # straddles z=3 boundary between cells
    puncta_mask[4:6, 10:12, 10:12] = 4
    puncta_mask[0:1, 12:14, 2:4] = 5

    cfg = {"voxel_size_z_um": 0.23, "voxel_size_xy_um": 0.10571}
    df = _assert_equivalent(
        puncta_mask, nucleolar_mask, cell_mask, threshold=1.0,
        cfg=cfg, extra_thresholds=[0.5, 2.0])
    assert int(df.loc[df.cell_id == 1, "n_puncta"].iloc[0]) == 2
    assert int(df.loc[df.cell_id == 2, "n_puncta"].iloc[0]) == 2


def test_equivalence_cell_with_no_puncta():
    """A cell with no puncta must still emit the n_puncta=0 NaN row."""
    cell_mask = np.zeros((10, 30), dtype=np.int32)
    cell_mask[:, :10] = 1
    cell_mask[:, 10:20] = 2   # this cell gets no puncta
    cell_mask[:, 20:] = 3

    nucleolar_mask = np.zeros((10, 30), dtype=np.int32)
    nucleolar_mask[4:6, 25:28] = 1

    puncta_mask = np.zeros((10, 30), dtype=np.int32)
    puncta_mask[2:4, 3:5] = 1
    puncta_mask[6:8, 24:26] = 2

    df = _assert_equivalent(puncta_mask, nucleolar_mask, cell_mask, threshold=4.0)
    empty = df.loc[df.cell_id == 2].iloc[0]
    assert empty["n_puncta"] == 0
    assert np.isnan(empty["mean_distance"])
