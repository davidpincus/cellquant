"""Structural regression tests for cellquant.py.

Runs the pipeline on the cropped example datasets and checks that output
files exist, have the expected columns, and contain a reasonable number
of rows. Does NOT check exact values (Cellpose is non-deterministic
across platforms).
"""

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
CELLQUANT = ROOT / "cellquant.py"
MAMMALIAN_DIR = ROOT / "example_data" / "mammalian_SGs"
YEAST_DIR = ROOT / "example_data" / "yeast_temperature"

# This timeout exists to stop a hung run from wedging the suite forever, not to
# assert anything about speed. Budget it for the slowest machine we expect to
# see, not a developer laptop: the whole suite takes ~40 s on Apple Silicon but
# Cellpose on a shared CPU CI runner is more than an order of magnitude slower,
# and a 600 s ceiling failed the mammalian run there intermittently.
TIMEOUT = int(os.environ.get("CELLQUANT_TEST_TIMEOUT", "1800"))


def _run(args: list[str], tmp_path: Path) -> Path:
    """Run cellquant.py and return the output directory."""
    out = tmp_path / "output"
    cmd = [sys.executable, str(CELLQUANT)] + args + ["--out", str(out)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
    assert result.returncode == 0, (
        f"cellquant.py failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    return out


# ── Mammalian SGs ──────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def mammalian_output(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("mammalian")
    return _run(
        [
            str(MAMMALIAN_DIR),
            "1:DAPI:nucleus",
            "2:G3BP1:quantify",
            "3:PABPC1:quantify",
            "--cell-type", "mammalian",
            "--filename-pattern", "MAX_{condition}_rep{replicate}",
            "--skip-plots",
            "--no-save-masks",
        ],
        tmp,
    )


def test_mammalian_output_files(mammalian_output):
    out = mammalian_output
    assert (out / "cells.csv").exists()
    assert (out / "images.csv").exists()
    assert (out / "config_used.yml").exists()
    qc_pngs = list((out / "qc").glob("*.png"))
    assert len(qc_pngs) == 2


def test_mammalian_cells_csv(mammalian_output):
    df = pd.read_csv(mammalian_output / "cells.csv")
    expected_cols = {
        "image", "condition", "replicate", "cell_id",
        "cell_area_px", "keep",
        "G3BP1_cell_mean", "PABPC1_cell_mean",
    }
    assert expected_cols.issubset(set(df.columns))
    assert len(df) > 0


def test_mammalian_images_csv(mammalian_output):
    df = pd.read_csv(mammalian_output / "images.csv")
    assert {"image", "condition", "replicate", "n_cells", "n_keep"}.issubset(
        set(df.columns)
    )
    assert len(df) == 2


def test_mammalian_config(mammalian_output):
    cfg = yaml.safe_load((mammalian_output / "config_used.yml").read_text())
    assert "cell_diameter" in cfg


# ── Yeast temperature ─────────────────────────────────────────────────


@pytest.fixture(scope="session")
def yeast_output(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("yeast")
    return _run(
        [
            str(YEAST_DIR),
            "1:Tif6:quantify",
            "2:Nsr1:nucleolus",
            "3:Sis1:quantify",
            "--cell-type", "yeast",
            "--filename-pattern", "MAX_{condition}_rep{replicate}",
            "--colocalization",
            "--allow-2d-colocalization",
            "--nucleolar-proximity", "Nsr1",
            "--puncta-channels", "Sis1", "Tif6",
            "--skip-plots",
            "--no-save-masks",
        ],
        tmp,
    )


def test_yeast_output_files(yeast_output):
    out = yeast_output
    assert (out / "cells.csv").exists()
    assert (out / "images.csv").exists()
    assert (out / "config_used.yml").exists()
    assert (out / "colocalization.csv").exists()
    assert (out / "nucleolar_proximity.csv").exists()
    assert (out / "nucleolar_morphology.csv").exists()
    qc_pngs = list((out / "qc").glob("*.png"))
    assert len(qc_pngs) == 2


def test_yeast_cells_csv(yeast_output):
    df = pd.read_csv(yeast_output / "cells.csv")
    expected_cols = {
        "image", "condition", "replicate", "cell_id",
        "cell_area_px", "keep",
        "Sis1_puncta_n", "Tif6_puncta_n",
        "pearson_r_Tif6_vs_Sis1",
        "Sis1_mean_distance", "Sis1_fraction_proximal",
        "nucleolar_area", "nucleolar_solidity",
    }
    assert expected_cols.issubset(set(df.columns))
    assert len(df) > 0


def test_yeast_colocalization_csv(yeast_output):
    df = pd.read_csv(yeast_output / "colocalization.csv")
    expected_cols = {
        "image", "cell_id", "pair",
        "pearson_r", "manders_m1", "manders_m2",
    }
    assert expected_cols.issubset(set(df.columns))
    assert len(df) > 0


def test_yeast_proximity_csv(yeast_output):
    df = pd.read_csv(yeast_output / "nucleolar_proximity.csv")
    expected_cols = {
        "image", "cell_id", "channel",
        "n_puncta", "mean_distance", "fraction_proximal",
    }
    assert expected_cols.issubset(set(df.columns))
    assert len(df) > 0


def test_yeast_morphology_csv(yeast_output):
    df = pd.read_csv(yeast_output / "nucleolar_morphology.csv")
    expected_cols = {
        "image", "cell_id",
        "nucleolar_area", "nucleolar_solidity",
        "nucleolar_circularity", "nucleolar_eccentricity",
    }
    assert expected_cols.issubset(set(df.columns))
    assert len(df) > 0


def test_yeast_config(yeast_output):
    cfg = yaml.safe_load((yeast_output / "config_used.yml").read_text())
    assert cfg["colocalization"] is True
