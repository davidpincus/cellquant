# Changelog

All notable changes to cellquant are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html) — where "breaking" means a change
that could make the same images produce different numbers.

## [1.1.1] — 2026-09-01

Repository and release infrastructure. No change to any measurement — the pipeline produces
identical output to 1.1.0.

### Added

- `--version`, so a bug report can name the version it came from.
- Continuous integration: the suite runs on Linux and macOS across Python 3.11 and 3.12.
- `CHANGELOG.md`, `CONTRIBUTING.md`, issue templates, and a pull-request template.
- `tools/` — the Zenodo uploader that archives the raw image data, its dataset metadata,
  and the release procedure for both deposits.
- A data-availability section in `README.md`, and the Zenodo DOI for the raw images.

### Changed

- `CITATION.cff` gained `date-released`, `url`, `abstract` and `keywords`; `.zenodo.json`
  gained `access_right`, `language`, and a cross-reference to the dataset record.
- The regression suite's subprocess timeout is budgeted for a slow CI runner rather than a
  developer laptop, and is overridable via `CELLQUANT_TEST_TIMEOUT`.

### Fixed

- The Zenodo uploader resolves the data root from the repository root rather than from its
  own directory, which broke when it moved into `tools/`.

## [1.1.0] — 2026-09-01

The release accompanying the manuscript. Adds a native 3D pipeline, user-defined
compartments, and a set of guardrails that refuse a run rather than produce a
quietly-wrong number.

### Added

- **Native 3D pipeline.** A folder of z-stacks now runs in 3D end to end — segmentation,
  puncta detection, colocalization, and proximity — with areas replaced by volumes. Mode is
  detected from image shape and reported at startup; `--mode {auto,2d,3d}` forces it, and
  `--project-z {max,sum,mean}` flattens stacks instead.
- **Voxel-size handling.** Read from OME-TIFF metadata, or supplied with
  `--voxel-size XY_UM Z_UM`. The size a run actually used is reported and recorded in
  `provenance.json`.
- **User-defined compartments by set algebra**, e.g.
  `--compartment "perinuc = cell - nucleolus"`, usable anywhere a region is accepted.
  Terms can be grown or shrunk in microns (`nucleolus~0.3`). Each region adds a size column
  and a per-channel mean column to `cells.csv`.
- **Multiformat input** via bioio — ND2, CZI, LIF, OME-TIFF alongside plain TIFF.
- **3D nucleolar morphology**: volume and equivalent diameter.
- **`--version` flag.**
- **Continuous integration**: the test suite runs on Linux and macOS across Python 3.11
  and 3.12.
- **`CONTRIBUTING.md`**, issue templates, and a pull-request template.
- **cellquant-vs-CellProfiler benchmark harness** (`benchmark/`) with centroid matching,
  correlation, and Bland-Altman agreement, plus a Dockerized headless run.
- **3D validation harness** (`validation_3d/`) — the reproducibility record for the
  2D-vs-3D comparison and the manuscript figure and table scripts, including the
  exploratory diagnostics that document dead ends.

### Changed

- **Colocalization on 2D projections is now refused by default.** Collapsing Z inflates
  apparent overlap. `--allow-2d-colocalization` forces it and stamps the output
  `projection_derived=True`. In 3D it runs natively with no gate.
- **A 3D run with no known voxel size now aborts** rather than assuming 1 µm cubes. Every
  volume, every micron distance, and the anisotropy correction on the puncta filter depend
  on it.
- **Micron-denominated compartment erosion aborts** on images with no known pixel size,
  rather than silently applying the distance in pixels.
- Yeast 3D segmentation defaults to stitched slices at threshold 0.65 — full 3D gave no
  roundness gain at 2.3× the cost.
- The 3D cell-shape warning was recalibrated to catch the oblate signature of a wrong
  voxel size.
- Masks are written as compressed TIFFs.
- Figure 2B and Tables S1–S3 were regenerated under a single statistical convention.
- `pytest` collection is scoped to `tests/`.

### Fixed

- Lateral pixel size is now read on the 2D path; micron-denominated parameters were
  previously interpreted as pixels.
- The user-supplied voxel-size rule is applied on the TIFF path as well.
- The `(1.0, 1.0)` voxel-size sentinel no longer masquerades as a real measurement.
- `compute_nucleolar_proximity` was vectorized into a single pass over puncta voxels.
- The validation harness errors on an empty input glob instead of reporting success.
- The numpy version check moved ahead of the cellpose import, so numpy 2.x is caught at
  startup rather than as a segfault.

### Removed

- 3D nucleolar sphericity and surface area. They are not reliably measurable at typical
  axial sampling; volume and equivalent diameter are reported instead. The 2D shape
  descriptors (solidity, circularity, eccentricity) are retained but describe the outline
  of a projected shadow and have no 3D counterpart.

## Earlier history (2026-02-22 – 2026-03-09)

v1.1.0 is the first tagged release; the work below was never cut as a separate version, and
is recorded here for context.

The original single-script 2D pipeline with Cellpose segmentation, puncta
detection, per-cell intensity and condensation metrics, colocalization, nucleolar proximity
and morphology, QC overlays, superplots with replicate-level statistics, and Prism-ready
CSV export — plus the docs, tutorials, and example datasets.

Notable work in this line: `--reuse-masks` to skip segmentation on re-runs,
`--reference-condition` for pairwise p-values against a control, replicate-level Wilcoxon
p-values on two-condition superplots, the `numpy<2.0` and `opencv-python-headless<4.10`
pins with startup compatibility checks, and a fix for the macOS duplicate-libomp crash.

[1.1.1]: https://github.com/davidpincus/cellquant/releases/tag/v1.1.1
[1.1.0]: https://github.com/davidpincus/cellquant/releases/tag/v1.1.0
