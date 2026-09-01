# Contributing to cellquant

Thanks for your interest. cellquant is a research tool that accompanies a published
manuscript, so contributions are welcome but the bar is "does this help a biologist measure
something correctly," not "does this make the code more elegant."

## Reporting a problem

Open an [issue](https://github.com/davidpincus/cellquant/issues). Please include the exact
command you ran, the full error as text, and the output of `python cellquant.py --version`.
The issue templates ask for these.

**Before filing a bug, check whether cellquant refused on purpose.** Several failures are
deliberate guardrails, not defects:

- A 3D run with no voxel size in the metadata and no `--voxel-size` aborts, rather than
  silently assuming 1 µm cubes.
- Colocalization on 2D projections is refused unless you pass `--allow-2d-colocalization`,
  because collapsing Z inflates apparent overlap.
- Micron-denominated compartment erosion aborts on images with no known pixel size.

[docs/PHILOSOPHY.md](docs/PHILOSOPHY.md) explains why refusing is the better default, and
[docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) covers the common cases.

## Asking how to analyze your data

Use the "Question about analyzing my data" issue template. cellquant is designed to be
driven with an AI assistant — describing your images to Claude or ChatGPT alongside
[docs/CLI_REFERENCE.md](docs/CLI_REFERENCE.md) is usually the fastest route to the right
command. Open an issue when that does not get you there, or when the answer should be
documented for everyone.

## Development setup

```bash
conda env create -f environment.yml
conda activate cellquant
pytest
```

The full suite takes about a minute. It runs the pipeline end to end on the cropped datasets
in `example_data/`.

## Making a change

1. **Keep it one file.** Everything lives in `cellquant.py`. That is deliberate — the target
   user shares and debugs a single file with an AI assistant. Please do not split it into a
   package.
2. **Match the surrounding style**: type hints, module-level functions, no classes for
   pipeline stages, and docstrings that explain *why* a threshold or default was chosen.
3. **Treat error messages as product surface.** Name the flag to set and the value that was
   seen. The user reading it is a biologist, not a developer.
4. **Document it.** A behavior change that is not reflected in `docs/CLI_REFERENCE.md` is
   incomplete. New flags need a `DEFAULTS` entry, a `parse_args` entry, and a reference
   entry — and consider whether any cell-type preset should override the default.
5. **Add a test** if you change what lands in an output CSV. The tests in `tests/` are
   structural: they assert that outputs exist with the expected columns and row counts, and
   deliberately do not assert exact numeric values, because Cellpose is not deterministic
   across platforms. Please keep new tests to that standard rather than pinning numbers.
6. **Run `pytest` before opening a PR.** CI runs the same suite on Linux and macOS across
   Python 3.11 and 3.12.

## Changing scientific behavior

If your change alters a measurement — a threshold, a default, a formula, a segmentation
parameter — say so explicitly in the PR and show the before/after on the example data. These
changes get more scrutiny than anything else, because published results depend on them. The
harness in `validation_3d/` is the model for how such a claim should be supported.

## Dependency pins

`numpy<2.0` and `opencv-python-headless<4.10` are load-bearing for Cellpose/PyTorch
compatibility; cellquant checks numpy at startup and exits if violated. Python is pinned to
3.11–3.12. Please do not relax these without demonstrating the whole suite passes.

## License

By contributing you agree that your contributions are licensed under the MIT License, the
same terms that cover the project.
