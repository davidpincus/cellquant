#!/usr/bin/env bash
# End-to-end cellquant vs CellProfiler agreement benchmark (HCT116 stress granules).
#
# Runs BOTH tools on the SAME images, then the comparison harness. Fails loud at
# every step (missing/empty inputs, no tool output, no matched cells). Nothing is
# tuned to improve agreement — parameters are fixed in build_cppipe.py and cellquant's
# mammalian preset.
#
# Usage:
#   ./run_benchmark.sh                      # uses the defaults below
#   HCT116_DIR=/path/to/MIPs ./run_benchmark.sh
#   OUT_ROOT=/tmp/bench ./run_benchmark.sh
#
# Configurable via environment:
#   HCT116_DIR    full-res 3-channel MAX_*.tif dir (1192x1200)   [default: <repo>/SG_zstacks/MIPs]
#   OUT_ROOT      where out_cq/ out_cp/ cp_input/ compare_out/ go [default: <this dir>]
#   CELLQUANT_PY  python with cellquant's deps (tifffile/numpy)   [default: miniforge cellquant env]
#   COMPARE_PY    python for compare_tools.py (numpy/scipy/pandas/mpl/skimage/tifffile) [default: CELLQUANT_PY]
#   CP_IMAGE      CellProfiler docker image                        [default: cellprofiler/cellprofiler:4.2.6]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

HCT116_DIR="${HCT116_DIR:-$REPO/SG_zstacks/MIPs}"
OUT_ROOT="${OUT_ROOT:-$HERE}"
CELLQUANT_PY="${CELLQUANT_PY:-/Users/dpincus/miniforge3/envs/cellquant/bin/python}"
COMPARE_PY="${COMPARE_PY:-$CELLQUANT_PY}"
CP_IMAGE="${CP_IMAGE:-cellprofiler/cellprofiler:4.2.6}"
CHANNELS=("DAPI" "G3BP1" "PABPC1")   # channel order 1,2,3 (paper's validated invocation)

OUT_CQ="$OUT_ROOT/out_cq"
OUT_CP="$OUT_ROOT/out_cp"
CP_IN="$OUT_ROOT/cp_input"
CMP="$OUT_ROOT/compare_out"
CPPIPE="$HERE/cellprofiler_hct116.cppipe"

die(){ echo "FATAL: $*" >&2; exit 1; }

echo "== preflight =="
command -v docker >/dev/null || die "docker not found"
docker image inspect "$CP_IMAGE" >/dev/null 2>&1 || die "docker image '$CP_IMAGE' not present (docker pull $CP_IMAGE)"
[ -x "$CELLQUANT_PY" ] || die "CELLQUANT_PY not executable: $CELLQUANT_PY"
[ -d "$HCT116_DIR" ]   || die "HCT116_DIR not a dir: $HCT116_DIR"
n_src=$(find "$HCT116_DIR" -maxdepth 1 -name 'MAX_*.tif' | wc -l | tr -d ' ')
[ "$n_src" -gt 0 ] || die "no MAX_*.tif in HCT116_DIR: $HCT116_DIR"
echo "  HCT116_DIR = $HCT116_DIR ($n_src images)"
echo "  OUT_ROOT   = $OUT_ROOT"

echo "== [1/5] cellquant (fresh) =="
rm -rf "$OUT_CQ"; mkdir -p "$OUT_CQ"
"$CELLQUANT_PY" "$REPO/cellquant.py" "$HCT116_DIR" \
    "1:DAPI:nucleus" "2:G3BP1:quantify" "3:PABPC1:quantify" \
    --cell-type mammalian --out "$OUT_CQ" \
    --puncta-channels G3BP1 \
    --filename-pattern "MAX_{condition}_rep{replicate}" \
    > "$OUT_CQ/cellquant_run.log" 2>&1 || die "cellquant failed (see $OUT_CQ/cellquant_run.log)"
[ -s "$OUT_CQ/cells.csv" ] || die "cellquant produced no cells.csv"
echo "  cells.csv rows: $(($(wc -l < "$OUT_CQ/cells.csv") - 1))"

echo "== [2/5] split channels for CellProfiler =="
rm -rf "$CP_IN"; mkdir -p "$CP_IN"
CP_IN="$CP_IN" HCT116_DIR="$HCT116_DIR" "$CELLQUANT_PY" - <<'PY'
import os, tifffile, numpy as np, glob
src=os.environ["HCT116_DIR"]; dst=os.environ["CP_IN"]; chans=["DAPI","G3BP1","PABPC1"]
fs=sorted(glob.glob(os.path.join(src,"MAX_*.tif")))
assert fs, "no MAX_*.tif to split"
for f in fs:
    a=tifffile.imread(f)
    assert a.ndim==3 and a.shape[0]==3, f"{f}: expected (3,Y,X), got {a.shape}"
    base=os.path.splitext(os.path.basename(f))[0]
    for ci,ch in enumerate(chans):
        out=os.path.join(dst,f"{base}_{ch}.tif")
        tifffile.imwrite(out, a[ci], photometric='minisblack')
        assert np.array_equal(tifffile.imread(out), a[ci]), out  # byte-identical slice
print(f"  split {len(fs)} images -> {len(fs)*3} single-channel TIFFs")
PY

echo "== [3/5] (re)build the CellProfiler pipeline in the container =="
docker run --rm -v "$HERE":/pipe --entrypoint python "$CP_IMAGE" \
    /pipe/build_cppipe.py -o /pipe/cellprofiler_hct116.cppipe 2>&1 | grep -E "saved|RELOAD" \
    || die "cppipe build failed"
[ -s "$CPPIPE" ] || die "no cppipe produced"

echo "== [4/5] CellProfiler (headless, docker) =="
rm -rf "$OUT_CP"; mkdir -p "$OUT_CP"   # ExportToSpreadsheet refuses to overwrite; start clean
docker run --rm \
    -v "$CP_IN":/in -v "$OUT_CP":/out -v "$HERE":/pipe \
    "$CP_IMAGE" -c -r -p /pipe/cellprofiler_hct116.cppipe -i /in -o /out \
    > "$OUT_CP/cellprofiler_run.log" 2>&1 || die "CellProfiler failed (see $OUT_CP/cellprofiler_run.log)"
for f in cp_Cells.csv cp_Puncta.csv cp_Image.csv; do
    [ -s "$OUT_CP/$f" ] || die "CellProfiler produced no $f"
done
echo "  cp_Cells rows: $(($(wc -l < "$OUT_CP/cp_Cells.csv") - 1)), cp_Puncta rows: $(($(wc -l < "$OUT_CP/cp_Puncta.csv") - 1))"

echo "== [5/5] compare =="
rm -rf "$CMP"; mkdir -p "$CMP"
"$COMPARE_PY" "$HERE/compare_tools.py" \
    --cellquant-cells "$OUT_CQ/cells.csv" \
    --cellquant-masks "$OUT_CQ/masks" \
    --cp-cells "$OUT_CP/cp_Cells.csv" \
    --cp-puncta "$OUT_CP/cp_Puncta.csv" \
    --cp-image "$OUT_CP/cp_Image.csv" \
    --out-dir "$CMP" || die "compare_tools.py failed"

echo "== DONE =="
echo "  results: $CMP/{benchmark_summary.csv,benchmark_summary.json,matches.csv,agreement_combined.pdf}"
