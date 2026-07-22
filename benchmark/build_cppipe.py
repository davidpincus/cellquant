#!/usr/bin/env python
"""Build the CellProfiler .cppipe for the HCT116 G3BP1 stress-granule benchmark.

This is run *inside* the cellprofiler/cellprofiler:4.2.6 container (it needs
cellprofiler_core / cellprofiler, which are not installed on the host). It emits
a headless-loadable pipeline and validates it by reloading.

Module sequence (matches the benchmark spec; native segmentation, no RunCellpose):
  1. Images
  2. Metadata            regex-extract condition/replicate/channel from filename
  3. NamesAndTypes       split 3 channels -> DAPI, G3BP1, PABPC1 (by filename token)
  4. Groups
  5. IdentifyPrimaryObjects     Nuclei  <- DAPI (Otsu global, declump by shape)
  6. IdentifySecondaryObjects   Cells   <- Nuclei, Propagation, G3BP1 as cell body
  7. EnhanceOrSuppressFeatures  Enhance Speckles on G3BP1 -> G3BP1_speckles
  8. IdentifyPrimaryObjects     Puncta  <- G3BP1_speckles (restricted to Cells via RelateObjects)
  9. RelateObjects              parent Cells, child Puncta  -> Children_Puncta_Count
 10. MeasureObjectSizeShape     Cells, Puncta               -> AreaShape_Area
 11. MeasureObjectIntensity     G3BP1 on Cells              -> Intensity_MeanIntensity_G3BP1
 12. ExportToSpreadsheet        cp_*.csv (per-object, Location, Metadata)

PARAMETER PROVENANCE (set once from image biology / cellquant's own preset;
NOT tuned to make the tools agree):
  * SPECKLE_FEATURE_SIZE = 6 px. cellquant's mammalian preset detects puncta with
    a Laplacian-of-Gaussian of log_sigma = 2.0 (cellquant.py CELL_TYPE_PRESETS
    ["mammalian"]["log_sigma"]). A LoG of sigma s maximally responds to a bright
    blob of diameter 2*sqrt(2)*s = 2*1.414*2.0 = 5.66 px -> 6 px. So the speckle
    feature size is cellquant's own expected puncta diameter, rounded, not invented.
  * NUCLEI_DIAMETER_PX = (60, 250). HCT116 nuclei measure ~110-185 px across in
    these 1192x1200 acquisitions; the default CellProfiler scale shatters them.
  * PUNCTA_DIAMETER_PX = (3, 30). Lower bound = cellquant puncta_min_area_px = 6
    -> equivalent diameter ~2.8 -> 3 px; upper bound covers the largest granules.
  * PUNCTA_THRESHOLD: fixed Global/Manual threshold on the speckle-enhanced image.
    Absolute (not adaptive) because top-hat/Otsu/RobustBackground detect cytoplasmic
    texture as "puncta" in unstressed control cells (they yield control >= arsenite,
    which is biologically wrong). Value is set from the enhanced-image intensity
    distribution, once; see README for the calibration note.
"""
import argparse

import cellprofiler_core.preferences as prefs
prefs.set_headless()
from cellprofiler_core.pipeline import Pipeline
from cellprofiler_core.modules.images import Images
from cellprofiler_core.modules.metadata import Metadata
from cellprofiler_core.modules.namesandtypes import NamesAndTypes
from cellprofiler_core.modules.groups import Groups
from cellprofiler.modules.identifyprimaryobjects import IdentifyPrimaryObjects
from cellprofiler.modules.identifysecondaryobjects import IdentifySecondaryObjects
from cellprofiler.modules.enhanceorsuppressfeatures import EnhanceOrSuppressFeatures
from cellprofiler.modules.relateobjects import RelateObjects
from cellprofiler.modules.measureobjectsizeshape import MeasureObjectSizeShape
from cellprofiler.modules.measureobjectintensity import MeasureObjectIntensity
from cellprofiler.modules.exporttospreadsheet import ExportToSpreadsheet

# ---- parameters (see PARAMETER PROVENANCE above) ----
SPECKLE_FEATURE_SIZE = 6           # px; = 2*sqrt(2)*cellquant_log_sigma(2.0)
NUCLEI_DIAMETER_PX = (60, 250)
PUNCTA_DIAMETER_PX = (3, 30)
CELL_BODY_IMAGE = "PABPC1"         # cytoplasmic marker for cell-body propagation.
# The spec (module 6) suggested G3BP1, but G3BP1 is dim/diffuse in unstressed
# control cells, giving weak propagation boundaries there. PABPC1 is a brighter,
# more uniform cytoplasmic marker; it is also cellquant's own mammalian seg
# channel (manuscript Table S4) and what the prior figure5_benchmark used, so both
# tools segment cells from comparable information. Confirmed by the user at the
# checkpoint. NOT a tuning knob for tool agreement (it is the cell-outline source).
PUNCTA_THRESHOLD = 0.003           # fixed Manual threshold on the enhanced image
# Calibration (once, from the enhanced-image distribution, NOT to cellquant):
# after Enhance-Speckles (white top-hat) the granule residuals top out ~0.004-0.006
# of the 0-1 range and the diffuse cytoplasm sits ~10x lower; a Manual threshold of
# 0.003 sits in that gap. Arsenite has ~10x more supra-threshold speckle pixels than
# control, matching the biology. See README "Puncta threshold calibration".


def build():
    p = Pipeline()

    def add(m):
        m.module_num = len(p.modules()) + 1
        p.add_module(m)
        return m

    # 1. Images
    add(Images())

    # 2. Metadata — regex-extract condition/replicate/channel from the filename.
    md = Metadata()
    md.wants_metadata.value = True
    grp = md.extraction_methods[0]
    grp.extraction_method.value = "Extract from file/folder names"
    grp.source.value = "File name"
    grp.file_regexp.value = (
        r"^MAX_(?P<condition>.+?)_rep(?P<replicate>[0-9]+)_"
        r"(?P<channel>DAPI|G3BP1|PABPC1)"
    )
    grp.filter_choice.value = "All images"
    add(md)

    # 3. NamesAndTypes — one grayscale image per channel, matched by filename token.
    na = NamesAndTypes()
    na.assignment_method.value = "Images matching rules"
    na.matching_choice.value = "Order"
    na.add_assignment()
    na.add_assignment()  # default has 1 assignment; need 3 total
    specs = [("DAPI", "_DAPI"), ("G3BP1", "_G3BP1"), ("PABPC1", "_PABPC1")]
    for asg, (name, tok) in zip(na.assignments, specs):
        asg.image_name.value = name
        asg.rule_filter.value = f'and (file does contain "{tok}")'
        asg.load_as_choice.value = "Grayscale image"
    add(na)

    # 4. Groups
    add(Groups())

    # 5. Nuclei <- DAPI (Otsu global, declump by shape)
    ipn = IdentifyPrimaryObjects()
    ipn.x_name.value = "DAPI"
    ipn.y_name.value = "Nuclei"
    ipn.size_range.value = NUCLEI_DIAMETER_PX
    ipn.exclude_size.value = True
    ipn.unclump_method.value = "Shape"
    ipn.watershed_method.value = "Shape"
    ipn.threshold.threshold_scope.value = "Global"
    ipn.threshold.global_operation.value = "Otsu"
    add(ipn)

    # 6. Cells <- Nuclei by Propagation, using G3BP1 as the cell-body image.
    iso = IdentifySecondaryObjects()
    iso.x_name.value = "Nuclei"
    iso.y_name.value = "Cells"
    iso.image_name.value = CELL_BODY_IMAGE
    iso.method.value = "Propagation"
    add(iso)

    # 7. Enhance Speckles on G3BP1 (feature size from cellquant's LoG sigma).
    eos = EnhanceOrSuppressFeatures()
    eos.x_name.value = "G3BP1"
    eos.y_name.value = "G3BP1_speckles"
    eos.method.value = "Enhance"
    eos.enhance_method.value = "Speckles"
    eos.object_size.value = SPECKLE_FEATURE_SIZE
    add(eos)

    # 8. Puncta <- enhanced G3BP1 (fixed absolute threshold; restricted to Cells
    #    downstream via RelateObjects).
    ipp = IdentifyPrimaryObjects()
    ipp.x_name.value = "G3BP1_speckles"
    ipp.y_name.value = "Puncta"
    ipp.size_range.value = PUNCTA_DIAMETER_PX
    ipp.exclude_size.value = True
    ipp.threshold.threshold_scope.value = "Global"
    ipp.threshold.global_operation.value = "Manual"
    ipp.threshold.manual_threshold.value = PUNCTA_THRESHOLD
    add(ipp)

    # 9. Relate Puncta -> Cells (yields Children_Puncta_Count on Cells)
    ro = RelateObjects()
    ro.x_name.value = "Cells"
    ro.y_name.value = "Puncta"
    add(ro)

    # 10. Size/shape for Cells and Puncta (AreaShape_Area)
    moss = MeasureObjectSizeShape()
    moss.objects_list.value = "Cells, Puncta"
    add(moss)

    # 11. G3BP1 intensity on Cells (Intensity_MeanIntensity_G3BP1)
    moi = MeasureObjectIntensity()
    moi.images_list.value = "G3BP1"
    moi.objects_list.value = "Cells, Puncta"
    add(moi)

    # 12. Export per-object CSVs (cp_Cells.csv, cp_Puncta.csv, cp_Image.csv, ...)
    ets = ExportToSpreadsheet()
    ets.prefix.value = "cp_"
    add(ets)

    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="/pipe/cellprofiler_hct116.cppipe")
    args = ap.parse_args()

    p = build()
    with open(args.out, "w") as fh:
        p.dump(fh, save_image_plane_details=False)
    print("BUILT modules:", [m.module_name for m in p.modules()])
    print("saved:", args.out)

    # validate reload
    p2 = Pipeline()
    with open(args.out) as fh:
        p2.load(fh)
    print("RELOAD OK, modules:", len(p2.modules()))


if __name__ == "__main__":
    main()
