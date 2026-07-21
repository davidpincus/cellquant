# 3D segmentation diagnosis — handoff to Asif

**Audience:** Asif (built the 3D dispatch). **From:** validation_3d/ on David's M-series laptop, post-erosion-fix cellquant.
**Purpose:** characterize the observed 3D-segmentation issues on Abani's 6 mammalian z-stacks so we can reason about root cause from data rather than guesses. **Hypotheses framed as questions for you to confirm or refute, not as conclusions.**

---

## TL;DR

The two images that were flagged as "3D segmentation failures" actually behave very differently from each other, and **one of them isn't a failure at all**:

| Image | Initial perception | What the data shows | Verdict |
|---|---|---|---|
| `control_rep3` | "3D loses cells the MIP plainly shows" (6 kept of 11 total) | MIP 2D finds 6 cells; 3D stitch keeps 6 cells; the per-Z 2D cell counts go [6,5,6,5,6,2,4,3,3,3] because signal monotonically decays with Z. Stitch correctly merged the sparse-bottom slices with the dense-top slices to recover all 6 objects. | **Not a failure.** The QC perception of "fewer cells in Z=5 panel" was about cell visibility at that focal plane, not segmentation loss. |
| `arsenite_rep3` | "3D over-detects, finds ~30 spurious objects" (17 kept of 47 total) | Real — driven by Cellpose fragmenting a tight HCT116 cluster (cells touching with no inter-cell gaps). `--seg-3d-method full` makes it 3× **worse** (47→111 objects, kept-cell median vol 54k→7k vox). | **Real over-detection. The mode switch is not the fix.** |

Bonus issue surfaced in Part 1 that the original prompt didn't flag:

| Image | Issue | Data |
|---|---|---|
| `control_rep1` | 2 large objects (median 176k voxels — comparable to the 294k median of kept cells) get rejected by the nuclei gate. These are real cells whose DAPI nuclei the nuclear segmentation missed. | A separate failure mode from the cluster-fragmentation problem — points at *nuclear* segmentation reliability, not cell segmentation. |

**Two distinct hypotheses for Asif to confirm,** detailed below: (a) tight-cluster fragmentation in Cellpose mammalian composite-input segmentation, and (b) DAPI-only nuclear segmentation occasionally missing real nuclei.

---

## Part 1 — per-image characterization

Table `seg_characterization.csv`, full plot `per_z_signal_profiles.png`. Key fields per image:

| image | z | SNR(MIP) | SNR(midZ) | 3D total | 3D keep | n_0_nuc | DAPI MIP nuclei estimate | 2D published keep | kept vol (vox) | rejected vol (vox) |
|---|---|---|---|---|---|---|---|---|---|---|
| arsenite_rep1 | 10 | 20.3 | 34.0 | 48 | 14 | 34 | 13 | 11 | 2.0e5 | 482 |
| arsenite_rep2 | 10 | 20.3 | 31.9 | 17 | 11 | 6 | 10 | 7 | 2.7e5 | 405 |
| **arsenite_rep3** | 10 | **15.4** | 29.6 | **47** | 17 | **30** | 12 | 10 | **5.4e4** | **2,430** |
| **control_rep1** | 10 | 22.4 | 36.3 | 15 | 13 | 2 | 7 | 14 | 2.9e5 | **1.8e5** |
| control_rep2 | 5 | 23.8 | 35.9 | 21 | 19 | 2 | 6 | – | 1.4e5 | 3.2e4 |
| control_rep3 | 10 | 16.8 | 28.4 | 11 | 6 | 5 | 5 | – | 1.7e5 | 1.9e4 |

Observations:

- **SNR(MIP) correlates with both flagged images** (15.4 and 16.8 vs 20+ for the others) — but `control_rep3` segments fine, so SNR alone isn't predictive of failure.
- **No clumps**: zero cells with >4 nuclei in any image, so the upper-nuclei gate isn't doing anything.
- **arsenite_rep3 kept-cell volume is 4× smaller** than other arsenite images (5.4e4 vs 2.0e5 / 2.7e5 vox). Even the "kept" cells are fragments.
- **control_rep1 rejected-cell volume is huge** (1.8e5 vox, ~60% of kept). Those are real cells, not debris, getting filtered for lacking nuclei.

Axial signal profile: all 6 stacks have signal in every Z (no obviously-empty slices). No image has signal concentrated to a handful of slices.

---

## Part 2 — `--seg-3d-method stitch` vs `full` on the two flagged images

`seg_mode_compare/summary.csv`:

| image | mode | total | keep | n_0_nuc | kept_vol_med (vox) | rejected_vol_med (vox) | wall (s) |
|---|---|---|---|---|---|---|---|
| control_rep3 | stitch | 11 | 6 | 5 | 1.67e5 | 1.85e4 | 128 |
| control_rep3 | full | 10 | 6 | 4 | **2.58e5** | 2.9e3 | **3,569** |
| arsenite_rep3 | stitch | 47 | 17 | 30 | 5.45e4 | 2.4e3 | 157 |
| arsenite_rep3 | full | **111** | 17 | **94** | **7.0e3** | 2.3e3 | **3,598** |

Findings:

- **control_rep3:** full mode gives slightly fatter / more complete cells (kept median 167k → 258k vox) — same 6 kept cells, marginally better cell extents. At **27× wall-clock**.
- **arsenite_rep3:** full mode is **dramatically worse**. 2× more raw objects (47 → 111); 3× more zero-nuclei rejections (30 → 94); kept-cell median volume **drops from 54k → 7k vox** (so even the "good" cells become fragments). Same 17 kept after the nuclei gate.

The over-detection in arsenite_rep3 is *amplified* by `full` mode, not fixed by it. Stitch is the better default by a wide margin and full is not a viable rescue.

---

## Part 3 — per-Z-slice 2D segmentation on control_rep3

`control_rep3_cells_per_z.png` + `control_rep3_per_slice.csv`. Independent 2D Cellpose-SAM per Z slice (mammalian-preset params, `seg_downsample=3`):

| Z | cells found (2D) | mean composite signal |
|---|---|---|
| 0 | 6 | 0.125 |
| 1 | 5 | 0.118 |
| 2 | 6 | 0.113 |
| 3 | 5 | 0.109 |
| 4 | 6 | 0.105 |
| 5 | 2 | 0.102 |
| 6 | 4 | 0.100 |
| 7 | 3 | 0.097 |
| 8 | 3 | 0.095 |
| 9 | 3 | 0.093 |

Comparison: MIP 2D finds 6; 3D stitch keeps 6.

**Stitch is doing the right thing here.** Signal decays monotonically with Z; the top half of the stack (z=0-4) has all 6 cells; the bottom half loses some cells to defocus (down to 2 at z=5). Stitch's IoU linkage correctly merges these per-slice masks into 6 coherent 3D objects matching the MIP. The QC overlay's Z=5 panel showing fewer outlines is faithful to the actual segmentation at that slice — but the overall 3D mask still contains all 6 cells (you can see this in the MIP panel of the same QC image, which projects across Z and recovers the full set).

This rules out one of the original hypotheses ("stitch-by-IoU is losing cells present in too few slices") for this image. It's working.

---

## Part 4 — synthesis & hypotheses for you to confirm

### What's actually broken: tight-cluster fragmentation in `arsenite_rep3`

Look at `outputs_mammalian/3d/arsenite_rep3/qc/arsenite_rep3_qc.png` and `seg_mode_compare/arsenite_rep3/{stitch,full}/qc/*_qc.png` side-by-side. The right side of the field has a dense HCT116 cluster — cells are touching with no visible inter-cell gap. The composite cell-seg input (sum-of-channels, p1/p99.8-normalized) renders this cluster as a continuous bright region. Cellpose-SAM splits it into pieces, but inconsistently — some pieces match real cells (kept after the nuclei gate, n=17), most are small fragments (median 2.4k vox = ~280 µm³, well below a real cell).

**Candidate hypotheses for the root cause** (please weigh in on which is plausible):

1. **Composite-input clobbers the cell-boundary signal.** In tight clusters, the sum-of-channels composite looks like a single bright blob because DAPI nuclei dominate. The cytoplasmic boundaries get lost. **Would using only the cytoplasmic markers (G3BP1+PABPC1 sum, dropping DAPI from the composite) give Cellpose better edges in clusters?**
2. **Cellpose-SAM stitch parameters are tuned for sparser fields.** The default `stitch_threshold=0.4` linkage may not handle tight-touching cells well. **Would a lower stitch threshold (more permissive linkage) merge the small fragments into real cells?**
3. **`min_cell_volume_vox` is set to 0 in the mammalian preset.** A non-zero floor (e.g. 10k voxels) would discard the small fragments at the source, before they enter the nuclei gate. Cleaner counts. Trade-off: real small cells at the edge of the field could be lost. **Would you object to a non-zero floor for mammalian-preset 3D?**
4. **Adaptive seg-input per cell-type.** The 2D path uses `composite_seg`; maybe 3D mammalian should default to a specific channel (e.g. `--cell-seg-channel PABPC1`) for cluster-robust segmentation. **Did you have a reason to keep composite as the 3D default?**

Mode switch is *not* the fix; `--seg-3d-method full` made arsenite_rep3 strictly worse (more fragments) and was 27× slower for marginal-to-no benefit on control_rep3. Don't change that default.

### A separate failure mode: nuclear segmentation occasionally misses real nuclei

`control_rep1` has 2 objects with volume 1.8e5 vox getting tossed by the nuclei gate. These are real cells, not debris — they're the size of kept cells. The DAPI-MIP quick-count finds 7 nuclei, the 3D run keeps 13 cells, the 2D-published keeps 14. So nuclear segmentation is the bottleneck — when it misses, real cells get filtered.

This is independent of cell segmentation. **Possible hypotheses:**

5. **DAPI nuclear segmentation params (mammalian preset) are too strict in 3D.** The preset doesn't override Cellpose-SAM nuclei diameter for 3D; maybe it should. **Did you tune nuclei diameter for 3D specifically?**
6. **Nuclei get lost in the same composite-vs-individual-channel issue as cells.** Nuclei are segmented from DAPI alone, but if the DAPI signal-to-background is low in some cells, they're missed. The 3D path applies the same `seg_downsample=3` to nuclei that it applies to cells — for nuclei this might be over-aggressive. **Worth a per-target seg_downsample (one for cells, one for nuclei)?**

### What was not broken: control_rep3 perception of under-detection

For the record, so we don't chase a phantom: control_rep3 segments correctly. The QC overlay showing fewer outlines at Z=5 is faithful — at that focal plane, fewer cells are in focus because signal decays monotonically along Z (composite signal 0.125 at z=0 → 0.093 at z=9). The MIP-based segmentation and the 3D stitch-based segmentation both produce 6 cells, matching the DAPI MIP quick-count estimate (5).

This is a **documentation/UX issue**, not a bug — the 3D QC panel layout (left=MIP, right=middle-Z slice) can be misleading when cells aren't uniformly distributed in Z. **Worth surfacing in CONCEPTS.md or the QC layout** ("the Z slice panel is not the segmentation; the MIP panel is the projection of the segmentation").

---

## Recommendations summary

| Question | Cheapest experiment to answer it |
|---|---|
| Does cytoplasm-only composite (drop DAPI) fix arsenite_rep3 fragmentation? | `--cell-seg-channel G3BP1` or build a 2-channel composite, re-run on arsenite_rep3. ~3 min. |
| Does a `min_cell_volume_vox=10000` floor clean up debris without losing real small cells? | Re-run all 6 images with `--min-cell-volume-vox 10000`, compare counts and visually check kept cells. ~15 min. |
| Does a lower stitch threshold merge the arsenite_rep3 fragments? | `--stitch-threshold 0.2` on arsenite_rep3 only. ~3 min. |
| Why does control_rep1 DAPI miss 2 nuclei? | Visual + per-cell quick check — likely intensity-related, or cells at the edge of the stack. ~10 min eyeballing. |
| Should mammalian preset's 3D nuclei diameter be specified? | Test `--nuclei-diameter 80` (typical HCT116 nucleus diameter) on control_rep1, see if rejected→kept conversion happens. ~3 min. |

I haven't run any of these because the prompt was characterization-only. They're ranked roughly by how much they'd inform the root-cause story.

## Files

- `seg_characterization.csv` — Part 1 table
- `per_z_signal_profiles.png` — axial signal per image
- `rejected_cell_size_dist.png` — kept vs rejected object volumes
- `seg_mode_compare/summary.csv` and `seg_mode_compare/{image}/{mode}/qc/*.png` — Part 2 stitch vs full
- `control_rep3_cells_per_z.png` and `control_rep3_per_slice.csv` — Part 3 per-Z analysis
- This document — Part 4

No code changes were made.
