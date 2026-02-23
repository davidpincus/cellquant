# cellquant Beta Testing Guide

**Version:** Pre-submission draft
**Date:** February 2026
**Contact:** David Pincus (pincus@uchicago.edu)

---

## What is this?

`cellquant` is a tool for quantifying fluorescence microscopy images by performing cell segmentation, puncta/foci detection, colocalization, and statistics all from a single command. It's designed for biologists, not programmers.

We're preparing this for a peer-reviewed tutorial and want real world feedback before submitting.

## What we need from you

We need you to try to break it. The most valuable feedback is "I got stuck here" or "I didn't understand this step." Specifically:

1. **Can you install it?** Follow the installation guide from scratch. If you hit a wall, tell us exactly where.
2. **Can you run the example?** The quickstart tutorial should produce results in about 15 minutes. Did it?
3. **Can you run your own data?** If you have fluorescence images (aggregates, foci, stress granules), try it. What happened?
4. **Are the docs clear?** Where did you get lost? What would have helped?
5. **Are the outputs useful?** Do the plots, CSVs, and QC overlays give you what you'd actually want?

## Getting started

Start with the installation guide (`docs/INSTALL.md`). It walks through everything from the beginning. No prior experience with Python or the command line is needed.

**If you get stuck during installation, that is useful feedback.** Please note exactly where you stopped and what happened (a screenshot of any error message is ideal). Then email me and we'll get you unstuck.

Once installed, follow the quickstart tutorial (`docs/QUICKSTART.md`) to run your first analysis on the included example data.

If you want to go further, `docs/TUTORIAL_1_mammalian_SGs.md` and `docs/TUTORIAL_2_yeast_temp.md` walk through more detailed analyses. These are closer to what you'd do with your own images.

## Try your own images

This is the most helpful thing you can do. The pipeline works with multi-channel fluorescence images (TIFFs, single-plane or max projections). It should handle the kinds of images common in chaperone and protein quality control labs — aggregation assays, stress-induced foci, protein relocalization, etc.

If you're not sure how to set up a command for your images, email me or bring it to the Zoom call and we'll work through it together.

## Giving feedback

Send feedback however is easiest:

- Email: pincus@uchicago.edu
- GitHub issues on the repository

A few specific questions:

- What computer did you use? (Mac/PC/Linux, and roughly how old)
- How far did you get through installation before needing help (if at all)?
- Did the quickstart tutorial run to completion?
- Did you try your own data? What kind of images?
- What's the single biggest thing that would make this more useful for your lab?

## Quick troubleshooting

| Problem | What to do |
|---------|------------|
| Installation failed | Note the error, email me — this is exactly the feedback we need |
| Segmentation looks wrong | Check the QC overlay images in the output folder; try adjusting `--cell-diameter` |
| Script crashes | Copy the full error message and send it over |
| Not sure what command to run | See `docs/CLI_REFERENCE.md` or email me |

For more details, see `docs/TROUBLESHOOTING.md`.
