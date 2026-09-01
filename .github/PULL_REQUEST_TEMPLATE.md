## What this changes

<!-- One or two sentences. What does a user get that they did not have before? -->

## Why

<!-- The problem this solves. Link an issue if there is one. -->

## Does this change a measurement?

<!--
Answer honestly: does this alter a threshold, default, formula, or segmentation
parameter such that the same images could produce different numbers?

If YES, show the before/after on the example data — published results depend on these.
If NO, say "No — docs/infrastructure/refactor only."
-->

## Checklist

- [ ] `pytest` passes locally
- [ ] `docs/CLI_REFERENCE.md` updated if flags or behavior changed
- [ ] New flags have a `DEFAULTS` entry and a `parse_args` entry
- [ ] Considered whether any cell-type preset should override the new default
- [ ] Error messages name the flag to set and the value seen
- [ ] Still a single file — no new modules or package structure
