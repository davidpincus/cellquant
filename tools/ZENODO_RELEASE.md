# Archiving to Zenodo

## Current state

| | Record | DOI | Status |
|---|---|---|---|
| **Data** — 99 raw z-stacks, 52.94 GB | [21843280](https://zenodo.org/records/21843280) | `10.5281/zenodo.21843280` (version)<br>`10.5281/zenodo.21829810` (concept) | **Published** 2026-08-07, version 3 of 3 |
| **Software** — a snapshot of this repository | none | none | **Does not exist yet** |

**The published record is mistyped.** Record 21843280 contains only raw image data — 99
TIFFs, verified byte-for-byte complete against the local directories — but its metadata
describes the *software*: resource type `Software`, title "cellquant: a single-script
pipeline for quantitative fluorescence microscopy", MIT license, and a description of a
pipeline that is not among its files. The code has never been archived at all.

Files on a published Zenodo record are frozen permanently. **Metadata is not** — the type,
title, description, license and related identifiers can all still be corrected in place,
without minting a new DOI or disturbing the existing one.

---

## 1. Correct the metadata on record 21843280

Edit at **<https://zenodo.org/uploads/21843280>** — that is the edit form for a published
record, reached from the record page's "Edit" button. (`/records/<id>/edit` is not a Zenodo
URL and returns 404.) Use [`zenodo_dataset_metadata.json`](zenodo_dataset_metadata.json) as
the source. The fields that are wrong today:

| Field | Currently | Should be |
|---|---|---|
| Resource type | Software | **Dataset** |
| Title | "cellquant: a single-script pipeline…" | The dataset title in the JSON |
| Description | Describes the pipeline | Describes the images |
| License | MIT | **CC-BY-4.0** — MIT is a software licence and does not fit image data |
| Keywords | Software keywords | The dataset keywords in the JSON |
| Related identifiers | GitHub URL only | Add the software DOI once step 2 mints it |

Do **not** re-upload anything. The file list is already correct and complete: 86 yeast
z-stacks (16/12/19/16/23 at 25/30/32/36/40 °C), 6 HCT116 z-stacks, and 7 MIPs.

Verify at any time with:

```bash
export ZENODO_TOKEN=...
python tools/zenodo_upload.py --deposition 21843280 --dry-run
```

## 2. Archive the software

The GitHub–Zenodo integration has never been enabled, which is why no code archive exists
and why `/.zenodo.json` has so far had no effect.

1. Sign in at <https://zenodo.org/account/settings/github/> with the GitHub account that
   owns the repository.
2. Find `davidpincus/cellquant` and flip the toggle **On**. Use "Sync now" if it is not
   listed — Zenodo only shows repositories you have admin rights to.
3. Zenodo archives on *release publication*, and only for releases published after the
   toggle was switched on. The existing `v1.1.0` release predates it, so it will never be
   archived. Publish a **new** release (`v1.1.1`) rather than deleting and re-creating
   `v1.1.0` — nothing is destroyed and no tag has to be moved.
4. Check the resulting record against `/.zenodo.json` — type Software, MIT, the four
   creators, and both related identifiers.

`/.zenodo.json` already points at the data via the **concept** DOI
`10.5281/zenodo.21829810`, which always resolves to the newest version, so that link stays
correct as the dataset gains versions.

## 3. Close the loop

- [ ] Add the software DOI to record 21843280's related identifiers
      (relation `isSupplementedBy`).
- [ ] Add the DOI badge to `README.md` and the `doi:` field to `CITATION.cff`.
- [ ] Replace "Submission and DOI pending" in the `README.md` citation block.
- [ ] Put the real DOI in the README "Data availability" section.
- [ ] Once the paper is out, add its DOI to both records with relation `isSupplementTo`.

## Notes

- Cite **concept** DOIs in the manuscript, not version DOIs — they resolve to the latest
  version and will not go stale.
- Versions 1 and 2 of the dataset record keep their own metadata; correcting version 3 does
  not rewrite them. The concept DOI resolves to version 3, so this is usually not worth
  chasing.
- `ZENODO_TOKEN` needs the `deposit:write` and `deposit:actions` scopes. Keep it in the
  environment; never on a command line, never committed.
