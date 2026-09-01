# Zenodo archives

## What exists

| | Record | DOI | Contents |
|---|---|---|---|
| **Software** | [22239479](https://zenodo.org/records/22239479) | `10.5281/zenodo.22239479` | `cellquant-v1.1.1.zip` — a snapshot of this repository |
| **Data** | [21843280](https://zenodo.org/records/21843280) | `10.5281/zenodo.21843280` | 99 raw z-stacks, 52.94 GB |

Both are published and correctly described: the software record is typed Software under MIT,
the data record is typed Dataset under CC-BY-4.0.

## Known wrinkle: the two records share a concept DOI

Zenodo groups versions of a record under a **concept DOI** that resolves to the newest
version. These two records are versions of the *same* lineage:

```
concept 10.5281/zenodo.21829810
  ├── … earlier versions
  ├── index 2 → 21843280   the raw image data   (Dataset, 99 files)
  └── index 3 → 22239479   the software v1.1.1  (Software, 1 zip)   ← newest
```

This happened because the image data was uploaded in August as a *new version* of a
pre-existing software deposit rather than as an independent record — which is also why it
arrived carrying software metadata that had to be corrected later.

**Consequence: never cite `10.5281/zenodo.21829810`.** It resolves to whichever of the two
was published most recently — currently the code, not the images. Every reference in this
repository deliberately uses a **version** DOI instead. Normally the advice is the reverse,
so this is worth remembering when writing the manuscript's data-availability statement.

Separating them would mean re-uploading 52.9 GB into a fresh record, roughly three hours,
and would mint a new data DOI. Citing version DOIs avoids the problem entirely, so this is
recorded rather than fixed.

Note that each future release published through the GitHub integration adds another version
to this same lineage, pushing the concept DOI further from the data.

## Cutting a new release

1. Bump `__version__` in `cellquant.py`, and `version` in `CITATION.cff` and `.zenodo.json`.
2. Add a `CHANGELOG.md` entry. Separate changes that alter a measurement from those that do
   not — that distinction is what a reader most needs.
3. Tag and push: `git tag -a vX.Y.Z -m "…" && git push origin vX.Y.Z`
4. Publish a GitHub release from that tag. The Zenodo integration archives on *publication*,
   so a draft does not trigger it.
5. Zenodo reads `/.zenodo.json` at that moment. Verify the new record afterwards.
6. Update the DOI badge in `README.md` and `doi:`/`identifiers:` in `CITATION.cff` to the
   new version DOI.

## Re-uploading data

`zenodo_upload.py` is resumable and checksum-verified: it lists what a draft already holds
and skips those files.

```bash
export ZENODO_TOKEN=...                                    # deposit:write, deposit:actions
python tools/zenodo_upload.py --deposition <id> --dry-run  # compare local vs stored
python tools/zenodo_upload.py --deposition <id>            # upload the difference
```

It cannot publish, deliberately — minting a DOI stays a decision made in the browser.

## Remaining

- [ ] Add `10.5281/zenodo.22239479` to record 21843280's related identifiers, relation
      `isSupplementedBy`. Metadata on a published record is editable at
      <https://zenodo.org/uploads/21843280>; only files are frozen.
- [ ] Once the manuscript is published, add its DOI to both records with relation
      `isSupplementTo`, and replace "Manuscript submission pending" in `README.md`.
