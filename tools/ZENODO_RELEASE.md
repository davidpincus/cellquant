# Archiving a release to Zenodo

Two separate deposits, and therefore two DOIs:

| | What | Metadata comes from | How it is created |
|---|---|---|---|
| **Data** | The raw z-stacks (99 files, ~52.9 GB) | `tools/zenodo_dataset_metadata.json`, applied by hand | `tools/zenodo_upload.py`, then published in the browser |
| **Software** | A snapshot of this repository | `/.zenodo.json`, applied automatically | The GitHub–Zenodo integration, fired by publishing a GitHub release |

**Do these in order.** `/.zenodo.json` already lists the data DOI
`10.5281/zenodo.21843280` as a related identifier, so publishing the data deposit first
means the software archive points at a DOI that resolves the moment it is minted.

---

## 1. The data deposit (Zenodo 21843280)

The files are already uploaded — the last run stored 89 files with 0 failures. What remains
is metadata and the decision to publish.

1. Open the draft: <https://zenodo.org/uploads/21843280>
2. Open `tools/zenodo_dataset_metadata.json` and copy each field into the matching box in
   the Zenodo edit form: title, upload type (Dataset), license (CC-BY-4.0), description,
   creators, keywords, related identifiers, notes. The `_comment` and `_note` keys are
   guidance for you and are not part of the metadata — do not paste them.
3. Confirm the file list is complete before publishing:
   ```bash
   export ZENODO_TOKEN=...
   python tools/zenodo_upload.py --deposition 21843280 --dry-run
   ```
   It prints what is stored versus what is local. Anything listed as needing upload means
   the archive is incomplete; re-run without `--dry-run` to finish.
4. **Publish in the browser.** The uploader deliberately cannot publish — minting a DOI is
   irreversible, and a published Zenodo record's files can never be changed.

Once published, note the DOI. It should be `10.5281/zenodo.21843280`.

## 2. The software deposit

The GitHub–Zenodo integration has **not** been enabled for this repository yet, so no code
archive exists. Until it is switched on, `/.zenodo.json` has no effect.

1. Sign in at <https://zenodo.org/account/settings/github/> with the GitHub account that
   owns the repository.
2. Find `davidpincus/cellquant` in the list and flip the toggle **On**. If it is missing,
   use "Sync now" — Zenodo only lists repositories you have admin rights to.
3. Zenodo archives on *release publication*, and only for releases published after the
   toggle was switched on. The existing `v1.1.0` release predates it, so re-trigger it:
   either delete and re-create the release from the same tag, or publish a new release.
4. Check that the resulting record matches `/.zenodo.json` — title, description, the four
   creators, MIT license, and both related identifiers.

## 3. Close the loop

After both DOIs exist:

- [ ] Add the software DOI to the data deposit's related identifiers
      (`isSupplementedBy` → the software DOI), replacing the placeholder GitHub URL.
      Zenodo allows metadata edits on a published record; only files are frozen.
- [ ] Add the DOI badge to `README.md` and the `doi:` field to `CITATION.cff`.
- [ ] Update the citation line in `README.md`, which currently reads
      "Submission and DOI pending."
- [ ] Update the "Data availability" section of `README.md` with the real DOI link.
- [ ] Once the paper is published, add its DOI to both deposits as a related identifier
      with relation `isSupplementTo`.

## Notes

- The `ZENODO_TOKEN` needs the `deposit:write` and `deposit:actions` scopes. Keep it in the
  environment; never pass it on a command line or commit it.
- A published record's **files are immutable**. Verify the file list before publishing, not
  after.
- Zenodo issues both a version DOI and a concept DOI that always resolves to the newest
  version. Cite the concept DOI in the manuscript so the reference does not go stale.
