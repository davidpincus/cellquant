# Maintainer tools

Scripts used to prepare a release. Not part of the analysis pipeline — nothing in `docs/`
or the tutorials depends on them, and users never need to run them.

## `zenodo_upload.py`

Resumable, checksum-verified uploader for a Zenodo draft deposit, used to archive the
example and validation image data that is too large to live in git. It lists what the draft
already holds and skips those files, so an interrupted run can simply be restarted.

Requires a Zenodo token with the `deposit:write` and `deposit:actions` scopes, supplied
through the `ZENODO_TOKEN` environment variable — never as a command-line argument, and
never committed:

```bash
export ZENODO_TOKEN=...
python tools/zenodo_upload.py --deposition <id> --probe   # cheap connectivity check
python tools/zenodo_upload.py --deposition <id>           # the real upload
```

The script does not publish the deposit. Review the draft in the browser and publish it
there, so minting a DOI is always a deliberate act.
