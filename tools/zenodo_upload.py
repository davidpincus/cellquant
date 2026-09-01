#!/usr/bin/env python3
"""Resumable, checksum-verified uploader for a Zenodo draft.

Re-runnable: it lists what the draft already holds and skips those files, so
interrupt it freely and start it again.

SETUP
  1. Token: https://zenodo.org/account/settings/applications/tokens/new/
     scopes  deposit:write  and  deposit:actions
  2. export ZENODO_TOKEN=...
  3. Draft id is the number in the URL: https://zenodo.org/uploads/21843280

FIRST, DIAGNOSE (cheap, uploads 8 bytes and deletes it):
  python tools/zenodo_upload.py --deposition 21843280 --probe

THEN RUN IT IN THE BACKGROUND (survives closing the terminal):
  nohup python tools/zenodo_upload.py --deposition 21843280 \
        --log zenodo_upload.log >/dev/null 2>&1 &
  tail -f zenodo_upload.log          # watch it
  pkill -f zenodo_upload.py          # stop it

Changes from the first version, made after 502s on every retry of a 609 MB PUT:
  - token goes in an Authorization header, not the query string
  - the response body is printed on failure; a bare status code says nothing
  - explicit Content-Length, so requests cannot fall back to chunked transfer
    encoding, which some gateways reject with exactly this 502
  - 5xx backoff is long and capped high (server-side problems need minutes,
    not seconds), and a file that keeps failing is skipped so one bad file
    cannot block the other 89
  - line-buffered logging suitable for nohup
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time
from pathlib import Path

import requests

API = "https://zenodo.org/api"
DATA_DIRS = ["Tif6_Nsr1_Sis1_6hr/z-stacks", "SG_zstacks"]
CHUNK = 8 * 1024 * 1024
TIMEOUT = (30, 3600)          # (connect, read); a 609 MB PUT on a slow uplink is long
ATTEMPTS = 8


def log(msg: str = "") -> None:
    print(msg, flush=True)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(CHUNK), b""):
            h.update(b)
    return h.hexdigest()


def human(n: float) -> str:
    return f"{n/1e9:.2f} GB" if n >= 1e9 else f"{n/1e6:.0f} MB"


def collect(root: Path) -> list[Path]:
    out: list[Path] = []
    for d in DATA_DIRS:
        base = root / d
        if not base.is_dir():
            sys.exit(f"missing directory: {base}")
        out += [p for p in sorted(base.rglob("*"))
                if p.is_file() and not p.is_symlink() and not p.name.startswith(".")]
    names = [p.name for p in out]
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        sys.exit(f"duplicate file names across directories: {sorted(dupes)}")
    return out


def draft_info(dep: str, sess: requests.Session):
    """Return (bucket_url, {name: md5}, raw_json)."""
    last = None
    for path in (f"{API}/deposit/depositions/{dep}", f"{API}/records/{dep}/draft"):
        r = sess.get(path, timeout=(30, 120))
        last = r
        if r.status_code == 200:
            j = r.json()
            bucket = (j.get("links") or {}).get("bucket")
            if not bucket:
                continue
            have = {}
            for f in j.get("files", []) or []:
                k = f.get("filename") or f.get("key")
                if k:
                    have[k] = (f.get("checksum") or "").replace("md5:", "")
            return bucket, have, j
    sys.exit(f"could not open draft {dep} (last HTTP {last.status_code if last else '?'}): "
             f"{(last.text[:400] if last else '')}")


def show_quota(j: dict) -> None:
    """Print whatever the draft says about size limits. Field names differ
    between Zenodo's legacy and InvenioRDM payloads, so probe several."""
    for k in ("files", "links", "metadata"):
        j.get(k)
    q = j.get("quota") or (j.get("files") or {}).get("quota") if isinstance(j.get("files"), dict) else None
    size = None
    if isinstance(j.get("files"), dict):
        size = j["files"].get("total_bytes") or j["files"].get("size")
    for key in ("size", "total_bytes"):
        if size is None and key in j:
            size = j[key]
    if q or size:
        log(f"  draft reports: used={human(size) if size else '?'}  quota={q or '?'}")


def upload(bucket: str, path: Path, sess: requests.Session) -> tuple[str | None, str]:
    """PUT one file. Returns (remote_md5, ''), or (None, reason)."""
    url = f"{bucket}/{path.name}"
    size = path.stat().st_size
    for attempt in range(1, ATTEMPTS + 1):
        try:
            with path.open("rb") as fh:
                r = sess.put(url, data=fh, timeout=TIMEOUT,
                             headers={"Content-Type": "application/octet-stream",
                                      "Content-Length": str(size)})
            if r.status_code in (200, 201):
                return (r.json().get("checksum") or "").replace("md5:", ""), ""
            body = (r.text or "")[:500].replace("\n", " ")
            if r.status_code in (400, 401, 403, 413):
                return None, f"HTTP {r.status_code} (not retryable): {body}"
            log(f"    HTTP {r.status_code}: {body}")
            raise RuntimeError(f"HTTP {r.status_code}")
        except Exception as exc:
            if attempt == ATTEMPTS:
                return None, f"gave up after {ATTEMPTS} attempts: {exc}"
            wait = min(300, 5 * 2 ** (attempt - 1))     # 5,10,20,40,80,160,300
            log(f"    attempt {attempt}/{ATTEMPTS} failed ({exc}); sleeping {wait}s")
            time.sleep(wait)
    return None, "unreachable"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deposition", required=True)
    # This script lives in tools/, so the repo root -- where DATA_DIRS are
    # resolved from -- is one level up.
    ap.add_argument("--root", type=Path,
                    default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--probe", action="store_true",
                    help="upload a tiny test file and print the full response")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--log", type=Path, help="append output here (for nohup)")
    ap.add_argument("--limit", type=int, default=0,
                    help="upload at most N files (use --limit 1 to test one real file)")
    ap.add_argument("--delete", metavar="NAME",
                    help="delete one file from the draft, e.g. _zenodo_probe.txt")
    args = ap.parse_args()

    if args.log:
        sys.stdout = open(args.log, "a", buffering=1)
        sys.stderr = sys.stdout
    log(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')} start =====")

    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token:
        sys.exit("set ZENODO_TOKEN first")
    sess = requests.Session()
    sess.headers.update({"Authorization": f"Bearer {token}"})

    files = collect(args.root)
    log(f"local: {len(files)} files, {human(sum(f.stat().st_size for f in files))}")

    bucket, have, raw = draft_info(args.deposition, sess)
    log(f"draft {args.deposition}: {len(have)} file(s) already stored")
    show_quota(raw)

    if args.delete:
        r = sess.delete(f"{bucket}/{args.delete}", timeout=(30, 300))
        log(f"delete {args.delete}: HTTP {r.status_code} {(r.text or '')[:200]}")
        if r.status_code not in (200, 204):
            log("  (if this keeps failing, remove it in the browser: draft -> Files)")
        return

    if args.probe:
        log("\n--- probe: uploading 8 bytes to test the endpoint ---")
        tmp = args.root / "_zenodo_probe.txt"
        tmp.write_bytes(b"probe123")
        try:
            r = sess.put(f"{bucket}/{tmp.name}", data=tmp.read_bytes(), timeout=(30, 120))
            log(f"  status : {r.status_code}")
            log(f"  headers: {dict(list(r.headers.items())[:6])}")
            log(f"  body   : {(r.text or '')[:600]}")
            if r.status_code in (200, 201):
                d = sess.delete(f"{bucket}/{tmp.name}", timeout=(30, 120))
                log(f"  cleanup: HTTP {d.status_code}")
                log("\n  endpoint works for small files -> the 502 is size/duration related.")
            else:
                log("\n  endpoint rejects even 8 bytes -> read the body above.")
        finally:
            tmp.unlink(missing_ok=True)
        return

    todo, skipped = [], 0
    for f in files:
        if f.name in have and not (args.verify and md5(f) != have[f.name]):
            skipped += 1
        else:
            todo.append(f)
    if args.limit:
        todo = todo[:args.limit]
        log(f"--limit {args.limit}: restricting this run to {len(todo)} file(s)")
    log(f"skip {skipped}, upload {len(todo)} ({human(sum(f.stat().st_size for f in todo))})\n")
    if args.dry_run:
        for f in todo[:10]:
            log(f"  would upload {f.name:34s} {human(f.stat().st_size)}")
        return
    if not todo:
        log("nothing to do. Publish in the browser.")
        return

    done, failed, t0 = 0, [], time.time()
    for i, f in enumerate(todo, 1):
        size = f.stat().st_size
        log(f"[{i}/{len(todo)}] {f.name}  {human(size)}")
        remote, err = upload(bucket, f, sess)
        if remote is None:
            failed.append((f.name, err))
            log(f"    SKIPPED: {err}")
            continue
        local = md5(f)
        if remote and remote != local:
            failed.append((f.name, f"checksum mismatch {local} != {remote}"))
            log("    CHECKSUM MISMATCH - will retry on the next run")
            continue
        done += size
        rate = done / max(1e-9, time.time() - t0)
        left = sum(x.stat().st_size for x in todo[i:])
        log(f"    ok  {rate/1e6:.1f} MB/s  eta {left/max(rate,1)/60:.0f} min")

    log(f"\nuploaded {human(done)} in {(time.time()-t0)/60:.0f} min; "
        f"{len(failed)} file(s) failed")
    for n, e in failed:
        log(f"  ! {n}: {e}")
    log("Draft NOT published. Re-run to retry failures, then Publish in the browser.")


if __name__ == "__main__":
    main()
