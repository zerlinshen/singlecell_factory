#!/usr/bin/env python3
"""Download missing GSE162170 GEO supplementary resources with resumable ranges.

Reads a Wave-5 raw/public resource manifest (created during the all-figures
reproduction lane), downloads entries marked ``missing_download`` or incomplete,
verifies gzip integrity where applicable, computes SHA256, and updates the
manifest in place.

This intentionally uses curl rather than aria2c because aria2c is not installed
in the current environment. For large files it performs per-file byte-range
parallel downloads and concatenates verified parts.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable


def _sha256(path: Path, chunk: int = 1024 * 1024 * 8) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _gzip_ok(path: Path) -> bool:
    if not str(path).endswith(".gz"):
        return True
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1024 * 1024):
                pass
        return True
    except Exception:
        return False


def _run(cmd: list[str], log_path: Path) -> None:
    with log_path.open("a") as log:
        log.write("+ " + " ".join(cmd) + "\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")


def _curl_base() -> list[str]:
    return [
        "curl", "-fL",
        "--retry", "10",
        "--retry-delay", "5",
        "--retry-all-errors",
        "--connect-timeout", "30",
        "--max-time", "0",
        # If a connection silently stalls, fail it so the next invocation can
        # resume from the already-written range-part bytes.
        "--speed-time", "120",
        "--speed-limit", "1024",
    ]


def _download_small(url: str, out: Path, log_path: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    _run([*_curl_base(), "-C", "-", "-o", str(out), url], log_path)


def _range_part(args: tuple[str, Path, int, int, int, Path, int]) -> Path:
    url, part_dir, i, start, end, log_path, max_subchunk = args
    part = part_dir / f"part_{i:03d}"
    expected = end - start + 1
    if part.exists() and part.stat().st_size == expected:
        return part
    if part.exists() and part.stat().st_size > expected:
        part.unlink()
    current = part.stat().st_size if part.exists() else 0
    # Download/append in bounded subranges.  Large single range requests were
    # observed to fail late on NCBI with curl code 18, losing the temporary
    # partial.  Small append-only subranges make retries cheap and preserve
    # progress in the canonical part file.
    subchunk = max(512 * 1024, min(8 * 1024 * 1024, max_subchunk))
    while current < expected:
        sub_start = start + current
        sub_end = min(end, sub_start + subchunk - 1)
        want = sub_end - sub_start + 1
        tmp = part.with_suffix(part.suffix + ".resume")
        if tmp.exists():
            tmp.unlink()
        _run([*_curl_base(), "--range", f"{sub_start}-{sub_end}", "-o", str(tmp), url], log_path)
        got_tmp = tmp.stat().st_size
        if got_tmp != want:
            raise RuntimeError(f"range subchunk {tmp} size mismatch: got {got_tmp}, expected {want}")
        with part.open("ab") as w, tmp.open("rb") as r:
            shutil.copyfileobj(r, w, length=1024 * 1024 * 8)
        tmp.unlink()
        current = part.stat().st_size
    got = part.stat().st_size
    if got != expected:
        raise RuntimeError(f"range part {part} size mismatch: got {got}, expected {expected}")
    return part


def _download_range(url: str, out: Path, expected_bytes: int, parts: int, log_path: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    part_dir = out.parent / (out.name + ".parts")
    part_dir.mkdir(parents=True, exist_ok=True)
    chunk = (expected_bytes + parts - 1) // parts
    tasks = []
    for i in range(parts):
        start = i * chunk
        if start >= expected_bytes:
            break
        end = min(expected_bytes - 1, (i + 1) * chunk - 1)
        tasks.append((url, part_dir, i, start, end, log_path, 512 * 1024))
    with cf.ThreadPoolExecutor(max_workers=min(parts, len(tasks))) as ex:
        futs = [ex.submit(_range_part, t) for t in tasks]
        for fut in cf.as_completed(futs):
            fut.result()
    tmp = out.with_suffix(out.suffix + ".tmp")
    with tmp.open("wb") as w:
        for i in range(len(tasks)):
            part = part_dir / f"part_{i:03d}"
            with part.open("rb") as r:
                shutil.copyfileobj(r, w, length=1024 * 1024 * 8)
    if tmp.stat().st_size != expected_bytes:
        raise RuntimeError(f"assembled size mismatch for {out}: got {tmp.stat().st_size}, expected {expected_bytes}")
    tmp.replace(out)
    shutil.rmtree(part_dir)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--log", type=Path)
    ap.add_argument("--parts", type=int, default=8)
    ap.add_argument("--range-threshold-mb", type=int, default=50)
    ap.add_argument("--limit", type=int, default=0, help="download at most N missing files (0=all)")
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    log_path = args.log or (Path(manifest["run_dir"]) / "download_gse162170_geo_parallel.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write(f"download_start_utc={dt.datetime.now(dt.timezone.utc).isoformat()}\n")

    threshold = args.range_threshold_mb * 1024 * 1024
    downloaded = 0
    for row in manifest["geo_supplementary_files"]:
        out = Path(row["project_input_path"])
        expected = int(row.get("remote_bytes") or 0)
        if out.exists() and expected and out.stat().st_size == expected and _gzip_ok(out):
            row["status"] = "present"
            row["actual_bytes"] = out.stat().st_size
            row["sha256"] = _sha256(out)
            continue
        if out.exists() and out.is_symlink() and _gzip_ok(out):
            row["status"] = "present_linked"
            row["actual_bytes"] = out.stat().st_size
            row["sha256"] = _sha256(out)
            continue
        if args.limit and downloaded >= args.limit:
            continue
        row["download_attempted_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
        try:
            if expected >= threshold:
                _download_range(row["url"], out, expected, args.parts, log_path)
            else:
                _download_small(row["url"], out, log_path)
            if expected and out.stat().st_size != expected:
                raise RuntimeError(f"final size mismatch: got {out.stat().st_size}, expected {expected}")
            if not _gzip_ok(out):
                raise RuntimeError("gzip integrity check failed")
            row["status"] = "present"
            row["actual_bytes"] = out.stat().st_size
            row["sha256"] = _sha256(out)
            row["download_error"] = None
            downloaded += 1
            print(f"OK {out.name} {row['actual_bytes']} {row['sha256']}", flush=True)
        except Exception as e:
            row["status"] = "download_failed"
            row["download_error"] = str(e)
            args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
            print(f"FAIL {out.name}: {e}", file=sys.stderr, flush=True)
            return 2
        args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")

    manifest["download_completed_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    manifest["disk_free_bytes_after_download"] = shutil.disk_usage(manifest["project_root"]).free
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    with log_path.open("a") as log:
        log.write(f"download_end_utc={dt.datetime.now(dt.timezone.utc).isoformat()}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
