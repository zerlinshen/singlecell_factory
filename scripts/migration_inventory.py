#!/usr/bin/env python3
"""Walk factory output/ trees and emit a migration inventory CSV.

Read-only — no files are moved or modified.
Output: /home/zerlinshen/.omc/state/migration_manifest.csv
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

SCAN_ROOTS = [
    Path("/home/zerlinshen/singlecell_factory/output"),
    Path("/home/zerlinshen/r_multiomics_factory/output"),
]
OUTPUT_CSV = Path("/home/zerlinshen/.omc/state/migration_manifest.csv")
PROJECTS_ROOT = Path("/home/zerlinshen/projects")

_PROJECT_ID_PATTERNS = [
    (re.compile(r"nc2024", re.I), "nc2024-nsclc"),
    (re.compile(r"nsclc", re.I), "nc2024-nsclc"),
    (re.compile(r"staircase[_-]nano", re.I), "staircase-nano"),
    (re.compile(r"staircase[_-]small", re.I), "staircase-small-real"),
    (re.compile(r"staircase[_-]medium", re.I), "staircase-medium-real"),
    (re.compile(r"staircase", re.I), "staircase-nano"),
]

_RUN_TS_RE = re.compile(r"(\d{8}T\d{6}|\d{4}-\d{2}-\d{2})")


def _guess_project_id(path: Path) -> str:
    text = str(path)
    for pattern, project_id in _PROJECT_ID_PATTERNS:
        if pattern.search(text):
            return project_id
    return "_unsorted"


def _guess_run_timestamp(path: Path) -> str:
    for part in reversed(path.parts):
        m = _RUN_TS_RE.search(part)
        if m:
            return m.group(1)
    return ""


def _load_manifest_hint(path: Path) -> dict:
    for candidate in (
        path.parent / "bundle_manifest.json",
        path.parent / "run_manifest.json",
        path.parent.parent / "bundle_manifest.json",
        path.parent.parent / "run_manifest.json",
    ):
        if candidate.is_file():
            try:
                data = json.loads(candidate.read_text(encoding="utf-8", errors="replace"))
                return data
            except Exception:
                pass
    return {}


def _target_relative_path(src: Path, project_id: str, run_ts: str) -> str:
    run_segment = f"runs/{run_ts or 'unknown-run'}/python"
    rel = src.name
    return f"{project_id}/{run_segment}/{rel}"


def scan() -> list[dict]:
    rows: list[dict] = []
    for root in SCAN_ROOTS:
        if not root.exists():
            continue
        for dirpath, _, filenames in os.walk(root):
            for fname in filenames:
                fpath = Path(dirpath) / fname
                try:
                    stat = fpath.stat()
                except OSError:
                    continue
                hint = _load_manifest_hint(fpath)
                project_id = hint.get("project_id") or _guess_project_id(fpath)
                run_ts = _guess_run_timestamp(fpath)
                rows.append({
                    "src_path": str(fpath),
                    "size_bytes": stat.st_size,
                    "mtime": int(stat.st_mtime),
                    "project_id_guess": project_id,
                    "run_timestamp_guess": run_ts,
                    "target_relative_path": _target_relative_path(fpath, project_id, run_ts),
                })
    return rows


def write_csv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "src_path", "size_bytes", "mtime",
        "project_id_guess", "run_timestamp_guess", "target_relative_path",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict]) -> None:
    total_files = len(rows)
    total_bytes = sum(r["size_bytes"] for r in rows)
    by_project: dict[str, dict] = {}
    for r in rows:
        pid = r["project_id_guess"]
        if pid not in by_project:
            by_project[pid] = {"files": 0, "bytes": 0}
        by_project[pid]["files"] += 1
        by_project[pid]["bytes"] += r["size_bytes"]

    print(f"Migration inventory — {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}")
    print(f"  Total files : {total_files}")
    print(f"  Total bytes : {total_bytes:,}  ({total_bytes / (1024**3):.2f} GiB)")
    print(f"  Per project:")
    for pid, counts in sorted(by_project.items()):
        gib = counts['bytes'] / (1024**3)
        print(f"    {pid:30s}  {counts['files']:6d} files  {gib:.2f} GiB")
    print(f"  CSV written : {OUTPUT_CSV}")


def main() -> int:
    rows = scan()
    write_csv(rows, OUTPUT_CSV)
    print_summary(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
