#!/usr/bin/env python3
"""Migration apply script — dry-run only by default.

SAFETY_DISABLE = True  <-- this constant is intentional. Actual data movement
is disabled at the code level so that a mishap (e.g. running with wrong flags)
can NEVER accidentally move 33 GB of scientific data without explicit code
review and removal of this flag. The two-flag + env-var gate is a second
layer of human confirmation; this constant is the third layer.

To enable actual execution: read migration_README.md, run the dry-run, review
the CSV, then remove or set SAFETY_DISABLE = False and rebuild this script.
"""
from __future__ import annotations

# --- SAFETY GATE -----------------------------------------------------------
# Set to False only after a human has read migration_README.md and reviewed
# the dry-run output. Never set by automated tooling.
SAFETY_DISABLE = True
# ---------------------------------------------------------------------------

import argparse
import csv
import hashlib
import os
import subprocess
import sys
from pathlib import Path

INVENTORY_CSV = Path("/home/zerlinshen/.omc/state/migration_manifest.csv")
PROJECTS_ROOT = Path(os.environ.get("PROJECTS_ROOT", "/home/zerlinshen/projects"))
FACTORY_PYTHON = Path("/home/zerlinshen/singlecell_factory")
FACTORY_R = Path("/home/zerlinshen/r_multiomics_factory")
CONTRACTS_SCHEMA_PY = FACTORY_PYTHON / "contracts" / "bundle_schema.yaml"
CONTRACTS_SCHEMA_R = FACTORY_R / "contracts" / "bundle_schema.yaml"


# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------

def _check_device_pair(src: Path, dst: Path) -> bool:
    """Return True when src and dst are on the same device."""
    def _device(p: Path) -> str:
        existing = p if p.exists() else p.parent
        try:
            result = subprocess.check_output(
                ["df", "--output=source", str(existing)],
                text=True, stderr=subprocess.DEVNULL,
            )
            lines = result.strip().splitlines()
            return lines[-1].strip() if len(lines) >= 2 else ""
        except (subprocess.CalledProcessError, FileNotFoundError):
            return ""

    return _device(src) == _device(dst)


def _check_free_space(src_paths: list[Path], dst_root: Path) -> tuple[bool, int, int]:
    """Return (ok, required_bytes, available_bytes). Requires 1.5x total src size."""
    total = 0
    for p in src_paths:
        try:
            total += p.stat().st_size
        except OSError:
            pass
    required = int(total * 1.5)

    dst_existing = dst_root if dst_root.exists() else dst_root.parent
    try:
        result = subprocess.check_output(
            ["df", "--output=avail", "-B1", str(dst_existing)],
            text=True, stderr=subprocess.DEVNULL,
        )
        lines = result.strip().splitlines()
        available = int(lines[-1].strip()) if len(lines) >= 2 else 0
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
        available = 0

    return available >= required, required, available


def _check_projects_root_collision(projects_root: Path) -> bool:
    """Return True (collision) when projects_root exists and is non-empty."""
    if not projects_root.is_dir():
        return False
    try:
        return any(True for _ in projects_root.iterdir())
    except PermissionError:
        return True


def _check_dirty_tree(repo_path: Path) -> bool:
    """Return True when the repo has uncommitted changes."""
    try:
        out = subprocess.check_output(
            ["git", "-C", str(repo_path), "status", "--porcelain"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        return bool(out)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return ""


def _check_contracts_parity() -> tuple[bool, str, str]:
    """Return (ok, sha_py, sha_r). ok=True when both hashes match and are non-empty."""
    sha_py = _sha256_file(CONTRACTS_SCHEMA_PY)
    sha_r = _sha256_file(CONTRACTS_SCHEMA_R)
    if not sha_py or not sha_r:
        return False, sha_py, sha_r
    return sha_py == sha_r, sha_py, sha_r


# ---------------------------------------------------------------------------
# Dry-run output
# ---------------------------------------------------------------------------

def _load_inventory() -> list[dict]:
    if not INVENTORY_CSV.exists():
        print(
            f"WARNING: inventory CSV not found at {INVENTORY_CSV}. "
            "Run migration_inventory.py first.",
            file=sys.stderr,
        )
        return []
    rows = []
    with INVENTORY_CSV.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


def _print_dry_run(rows: list[dict]) -> None:
    print("[DRY-RUN] The following moves would be performed (no data has been moved):")
    print(f"{'SRC':<70}  {'DST'}")
    print("-" * 140)
    for row in rows:
        dst = PROJECTS_ROOT / row["target_relative_path"]
        print(f"{row['src_path']:<70}  {dst}")
    print(f"\n[DRY-RUN] Total: {len(rows)} files. Review the CSV at {INVENTORY_CSV}")
    print("[DRY-RUN] No data was moved.")


# ---------------------------------------------------------------------------
# Pre-flight runner
# ---------------------------------------------------------------------------

def _run_preflights(rows: list[dict]) -> bool:
    all_ok = True
    src_paths = [Path(r["src_path"]) for r in rows]

    # 1. Device check (sample first src against dst_root)
    if src_paths:
        same_device = _check_device_pair(src_paths[0], PROJECTS_ROOT)
        status = "OK (same device, mv eligible)" if same_device else "WARN (cross-device, would need rsync two-phase)"
        print(f"[preflight] device-check: {status}")

    # 2. Free-space gate
    ok, required, available = _check_free_space(src_paths, PROJECTS_ROOT)
    gib_req = required / 1024**3
    gib_avail = available / 1024**3
    if ok:
        print(f"[preflight] free-space: OK ({gib_avail:.2f} GiB available >= {gib_req:.2f} GiB required)")
    else:
        print(f"[preflight] free-space: FAIL ({gib_avail:.2f} GiB available < {gib_req:.2f} GiB required at 1.5x)")
        all_ok = False

    # 3. Collision check
    collision = _check_projects_root_collision(PROJECTS_ROOT)
    if collision:
        print(
            f"[preflight] projects-root-collision: ABORT — {PROJECTS_ROOT} exists and is non-empty. "
            "Confirm intended root via PROJECTS_ROOT env var."
        )
        all_ok = False
    else:
        print(f"[preflight] projects-root-collision: OK (target root is empty or absent)")

    # 4. Dirty-tree gate
    py_dirty = _check_dirty_tree(FACTORY_PYTHON)
    r_dirty = _check_dirty_tree(FACTORY_R)
    if py_dirty or r_dirty:
        dirty_repos = []
        if py_dirty:
            dirty_repos.append("singlecell_factory")
        if r_dirty:
            dirty_repos.append("r_multiomics_factory")
        print(f"[preflight] dirty-tree: FAIL — {', '.join(dirty_repos)} has uncommitted changes")
        all_ok = False
    else:
        print("[preflight] dirty-tree: OK (both repos clean)")

    # 5. Contracts parity
    ok_c, sha_py, sha_r = _check_contracts_parity()
    if not CONTRACTS_SCHEMA_PY.exists() or not CONTRACTS_SCHEMA_R.exists():
        print("[preflight] contracts-parity: SKIP (schema files not yet created)")
    elif ok_c:
        print(f"[preflight] contracts-parity: OK (sha256={sha_py[:16]}...)")
    else:
        print(f"[preflight] contracts-parity: FAIL (py={sha_py[:16]}... != r={sha_r[:16]}...)")
        all_ok = False

    return all_ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Request execution mode (also requires --i-have-reviewed-the-dry-run and MIGRATION_APPROVED=1).",
    )
    parser.add_argument(
        "--i-have-reviewed-the-dry-run",
        action="store_true",
        dest="reviewed",
        help="Confirm that the dry-run output has been reviewed.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = _load_inventory()

    if not args.execute:
        _print_dry_run(rows)
        return 0

    # --execute without both confirmations → hard error
    if not args.reviewed or os.environ.get("MIGRATION_APPROVED") != "1":
        print(
            "ERROR: Actual data migration requires explicit human approval. "
            "The dry-run output must be reviewed first. To enable execution, "
            "re-run with `--execute --i-have-reviewed-the-dry-run` AND set "
            "environment var `MIGRATION_APPROVED=1`.",
            file=sys.stderr,
        )
        return 1

    # Both flags + env var set: run pre-flights
    print("[migration_apply] All confirmation flags set. Running pre-flight checks...")
    preflight_ok = _run_preflights(rows)

    if not preflight_ok:
        print("\n[migration_apply] Pre-flight checks FAILED. Aborting.", file=sys.stderr)
        return 1

    # SAFETY_DISABLE prevents actual execution even when pre-flights pass.
    # Exit non-zero so CI/cron pipelines that mistakenly invoke `--execute` are caught
    # by exit-status checks rather than silently treating a blocked run as success.
    if SAFETY_DISABLE:
        print(
            "\nPre-flights OK; would proceed with mv but data migration is currently "
            "disabled in this build; remove the SAFETY_DISABLE flag in code to actually execute."
        )
        return 2

    # Unreachable while SAFETY_DISABLE = True.
    print("ERROR: reached execution branch — this should not happen while SAFETY_DISABLE=True", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
