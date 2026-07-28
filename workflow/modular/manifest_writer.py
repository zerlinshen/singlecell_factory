from __future__ import annotations

import hashlib
import json
import os
import subprocess
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from workflow.factory_paths import R_MULTIOMICS_FACTORY_ROOT, SINGLECELL_FACTORY_ROOT
from workflow.modular.software_provenance import collect_software_versions

_MAX_UNTRACKED_HASH_BYTES = 16 * 1024 * 1024
_MAX_UNTRACKED_FILE_HASH_BYTES = 1024 * 1024


def factory_git_state(repo_path: Path) -> dict:
    """Return a complete, machine-auditable source state for a git repo.

    sha is the 7-char short commit hash.
    ``diff_sha256`` is retained as the legacy unstaged tracked-diff hash.
    Staged changes and untracked paths/content are fingerprinted separately so
    a dirty tree cannot be represented by the empty-diff hash alone.
    Returns {"sha": "", "dirty": False, "diff_sha256": ""} for non-git paths.
    """
    repo = Path(repo_path)
    try:
        sha = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--short=7", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"sha": "", "dirty": False, "diff_sha256": ""}

    try:
        status = subprocess.check_output(
            ["git", "-C", str(repo), "status", "--porcelain"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        dirty = bool(status)
    except subprocess.CalledProcessError:
        dirty = False

    def _git_bytes(*args: str) -> bytes:
        return subprocess.check_output(
            ["git", "-C", str(repo), *args], stderr=subprocess.DEVNULL
        )

    try:
        tracked_diff = _git_bytes("diff", "--no-ext-diff", "--binary")
        staged_diff = _git_bytes("diff", "--cached", "--no-ext-diff", "--binary")
        untracked_raw = _git_bytes("ls-files", "--others", "--exclude-standard", "-z")
    except subprocess.CalledProcessError:
        tracked_diff = staged_diff = untracked_raw = b""

    untracked_files = sorted(
        path.decode("utf-8", errors="surrogateescape")
        for path in untracked_raw.split(b"\0")
        if path
    )
    content_hash = hashlib.sha256()
    content_complete = True
    content_hashed_bytes = 0
    for relative in untracked_files:
        path = repo / relative
        encoded_path = os.fsencode(relative)
        content_hash.update(len(encoded_path).to_bytes(8, "big"))
        content_hash.update(encoded_path)
        try:
            if path.is_symlink():
                payload = os.fsencode(os.readlink(path))
            else:
                remaining = max(0, _MAX_UNTRACKED_HASH_BYTES - content_hashed_bytes)
                limit = min(_MAX_UNTRACKED_FILE_HASH_BYTES, remaining)
                with path.open("rb") as handle:
                    payload = handle.read(limit + 1)
                if len(payload) > limit:
                    payload = payload[:limit]
                    content_complete = False
                if path.stat().st_size > len(payload):
                    content_complete = False
        except (OSError, UnicodeError):
            content_complete = False
            payload = b""
        content_hash.update(len(payload).to_bytes(8, "big"))
        content_hash.update(payload)
        content_hashed_bytes += len(payload)

    tracked_sha = hashlib.sha256(tracked_diff).hexdigest() if tracked_diff else ""
    staged_sha = hashlib.sha256(staged_diff).hexdigest() if staged_diff else ""
    inventory_bytes = "\0".join(untracked_files).encode("utf-8", errors="surrogateescape")
    return {
        "sha": sha,
        "dirty": dirty,
        "diff_sha256": tracked_sha,
        "tracked_diff_sha256": tracked_sha,
        "staged_diff_sha256": staged_sha,
        "untracked_files": untracked_files,
        "untracked_inventory_sha256": (
            hashlib.sha256(inventory_bytes).hexdigest() if untracked_files else ""
        ),
        "untracked_content_sha256": content_hash.hexdigest() if untracked_files else "",
        "untracked_content_hash_complete": content_complete,
        "untracked_content_hashed_bytes": content_hashed_bytes,
    }


def write_manifest(
    run_dir: Path,
    *,
    project_id: str,
    run_id: str,
    modules_run: list[str],
    bundle_sha256: Optional[str] = None,
    produced_on: str = "",
    factory_python_path: Path = SINGLECELL_FACTORY_ROOT,
    factory_r_path: Path = R_MULTIOMICS_FACTORY_ROOT,
    extra: Optional[dict] = None,
    requested_modules: Optional[list[str]] = None,
    planned_modules: Optional[list[str]] = None,
    executed_modules: Optional[list[str]] = None,
    completed_modules: Optional[list[str]] = None,
    skipped_modules: Optional[list[str]] = None,
    failed_modules: Optional[list[str]] = None,
    overall_status: str = "",
    producer_manifest: str = "",
    factory_python_state: Optional[dict] = None,
) -> Path:
    """Write <run-dir>/manifest.json with the canonical schema. Returns the path.

    Dual R-factory SHA contract (PREC-1):
      factory_r.sha  — the R factory short SHA at manifest-write time (legacy field, kept).
      r_factory_sha_at_manifest_write — explicit top-level copy for unambiguous consumption.
      The bundle exporter writes r_factory_sha_at_export to bundle/provenance.json at export
      time (may differ from r_factory_sha_at_manifest_write if a pull happened between the
      two steps). The R loader warns on mismatch between the two fields.
    """
    created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    if "multiomics_r_factory" in str(factory_r_path):  # legacy name accepted during transition
        warnings.warn(
            "factory_r_path uses legacy name 'multiomics_r_factory';"  # legacy name accepted during transition
            " use 'r_multiomics_factory' instead.",
            DeprecationWarning,
            stacklevel=2,
        )

    r_state = factory_git_state(Path(factory_r_path))

    manifest = {
        "project_id": project_id,
        "run_id": run_id,
        "created_at": created_at,
        "produced_on": produced_on,
        "factory_python": (
            dict(factory_python_state)
            if factory_python_state is not None
            else factory_git_state(Path(factory_python_path))
        ),
        "factory_r": r_state,
        "r_factory_sha_at_manifest_write": r_state["sha"],
        # The git SHAs above pin the *code*; software_versions pins the stack it
        # ran against. Both are needed: the environment files intentionally leave
        # most of the scientific stack unpinned, so the same SHA can produce
        # different numbers on two hosts. See software_provenance.py.
        "software_versions": collect_software_versions(),
        "modules_run": list(modules_run),
        "requested_modules": list(requested_modules or modules_run),
        "planned_modules": list(planned_modules or requested_modules or modules_run),
        "executed_modules": list(executed_modules or modules_run),
        "completed_modules": list(completed_modules or modules_run),
        "skipped_modules": list(skipped_modules or []),
        "failed_modules": list(failed_modules or []),
        "overall_status": overall_status or ("failed" if failed_modules else "complete"),
        "bundle_sha256": bundle_sha256 or "",
        "producer_manifest": producer_manifest,
        "claim_guard": "not_for_de_or_new_quantitative_claims_without_full_object_validation",
        "extra": extra or {},
    }

    out_path = Path(run_dir) / "manifest.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return out_path


if __name__ == "__main__":
    import tempfile
    from pathlib import Path

    state = factory_git_state(SINGLECELL_FACTORY_ROOT)
    print(f"factory_git_state: {state}")

    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp) / "runs" / "2026-05-14T1530Z-a3f4c2b"
        run_dir.mkdir(parents=True)
        out = write_manifest(
            run_dir,
            project_id="test-project",
            run_id="2026-05-14T1530Z-a3f4c2b",
            modules_run=["qc", "clustering"],
            produced_on="ubuntu-tail",
        )
        print(f"manifest written to: {out}")
        print(out.read_text())
