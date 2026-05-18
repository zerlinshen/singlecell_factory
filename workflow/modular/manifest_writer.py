from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def factory_git_state(repo_path: Path) -> dict:
    """Return {sha, dirty, diff_sha256} for the git repo at repo_path.

    sha is the 7-char short commit hash.
    diff_sha256 is sha256 of `git diff` output when dirty, else "".
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

    diff_sha256 = ""
    if dirty:
        try:
            diff_bytes = subprocess.check_output(
                ["git", "-C", str(repo), "diff"],
                stderr=subprocess.DEVNULL,
            )
            diff_sha256 = hashlib.sha256(diff_bytes).hexdigest()
        except subprocess.CalledProcessError:
            diff_sha256 = ""

    return {"sha": sha, "dirty": dirty, "diff_sha256": diff_sha256}


def write_manifest(
    run_dir: Path,
    *,
    project_id: str,
    run_id: str,
    modules_run: list[str],
    bundle_sha256: Optional[str] = None,
    produced_on: str = "",
    factory_python_path: Path = Path("/home/zerlinshen/singlecell_factory"),
    factory_r_path: Path = Path("/home/zerlinshen/multiomics_r_factory"),
    extra: Optional[dict] = None,
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

    r_state = factory_git_state(Path(factory_r_path))

    manifest = {
        "project_id": project_id,
        "run_id": run_id,
        "created_at": created_at,
        "produced_on": produced_on,
        "factory_python": factory_git_state(Path(factory_python_path)),
        "factory_r": r_state,
        "r_factory_sha_at_manifest_write": r_state["sha"],
        "modules_run": list(modules_run),
        "bundle_sha256": bundle_sha256 or "",
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

    state = factory_git_state(Path("/home/zerlinshen/singlecell_factory"))
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
