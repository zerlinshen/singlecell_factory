from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

_RUN_ID_RE = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{4}Z-[0-9a-f]{7}$")


def generate_run_id(factory_sha: Optional[str] = None) -> str:
    """Format: 2026-05-14T1530Z-a3f4c2b.

    Uses a valid all-zero SHA placeholder when the factory SHA is unavailable so
    auto-generated IDs still satisfy the public run-id contract.
    """
    now = datetime.now(timezone.utc)
    ts = now.strftime("%Y-%m-%dT%H%MZ")
    sha = (factory_sha or "0000000")[:7].lower()
    return f"{ts}-{sha}"


def validate_run_id(run_id: str) -> bool:
    """Return True iff run_id matches ^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{4}Z-[0-9a-f]{7}$."""
    return bool(_RUN_ID_RE.match(run_id))


def resolve_run_dir(
    project_root: Path,
    run_id: Optional[str] = None,
    factory_sha: Optional[str] = None,
) -> Path:
    """Resolve <project-root>/runs/<run-id>/ and create subdirectories.

    Auto-generates run_id as <UTC-ts>-<short-sha> when absent.
    Creates: python/, python/bundle/, r/, logs/.
    Returns the run_dir Path.
    """
    if run_id is None:
        run_id = generate_run_id(factory_sha)
    if not validate_run_id(run_id):
        raise ValueError(
            "invalid run_id. Expected format YYYY-MM-DDTHHMMZ-<7hex>; "
            f"got {run_id!r}"
        )
    run_dir = Path(project_root) / "runs" / run_id
    for subdir in [
        run_dir / "python" / "bundle",
        run_dir / "r",
        run_dir / "logs",
    ]:
        subdir.mkdir(parents=True, exist_ok=True)
    return run_dir


def python_dir(run_dir: Path) -> Path:
    return Path(run_dir) / "python"


def bundle_dir(run_dir: Path) -> Path:
    return Path(run_dir) / "python" / "bundle"


def r_dir(run_dir: Path) -> Path:
    return Path(run_dir) / "r"


def logs_dir(run_dir: Path) -> Path:
    return Path(run_dir) / "logs"


def manifest_path(run_dir: Path) -> Path:
    return Path(run_dir) / "manifest.json"


if __name__ == "__main__":
    import tempfile

    sha = "a3f4c2b"
    rid = generate_run_id(sha)
    print(f"run_id: {rid}")
    print(f"valid:  {validate_run_id(rid)}")

    with tempfile.TemporaryDirectory() as tmp:
        run_dir = resolve_run_dir(Path(tmp) / "my-project", run_id=rid)
        print(f"run_dir:    {run_dir}")
        print(f"python_dir: {python_dir(run_dir)}")
        print(f"bundle_dir: {bundle_dir(run_dir)}")
        print(f"r_dir:      {r_dir(run_dir)}")
        print(f"logs_dir:   {logs_dir(run_dir)}")
        print(f"manifest:   {manifest_path(run_dir)}")
        assert python_dir(run_dir).exists()
        assert bundle_dir(run_dir).exists()
        assert r_dir(run_dir).exists()
        assert logs_dir(run_dir).exists()
        print("all subdirs created OK")
