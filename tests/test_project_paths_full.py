"""Tests for workflow.modular.project_paths — coverage for uncovered paths."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from workflow.modular.project_paths import (
    bundle_dir,
    generate_run_id,
    logs_dir,
    manifest_path,
    python_dir,
    r_dir,
    resolve_run_dir,
    validate_run_id,
)

_VALID_RUN_IDS = [
    "2026-05-14T1530Z-a3f4c2b",
    "2024-01-01T0000Z-0000000",
    "2099-12-31T2359Z-fffffff",
]

_INVALID_RUN_IDS = [
    "",
    "2026-05-14T1530Z-a3f4c2",       # sha too short (6 chars)
    "2026-05-14T1530Z-a3f4c2bb",      # sha too long (8 chars)
    "2026-05-14 1530Z-a3f4c2b",       # space instead of T
    "2026-05-14T1530Z-A3F4C2B",       # uppercase sha
    "26-05-14T1530Z-a3f4c2b",         # 2-digit year
    "2026-05-14T1530-a3f4c2b",        # missing Z
    "not-a-run-id",
    "2026-05-14T1530Za3f4c2b",        # missing hyphen before sha
]


def test_validate_run_id_valid_samples():
    for rid in _VALID_RUN_IDS:
        assert validate_run_id(rid), f"Expected valid: {rid}"


def test_validate_run_id_invalid_samples():
    for rid in _INVALID_RUN_IDS:
        assert not validate_run_id(rid), f"Expected invalid: {rid}"


def test_generate_run_id_with_sha_matches_pattern():
    rid = generate_run_id(factory_sha="a3f4c2b")
    assert validate_run_id(rid), f"generate_run_id(sha) produced invalid id: {rid}"


def test_generate_run_id_uses_provided_sha():
    rid = generate_run_id(factory_sha="a3f4c2b")
    assert rid.endswith("-a3f4c2b"), f"Expected sha suffix: {rid}"


def test_generate_run_id_truncates_long_sha():
    rid = generate_run_id(factory_sha="a3f4c2bdeadbeef")
    # Only first 7 chars of sha should appear
    assert re.search(r"-[0-9a-f]{7}$", rid), f"sha not truncated: {rid}"
    assert rid.endswith("-a3f4c2b"), f"wrong sha in id: {rid}"


def test_generate_run_id_unknown_placeholder_when_no_sha():
    rid = generate_run_id(factory_sha=None)
    assert rid.endswith("-0000000"), f"Expected valid zero placeholder: {rid}"
    assert validate_run_id(rid), f"Expected valid generated run_id without sha: {rid}"
    assert "T" in rid and "Z-" in rid, f"Timestamp format wrong: {rid}"


def test_resolve_run_dir_creates_required_subdirs(tmp_path):
    project_root = tmp_path / "my-project"
    run_dir = resolve_run_dir(project_root, run_id="2026-05-14T1530Z-a3f4c2b")
    assert run_dir == project_root / "runs" / "2026-05-14T1530Z-a3f4c2b"
    assert (run_dir / "python" / "bundle").is_dir()
    assert (run_dir / "r").is_dir()
    assert (run_dir / "logs").is_dir()


def test_resolve_run_dir_rejects_invalid_run_id_before_mkdir(tmp_path):
    project_root = tmp_path / "my-project"
    with pytest.raises(ValueError, match="invalid run_id"):
        resolve_run_dir(project_root, run_id="../../outside")
    assert not (tmp_path / "outside").exists()
    assert not (project_root / "runs").exists()


def test_resolve_run_dir_auto_generates_run_id(tmp_path):
    project_root = tmp_path / "proj"
    run_dir = resolve_run_dir(project_root, factory_sha="a3f4c2b")
    # With a real hex sha, auto-generated run_id should be valid
    run_id = run_dir.name
    assert validate_run_id(run_id), f"Auto-generated run_id invalid: {run_id}"
    assert (run_dir / "python" / "bundle").is_dir()


def test_resolve_run_dir_is_idempotent(tmp_path):
    project_root = tmp_path / "proj"
    rid = "2026-05-14T0001Z-aaaaaaa"
    run_dir1 = resolve_run_dir(project_root, run_id=rid)
    run_dir2 = resolve_run_dir(project_root, run_id=rid)
    assert run_dir1 == run_dir2
    assert run_dir2.is_dir()


def test_path_helpers_return_correct_subpaths(tmp_path):
    run_dir = tmp_path / "runs" / "2026-05-14T0001Z-aaaaaaa"
    assert python_dir(run_dir) == run_dir / "python"
    assert bundle_dir(run_dir) == run_dir / "python" / "bundle"
    assert r_dir(run_dir) == run_dir / "r"
    assert logs_dir(run_dir) == run_dir / "logs"
    assert manifest_path(run_dir) == run_dir / "manifest.json"
