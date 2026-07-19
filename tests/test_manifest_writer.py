"""Tests for workflow.modular.manifest_writer."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from unittest.mock import patch

from workflow.modular.manifest_writer import factory_git_state, write_manifest


# ---------------------------------------------------------------------------
# factory_git_state tests
# ---------------------------------------------------------------------------

def test_factory_git_state_real_repo():
    """Happy path: singlecell_factory is a git repo — sha must be 7 hex chars."""
    repo = Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory")
    state = factory_git_state(repo)
    assert isinstance(state["sha"], str)
    assert len(state["sha"]) == 7
    assert all(c in "0123456789abcdef" for c in state["sha"])
    assert isinstance(state["dirty"], bool)
    assert isinstance(state["diff_sha256"], str)


def test_factory_git_state_not_git_path(tmp_path):
    """Non-git directory returns sentinel empty values."""
    state = factory_git_state(tmp_path)
    assert state == {"sha": "", "dirty": False, "diff_sha256": ""}


def test_factory_git_state_nonexistent_path():
    """Nonexistent path returns sentinel empty values."""
    state = factory_git_state(Path("/nonexistent/path/xyz"))
    assert state == {"sha": "", "dirty": False, "diff_sha256": ""}


def test_factory_git_state_dirty_tree_returns_nonempty_diff_sha256(tmp_path):
    """When repo is dirty, diff_sha256 should be a 64-char hex string."""
    sha_mock = "abc1234\n"
    status_mock = "M  workflow/modular/cli.py\n"
    diff_mock = b"fake diff content for test"

    with patch("subprocess.check_output") as mock_cmd:
        def side_effect(args, **kwargs):
            if "--short=7" in args:
                return sha_mock
            if "status" in args:
                return status_mock
            if "diff" in args and "--cached" not in args:
                return diff_mock
            return b""
        mock_cmd.side_effect = side_effect

        state = factory_git_state(tmp_path)

    assert state["sha"] == "abc1234"
    assert state["dirty"] is True
    assert len(state["diff_sha256"]) == 64
    assert all(c in "0123456789abcdef" for c in state["diff_sha256"])


def test_factory_git_state_records_staged_and_untracked_sources(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(
        ["git", "-C", str(tmp_path), "config", "user.email", "test@example.com"],
        check=True,
    )
    subprocess.run(["git", "-C", str(tmp_path), "config", "user.name", "Test"], check=True)
    tracked = tmp_path / "tracked.txt"
    tracked.write_text("base\n", encoding="utf-8")
    unstaged = tmp_path / "unstaged.txt"
    unstaged.write_text("base\n", encoding="utf-8")
    subprocess.run(
        ["git", "-C", str(tmp_path), "add", "tracked.txt", "unstaged.txt"],
        check=True,
    )
    subprocess.run(["git", "-C", str(tmp_path), "commit", "-qm", "base"], check=True)

    tracked.write_text("staged\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(tmp_path), "add", "tracked.txt"], check=True)
    unstaged.write_text("unstaged\n", encoding="utf-8")
    untracked = tmp_path / "untracked.txt"
    untracked.write_text("untracked payload\n", encoding="utf-8")

    state = factory_git_state(tmp_path)
    assert state["dirty"] is True
    assert len(state["tracked_diff_sha256"]) == 64
    assert len(state["staged_diff_sha256"]) == 64
    assert state["untracked_files"] == ["untracked.txt"]
    assert len(state["untracked_inventory_sha256"]) == 64
    assert len(state["untracked_content_sha256"]) == 64
    assert state["untracked_content_hash_complete"] is True


def test_factory_git_state_clean_tree_has_empty_diff_sha256(tmp_path):
    """When repo is clean, diff_sha256 must be empty string."""
    sha_mock = "abc1234\n"
    status_mock = ""  # clean

    with patch("subprocess.check_output") as mock_cmd:
        def side_effect(args, **kwargs):
            if "--short=7" in args:
                return sha_mock
            if "status" in args:
                return status_mock
            return b""
        mock_cmd.side_effect = side_effect

        state = factory_git_state(tmp_path)

    assert state["dirty"] is False
    assert state["diff_sha256"] == ""


# ---------------------------------------------------------------------------
# write_manifest tests
# ---------------------------------------------------------------------------

def test_write_manifest_produces_valid_json(tmp_path):
    run_dir = tmp_path / "runs" / "2026-05-14T1530Z-a3f4c2b"
    run_dir.mkdir(parents=True)
    out = write_manifest(
        run_dir,
        project_id="test-project",
        run_id="2026-05-14T1530Z-a3f4c2b",
        modules_run=["cellranger", "qc", "doublet_detection"],
        produced_on="ubuntu-tail",
        factory_python_path=Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory"),
        factory_r_path=Path("/tmp/nonexistent_r_factory"),
    )
    assert out.exists()
    data = json.loads(out.read_text())
    # Required top-level fields
    assert data["project_id"] == "test-project"
    assert data["run_id"] == "2026-05-14T1530Z-a3f4c2b"
    assert data["modules_run"] == ["cellranger", "qc", "doublet_detection"]
    assert "created_at" in data
    assert "factory_python" in data
    assert "factory_r" in data
    assert "claim_guard" in data


def test_write_manifest_has_r_factory_sha_at_manifest_write(tmp_path):
    """PREC-1 contract: r_factory_sha_at_manifest_write must be a top-level field."""
    run_dir = tmp_path / "runs" / "2026-05-14T0001Z-aaaaaaa"
    run_dir.mkdir(parents=True)
    out = write_manifest(
        run_dir,
        project_id="proj",
        run_id="2026-05-14T0001Z-aaaaaaa",
        modules_run=["qc"],
        factory_r_path=Path("/tmp/nonexistent_r_factory"),
    )
    data = json.loads(out.read_text())
    assert "r_factory_sha_at_manifest_write" in data
    # Must match factory_r.sha (the r_state sha)
    assert data["r_factory_sha_at_manifest_write"] == data["factory_r"]["sha"]


def test_write_manifest_bundle_sha256_defaults_empty(tmp_path):
    run_dir = tmp_path / "runs" / "2026-05-14T0001Z-aaaaaaa"
    run_dir.mkdir(parents=True)
    out = write_manifest(
        run_dir,
        project_id="proj",
        run_id="2026-05-14T0001Z-aaaaaaa",
        modules_run=[],
    )
    data = json.loads(out.read_text())
    assert data["bundle_sha256"] == ""


def test_write_manifest_accepts_extra_fields(tmp_path):
    run_dir = tmp_path / "runs" / "2026-05-14T0001Z-bbbbbbb"
    run_dir.mkdir(parents=True)
    extra = {"sc_require_project_root": True, "custom_key": "hello"}
    out = write_manifest(
        run_dir,
        project_id="proj",
        run_id="2026-05-14T0001Z-bbbbbbb",
        modules_run=["qc"],
        extra=extra,
    )
    data = json.loads(out.read_text())
    assert data["extra"] == extra


def test_write_manifest_creates_parent_dirs(tmp_path):
    run_dir = tmp_path / "deep" / "nested" / "runs" / "2026-05-14T0001Z-ccccccc"
    # Do NOT create run_dir — write_manifest must handle it
    out = write_manifest(
        run_dir,
        project_id="proj",
        run_id="2026-05-14T0001Z-ccccccc",
        modules_run=["qc"],
    )
    assert out.exists()
