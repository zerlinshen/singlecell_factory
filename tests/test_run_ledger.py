from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from workflow.modular._run_ledger import RunLedger, _filter_env


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ctx(tmp_path: Path) -> MagicMock:
    cfg = MagicMock()
    cfg.cellranger.sample_root = tmp_path / "data" / "raw" / "sample"
    cfg.optional_modules = ["qc", "clustering"]
    ctx = MagicMock()
    ctx.cfg = cfg
    return ctx


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_record_start_captures_required_fields(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "test_project", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()

    r = ledger._record
    assert r["schema_version"] == "run_ledger_v1"
    assert r["project"] == "test_project"
    assert isinstance(r["timestamp_utc"], str) and r["timestamp_utc"].endswith("Z")
    assert isinstance(r["cli_args"], list)
    assert isinstance(r["env"], dict)
    assert "git_sha" in r
    assert "git_branch" in r
    assert isinstance(r["git_dirty"], bool)
    assert "conda_env" in r
    assert "python_version" in r
    assert "hostname" in r
    assert "sample_root" in r
    assert "optional_modules" in r


def test_record_end_writes_json_with_sha256(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "test_project", tmp_path)
    ledger.record_start()

    # Write a small fake h5ad file
    fake_h5ad = tmp_path / "final_adata.h5ad"
    fake_h5ad.write_bytes(b"FAKE_ADATA_CONTENT")

    ledger.record_end(fake_h5ad)
    out_path = ledger.write()

    assert out_path is not None and out_path.exists()
    data = json.loads(out_path.read_text())
    assert data["schema_version"] == "run_ledger_v1"
    assert len(data["final_adata_sha256"]) == 64  # SHA-256 hex
    assert "end_timestamp_utc" in data
    assert "peak_rss_bytes" in data


def test_record_module_appends(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "test_project", tmp_path)
    ledger.record_start()

    ledger.record_module("qc", "ok", "completed", 10.5, 1024 * 1024)
    ledger.record_module("clustering", "ok", "completed", 45.2, 2 * 1024 * 1024)
    ledger.record_module("differential_expression", "skipped", "missing dep", 0.0, 0)

    assert len(ledger._module_results) == 3
    names = [m["name"] for m in ledger._module_results]
    assert names == ["qc", "clustering", "differential_expression"]
    assert ledger._module_results[1]["wall_seconds"] == 45.2


def test_ledger_failure_does_not_raise(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "test_project", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()

    # Patch open to raise on write
    with patch("builtins.open", side_effect=OSError("disk full")):
        # record_end touches the rss poller but not open, so just call write
        try:
            ledger.record_end(None)
            ledger.write()
        except Exception as exc:
            pytest.fail(f"Ledger failure raised unexpectedly: {exc}")


def test_env_filter_only_sc_and_safe_keys(tmp_path):
    injection = {
        "SC_DE_ENGINE": "sparse",
        "SCF_MASSIVE_CHECKPOINT_POLICY": "metadata_only",
        "CUDA_VISIBLE_DEVICES": "0",
        "CONDA_DEFAULT_ENV": "sc_gpu",
        "PATH": "/usr/bin:/bin",        # must NOT appear
        "PYTHONSTARTUP": "/evil.py",    # must NOT appear
        "HOME": "/home/testuser",       # in allowlist, should appear
    }
    with patch.dict(os.environ, injection, clear=True):
        result = _filter_env()

    assert "SC_DE_ENGINE" in result
    assert "SCF_MASSIVE_CHECKPOINT_POLICY" in result
    assert "CUDA_VISIBLE_DEVICES" in result
    assert "CONDA_DEFAULT_ENV" in result
    assert "HOME" in result
    assert "PATH" not in result
    assert "PYTHONSTARTUP" not in result
