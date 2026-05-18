"""Extended tests for RunLedger covering paths not in test_run_ledger.py."""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from workflow.modular._run_ledger import RunLedger, _filter_env


def _make_ctx(tmp_path: Path) -> MagicMock:
    cfg = MagicMock()
    cfg.cellranger.sample_root = tmp_path / "data" / "raw" / "sample"
    cfg.optional_modules = ["qc"]
    ctx = MagicMock()
    ctx.cfg = cfg
    return ctx


def test_write_returns_path_with_correct_filename(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "myproject", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()
    ledger.record_end(None)
    out = ledger.write()
    assert out is not None
    assert out.exists()
    assert out.name.startswith("myproject_")
    assert out.suffix == ".json"


def test_write_produces_parseable_json(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "myproject", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()
    ledger.record_end(None)
    out = ledger.write()
    data = json.loads(out.read_text())
    assert data["schema_version"] == "run_ledger_v1"


def test_record_end_with_missing_h5ad_has_empty_sha(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "proj", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()
    missing = tmp_path / "nonexistent_final_adata.h5ad"
    ledger.record_end(missing)
    assert ledger._record.get("final_adata_sha256") == ""


def test_record_end_with_none_has_empty_sha(tmp_path):
    ctx = _make_ctx(tmp_path)
    ledger = RunLedger(ctx, "proj", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()
    ledger.record_end(None)
    assert ledger._record.get("final_adata_sha256") == ""


def test_ctx_without_cfg_does_not_raise(tmp_path):
    ctx = MagicMock()
    del ctx.cfg  # no cfg attribute
    ledger = RunLedger(ctx, "proj", tmp_path)
    ledger.record_start()
    ledger._rss_poller.stop()


def test_filter_env_excludes_path(tmp_path):
    with patch.dict("os.environ", {"PATH": "/usr/bin", "SC_FOO": "bar"}, clear=True):
        result = _filter_env()
    assert "PATH" not in result
    assert "SC_FOO" in result
