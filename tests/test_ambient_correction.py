"""Unit tests for ambient_correction.py (AmbientCorrectionModule).

All tests mock the subprocess boundary -- no R or conda is invoked.
Covers: T1 top50 trigger, T2 mt trigger, T4 count-correlation trigger,
T5 cross-sample CV trigger, dry-run mode, env-disable flag, R-failure mock,
and the fail-loud missing-QC path (Principle 9).

NOTE: T3 (doublet_excess) was removed 2026-05-20 -- DAG ordering makes
obs["predicted_doublet"] unavailable when triggers are evaluated.
No T3 test exists because T3 no longer exists.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

import anndata as ad
from scipy.sparse import csr_matrix


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_adata(
    n_obs: int = 50,
    n_vars: int = 20,
    top50: float = 30.0,
    mt_pct: float = 5.0,
    total_counts_factor: float = 1.0,
    sample_col: str | None = None,
    n_samples: int = 2,
) -> ad.AnnData:
    """Build a minimal AnnData with QC metrics pre-populated."""
    rng = np.random.default_rng(42)
    X = csr_matrix(rng.poisson(2.0, (n_obs, n_vars)).astype(np.float32))
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    obs["pct_counts_in_top_50"] = top50
    obs["pct_counts_mt"] = mt_pct
    # total_counts and n_genes_by_counts are tightly correlated by default
    base = rng.poisson(100, n_obs).astype(float) * total_counts_factor
    obs["total_counts"] = base
    obs["n_genes_by_counts"] = (base * 0.9 + rng.normal(0, 1, n_obs)).clip(1)
    if sample_col is not None:
        labels = [f"sample_{i % n_samples}" for i in range(n_obs)]
        obs[sample_col] = labels
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    return ad.AnnData(X=X, obs=obs, var=var)


def _make_ctx(adata: ad.AnnData, cfg_overrides: dict | None = None) -> SimpleNamespace:
    """Build a minimal PipelineContext-like namespace."""
    from workflow.modular.config import AmbientCorrectionConfig
    cfg_kwargs = cfg_overrides or {}
    ambient_cfg = AmbientCorrectionConfig(**cfg_kwargs)
    pipeline_cfg = SimpleNamespace(ambient=ambient_cfg)
    ctx = SimpleNamespace(
        adata=adata,
        cfg=pipeline_cfg,
        metadata={},
        module_status=[],
    )
    return ctx


# ---------------------------------------------------------------------------
# T1: high top50 percentile fires
# ---------------------------------------------------------------------------

def test_t1_top50_triggers(tmp_path):
    """T1 fires when median pct_counts_in_top_50 > trigger_top50 (default 50.0)."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    adata = _make_adata(top50=65.0, mt_pct=3.0)  # T1 fires, T2 does not
    ctx = _make_ctx(adata, {"dry_run_triggers_only": True})

    mod = AmbientCorrectionModule()
    mod.run(ctx)

    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "dry_run_triggers_only"
    assert "T1_top50_high" in rec["triggers_fired"]
    assert "T2_mt_excess" not in rec["triggers_fired"]


# ---------------------------------------------------------------------------
# T2: high mt percent fires
# ---------------------------------------------------------------------------

def test_t2_mt_triggers(tmp_path):
    """T2 fires when median pct_counts_mt > trigger_mt (default 15.0)."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    adata = _make_adata(top50=30.0, mt_pct=20.0)  # T2 fires, T1 does not
    ctx = _make_ctx(adata, {"dry_run_triggers_only": True})

    mod = AmbientCorrectionModule()
    mod.run(ctx)

    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "dry_run_triggers_only"
    assert "T2_mt_excess" in rec["triggers_fired"]
    assert "T1_top50_high" not in rec["triggers_fired"]


# ---------------------------------------------------------------------------
# T4: low count correlation fires
# ---------------------------------------------------------------------------

def test_t4_low_count_correlation_triggers(tmp_path):
    """T4 fires when Spearman(total_counts, n_genes) < trigger_count_correlation (0.85)."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    rng = np.random.default_rng(0)
    adata = _make_adata(top50=30.0, mt_pct=5.0)
    # Destroy total_counts/n_genes correlation by randomising one column.
    adata.obs["total_counts"] = rng.uniform(10, 1000, len(adata.obs))
    adata.obs["n_genes_by_counts"] = rng.uniform(10, 1000, len(adata.obs))

    ctx = _make_ctx(adata, {"dry_run_triggers_only": True, "trigger_count_correlation": 0.85})
    mod = AmbientCorrectionModule()
    mod.run(ctx)

    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "dry_run_triggers_only"
    assert "T4_low_count_correlation" in rec["triggers_fired"]


# ---------------------------------------------------------------------------
# T5: high cohort CV fires
# ---------------------------------------------------------------------------

def test_t5_cross_sample_variance_triggers(tmp_path):
    """T5 fires when housekeeping gene CV across samples > trigger_cross_sample_cv."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    rng = np.random.default_rng(7)
    n_obs = 40
    n_vars = 5
    # Gene names include housekeeping genes T5 checks for.
    var_names = ["ACTB", "GAPDH", "B2M", "gene_3", "gene_4"]
    # Sample 0 has very high expression, sample 1 has very low.
    X_data = np.zeros((n_obs, n_vars), dtype=np.float32)
    X_data[:20, :3] = 100.0   # sample_0: high housekeeping
    X_data[20:, :3] = 1.0    # sample_1: low housekeeping
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    obs["pct_counts_in_top_50"] = 30.0
    obs["pct_counts_mt"] = 5.0
    obs["total_counts"] = rng.poisson(100, n_obs).astype(float)
    obs["n_genes_by_counts"] = rng.poisson(50, n_obs).astype(float)
    obs["sample_id"] = ["sample_0"] * 20 + ["sample_1"] * 20
    var = pd.DataFrame(index=var_names)
    adata = ad.AnnData(X=csr_matrix(X_data), obs=obs, var=var)

    ctx = _make_ctx(adata, {"dry_run_triggers_only": True, "trigger_cross_sample_cv": 0.10})
    mod = AmbientCorrectionModule()
    mod.run(ctx)

    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "dry_run_triggers_only"
    assert "T5_cross_sample_variance" in rec["triggers_fired"]


# ---------------------------------------------------------------------------
# Dry-run mode: triggers fire but R not invoked
# ---------------------------------------------------------------------------

def test_dry_run_no_r_invocation(tmp_path):
    """With dry_run_triggers_only=True, triggers are recorded but R is never called."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule
    import workflow.modular.modules.ambient_correction as ac_mod

    adata = _make_adata(top50=80.0, mt_pct=5.0)
    ctx = _make_ctx(adata, {"dry_run_triggers_only": True})

    invoked = []
    original_invoke = AmbientCorrectionModule._invoke_decontx

    def mock_invoke(self, *a, **kw):
        invoked.append(True)
        return original_invoke(self, *a, **kw)

    mod = AmbientCorrectionModule()
    with patch.object(AmbientCorrectionModule, "_invoke_decontx", mock_invoke):
        mod.run(ctx)

    assert len(invoked) == 0, "_invoke_decontx must not be called in dry-run mode"
    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "dry_run_triggers_only"
    assert "T1_top50_high" in rec["triggers_fired"]


# ---------------------------------------------------------------------------
# Env-disable: SC_AMBIENT_TRIGGERS_DISABLE=1 skips everything
# ---------------------------------------------------------------------------

def test_env_disable_skips_all_triggers(monkeypatch, tmp_path):
    """SC_AMBIENT_TRIGGERS_DISABLE=1 returns disabled_by_user without evaluating triggers."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    monkeypatch.setenv("SC_AMBIENT_TRIGGERS_DISABLE", "1")
    adata = _make_adata(top50=90.0, mt_pct=30.0)  # both T1 and T2 would fire
    ctx = _make_ctx(adata)

    mod = AmbientCorrectionModule()
    mod.run(ctx)

    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "disabled_by_user"
    assert rec["triggers_fired"] == []
    # X must not be mutated
    assert "decontx_contamination" not in adata.obs.columns


# ---------------------------------------------------------------------------
# R-failure mock: non-zero exit records decontx_failed decision
# ---------------------------------------------------------------------------

def test_r_failure_raises_runtime_error(monkeypatch, tmp_path):
    """When R subprocess returns non-zero, RuntimeError is raised and adata.X is not mutated."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule
    import shutil as _shutil

    # T1 fires so the module will attempt to call R.
    adata = _make_adata(top50=80.0, mt_pct=5.0)
    original_X_sum = float(np.asarray(adata.X.todense()).sum())

    ctx = _make_ctx(adata)

    # Mock conda on PATH so preflight passes.
    monkeypatch.setattr(_shutil, "which", lambda name: "/usr/bin/conda" if name == "conda" else None)

    # Mock _resolve_r_driver_path to return a dummy path that "exists".
    fake_r_script = tmp_path / "decontx_run.R"
    fake_r_script.write_text("# fake\n")
    import workflow.modular.modules.ambient_correction as ac_mod
    monkeypatch.setattr(ac_mod, "_resolve_r_driver_path", lambda: fake_r_script)

    # Mock Popen: return exit code 1, stderr "R error".
    class FakePopen:
        pid = 99999
        returncode = 1

        def communicate(self, timeout=None):
            return ("", "R subprocess error: something went wrong")

        def poll(self):
            return self.returncode

    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: FakePopen())

    import os as _os
    monkeypatch.setattr(_os, "killpg", lambda pgid, sig: None)
    monkeypatch.setattr(_os, "getpgid", lambda pid: pid)

    mod = AmbientCorrectionModule()
    with pytest.raises(RuntimeError, match="DecontX R driver exited 1"):
        mod.run(ctx)

    # adata.X must not have been modified.
    assert float(np.asarray(adata.X.todense()).sum()) == pytest.approx(original_X_sum)
    assert "decontx_contamination" not in adata.obs.columns


# ---------------------------------------------------------------------------
# Missing-QC fail-loud (Principle 9)
# ---------------------------------------------------------------------------

def test_missing_qc_metrics_raises_fail_loud(tmp_path):
    """When both pct_counts_in_top_50 and pct_counts_mt are absent, raise RuntimeError
    with decision='qc_metrics_unavailable' recorded in adata.uns. Principle 9:
    silent skip forbidden."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    # Build AnnData WITHOUT the key QC columns.
    rng = np.random.default_rng(1)
    X = csr_matrix(rng.poisson(2, (30, 10)).astype(np.float32))
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(30)])
    # Deliberately omit pct_counts_in_top_50 and pct_counts_mt.
    obs["total_counts"] = rng.poisson(100, 30).astype(float)
    obs["n_genes_by_counts"] = rng.poisson(50, 30).astype(float)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    ctx = _make_ctx(adata)
    mod = AmbientCorrectionModule()

    with pytest.raises(RuntimeError, match="qc_metrics_unavailable|required QC metrics"):
        mod.run(ctx)

    # Provenance must be recorded before the raise.
    assert "ambient_correction" in adata.uns
    assert adata.uns["ambient_correction"]["decision"] == "qc_metrics_unavailable"
    assert adata.uns["ambient_correction"]["triggers_fired"] == []


# ---------------------------------------------------------------------------
# No-trigger: both metrics present but below thresholds -> skipped_no_trigger
# ---------------------------------------------------------------------------

def test_no_trigger_records_skipped(tmp_path):
    """When QC metrics are present but no trigger fires, decision is skipped_no_trigger."""
    from workflow.modular.modules.ambient_correction import AmbientCorrectionModule

    adata = _make_adata(top50=20.0, mt_pct=3.0)  # both well below thresholds
    # Ensure T4 does not fire either: high correlation by construction in _make_adata.
    ctx = _make_ctx(adata)

    mod = AmbientCorrectionModule()
    mod.run(ctx)

    rec = adata.uns["ambient_correction"]
    assert rec["decision"] == "skipped_no_trigger"
    assert rec["triggers_fired"] == []
    assert "decontx_contamination" not in adata.obs.columns
