"""Wave 3 / US-W3-2 — ClusteringContractViolation + orchestrator poison guard.

Covers:
- ClusteringContractViolation is an Exception subclass (not BaseException-only).
- poison_adata sets X=None + uns["__corrupted__"]=True.
- assert_not_corrupted raises on poisoned adata.
- resolve_gpu_failure_policy honors env var > cfg > default.
- pipeline._run_module fires the guard at module entry.
"""
from __future__ import annotations

import os
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


def _make_adata():
    X = sp.random(20, 10, density=0.3, dtype=np.float32, random_state=0).tocsr()
    obs = pd.DataFrame({"s": [f"S{i}" for i in range(20)]}, index=[f"c{i}" for i in range(20)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_contract_violation_is_exception_subclass():
    from workflow.modular._contract_violation import ClusteringContractViolation
    assert issubclass(ClusteringContractViolation, Exception)
    assert ClusteringContractViolation is not Exception
    assert ClusteringContractViolation.__mro__[1] is Exception, (
        "ClusteringContractViolation must inherit directly from Exception, "
        "not BaseException — legitimate finally cleanup must still run."
    )


def test_poison_adata_sets_X_none_and_flag():
    from workflow.modular._contract_violation import poison_adata, CORRUPTED_FLAG
    adata = _make_adata()
    assert adata.X is not None
    poison_adata(adata, reason="test")
    assert adata.X is None
    assert adata.uns[CORRUPTED_FLAG] is True
    assert adata.uns[CORRUPTED_FLAG + "_reason"] == "test"


def test_assert_not_corrupted_raises_on_poisoned_adata():
    from workflow.modular._contract_violation import (
        ClusteringContractViolation,
        assert_not_corrupted,
        poison_adata,
    )
    adata = _make_adata()
    assert_not_corrupted(adata, "some_module")  # no raise on healthy
    poison_adata(adata, reason="gpu failure")
    with pytest.raises(ClusteringContractViolation) as exc:
        assert_not_corrupted(adata, "next_module")
    assert "next_module" in str(exc.value)
    assert "gpu failure" in str(exc.value)


def test_assert_not_corrupted_silent_on_healthy_or_none():
    from workflow.modular._contract_violation import assert_not_corrupted
    assert_not_corrupted(None, "any_module")  # tolerate None
    assert_not_corrupted(_make_adata(), "any_module")  # tolerate healthy


def test_resolve_policy_default_is_derived_from_gpu_mode(monkeypatch):
    """An opportunistic GPU choice must be retractable; a forced one must not be.

    `--gpu-mode auto` lets the pipeline elect the GPU on the operator's behalf,
    so a GPU-only fault must fall back to CPU rather than destroying the run
    (rsc.pp.pca is a documented cuSOLVER incompatibility on this workstation —
    docs/HOTSPOT1_DIAGNOSIS.md). `--gpu-mode force` is an explicit operator
    demand, so it still raises rather than silently delivering a CPU result.
    """
    from workflow.modular._contract_violation import resolve_gpu_failure_policy
    monkeypatch.delenv("SC_GPU_FAILURE_POLICY", raising=False)

    assert resolve_gpu_failure_policy(SimpleNamespace()) == "restore-cpu"
    assert resolve_gpu_failure_policy(SimpleNamespace(gpu_mode="auto")) == "restore-cpu"
    assert resolve_gpu_failure_policy(SimpleNamespace(gpu_mode="off")) == "restore-cpu"
    assert resolve_gpu_failure_policy(SimpleNamespace(gpu_mode="force")) == "raise"


def test_resolve_policy_explicit_cfg_overrides_derived_default(monkeypatch):
    from workflow.modular._contract_violation import resolve_gpu_failure_policy
    monkeypatch.delenv("SC_GPU_FAILURE_POLICY", raising=False)
    cfg = SimpleNamespace(gpu_mode="auto", gpu_failure_policy="raise")
    assert resolve_gpu_failure_policy(cfg) == "raise"


def test_resolve_policy_env_overrides_cfg(monkeypatch):
    from workflow.modular._contract_violation import resolve_gpu_failure_policy
    monkeypatch.setenv("SC_GPU_FAILURE_POLICY", "restore-cpu")
    cfg = SimpleNamespace(gpu_failure_policy="reload-checkpoint")
    assert resolve_gpu_failure_policy(cfg) == "restore-cpu"


def test_resolve_policy_cfg_when_env_absent(monkeypatch):
    from workflow.modular._contract_violation import resolve_gpu_failure_policy
    monkeypatch.delenv("SC_GPU_FAILURE_POLICY", raising=False)
    cfg = SimpleNamespace(gpu_failure_policy="reload-checkpoint")
    assert resolve_gpu_failure_policy(cfg) == "reload-checkpoint"


def test_resolve_policy_rejects_invalid(monkeypatch):
    from workflow.modular._contract_violation import resolve_gpu_failure_policy
    monkeypatch.setenv("SC_GPU_FAILURE_POLICY", "bogus")
    cfg = SimpleNamespace()
    with pytest.raises(ValueError):
        resolve_gpu_failure_policy(cfg)


def test_pipeline_run_module_fires_guard_on_poisoned_adata():
    """_run_module aborts before delegating .run() when adata is poisoned."""
    from workflow.modular._contract_violation import (
        ClusteringContractViolation,
        poison_adata,
    )
    from workflow.modular import pipeline as pipe

    adata = _make_adata()
    poison_adata(adata, reason="upstream contract violation")

    called = {"yes": False}

    class _Mod:
        name = "downstream_mod"
        def run(self, ctx):
            called["yes"] = True

    ctx = SimpleNamespace(adata=adata, metadata={}, status=lambda *a, **k: None)
    with pytest.raises(ClusteringContractViolation):
        pipe._run_module(_Mod(), ctx, mandatory=False)
    assert called["yes"] is False, "guard must abort BEFORE delegating to .run()"


def test_config_carries_gpu_failure_policy_default():
    """PipelineConfig exposes gpu_failure_policy defaulting to "" (derive from gpu_mode).

    The empty default is what lets resolve_gpu_failure_policy() distinguish
    "operator explicitly chose a policy" from "nobody chose", which is required
    to derive raise-vs-restore-cpu from gpu_mode.
    """
    from workflow.modular.config import PipelineConfig
    # Verify the field exists with the expected default via dataclass fields.
    import dataclasses
    fields = {f.name: f for f in dataclasses.fields(PipelineConfig)}
    assert "gpu_failure_policy" in fields, "PipelineConfig must expose gpu_failure_policy"
    assert fields["gpu_failure_policy"].default == ""
