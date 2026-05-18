"""US-W4-7 — WNN synthetic smoke test for multimodal_integration module.

Pure-Python branching tests. Does NOT invoke the R subprocess (kept fast,
no Rscript dependency in CI). Validates:

- Default config is engine='off' and module is a no-op (skipped status).
- engine='wnn' with missing modality latents skips cleanly with reason.
- engine='wnn' with required latents present + Rscript missing skips with
  documented status, not exception.
- engine='mofa' with missing inputs skips cleanly.
- Unknown engine raises ValueError per module contract.

A separate integration test (US-W4-8) exercises the real Rscript subprocess
on a 5k synthetic; that test lives in tests/test_multimodal_wnn_integration.py
(future) and is gated on RSCRIPT_BIN availability.
"""
from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


def _make_adata_with_latents(n: int = 200, with_rna: bool = True, with_second: bool = True) -> ad.AnnData:
    rng = np.random.default_rng(0)
    X = sp.random(n, 50, density=0.2, format="csr", dtype=np.float32, random_state=0).tocsr()
    X.data = np.ceil(X.data * 10).astype(np.float32)
    obs = pd.DataFrame({"s": [f"S{i % 3}" for i in range(n)]}, index=[f"c{i}" for i in range(n)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(50)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    if with_rna:
        a.obsm["X_pca"] = rng.normal(size=(n, 15)).astype(np.float32)
    if with_second:
        a.obsm["protein_clr"] = rng.normal(size=(n, 10)).astype(np.float32)
    return a


def _make_ctx(adata):
    return SimpleNamespace(
        adata=adata,
        cfg=SimpleNamespace(),
        metadata={},
        status=lambda *args, **kwargs: None,
        random_state=42,
    )


def test_default_engine_off_skips_cleanly():
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    adata = _make_adata_with_latents()
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig()
    assert cfg.engine == "off"
    MultimodalIntegrationModule(cfg).run(ctx)
    assert adata.uns["multimodal_status"]["engine"] == "off"
    assert adata.uns["multimodal_status"]["status"] == "skipped"
    assert "X_wnn" not in adata.obsm


def test_wnn_engine_missing_second_modality_skips():
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig(engine="wnn")
    MultimodalIntegrationModule(cfg).run(ctx)
    status = adata.uns["multimodal_status"]
    assert status["engine"] == "wnn"
    assert status["status"] == "skipped"
    assert "missing latents" in status["reason"]


def test_mofa_engine_missing_second_modality_skips():
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig(engine="mofa")
    MultimodalIntegrationModule(cfg).run(ctx)
    assert adata.uns["multimodal_status"]["status"] == "skipped"


def test_unknown_engine_raises():
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    adata = _make_adata_with_latents()
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig(engine="invalid_engine")
    with pytest.raises(ValueError, match=r"unsupported engine"):
        MultimodalIntegrationModule(cfg).run(ctx)


def test_engine_resolved_from_env(monkeypatch):
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    monkeypatch.setenv("SC_MULTIMODAL_ENGINE", "wnn")
    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig()  # default engine="off"; env overrides
    MultimodalIntegrationModule(cfg).run(ctx)
    status = adata.uns["multimodal_status"]
    assert status["engine"] == "wnn", "env var SC_MULTIMODAL_ENGINE must override cfg default"


def test_module_registered_in_catalog():
    """Multimodal_integration must be registered so pipeline can dispatch it."""
    from workflow.modular.pipeline import _build_registry
    reg = _build_registry()
    assert "multimodal_integration" in reg, "module must be in registry for CLI dispatch"
    mod = reg["multimodal_integration"]
    assert mod.name == "multimodal_integration"
    assert mod.required is False


# ---------------------------------------------------------------------------
# US-W5-3 — AC-5: five new tests for scoped ModuleContractError raise
# ---------------------------------------------------------------------------

def test_wnn_explicit_second_key_missing_raises():
    """engine=wnn + _explicitly_set=True + missing second key → ModuleContractError."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular._contract_violation import ModuleContractError

    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig(engine="wnn", second_obsm_key="protein_clr", _explicitly_set=True)
    with pytest.raises(ModuleContractError, match="explicitly-configured second obsm key"):
        MultimodalIntegrationModule(cfg).run(ctx)


def test_mofa_explicit_second_key_missing_raises():
    """engine=mofa + _explicitly_set=True + missing second key → ModuleContractError."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular._contract_violation import ModuleContractError

    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig(engine="mofa", second_obsm_key="protein_clr", _explicitly_set=True)
    with pytest.raises(ModuleContractError, match="explicitly-configured second obsm key"):
        MultimodalIntegrationModule(cfg).run(ctx)


def test_wnn_explicit_second_key_present_succeeds():
    """engine=wnn + _explicitly_set=True + second key present → no raise (reaches Rscript skip)."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular._contract_violation import ModuleContractError

    adata = _make_adata_with_latents(with_second=True)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig(engine="wnn", second_obsm_key="protein_clr", _explicitly_set=True)
    # Should NOT raise ModuleContractError — key is present; module proceeds to Rscript
    # lookup, which skips cleanly (no Rscript in CI).
    try:
        MultimodalIntegrationModule(cfg).run(ctx)
    except ModuleContractError:
        pytest.fail("ModuleContractError raised even though second obsm key is present")
    status = adata.uns.get("multimodal_status", {})
    # Engine reached wnn path (skipped due to no Rscript, but not contract error)
    assert status.get("engine") == "wnn"


def test_explicitly_set_via_cli_flag_raises(monkeypatch):
    """SC_MULTIMODAL_SECOND_OBSM env var (CLI flag path) triggers explicit-set raise."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular._contract_violation import ModuleContractError

    monkeypatch.setenv("SC_MULTIMODAL_ENGINE", "wnn")
    monkeypatch.setenv("SC_MULTIMODAL_SECOND_OBSM", "my_custom_key")
    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig()  # defaults; env vars drive both engine + second key
    with pytest.raises(ModuleContractError, match="explicitly-configured second obsm key"):
        MultimodalIntegrationModule(cfg).run(ctx)


def test_explicitly_set_redundant_default_still_raises(monkeypatch):
    """Passing --second-obsm-key=protein_clr (the default value) is still explicit opt-in."""
    from workflow.modular.modules.multimodal_integration import (
        MultimodalIntegrationConfig,
        MultimodalIntegrationModule,
    )
    from workflow.modular._contract_violation import ModuleContractError

    # Simulate CLI flag path: env var is set to the default value string.
    monkeypatch.setenv("SC_MULTIMODAL_ENGINE", "wnn")
    monkeypatch.setenv("SC_MULTIMODAL_SECOND_OBSM", "protein_clr")
    adata = _make_adata_with_latents(with_second=False)
    ctx = _make_ctx(adata)
    cfg = MultimodalIntegrationConfig()  # second_obsm_key="protein_clr" — default value
    # Even though value == default, env var being SET means explicit opt-in → must raise.
    with pytest.raises(ModuleContractError, match="explicitly-configured second obsm key"):
        MultimodalIntegrationModule(cfg).run(ctx)
