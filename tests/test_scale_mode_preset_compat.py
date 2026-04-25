"""Tests that scale_mode preset bundles expand to exact capability flag values.

Guarantees that --scale-mode massive still produces identical behavior after
the scale_mode → capability-flag decomposition (Phase 7A.1).

Run with:
    pytest -q --no-cov tests/test_scale_mode_preset_compat.py
"""
from __future__ import annotations

import os
import types

import pytest

from workflow.modular.config import PipelineConfig, scale_mode_to_capabilities, _SCALE_MODE_PRESETS
from workflow.modular.cli import _apply_scale_mode, DEFAULT_OPTIONAL_MODULES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _minimal_args(**overrides):
    """Return a SimpleNamespace that mimics the default argparse.Namespace."""
    defaults = dict(
        scale_mode="standard",
        lazy_read="",
        doublet_strategy="",
        clustering_engine="",
        checkpoint_policy="",
        optional_modules=DEFAULT_OPTIONAL_MODULES,
        n_top_genes=3000,
        n_pcs=40,
        n_neighbors=15,
        leiden_resolution=0.8,
        de_n_genes=300,
        parallel_workers=1,
    )
    defaults.update(overrides)
    return types.SimpleNamespace(**defaults)


# ---------------------------------------------------------------------------
# Test 1: massive preset bundle exact match
# ---------------------------------------------------------------------------

def test_massive_preset_bundle_exact():
    """scale_mode_to_capabilities('massive') must equal the documented preset."""
    caps = scale_mode_to_capabilities("massive")
    assert caps["lazy_read"] == "true"
    assert caps["doublet_strategy"] == "grouped"
    assert caps["clustering_engine"] == "css"
    assert caps["checkpoint_policy"] == "mandatory_only"


# ---------------------------------------------------------------------------
# Test 2: _apply_scale_mode expands massive into all 4 flags
# ---------------------------------------------------------------------------

def test_apply_scale_mode_massive_expands_flags():
    """_apply_scale_mode with massive must set all 4 capability flags."""
    args = _minimal_args(scale_mode="massive")
    result = _apply_scale_mode(args)
    assert result.lazy_read == "true"
    assert result.doublet_strategy == "grouped"
    assert result.clustering_engine == "css"
    assert result.checkpoint_policy == "mandatory_only"


# ---------------------------------------------------------------------------
# Test 3: explicit flag overrides the preset bundle
# ---------------------------------------------------------------------------

def test_explicit_flag_overrides_preset():
    """An explicit capability flag must win over the scale_mode preset."""
    # Override clustering_engine while keeping massive for everything else.
    args = _minimal_args(scale_mode="massive", clustering_engine="sparse_exact")
    result = _apply_scale_mode(args)
    assert result.clustering_engine == "sparse_exact", (
        "explicit --clustering-engine sparse_exact must override massive preset css"
    )
    # Other flags should still come from the massive preset.
    assert result.lazy_read == "true"
    assert result.doublet_strategy == "grouped"
    assert result.checkpoint_policy == "mandatory_only"


# ---------------------------------------------------------------------------
# Test 4: standard preset leaves all capability flags at safe defaults
# ---------------------------------------------------------------------------

def test_standard_preset_safe_defaults():
    """scale_mode=standard must leave capability flags at their safe defaults."""
    caps = scale_mode_to_capabilities("standard")
    assert caps["lazy_read"] == "auto"
    assert caps["doublet_strategy"] == "auto"
    assert caps["clustering_engine"] == "auto"
    assert caps["checkpoint_policy"] == "full"


# ---------------------------------------------------------------------------
# Test 5: SC_DOUBLET_STRATEGY env var overrides cfg.doublet_strategy
# ---------------------------------------------------------------------------

def test_env_var_overrides_cfg_doublet_strategy(monkeypatch):
    """SC_DOUBLET_STRATEGY env var must take priority over cfg.doublet_strategy."""
    import numpy as np
    import pandas as pd
    import anndata as ad
    from unittest.mock import MagicMock
    from workflow.modular.modules.doublet_detection import DoubletDetectionModule

    monkeypatch.setenv("SC_DOUBLET_STRATEGY", "whole")

    # Build a tiny mock context with doublet_strategy=grouped in cfg.
    cfg = MagicMock()
    cfg.doublet_strategy = "grouped"
    ctx = MagicMock()
    ctx.cfg = cfg

    # Create a minimal AnnData with sample column and n_obs >= 100k to trigger
    # the "auto" grouped path — but env var says "whole".
    n = 120000
    obs = pd.DataFrame({"sample": ["A"] * (n // 2) + ["B"] * (n // 2)})
    adata = ad.AnnData(X=np.zeros((n, 10), dtype="float32"), obs=obs)

    strategy = DoubletDetectionModule._resolve_doublet_strategy(ctx, adata)
    assert strategy == "whole", (
        "SC_DOUBLET_STRATEGY=whole must override cfg.doublet_strategy=grouped"
    )
