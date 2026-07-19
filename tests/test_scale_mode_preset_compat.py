"""Tests that scale_mode preset bundles expand to exact capability flag values.

Guarantees that scale_mode preset → capability flag bundles match the
documented contract. Updated 2026-05-19 for F-4 (Plan
~/.omc/plans/nc-cell-clustering-final-strategy-plan.md, Principle 2): the
`massive` preset's clustering_engine slot was remapped from "css" (removed
from production science path per F-1) to "auto". Resource settings preserve
lazy reading/checkpoint behavior while doublet strategy stays canonical.

Run with:
    pytest -q --no-cov tests/test_scale_mode_preset_compat.py
"""
from __future__ import annotations

from contextlib import nullcontext
import os
import types

import pytest

from workflow.modular.config import PipelineConfig, scale_mode_to_capabilities, _SCALE_MODE_PRESETS
from workflow.modular.cli import (
    DEFAULT_OPTIONAL_MODULES,
    _apply_scale_mode,
    _apply_scientific_profile,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _minimal_args(**overrides):
    """Return a SimpleNamespace that mimics the default argparse.Namespace."""
    defaults = dict(
        scale_mode="standard",
        scientific_profile="canonical",
        acknowledge_scientific_non_equivalence=False,
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
    """scale_mode_to_capabilities('massive') must match the F-4-remapped preset.

    F-4 (2026-05-19, Plan ~/.omc/plans/nc-cell-clustering-final-strategy-plan.md,
    Principle 2): clustering_engine slot was remapped from "css" → "auto".
    Resource flags (lazy_read, checkpoint_policy) are preserved. Scientific
    choices (doublet strategy and clustering method) remain canonical.
    """
    caps = scale_mode_to_capabilities("massive")
    assert caps["lazy_read"] == "true"
    assert "doublet_strategy" not in caps
    assert "clustering_engine" not in caps
    assert caps["checkpoint_policy"] == "full"


# ---------------------------------------------------------------------------
# Test 1b: massive preset must NEVER route to CSS post-F-4
# ---------------------------------------------------------------------------

def test_massive_preset_does_not_route_to_css():
    """F-1 + F-4 invariant: scale_mode=massive must not activate CSS by any route."""
    caps = scale_mode_to_capabilities("massive")
    assert "clustering_engine" not in caps


# ---------------------------------------------------------------------------
# Test 2: _apply_scale_mode expands only resource flags
# ---------------------------------------------------------------------------

def test_apply_scale_mode_massive_expands_only_resource_flags():
    """Massive scale mode must leave scientific/method fields unresolved."""
    args = _minimal_args(scale_mode="massive")
    # F-4 emits a DeprecationWarning when --scale-mode=massive is used.
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore", DeprecationWarning)
        result = _apply_scale_mode(args)
    assert result.lazy_read == "true"
    assert result.checkpoint_policy == "full"
    assert result.doublet_strategy == ""
    assert result.clustering_engine == ""
    assert result.optional_modules == DEFAULT_OPTIONAL_MODULES
    assert result.n_top_genes == 3000
    assert result.n_pcs == 40
    assert result.n_neighbors == 15
    assert result.leiden_resolution == 0.8
    assert result.de_n_genes == 300


@pytest.mark.parametrize("scale_mode", ["standard", "large", "massive"])
def test_resource_strategy_preserves_all_scientific_parameters(scale_mode):
    args = _minimal_args(scale_mode=scale_mode)
    with pytest.warns(DeprecationWarning) if scale_mode == "massive" else nullcontext():
        result = _apply_scientific_profile(_apply_scale_mode(args))

    assert {
        "optional_modules": result.optional_modules,
        "n_top_genes": result.n_top_genes,
        "n_pcs": result.n_pcs,
        "n_neighbors": result.n_neighbors,
        "leiden_resolution": result.leiden_resolution,
        "de_n_genes": result.de_n_genes,
        "doublet_strategy": result.doublet_strategy,
        "clustering_engine": result.clustering_engine,
    } == {
        "optional_modules": DEFAULT_OPTIONAL_MODULES,
        "n_top_genes": 3000,
        "n_pcs": 40,
        "n_neighbors": 15,
        "leiden_resolution": 0.8,
        "de_n_genes": 300,
        "doublet_strategy": "auto",
        "clustering_engine": "auto",
    }
    assert result.resolved_scientific_parameter_diff == {}


def test_scientific_change_requires_non_equivalence_acknowledgement():
    args = _minimal_args(scientific_profile="legacy-large")
    with pytest.raises(SystemExit, match="not scientifically equivalent"):
        _apply_scientific_profile(args)


def test_acknowledged_scientific_profile_records_resolved_diff():
    args = _minimal_args(
        scientific_profile="legacy-massive",
        acknowledge_scientific_non_equivalence=True,
    )

    result = _apply_scientific_profile(args)

    assert result.scientific_non_equivalence_acknowledged is True
    assert result.resolved_scientific_parameter_diff["n_top_genes"] == {
        "canonical": 3000,
        "resolved": 1000,
    }
    assert result.resolved_scientific_parameter_diff["optional_modules"] == {
        "canonical": DEFAULT_OPTIONAL_MODULES,
        "resolved": "clustering",
    }
    assert result.resolved_scientific_parameter_diff["doublet_strategy"] == {
        "canonical": "auto",
        "resolved": "grouped",
    }


# ---------------------------------------------------------------------------
# Test 2b: --scale-mode=massive emits a DeprecationWarning (F-4)
# ---------------------------------------------------------------------------

def test_apply_scale_mode_massive_emits_deprecation_warning():
    """F-4: --scale-mode=massive must emit a DeprecationWarning explaining the css → auto remap."""
    import warnings as _warnings
    args = _minimal_args(scale_mode="massive")
    with _warnings.catch_warnings(record=True) as w:
        _warnings.simplefilter("always")
        _apply_scale_mode(args)
        dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
        assert dep_warnings, (
            "F-4: scale_mode=massive must emit DeprecationWarning explaining "
            "the css → auto remap. None emitted."
        )
        assert any("css" in str(x.message).lower() or "F-1" in str(x.message)
                   for x in dep_warnings), (
            "F-4: DeprecationWarning must explain the css removal context."
        )


# ---------------------------------------------------------------------------
# Test 3: explicit flag overrides the preset bundle
# ---------------------------------------------------------------------------

def test_scale_mode_does_not_overwrite_explicit_scientific_choice():
    """Resource resolution must leave an explicit scientific choice untouched."""
    args = _minimal_args(scale_mode="massive", clustering_engine="sparse_exact")
    import warnings as _warnings
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore", DeprecationWarning)
        result = _apply_scale_mode(args)
    assert result.clustering_engine == "sparse_exact", (
        "explicit --clustering-engine sparse_exact must be independent of scale mode."
    )
    assert result.lazy_read == "true"
    assert result.doublet_strategy == ""
    assert result.checkpoint_policy == "full"


# ---------------------------------------------------------------------------
# Test 4: standard resource preset contains no scientific/method choices
# ---------------------------------------------------------------------------

def test_standard_preset_safe_defaults():
    """scale_mode=standard contains only resource execution settings."""
    caps = scale_mode_to_capabilities("standard")
    assert caps["lazy_read"] == "auto"
    assert caps["checkpoint_policy"] == "full"
    assert set(caps) == {"lazy_read", "checkpoint_policy"}


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
