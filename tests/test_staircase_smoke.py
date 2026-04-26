"""Phase 7C smoke tests: real-data tier sanity checks (opt-in via markers).

These tests are SKIPPED by default. They require the staircase fixtures to be
generated first via:

    python scripts/build_staircase_fixtures.py --tiers small_real medium_real

To run:
    pytest -m small_real tests/test_staircase_smoke.py
    pytest -m medium_real tests/test_staircase_smoke.py
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

ad = pytest.importorskip("anndata")

ROOT = Path(__file__).resolve().parents[1]
STAIRCASE_DIR = ROOT / "tests" / "data" / "staircase"


def _load_tier(tier: str):
    # build_staircase_fixtures.py writes <tier>.h5ad with anndata-native encoding.
    h5ad_canonical = STAIRCASE_DIR / f"{tier}.h5ad"
    h5ad_legacy = STAIRCASE_DIR / f"{tier}_anndata.h5ad"
    zarr_named = STAIRCASE_DIR / f"{tier}.zarr"
    if h5ad_canonical.exists():
        return ad.read_h5ad(h5ad_canonical)
    if h5ad_legacy.exists():
        return ad.read_h5ad(h5ad_legacy)
    if zarr_named.exists():
        return ad.read_zarr(zarr_named)
    pytest.skip(
        f"{tier} fixture not found (looked for {h5ad_canonical}, {h5ad_legacy}, {zarr_named}). "
        f"Run: python scripts/build_staircase_fixtures.py --tiers {tier}"
    )


@pytest.fixture(scope="module")
def small_real_adata():
    return _load_tier("small_real")


@pytest.fixture(scope="module")
def medium_real_adata():
    return _load_tier("medium_real")


# ---------------------------------------------------------------------------
# small_real (~100k cells, ~10 samples)
# ---------------------------------------------------------------------------


@pytest.mark.small_real
def test_small_real_fixture_has_required_obs_columns(small_real_adata):
    required = {"sample"}
    missing = required - set(small_real_adata.obs.columns)
    assert not missing, f"small_real missing obs columns: {missing}"


@pytest.mark.small_real
def test_small_real_fixture_size_at_least_50k(small_real_adata):
    assert small_real_adata.n_obs >= 50_000, (
        f"small_real has only {small_real_adata.n_obs} cells; "
        "expected at least 50k to exercise grouped Scrublet path."
    )


@pytest.mark.small_real
def test_doublet_strategy_auto_picks_grouped_at_100k(small_real_adata):
    """At >= 100k cells with sample column, auto strategy resolves to grouped."""
    from types import SimpleNamespace

    from workflow.modular.modules.doublet_detection import DoubletDetectionModule

    cfg = SimpleNamespace(doublet_strategy="auto")
    ctx = SimpleNamespace(cfg=cfg)
    strategy = DoubletDetectionModule._resolve_doublet_strategy(ctx, small_real_adata)
    if small_real_adata.n_obs >= 100_000:
        assert strategy == "grouped", (
            f"auto with n_obs={small_real_adata.n_obs} "
            f"and sample col present should pick 'grouped', got {strategy!r}"
        )


# ---------------------------------------------------------------------------
# medium_real (~500k cells, ~50 samples)
# ---------------------------------------------------------------------------


@pytest.mark.medium_real
def test_medium_real_fixture_size_at_least_300k(medium_real_adata):
    assert medium_real_adata.n_obs >= 300_000, (
        f"medium_real has only {medium_real_adata.n_obs} cells; "
        "expected at least 300k."
    )


@pytest.mark.medium_real
def test_medium_real_fixture_has_multiple_diseases(medium_real_adata):
    """medium_real should span multiple disease categories for cohort_subset tests."""
    if "disease" not in medium_real_adata.obs.columns:
        pytest.skip("disease column not present in medium_real fixture")
    n_disease = medium_real_adata.obs["disease"].nunique()
    assert n_disease >= 2, f"medium_real has only {n_disease} disease categories"
