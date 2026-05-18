"""Unit tests for _check_requires obsm_any_of primitive (US-W5-8-CONTRACT)."""
import numpy as np
import pytest
import anndata as ad

from workflow.modular.pipeline import _check_requires


def _make_ctx(obsm_keys):
    adata = ad.AnnData(np.zeros((3, 2)))
    for key in obsm_keys:
        adata.obsm[key] = np.zeros((3, 2))

    class _Ctx:
        pass

    ctx = _Ctx()
    ctx.adata = adata
    return ctx


class _ModObsmAnyOf:
    requires_keys = {"obsm_any_of": ["X_wnn", "X_pca"]}


class _ModObsmAll:
    requires_keys = {"obsm": ["X_pca", "X_umap"]}


class _ModMixed:
    requires_keys = {"obsm": ["X_pca"], "obsm_any_of": ["X_wnn", "X_lsi"]}


def test_obsm_any_of_passes_when_first_key_present():
    ctx = _make_ctx(["X_wnn"])
    assert _check_requires(_ModObsmAnyOf(), ctx) == []


def test_obsm_any_of_passes_when_second_key_present():
    ctx = _make_ctx(["X_pca"])
    assert _check_requires(_ModObsmAnyOf(), ctx) == []


def test_obsm_any_of_passes_when_both_keys_present():
    ctx = _make_ctx(["X_wnn", "X_pca"])
    assert _check_requires(_ModObsmAnyOf(), ctx) == []


def test_obsm_any_of_fails_when_none_present():
    ctx = _make_ctx([])
    missing = _check_requires(_ModObsmAnyOf(), ctx)
    assert len(missing) == 1
    assert "obsm_any_of" in missing[0]
    assert "X_wnn" in missing[0]
    assert "X_pca" in missing[0]


def test_obsm_all_fails_when_one_missing():
    ctx = _make_ctx(["X_pca"])
    missing = _check_requires(_ModObsmAll(), ctx)
    assert missing == ["obsm.X_umap"]


def test_obsm_all_passes_when_all_present():
    ctx = _make_ctx(["X_pca", "X_umap"])
    assert _check_requires(_ModObsmAll(), ctx) == []


def test_mixed_obsm_and_obsm_any_of():
    ctx = _make_ctx(["X_pca", "X_lsi"])
    assert _check_requires(_ModMixed(), ctx) == []


def test_mixed_fails_when_obsm_all_missing():
    ctx = _make_ctx(["X_lsi"])
    missing = _check_requires(_ModMixed(), ctx)
    assert "obsm.X_pca" in missing


def test_mixed_fails_when_obsm_any_of_missing():
    ctx = _make_ctx(["X_pca"])
    missing = _check_requires(_ModMixed(), ctx)
    assert any("obsm_any_of" in m for m in missing)


def test_none_adata_returns_empty():
    class _Ctx:
        adata = None

    assert _check_requires(_ModObsmAnyOf(), _Ctx()) == []
