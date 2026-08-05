"""Annotation claim honesty for DEFAULT_MARKERS starter pack (audit W1.2)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import pytest

from workflow.modular.modules.annotation import (
    AnnotationModule,
    DEFAULT_MARKER_SOURCE,
    DEFAULT_MARKERS,
)


class _Cfg:
    markers = {}
    tissue = "lung"
    annotation_strategy = "cluster_voting"
    annotation_confidence_threshold = -1e9  # keep all assignments
    reference_adata = None
    reference_label_key = "cell_type"
    reference_k = 15
    reference_min_confidence = 0.0
    reference_override_mode = "conservative"


class _Ctx:
    def __init__(self, adata, cfg=None):
        self.adata = adata
        self.cfg = cfg or _Cfg()
        self.metadata = {}
        self.random_state = 0
        import tempfile
        from pathlib import Path

        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        self.figure_dir = base / "figures"
        self.table_dir = base / "tables"
        self.figure_dir.mkdir()
        self.table_dir.mkdir()

    def __del__(self):
        try:
            self._tmp.cleanup()
        except Exception:
            pass


def _brainish_adata():
    # Genes include default pack markers so scoring can run.
    genes = sorted({g for gs in DEFAULT_MARKERS.values() for g in gs})
    # pad a few neural-ish symbols that are NOT the default pack
    genes = genes + ["MAP2", "DCX", "SOX2"]
    n_obs = 30
    rng = np.random.default_rng(0)
    X = rng.poisson(1.0, size=(n_obs, len(genes))).astype(np.float32)
    adata = ad.AnnData(X)
    adata.var_names = genes
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    adata.obs["leiden"] = (["0"] * 10 + ["1"] * 10 + ["2"] * 10)
    adata.obsm["X_umap"] = rng.normal(size=(n_obs, 2)).astype(np.float32)
    return adata


def test_default_markers_are_nonclaimable_even_on_lung_tissue(monkeypatch):
    """Starter pack must never be confirmatory without operator markers."""
    import scanpy as sc

    monkeypatch.setattr(
        sc.tl,
        "score_genes",
        lambda adata, gene_list, score_name=None, **kw: adata.obs.__setitem__(
            score_name, np.linspace(0, 1, adata.n_obs)
        ),
        raising=False,
    )
    # Skip heavy plots
    monkeypatch.setattr(
        AnnotationModule, "_plot_label_umap", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_score_umap", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_marker_umaps", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_write_epithelial_marker_qc", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_try_reference_mapping", staticmethod(lambda *a, **k: None)
    )

    adata = _brainish_adata()
    ctx = _Ctx(adata)
    AnnotationModule().run(ctx)
    assert ctx.metadata["annotation_marker_source"] == DEFAULT_MARKER_SOURCE
    assert ctx.metadata["annotation_claimable"] is False
    assert "starter" in ctx.metadata["annotation_claim_reason"]
    assert adata.uns["annotation"]["claimable"] is False


def test_operator_markers_are_claimable(monkeypatch):
    import scanpy as sc

    monkeypatch.setattr(
        sc.tl,
        "score_genes",
        lambda adata, gene_list, score_name=None, **kw: adata.obs.__setitem__(
            score_name, np.linspace(0, 1, adata.n_obs)
        ),
        raising=False,
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_label_umap", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_score_umap", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_marker_umaps", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_write_epithelial_marker_qc", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_try_reference_mapping", staticmethod(lambda *a, **k: None)
    )

    adata = _brainish_adata()
    cfg = _Cfg()
    cfg.markers = {"Neuron": ["MAP2", "DCX"], "Progenitor": ["SOX2"]}
    cfg.tissue = "brain"
    ctx = _Ctx(adata, cfg=cfg)
    AnnotationModule().run(ctx)
    assert ctx.metadata["annotation_marker_source"] == "operator_markers_json"
    assert ctx.metadata["annotation_claimable"] is True


def test_default_markers_tissue_mismatch_reason(monkeypatch):
    import scanpy as sc

    monkeypatch.setattr(
        sc.tl,
        "score_genes",
        lambda adata, gene_list, score_name=None, **kw: adata.obs.__setitem__(
            score_name, np.linspace(0, 1, adata.n_obs)
        ),
        raising=False,
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_label_umap", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_score_umap", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_plot_marker_umaps", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_write_epithelial_marker_qc", staticmethod(lambda *a, **k: None)
    )
    monkeypatch.setattr(
        AnnotationModule, "_try_reference_mapping", staticmethod(lambda *a, **k: None)
    )

    adata = _brainish_adata()
    cfg = _Cfg()
    cfg.tissue = "brain"
    ctx = _Ctx(adata, cfg=cfg)
    AnnotationModule().run(ctx)
    assert ctx.metadata["annotation_claimable"] is False
    assert "tissue_mismatch" in ctx.metadata["annotation_claim_reason"]
