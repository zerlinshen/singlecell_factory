from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import pytest


def _make_synthetic_adata() -> ad.AnnData:
    """5 clusters, each with a dominant canonical marker set."""
    n_cells_per_cluster = 80
    n_clusters = 5
    n_cells = n_cells_per_cluster * n_clusters

    # Gene panel: marker genes + enough background genes for sc.tl.score_genes binning (needs >50 control genes)
    marker_genes = [
        "MS4A1", "CD79A", "CD79B", "CD19",           # B cell (cluster 0)
        "CD3D", "CD3E", "TRAC", "CD4",               # T cell (cluster 1)
        "TPSAB1", "TPSB2", "CPA3", "KIT",            # Mast (cluster 2)
        "LYZ", "CD68", "CST3", "FCGR3A",             # Myeloid (cluster 3)
        "GNLY", "NKG7", "KLRD1", "NCAM1",            # NK (cluster 4)
        "EPCAM", "KRT7", "KRT8", "KRT18",            # Tumor
        "DCN", "LUM", "COL1A1", "COL3A1",            # Fibroblast
        "PECAM1", "VWF", "CDH5", "CLDN5",            # Endothelial
        "IGHG1", "IGKC", "MZB1", "XBP1",            # Plasma
        "CD1C", "CLEC9A", "FCER1A", "IRF8",          # DC
        "AIF1", "CD14", "IL7R", "CD8A",              # misc
    ]
    # Pad with unique background genes so sc.tl.score_genes has enough pool for binning
    background_genes = [f"GENE{i:04d}" for i in range(160)]
    genes = list(dict.fromkeys(marker_genes + background_genes))  # dedup preserving order
    n_genes = len(genes)
    gene_idx = {g: i for i, g in enumerate(genes)}

    rng = np.random.default_rng(42)
    X = rng.uniform(0.0, 0.1, size=(n_cells, n_genes)).astype(np.float32)

    # Each cluster strongly expresses its canonical markers
    cluster_markers = [
        ["MS4A1", "CD79A", "CD79B", "CD19"],
        ["CD3D", "CD3E", "TRAC", "CD4"],
        ["TPSAB1", "TPSB2", "CPA3", "KIT"],
        ["LYZ", "CD68", "CST3", "FCGR3A"],
        ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    ]
    leiden_labels = []
    for cluster_id in range(n_clusters):
        start = cluster_id * n_cells_per_cluster
        end = start + n_cells_per_cluster
        for gene in cluster_markers[cluster_id]:
            idx = gene_idx[gene]
            X[start:end, idx] = rng.uniform(3.0, 5.0, size=n_cells_per_cluster).astype(np.float32)
        leiden_labels.extend([str(cluster_id)] * n_cells_per_cluster)

    obs = pd.DataFrame({"leiden": pd.Categorical(leiden_labels)}, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


class _MinimalCfg:
    markers: dict = {}
    annotation_confidence_threshold: float = -999.0  # accept all
    annotation_strategy: str = "cluster_voting"
    reference_adata = None


class _MinimalCtx:
    def __init__(self, strategy: str, tmp_path):
        self.adata = _make_synthetic_adata()
        self.cfg = _MinimalCfg()
        self.cfg.annotation_strategy = strategy
        self.metadata: dict = {}
        self.table_dir = tmp_path
        self.figure_dir = tmp_path

    @property
    def _skip_viz(self):
        return True


def _run_annotation(strategy: str, tmp_path):
    import unittest.mock as mock
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc_real
    from workflow.modular.modules.annotation import AnnotationModule
    import workflow.modular.modules.annotation as ann_mod

    ctx = _MinimalCtx(strategy, tmp_path)
    mod = AnnotationModule()

    _orig_viz = AnnotationModule.__dict__["_plot_composition"]
    _orig_ref = AnnotationModule.__dict__["_try_reference_mapping"]
    try:
        AnnotationModule._plot_composition = staticmethod(lambda adata, ctx: None)
        AnnotationModule._try_reference_mapping = staticmethod(lambda adata, ctx: None)
        with mock.patch.object(sc_real.pl, "umap", return_value=None), \
             mock.patch.object(ann_mod.plt, "savefig", return_value=None), \
             mock.patch.object(ann_mod.plt, "close", return_value=None):
            mod.run(ctx)
    finally:
        AnnotationModule._plot_composition = _orig_viz
        AnnotationModule._try_reference_mapping = _orig_ref
    return ctx


def test_cluster_voting_b_cell(tmp_path):
    """cluster 0 (MS4A1+CD79A+) must map to B cell with cluster_voting."""
    ctx = _run_annotation("cluster_voting", tmp_path)
    majority = pd.read_csv(tmp_path / "cluster_majority_cell_type.csv", index_col=0)
    majority.index = majority.index.astype(str)
    assert majority.loc["0", "majority_cell_type"] == "B cell", (
        f"cluster 0 expected 'B cell', got '{majority.loc['0', 'majority_cell_type']}'"
    )


def test_cluster_voting_t_cell(tmp_path):
    """cluster 1 (CD3D+CD3E+) must map to T cell with cluster_voting."""
    ctx = _run_annotation("cluster_voting", tmp_path)
    majority = pd.read_csv(tmp_path / "cluster_majority_cell_type.csv", index_col=0)
    majority.index = majority.index.astype(str)
    assert majority.loc["1", "majority_cell_type"] == "T cell", (
        f"cluster 1 expected 'T cell', got '{majority.loc['1', 'majority_cell_type']}'"
    )


def test_cluster_voting_mast_cell(tmp_path):
    """cluster 2 (TPSAB1+CPA3+) must map to Mast cell with cluster_voting."""
    ctx = _run_annotation("cluster_voting", tmp_path)
    majority = pd.read_csv(tmp_path / "cluster_majority_cell_type.csv", index_col=0)
    majority.index = majority.index.astype(str)
    assert majority.loc["2", "majority_cell_type"] == "Mast cell", (
        f"cluster 2 expected 'Mast cell', got '{majority.loc['2', 'majority_cell_type']}'"
    )


def test_cluster_voting_writes_score_matrix(tmp_path):
    """cluster_voting must write cluster_score_matrix.csv."""
    _run_annotation("cluster_voting", tmp_path)
    score_matrix_path = tmp_path / "cluster_score_matrix.csv"
    assert score_matrix_path.exists(), "cluster_score_matrix.csv not written"
    df = pd.read_csv(score_matrix_path, index_col=0)
    assert "B cell" in df.columns
    assert len(df) == 5  # 5 clusters


def test_cluster_voting_strategy_recorded(tmp_path):
    """annotation_strategy must be recorded in ctx.metadata."""
    ctx = _run_annotation("cluster_voting", tmp_path)
    assert ctx.metadata.get("annotation_strategy") == "cluster_voting"


def test_cell_argmax_runs(tmp_path):
    """cell_argmax path must complete without error (no label accuracy gate)."""
    ctx = _run_annotation("cell_argmax", tmp_path)
    majority = pd.read_csv(tmp_path / "cluster_majority_cell_type.csv", index_col=0)
    assert len(majority) == 5
    assert ctx.metadata.get("annotation_strategy") == "cell_argmax"
