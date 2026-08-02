"""Tests for trajectory module contract gate (US-W5-8-CONTRACT)."""
import numpy as np
import pandas as pd
import pytest
import anndata as ad

from workflow.modular._contract_violation import ModuleContractError
from workflow.modular.modules.trajectory import TrajectoryModule


def _make_adata(**obsm_kwargs):
    adata = ad.AnnData(np.zeros((5, 3)))
    adata.obs["leiden"] = ["0", "0", "1", "1", "2"]
    for key, val in obsm_kwargs.items():
        adata.obsm[key] = val
    return adata


class _MinimalCtx:
    def __init__(self, adata, *, root_cluster=None, root_justification=None):
        self.adata = adata
        self.metadata = {}

        class _Cfg:
            trajectory_root_cluster = root_cluster
            trajectory_root_justification = root_justification

        self.cfg = _Cfg()

        import tempfile, pathlib
        self._tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self._tmp.name)
        self.figure_dir = base / "figures"
        self.table_dir = base / "tables"
        self.figure_dir.mkdir()
        self.table_dir.mkdir()

    def __del__(self):
        try:
            self._tmp.cleanup()
        except Exception:
            pass


def test_requires_keys_uses_obsm_any_of():
    mod = TrajectoryModule()
    assert "obsm_any_of" in mod.requires_keys
    assert "X_wnn" in mod.requires_keys["obsm_any_of"]
    assert "X_pca" in mod.requires_keys["obsm_any_of"]
    assert "obsm" not in mod.requires_keys or "X_umap" not in mod.requires_keys.get("obsm", [])


def test_raises_module_contract_error_when_both_absent():
    adata = _make_adata()
    ctx = _MinimalCtx(adata)
    with pytest.raises(ModuleContractError, match="X_wnn"):
        TrajectoryModule().run(ctx)


def _patch_scanpy_tl(monkeypatch):
    import scanpy as sc

    def _stub_diffmap(a, **kw):
        a.obsm["X_diffmap"] = np.zeros((a.n_obs, 2), dtype=np.float32)

    def _stub_dpt(a, **kw):
        a.obs["dpt_pseudotime"] = np.linspace(0, 1, a.n_obs)

    monkeypatch.setattr(sc.tl, "diffmap", _stub_diffmap, raising=False)
    monkeypatch.setattr(sc.tl, "dpt", _stub_dpt, raising=False)
    monkeypatch.setattr(sc.tl, "paga", lambda a, **kw: None, raising=False)


def test_x_pca_used_when_only_pca_present(monkeypatch):
    pca_vals = np.random.default_rng(0).random((5, 2)).astype(np.float32)
    adata = _make_adata(X_pca=pca_vals)
    ctx = _MinimalCtx(adata)
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)
    assert ctx.metadata["trajectory_embedding_key"] == "X_pca"
    assert ctx.metadata["trajectory_claimable"] is False
    assert adata.uns["trajectory"]["inference_class"] == "exploratory_unrooted_dpt"


def test_x_wnn_preferred_over_x_pca(monkeypatch):
    wnn_vals = np.ones((5, 2), dtype=np.float32) * 99.0
    pca_vals = np.zeros((5, 2), dtype=np.float32)
    adata = _make_adata(X_wnn=wnn_vals, X_pca=pca_vals)
    ctx = _MinimalCtx(adata)
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)
    assert ctx.metadata["trajectory_embedding_key"] == "X_wnn"


def test_explicit_biologically_justified_root_is_claimable(monkeypatch):
    pca_vals = np.random.default_rng(1).random((5, 2)).astype(np.float32)
    adata = _make_adata(X_pca=pca_vals)
    ctx = _MinimalCtx(
        adata,
        root_cluster="1",
        root_justification="Cluster 1 expresses validated early progenitor markers.",
    )
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)

    provenance = adata.uns["trajectory"]
    assert provenance["claimable"] is True
    assert provenance["root_cluster"] == "1"
    assert provenance["root_cluster_found"] is True
    assert provenance["root_justification"].startswith("Cluster 1")
    assert provenance["inference_class"] == "biologically_rooted_dpt"


@pytest.mark.parametrize(
    ("root_cluster", "root_justification"),
    [
        (None, None),
        ("1", None),
        ("missing", "Expected early state"),
    ],
)
def test_missing_or_invalid_biological_root_is_non_claimable(
    monkeypatch, root_cluster, root_justification
):
    pca_vals = np.random.default_rng(2).random((5, 2)).astype(np.float32)
    adata = _make_adata(X_pca=pca_vals)
    ctx = _MinimalCtx(
        adata,
        root_cluster=root_cluster,
        root_justification=root_justification,
    )
    _patch_scanpy_tl(monkeypatch)
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_umap", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_paga", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_gene_trends", staticmethod(lambda a, c: None))
    monkeypatch.setattr(TrajectoryModule, "_plot_pseudotime_density", staticmethod(lambda a, c: None))

    TrajectoryModule().run(ctx)

    provenance = adata.uns["trajectory"]
    assert provenance["claimable"] is False
    assert provenance["claim_status"] == "non_claimable_missing_biologically_justified_root"
    assert provenance["root_selection_method"] == "exploratory_diffmap_minimum_fallback"
    pseudotime_table = pd.read_csv(ctx.table_dir / "dpt_pseudotime.csv")
    assert not pseudotime_table["trajectory_claimable"].any()
    assert set(pseudotime_table["trajectory_claim_status"]) == {
        "non_claimable_missing_biologically_justified_root"
    }


# ------------------------------------------------- HVG axis alignment (2026-08-02)
# Real-data defect: gene trends are correlated over `adata.raw`, but the HVG flag was
# reindexed straight off `adata.var`. After the ingest boundary converts var_names
# Ensembl->symbol and leaves `.raw` alone, those axes are different namespaces. The
# reindex returned a NON-empty mask of the few features whose symbol was itself an
# Ensembl ID (36/17,764 on the real LUSC file), so the `size == 0` fallback never fired
# and the reported "top pseudotime-correlated genes" came from that handful.
# See governance/raw_axis_namespace_divergence_2026-08-02.md.

_ENS = [f"ENSG{i:011d}" for i in range(8)]
# "ENSG00000000007" is a var_name that is ITSELF an Ensembl ID -- the real file has 36
# such unnamed features -- and it is inside the HVG-flagged range on purpose, because
# that is what makes the naive reindex return a non-empty-but-meaningless mask.
_SYM = ["GZMA", "ENSG00000000007", "CD8A", "MKI67", "EPCAM", "PTPRC", "COL1A1", "PRF1"]


def _diverged_pair(*, keep_ensembl_id=True, n_hvg=6):
    """adata in symbol space, expression matrix still in Ensembl space."""
    var = pd.DataFrame(index=pd.Index(_SYM, dtype=object))
    var["highly_variable"] = [True] * n_hvg + [False] * (len(_SYM) - n_hvg)
    if keep_ensembl_id:
        var["ensembl_id"] = _ENS
    adata = ad.AnnData(np.zeros((4, len(_SYM))), var=var)
    expr_sub = ad.AnnData(
        np.zeros((4, len(_ENS))), var=pd.DataFrame(index=pd.Index(_ENS, dtype=object))
    )
    return adata, expr_sub


def test_hvg_mask_is_rekeyed_through_ensembl_id_when_axes_diverge():
    adata, expr_sub = _diverged_pair()
    ctx = _MinimalCtx(adata)
    mask = TrajectoryModule._hvg_mask_on_expression_axis(adata, expr_sub, ctx)
    assert mask is not None
    assert int(mask.sum()) == 6, "all flagged HVGs must survive alignment"
    assert ctx.metadata["trajectory_hvg_axis_alignment"] == "recovered_via_ensembl_id"


def test_partial_axis_overlap_is_refused_rather_than_silently_used():
    """The exact failure shape: non-empty but meaningless overlap must NOT be trusted."""
    adata, expr_sub = _diverged_pair(keep_ensembl_id=False)
    # 'ENSG00000000007' is both a var_name and an expr var_name, so the naive reindex
    # returns a non-empty mask -- which is precisely what defeated the size==0 guard.
    naive = adata.var["highly_variable"].reindex(expr_sub.var_names, fill_value=False)
    assert 0 < int(naive.to_numpy().sum()) < 6, "fixture must reproduce partial overlap"

    ctx = _MinimalCtx(adata)
    mask = TrajectoryModule._hvg_mask_on_expression_axis(adata, expr_sub, ctx)
    assert mask is None, "an unalignable mask must fall back to the full axis"
    assert ctx.metadata["trajectory_hvg_axis_alignment"] == "unavailable_axis_mismatch"


def test_matching_axes_need_no_realignment_and_record_nothing():
    adata, _ = _diverged_pair()
    ctx = _MinimalCtx(adata)
    mask = TrajectoryModule._hvg_mask_on_expression_axis(adata, adata, ctx)
    assert int(mask.sum()) == 6
    assert "trajectory_hvg_axis_alignment" not in ctx.metadata


def test_absent_or_empty_hvg_flag_yields_no_mask():
    adata, expr_sub = _diverged_pair(n_hvg=0)
    ctx = _MinimalCtx(adata)
    assert TrajectoryModule._hvg_mask_on_expression_axis(adata, expr_sub, ctx) is None
    del adata.var["highly_variable"]
    assert TrajectoryModule._hvg_mask_on_expression_axis(adata, expr_sub, ctx) is None
