"""Tests for the refactored modular RNA velocity module."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from workflow.modular.config import VelocityConfig


# ---------------------------------------------------------------------------
# VelocityConfig defaults
# ---------------------------------------------------------------------------


def test_velocity_config_defaults():
    cfg = VelocityConfig()
    assert cfg.min_shared_counts == 20
    assert cfg.n_pcs == 30
    assert cfg.n_neighbors == 30
    assert cfg.n_jobs == 4
    assert cfg.mode == "stochastic"


# ---------------------------------------------------------------------------
# MUTATING_MODULES
# ---------------------------------------------------------------------------


def test_rna_velocity_not_mutating():
    from workflow.modular.pipeline import MUTATING_MODULES

    assert "rna_velocity" not in MUTATING_MODULES
    assert "batch_correction" in MUTATING_MODULES


# ---------------------------------------------------------------------------
# _transfer_velocity_results
# ---------------------------------------------------------------------------


def test_transfer_velocity_results():
    from workflow.modular.modules.rna_velocity import _transfer_velocity_results

    n_cells, n_genes_orig, n_genes_filtered = 100, 500, 200

    adata_main = AnnData(
        sparse.random(n_cells, n_genes_orig, format="csr", dtype=np.float32)
    )
    adata_main.obs_names = [f"cell_{i}" for i in range(n_cells)]

    adata_v = AnnData(
        sparse.random(n_cells, n_genes_filtered, format="csr", dtype=np.float32)
    )
    adata_v.obs_names = adata_main.obs_names.copy()
    adata_v.obs["velocity_confidence"] = np.random.rand(n_cells)
    adata_v.obs["velocity_length"] = np.random.rand(n_cells)
    adata_v.obs["latent_time"] = np.random.rand(n_cells)
    adata_v.obsm["velocity_umap"] = np.random.rand(n_cells, 2)
    adata_v.uns["velocity_params"] = {"mode": "dynamical"}
    # This should NOT be transferred (gene-indexed)
    adata_v.layers["velocity"] = sparse.random(
        n_cells, n_genes_filtered, format="csr", dtype=np.float32
    )

    _transfer_velocity_results(adata_v, adata_main)

    # Cell-level obs transferred
    assert "velocity_confidence" in adata_main.obs
    assert "velocity_length" in adata_main.obs
    assert "latent_time" in adata_main.obs
    np.testing.assert_array_equal(
        adata_main.obs["velocity_confidence"].values,
        adata_v.obs["velocity_confidence"].values,
    )

    # obsm transferred
    assert "velocity_umap" in adata_main.obsm

    # uns transferred
    assert "velocity_params" in adata_main.uns

    # Gene index unchanged
    assert adata_main.n_vars == n_genes_orig

    # Layers NOT transferred
    assert "velocity" not in adata_main.layers


# ---------------------------------------------------------------------------
# numpy 2.x patch version guard
# ---------------------------------------------------------------------------


def test_numpy2_patch_version_guard(monkeypatch):
    from workflow.modular.modules import rna_velocity as mod

    # Simulate numpy 1.x -- patch should not be applied
    monkeypatch.setattr(np, "__version__", "1.26.4")
    result = mod._patch_scvelo_numpy2()
    assert result is False


def test_version_compare_numeric_semantics():
    from workflow.modular.modules.rna_velocity import _version_at_least

    assert _version_at_least("0.10.1", (0, 4))
    assert _version_at_least("1.0.0", (0, 4))
    assert not _version_at_least("0.3.9", (0, 4))


def test_numpy2_patch_skipped_for_scvelo_0_5_plus(monkeypatch):
    from workflow.modular.modules import rna_velocity as mod
    import types
    import sys

    monkeypatch.setattr(np, "__version__", "2.0.0")
    monkeypatch.setattr(mod, "_SCVELO_NUMPY2_PATCH_DONE", False)

    mock_scv = types.ModuleType("scvelo")
    mock_scv.__version__ = "0.5.0"
    monkeypatch.setitem(sys.modules, "scvelo", mock_scv)

    assert mod._patch_scvelo_numpy2() is False


# ---------------------------------------------------------------------------
# COO sparse accumulation
# ---------------------------------------------------------------------------


def test_coo_accumulation():
    """Duplicate (row, col) entries should be summed by COO -> CSR conversion."""
    rows = [0, 0, 1, 1, 1]
    cols = [0, 0, 2, 2, 3]
    data = np.ones(len(rows), dtype=np.float32)
    mat = sparse.coo_matrix((data, (rows, cols)), shape=(3, 5)).tocsr()
    assert mat[0, 0] == 2.0  # two entries summed
    assert mat[1, 2] == 2.0
    assert mat[1, 3] == 1.0
    assert mat[2, 0] == 0.0


# ---------------------------------------------------------------------------
# CLI new args
# ---------------------------------------------------------------------------


def test_velocity_cli_new_args(monkeypatch):
    import workflow.modular.cli as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project", "x",
            "--sample-root", "/tmp/s",
            "--velocity-n-jobs", "8",
            "--velocity-min-shared-counts", "10",
            "--velocity-n-pcs", "50",
            "--velocity-n-neighbors", "20",
        ],
    )
    args = mod.parse_args()
    assert args.velocity_n_jobs == 8
    assert args.velocity_min_shared_counts == 10
    assert args.velocity_n_pcs == 50
    assert args.velocity_n_neighbors == 20


# ---------------------------------------------------------------------------
# GTF auto-discovery from transcriptome_dir
# ---------------------------------------------------------------------------


def test_velocity_gtf_autodiscovery_from_transcriptome_dir(tmp_path):
    from workflow.modular.config import CellRangerConfig, PipelineConfig, VelocityConfig
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.rna_velocity import RNAVelocityModule

    ref = tmp_path / "ref"
    (ref / "genes").mkdir(parents=True)
    gtf = ref / "genes" / "genes.gtf.gz"
    gtf.write_text("dummy", encoding="utf-8")

    cfg = PipelineConfig(
        project="x",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "dataset",
            outs_dir=tmp_path / "outs",
            transcriptome_dir=ref,
        ),
        velocity=VelocityConfig(bam_path=tmp_path / "outs" / "possorted_genome_bam.bam"),
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=tmp_path, table_dir=tmp_path)

    resolved = RNAVelocityModule._resolve_gtf_path(ctx)
    assert resolved == gtf


def test_velocity_gtf_explicit_has_priority(tmp_path):
    from workflow.modular.config import CellRangerConfig, PipelineConfig, VelocityConfig
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.rna_velocity import RNAVelocityModule

    ref = tmp_path / "ref"
    (ref / "genes").mkdir(parents=True)
    auto = ref / "genes" / "genes.gtf.gz"
    auto.write_text("auto", encoding="utf-8")

    explicit = tmp_path / "custom.gtf"
    explicit.write_text("explicit", encoding="utf-8")

    cfg = PipelineConfig(
        project="x",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "dataset",
            outs_dir=tmp_path / "outs",
            transcriptome_dir=ref,
        ),
        velocity=VelocityConfig(
            bam_path=tmp_path / "outs" / "possorted_genome_bam.bam",
            gtf_path=explicit,
        ),
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=tmp_path, table_dir=tmp_path)

    resolved = RNAVelocityModule._resolve_gtf_path(ctx)
    assert resolved == explicit


# ---------------------------------------------------------------------------
# Integration test with monkeypatched scVelo
# ---------------------------------------------------------------------------


def test_rna_velocity_module_integration(monkeypatch, tmp_path):
    """Full module run with mocked scVelo, verify main adata gene index unchanged."""
    from workflow.modular.config import (
        CellRangerConfig,
        PipelineConfig,
        VelocityConfig,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.rna_velocity import RNAVelocityModule

    n_cells, n_genes = 50, 100
    X = sparse.random(n_cells, n_genes, format="csr", dtype=np.float32)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm["X_umap"] = np.random.rand(n_cells, 2).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical(
        [str(i % 3) for i in range(n_cells)]
    )
    # Provide spliced/unspliced directly
    adata.layers["spliced"] = sparse.random(
        n_cells, n_genes, format="csr", dtype=np.float32
    )
    adata.layers["unspliced"] = sparse.random(
        n_cells, n_genes, format="csr", dtype=np.float32
    )

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    mod_dir = run_dir / "rna_velocity"
    mod_dir.mkdir()

    cfg = PipelineConfig(
        project="test",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "dataset",
            outs_dir=tmp_path / "outs",
        ),
        velocity=VelocityConfig(mode="stochastic", min_shared_counts=1),
        optional_modules=["rna_velocity"],
    )

    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=mod_dir,
        table_dir=mod_dir,
        adata=adata,
    )

    # Mock scVelo functions
    import types

    mock_scv = types.ModuleType("scvelo")
    mock_scv.pp = types.SimpleNamespace()
    mock_scv.tl = types.SimpleNamespace()
    mock_scv.pl = types.SimpleNamespace()
    mock_scv.utils = types.SimpleNamespace()

    def fake_filter_and_normalize(adata, min_shared_counts=20):
        # Simulate filtering: remove half the genes
        keep = list(range(0, adata.n_vars, 2))
        subset = adata[:, keep].copy()
        adata.__dict__.update(subset.__dict__)

    def fake_moments(adata, n_pcs=30, n_neighbors=30):
        adata.layers["Ms"] = sparse.random(adata.n_obs, adata.n_vars, format="csr")
        adata.layers["Mu"] = sparse.random(adata.n_obs, adata.n_vars, format="csr")

    def fake_velocity(adata, mode="stochastic"):
        adata.layers["velocity"] = sparse.random(adata.n_obs, adata.n_vars, format="csr")

    def fake_velocity_graph(adata):
        adata.uns["velocity_graph"] = True

    def fake_velocity_confidence(adata):
        adata.obs["velocity_confidence"] = np.random.rand(adata.n_obs)
        adata.obs["velocity_length"] = np.random.rand(adata.n_obs)

    def fake_rank_velocity_genes(adata, groupby=None, n_genes=10):
        adata.uns["rank_velocity_genes"] = {
            "names": np.array([adata.var_names[:n_genes].tolist()] * 3).T
        }

    mock_scv.pp.filter_and_normalize = fake_filter_and_normalize
    mock_scv.pp.moments = fake_moments
    mock_scv.tl.velocity = fake_velocity
    mock_scv.tl.velocity_graph = fake_velocity_graph
    mock_scv.tl.velocity_confidence = fake_velocity_confidence
    mock_scv.tl.rank_velocity_genes = fake_rank_velocity_genes
    mock_scv.DataFrame = pd.DataFrame

    # Noop plotting
    def fake_plot(*args, **kwargs):
        pass

    mock_scv.pl.velocity_embedding_stream = fake_plot
    mock_scv.pl.velocity_embedding_grid = fake_plot
    mock_scv.pl.velocity_embedding = fake_plot
    mock_scv.pl.scatter = fake_plot

    import sys

    monkeypatch.setitem(sys.modules, "scvelo", mock_scv)
    monkeypatch.setitem(sys.modules, "scvelo.pp", mock_scv.pp)
    monkeypatch.setitem(sys.modules, "scvelo.tl", mock_scv.tl)
    monkeypatch.setitem(sys.modules, "scvelo.pl", mock_scv.pl)

    # Patch the numpy2 patch to be a noop
    import workflow.modular.modules.rna_velocity as rna_mod
    monkeypatch.setattr(rna_mod, "_patch_scvelo_numpy2", lambda: False)

    module = RNAVelocityModule()
    module.run(ctx)

    # Main adata gene index should be UNCHANGED
    assert ctx.adata.n_vars == n_genes
    assert list(ctx.adata.var_names) == [f"gene_{i}" for i in range(n_genes)]

    # Cell-level results should be transferred
    assert "velocity_confidence" in ctx.adata.obs
    assert "velocity_length" in ctx.adata.obs

    # Metadata
    assert ctx.metadata["velocity_mode"] == "stochastic"
    assert "velocity_mean_confidence" in ctx.metadata


# ---------------------------------------------------------------------------
# Gene count validation
# ---------------------------------------------------------------------------


def test_gene_count_validation(monkeypatch, tmp_path):
    """If filter_and_normalize removes all genes, raise ValueError."""
    from workflow.modular.config import (
        CellRangerConfig,
        PipelineConfig,
        VelocityConfig,
    )
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.rna_velocity import RNAVelocityModule

    n_cells, n_genes = 20, 50
    adata = AnnData(sparse.random(n_cells, n_genes, format="csr", dtype=np.float32))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm["X_umap"] = np.random.rand(n_cells, 2).astype(np.float32)
    adata.layers["spliced"] = sparse.random(n_cells, n_genes, format="csr", dtype=np.float32)
    adata.layers["unspliced"] = sparse.random(n_cells, n_genes, format="csr", dtype=np.float32)

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    mod_dir = run_dir / "rna_velocity"
    mod_dir.mkdir()

    cfg = PipelineConfig(
        project="test",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        velocity=VelocityConfig(min_shared_counts=20),
        optional_modules=["rna_velocity"],
    )
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=mod_dir, table_dir=mod_dir, adata=adata)

    import types
    import sys

    mock_scv = types.ModuleType("scvelo")
    mock_scv.pp = types.SimpleNamespace()
    mock_scv.tl = types.SimpleNamespace()
    mock_scv.pl = types.SimpleNamespace()
    mock_scv.utils = types.SimpleNamespace()

    def filter_remove_all(adata, min_shared_counts=20):
        # Remove ALL genes
        subset = adata[:, []].copy()
        adata.__dict__.update(subset.__dict__)

    mock_scv.pp.filter_and_normalize = filter_remove_all

    monkeypatch.setitem(sys.modules, "scvelo", mock_scv)

    import workflow.modular.modules.rna_velocity as rna_mod
    monkeypatch.setattr(rna_mod, "_patch_scvelo_numpy2", lambda: False)

    module = RNAVelocityModule()
    with pytest.raises(ValueError, match="filter_and_normalize removed all genes"):
        module.run(ctx)
