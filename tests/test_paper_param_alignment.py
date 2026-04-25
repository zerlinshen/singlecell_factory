"""Tests verifying NC2024 paper-aligned parameter defaults and cohort subset behavior."""
from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. DoubletConfig: n_prin_comps=30, expected_doublet_rate=0.06
# ---------------------------------------------------------------------------
def test_doublet_config_paper_defaults():
    from workflow.modular.config import DoubletConfig
    cfg = DoubletConfig()
    assert cfg.n_prin_comps == 30, "Paper uses n_prin_comps=30 for Scrublet"
    assert cfg.expected_doublet_rate == 0.06, "Paper uses expected_doublet_rate=0.06"


# ---------------------------------------------------------------------------
# 2. ClusteringConfig: n_pcs=15, leiden_resolution=1.0
# ---------------------------------------------------------------------------
def test_clustering_config_paper_defaults():
    from workflow.modular.config import ClusteringConfig
    cfg = ClusteringConfig()
    assert cfg.n_pcs == 15, "Paper clusters on 15-PC Harmony space"
    assert cfg.leiden_resolution == 1.0, "Paper uses Leiden resolution=1.0"


# ---------------------------------------------------------------------------
# 3. DEConfig: marker_correction default is 'benjamini-hochberg' (backward compat)
# ---------------------------------------------------------------------------
def test_de_config_default_correction():
    from workflow.modular.config import DEConfig
    cfg = DEConfig()
    assert cfg.marker_correction == "benjamini-hochberg"


# ---------------------------------------------------------------------------
# 4. PseudobulkDEConfig: padj_threshold=0.05, abs_logfc_threshold=1.0
# ---------------------------------------------------------------------------
def test_pseudobulk_de_config_paper_defaults():
    from workflow.modular.config import PseudobulkDEConfig
    cfg = PseudobulkDEConfig()
    assert cfg.padj_threshold == 0.05
    assert cfg.abs_logfc_threshold == 1.0


# ---------------------------------------------------------------------------
# 5. PipelineConfig.de_correction default is 'benjamini-hochberg'
# ---------------------------------------------------------------------------
def test_pipeline_config_de_correction_default():
    from workflow.modular.config import PipelineConfig
    cfg = PipelineConfig(project="x", output_dir=Path("/tmp"), cellranger=None)
    assert cfg.de_correction == "benjamini-hochberg"


# ---------------------------------------------------------------------------
# 6. CLI --cohort-subset parses correctly into cfg.cohort_subset
# ---------------------------------------------------------------------------
def test_cli_cohort_subset_parsed(monkeypatch):
    import workflow.modular.cli as cli_mod
    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project", "test",
            "--sample-root", "/tmp",
            "--cohort-subset", "disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma",
        ],
    )
    args = cli_mod.parse_args()
    assert args.cohort_subset == ["disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma"]


# ---------------------------------------------------------------------------
# 7. cohort_subset filter applied correctly in CellRangerModule
# ---------------------------------------------------------------------------
def test_cohort_subset_filter_applied(tmp_path):
    from workflow.modular.modules.cellranger import CellRangerModule
    from workflow.modular.config import PipelineConfig, CellRangerConfig
    from workflow.modular.context import PipelineContext

    adata = AnnData(np.ones((6, 3), dtype=float))
    adata.obs["disease"] = ["lung_adenocarcinoma", "lung_squamous_cell_carcinoma",
                            "normal", "normal", "lung_adenocarcinoma", "other"]
    adata.obs_names = [f"C{i}" for i in range(6)]
    adata.var_names = ["G1", "G2", "G3"]

    cfg = PipelineConfig(
        project="t",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        cohort_subset=["disease=lung_adenocarcinoma,lung_squamous_cell_carcinoma"],
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=tmp_path, table_dir=tmp_path)
    ctx.adata = adata

    mod = CellRangerModule()
    result = mod._apply_cohort_subset(adata, ctx)

    assert result.n_obs == 3
    assert set(result.obs["disease"].unique()) == {"lung_adenocarcinoma", "lung_squamous_cell_carcinoma"}
    assert ctx.metadata["cohort_subset_n_obs_before"] == 6
    assert ctx.metadata["cohort_subset_n_obs_after"] == 3


# ---------------------------------------------------------------------------
# 8. cohort_subset=None is a no-op
# ---------------------------------------------------------------------------
def test_cohort_subset_none_is_noop(tmp_path):
    from workflow.modular.modules.cellranger import CellRangerModule
    from workflow.modular.config import PipelineConfig, CellRangerConfig
    from workflow.modular.context import PipelineContext

    adata = AnnData(np.ones((4, 2), dtype=float))
    adata.obs_names = [f"C{i}" for i in range(4)]
    adata.var_names = ["G1", "G2"]

    cfg = PipelineConfig(
        project="t",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        cohort_subset=None,
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=tmp_path, table_dir=tmp_path)

    mod = CellRangerModule()
    result = mod._apply_cohort_subset(adata, ctx)

    assert result.n_obs == 4
    assert "cohort_subset_applied" not in ctx.metadata
