"""Unit tests for --sample-id post-load filter in the prepared-zarr branch."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from workflow.modular.config import CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.cellranger import CellRangerModule


def _make_ctx(tmp_path, sample_id: str = "lusc") -> PipelineContext:
    dataset = tmp_path / "dataset"
    cfg = PipelineConfig(
        project="test",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=dataset,
            outs_dir=dataset / "outs" / "filtered_feature_bc_matrix",
            sample_id=sample_id,
        ),
        optional_modules=[],
    )
    run_dir = tmp_path / "run"
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    return ctx


def _make_prepared_zarr(tmp_path, samples: list[str]) -> tuple[AnnData, object]:
    """Write a synthetic prepared_input.zarr with a 'sample' obs column."""
    import anndata as ad

    n = len(samples)
    adata = AnnData(np.ones((n, 2), dtype=float))
    adata.obs_names = [f"C{i}" for i in range(n)]
    adata.var_names = ["GeneA", "GeneB"]
    adata.obs["sample"] = pd.Categorical(samples)

    zarr_path = tmp_path / "dataset" / "prepared_input.zarr"
    zarr_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_zarr(zarr_path)
    return adata, zarr_path


def _run_module_with_prepared_zarr(tmp_path, sample_id: str) -> PipelineContext:
    samples = ["P15_T1", "P15_T1", "P15_T2", "P15_T3"]
    _make_prepared_zarr(tmp_path, samples)
    ctx = _make_ctx(tmp_path, sample_id=sample_id)
    mod = CellRangerModule()
    mod.run(ctx)
    return ctx


def test_default_lusc_no_filter(tmp_path):
    """Default sample_id='lusc' must not filter any cells."""
    ctx = _run_module_with_prepared_zarr(tmp_path, sample_id="lusc")
    assert ctx.adata.n_obs == 4, "All rows should be kept with default sentinel"
    assert "sample_id_filter_applied" not in ctx.metadata


def test_single_sample_id_filter(tmp_path):
    """--sample-id P15_T1 keeps only the two P15_T1 rows."""
    ctx = _run_module_with_prepared_zarr(tmp_path, sample_id="P15_T1")
    assert ctx.adata.n_obs == 2
    assert ctx.metadata["sample_id_filter_applied"] == ["P15_T1"]
    assert ctx.metadata["sample_id_filter_n_before"] == 4
    assert ctx.metadata["sample_id_filter_n_after"] == 2
    assert set(ctx.adata.obs["sample"].astype(str).unique()) == {"P15_T1"}


def test_comma_separated_sample_ids(tmp_path):
    """--sample-id P15_T1,P15_T2 keeps rows from both samples."""
    ctx = _run_module_with_prepared_zarr(tmp_path, sample_id="P15_T1,P15_T2")
    assert ctx.adata.n_obs == 3
    assert ctx.metadata["sample_id_filter_applied"] == ["P15_T1", "P15_T2"]
    assert ctx.metadata["sample_id_filter_n_before"] == 4
    assert ctx.metadata["sample_id_filter_n_after"] == 3


def test_missing_sample_column_raises(tmp_path):
    """When sample_id is set but obs lacks 'sample' column, raise ValueError."""
    import anndata as ad

    adata = AnnData(np.ones((3, 2), dtype=float))
    adata.obs_names = ["C0", "C1", "C2"]
    adata.var_names = ["GeneA", "GeneB"]
    # No 'sample' column

    zarr_path = tmp_path / "dataset" / "prepared_input.zarr"
    zarr_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_zarr(zarr_path)

    ctx = _make_ctx(tmp_path, sample_id="P15_T1")
    mod = CellRangerModule()
    with pytest.raises(ValueError, match="adata.obs has no 'sample' column"):
        mod.run(ctx)
