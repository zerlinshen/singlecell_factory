"""Tests that load_checkpoint raises CheckpointSchemaMismatch on stale v1 sidecars."""
from __future__ import annotations

import json

import numpy as np
import pytest
from anndata import AnnData

from workflow.modular.config import BatchConfig, CellRangerConfig, PipelineConfig
from workflow.modular.context import (
    SC_CHECKPOINT_SCHEMA_VERSION,
    CheckpointSchemaMismatch,
    PipelineContext,
)


def _make_ctx(tmp_path):
    cfg = PipelineConfig(
        project="test_version",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "samples",
            outs_dir=tmp_path / "samples" / "outs",
        ),
        optional_modules=[],
        checkpoint=True,
    )
    run_dir = tmp_path / "out" / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    return PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)


def _write_minimal_adata_h5ad(cp_dir: object, module_name: str) -> None:
    """Write a tiny h5ad so the checkpoint file path resolves."""
    from pathlib import Path

    cp_dir = Path(cp_dir)
    cp_dir.mkdir(parents=True, exist_ok=True)
    adata = AnnData(np.array([[1.0, 2.0]], dtype=np.float32))
    adata.obs_names = ["cell_0"]
    adata.var_names = ["gene_0", "gene_1"]
    adata.write(cp_dir / f"after_{module_name}.h5ad")


def test_load_checkpoint_raises_on_v1_sidecar(tmp_path):
    """load_checkpoint must raise CheckpointSchemaMismatch when sidecar schema_version < SC_CHECKPOINT_SCHEMA_VERSION."""
    ctx = _make_ctx(tmp_path)
    cp_dir = ctx._checkpoint_dir
    module_name = "qc"

    _write_minimal_adata_h5ad(cp_dir, module_name)

    # Write a v1 sidecar (no schema_version key -> defaults to 1)
    sidecar_v1 = {
        "module": module_name,
        "adata_checkpoint_policy": "",
        "metadata": {},
        "module_status": [],
        "module_dirs": {},
    }
    (cp_dir / f"after_{module_name}.json").write_text(
        json.dumps(sidecar_v1, indent=2), encoding="utf-8"
    )

    with pytest.raises(CheckpointSchemaMismatch) as exc_info:
        ctx.load_checkpoint(module_name)

    msg = str(exc_info.value)
    assert f"v1" in msg
    assert f"v{SC_CHECKPOINT_SCHEMA_VERSION}" in msg
    assert "regenerate_prepared_zarr.py" in msg


def test_load_checkpoint_raises_on_explicit_old_version(tmp_path):
    """load_checkpoint must raise when schema_version is explicitly set below current."""
    ctx = _make_ctx(tmp_path)
    cp_dir = ctx._checkpoint_dir
    module_name = "doublet_detection"

    _write_minimal_adata_h5ad(cp_dir, module_name)

    sidecar_old = {
        "schema_version": SC_CHECKPOINT_SCHEMA_VERSION - 1,
        "module": module_name,
        "adata_checkpoint_policy": "",
        "metadata": {},
        "module_status": [],
        "module_dirs": {},
    }
    (cp_dir / f"after_{module_name}.json").write_text(
        json.dumps(sidecar_old, indent=2), encoding="utf-8"
    )

    with pytest.raises(CheckpointSchemaMismatch):
        ctx.load_checkpoint(module_name)


def test_load_checkpoint_succeeds_on_current_version(tmp_path):
    """load_checkpoint must succeed when schema_version matches current."""
    ctx = _make_ctx(tmp_path)
    cp_dir = ctx._checkpoint_dir
    module_name = "cellranger"

    _write_minimal_adata_h5ad(cp_dir, module_name)

    sidecar_current = {
        "schema_version": SC_CHECKPOINT_SCHEMA_VERSION,
        "module": module_name,
        "adata_checkpoint_policy": "",
        "metadata": {"foo": "bar"},
        "module_status": [],
        "module_dirs": {},
    }
    (cp_dir / f"after_{module_name}.json").write_text(
        json.dumps(sidecar_current, indent=2), encoding="utf-8"
    )

    result = ctx.load_checkpoint(module_name)
    assert result is True
    assert ctx.metadata == {"foo": "bar"}
