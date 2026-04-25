from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from workflow.modular.config import CellRangerConfig, PipelineConfig, PseudobulkConfig
from workflow.modular.context import PipelineContext


def test_prepare_output_resume_uses_latest_checkpoint_run(tmp_path):
    import workflow.modular.pipeline as pipe

    out = tmp_path / "results"
    old_run = out / "proj_20240101_000000"
    new_run = out / "proj_20240102_000000"
    (old_run / ".checkpoints").mkdir(parents=True)
    (new_run / ".checkpoints").mkdir(parents=True)

    # Ensure deterministic ordering by mtime
    old_ts = 1_700_000_000
    new_ts = 1_800_000_000
    import os

    os.utime(old_run, (old_ts, old_ts))
    os.utime(new_run, (new_ts, new_ts))

    cfg = PipelineConfig(
        project="proj",
        output_dir=out,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        checkpoint=True,
        resume_from="qc",
    )
    ctx = pipe._prepare_output(cfg)
    assert ctx.run_dir == new_run


def test_checkpoint_roundtrip_restores_module_dirs(monkeypatch, tmp_path):
    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        checkpoint=True,
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    de_dir = run_dir / "differential_expression"
    de_dir.mkdir(parents=True, exist_ok=True)

    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=AnnData(np.ones((2, 2), dtype=np.float32)),
    )
    monkeypatch.setattr(PipelineContext, "_use_zarr_checkpoints", lambda self: False)
    ctx._module_dirs["differential_expression"] = de_dir
    ctx.save_checkpoint("qc")

    ctx2 = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
    )
    assert ctx2.load_checkpoint("qc")
    assert ctx2.module_output_dir("differential_expression") == de_dir


def test_pseudobulk_aggregate_materializes_lazy_counts_to_csr():
    from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule

    class FakeLazyMatrix:
        def __init__(self, arr):
            self._arr = sparse.csr_matrix(arr)

        def compute(self):
            return self._arr

    adata = AnnData(np.ones((2, 2), dtype=np.float32))
    adata.layers["counts"] = FakeLazyMatrix([[1, 0], [0, 2]])
    adata.obs["sample"] = ["S1", "S1"]
    adata.obs["cell_type"] = ["Tumor", "Tumor"]

    counts, meta = PseudobulkDEModule._aggregate(adata, "sample", "cell_type", min_cells_per_sample=1)
    assert sparse.isspmatrix_csr(adata.layers["counts"].compute())
    np.testing.assert_array_equal(counts.to_numpy(), np.array([[1, 2]]))


def test_pseudobulk_aggregate_uses_counts_layer():
    from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule

    counts = np.array(
        [
            [10, 1],
            [20, 2],
            [30, 3],
            [1, 10],
            [2, 20],
            [3, 30],
        ],
        dtype=np.int32,
    )
    adata = AnnData(np.log1p(counts.astype(np.float32)))
    adata.layers["counts"] = sparse.csr_matrix(counts)
    adata.var_names = ["g1", "g2"]
    adata.obs["sample"] = ["S1"] * 3 + ["S2"] * 3
    adata.obs["cell_type"] = ["Tumor"] * 6

    pb_counts, pb_meta = PseudobulkDEModule._aggregate(adata, "sample", "cell_type")
    idx_s1 = pb_meta.index[pb_meta["sample"] == "S1"][0]
    idx_s2 = pb_meta.index[pb_meta["sample"] == "S2"][0]
    np.testing.assert_array_equal(pb_counts.loc[idx_s1].to_numpy(), np.array([60, 6]))
    np.testing.assert_array_equal(pb_counts.loc[idx_s2].to_numpy(), np.array([6, 60]))


def test_pseudobulk_skips_without_counts_layer(tmp_path):
    from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule

    adata = AnnData(np.ones((8, 4), dtype=np.float32))
    adata.obs["sample"] = ["A"] * 4 + ["B"] * 4
    adata.obs["cell_type"] = ["Tumor"] * 8

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )

    PseudobulkDEModule().run(ctx)
    assert ctx.metadata["pseudobulk_de_status"] == "skipped_missing_counts_layer"
    assert ctx.module_status[-1]["module"] == "pseudobulk_de"
    assert ctx.module_status[-1]["status"] == "skipped"


def test_pseudobulk_skips_without_explicit_contract(tmp_path):
    from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule

    counts = np.ones((8, 4), dtype=np.int32)
    adata = AnnData(np.log1p(counts.astype(np.float32)))
    adata.layers["counts"] = sparse.csr_matrix(counts)
    adata.obs["sample"] = ["A"] * 4 + ["B"] * 4
    adata.obs["cell_type"] = ["Tumor"] * 4 + ["Immune"] * 4

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )

    PseudobulkDEModule().run(ctx)
    assert ctx.metadata["pseudobulk_de_status"] == "skipped_missing_contrast_contract"
    assert ctx.metadata["pseudobulk_de_mode"] == "skipped"
    assert ctx.module_status[-1]["status"] == "skipped"


def test_pseudobulk_confirmatory_mode_records_contract(tmp_path, monkeypatch):
    from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule

    counts = np.tile(np.array([[10, 1], [8, 2], [9, 1]], dtype=np.int32), (4, 1))
    adata = AnnData(np.log1p(counts.astype(np.float32)))
    adata.layers["counts"] = sparse.csr_matrix(counts)
    adata.var_names = ["g1", "g2"]
    adata.obs["sample"] = ["S1"] * 3 + ["S2"] * 3 + ["S3"] * 3 + ["S4"] * 3
    adata.obs["cell_type"] = ["Tumor"] * 12
    adata.obs["condition"] = ["A"] * 6 + ["B"] * 6

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        pseudobulk=PseudobulkConfig(
            contrast_col="condition",
            contrast_a="A",
            contrast_b="B",
        ),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )

    def fake_run_de(self, counts, meta, cond_col, ca, cb, min_samples_per_condition=2):
        return pd.DataFrame(
            [{"gene": "g1", "log2FC": 1.5, "pvalue": 0.01, "padj": 0.02}]
        )

    monkeypatch.setattr(PseudobulkDEModule, "_run_de", fake_run_de)
    monkeypatch.setattr(PseudobulkDEModule, "_plot_volcano", lambda *a, **k: None)
    monkeypatch.setattr(PseudobulkDEModule, "_plot_heatmap", lambda *a, **k: None)

    PseudobulkDEModule().run(ctx)
    out = pd.read_csv(run_dir / "pseudobulk_de_results.csv")
    assert list(out["group"]) == ["Tumor|A_vs_B"]
    assert ctx.metadata["pseudobulk_de_mode"] == "confirmatory"
    assert ctx.metadata["pseudobulk_de_status"] == "completed"
    assert ctx.metadata["pseudobulk_de_contrast_contract"]["contrast_col"] == "condition"


def test_rna_velocity_skips_without_splicing_inputs(tmp_path):
    from workflow.modular.modules.rna_velocity import RNAVelocityModule

    adata = AnnData(np.ones((6, 4), dtype=np.float32))
    adata.obsm["X_umap"] = np.zeros((6, 2), dtype=np.float32)

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        optional_modules=["rna_velocity"],
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )

    RNAVelocityModule().run(ctx)
    assert ctx.metadata["rna_velocity_status"] == "skipped_missing_splicing_modality"
    assert ctx.module_status[-1]["module"] == "rna_velocity"
    assert ctx.module_status[-1]["status"] == "skipped"


def test_annotation_handles_single_available_marker_set(monkeypatch, tmp_path):
    import sys
    import types

    fake_scanpy = types.SimpleNamespace(
        tl=types.SimpleNamespace(score_genes=lambda *a, **k: None),
        pl=types.SimpleNamespace(umap=lambda *a, **k: None),
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    import workflow.modular.modules.annotation as ann_mod

    adata = AnnData(np.array([[5.0], [4.0], [6.0]], dtype=np.float32))
    adata.var_names = ["EPCAM"]
    adata.obs["leiden"] = ["0", "0", "1"]
    adata.obsm["X_umap"] = np.zeros((3, 2), dtype=np.float32)

    def fake_score_genes(a, genes, score_name, use_raw=False):
        a.obs[score_name] = np.array([0.5, 0.4, 0.6], dtype=np.float32)

    monkeypatch.setattr(ann_mod.sc.tl, "score_genes", fake_score_genes)
    monkeypatch.setattr(ann_mod.sc.pl, "umap", lambda *a, **k: None)
    monkeypatch.setattr(ann_mod.plt, "savefig", lambda *a, **k: None)
    monkeypatch.setattr(ann_mod.plt, "close", lambda *a, **k: None)
    monkeypatch.setattr(
        ann_mod.AnnotationModule,
        "_plot_composition",
        staticmethod(lambda _adata, _ctx: None),
    )

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )
    ann_mod.AnnotationModule().run(ctx)

    np.testing.assert_allclose(
        adata.obs["annotation_confidence"].to_numpy(),
        adata.obs["score_Tumor epithelial"].to_numpy(),
    )


def test_batch_correction_records_skipped_status_when_batch_key_missing(monkeypatch, tmp_path):
    import sys
    import types

    fake_scanpy = types.SimpleNamespace(
        pp=types.SimpleNamespace(),
        pl=types.SimpleNamespace(umap=lambda *a, **k: None),
        tl=types.SimpleNamespace(),
        external=types.SimpleNamespace(pp=types.SimpleNamespace()),
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)

    from workflow.modular.modules.batch_correction import BatchCorrectionModule

    adata = AnnData(np.ones((6, 4), dtype=np.float32))
    adata.obs["leiden"] = ["0", "0", "1", "1", "2", "2"]

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert ctx.metadata["batch_correction_status"] == "skipped_missing_batch_key"
    assert ctx.module_status[-1]["module"] == "batch_correction"
    assert ctx.module_status[-1]["status"] == "skipped"
