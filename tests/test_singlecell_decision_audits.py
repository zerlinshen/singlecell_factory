from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd

from workflow.modular.config import (
    BatchConfig,
    CellRangerConfig,
    ClusteringConfig,
    PipelineConfig,
    QCConfig,
)
from workflow.modular.context import PipelineContext


def _cfg(tmp_path, **kwargs) -> PipelineConfig:
    return PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        **kwargs,
    )


def _ctx(tmp_path, adata, **cfg_kwargs) -> PipelineContext:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    return PipelineContext(
        cfg=_cfg(tmp_path, **cfg_kwargs),
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )


def test_qc_writes_threshold_audit(tmp_path):
    from workflow.modular.modules.qc import QCModule

    x = np.array(
        [
            [5, 0, 0, 0],
            [3, 2, 0, 0],
            [0, 0, 10, 1],
            [0, 4, 0, 0],
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(x)
    adata.var_names = ["MT-CO1", "GENE1", "GENE2", "RPL1"]
    ctx = _ctx(
        tmp_path,
        adata,
        qc=QCConfig(
            min_genes=1,
            max_genes=10,
            min_counts=1,
            max_counts=100,
            max_mito_pct=90.0,
            max_ribo_pct=100.0,
            min_cells=1,
        ),
    )

    QCModule().run(ctx)

    audit = json.loads((tmp_path / "run" / "qc_threshold_audit.json").read_text())
    assert audit["threshold_mode"] == "fixed_configurable_audit_only"
    assert audit["thresholds"]["max_mito_pct"] == 90.0
    assert audit["cells_after_qc"] == ctx.adata.n_obs
    assert "pct_counts_mt" in audit["metric_quantiles"]


def test_batch_mixing_metric_improves_when_batches_are_mixed():
    from workflow.modular.modules.batch_correction import BatchCorrectionModule

    adata = ad.AnnData(np.ones((8, 3), dtype=np.float32))
    adata.obs["sample"] = ["A"] * 4 + ["B"] * 4
    adata.obsm["separated"] = np.array(
        [[0, 0], [0, 1], [1, 0], [1, 1], [10, 0], [10, 1], [11, 0], [11, 1]],
        dtype=np.float32,
    )
    adata.obsm["mixed"] = np.array(
        [[0, 0], [0, 1], [1, 0], [1, 1], [0.1, 0], [0.1, 1], [1.1, 0], [1.1, 1]],
        dtype=np.float32,
    )

    before = BatchCorrectionModule._compute_batch_mixing_metrics(
        adata,
        "sample",
        "separated",
        n_neighbors=3,
    )
    after = BatchCorrectionModule._compute_batch_mixing_metrics(
        adata,
        "sample",
        "mixed",
        n_neighbors=3,
    )

    assert before["status"] == "ok"
    assert after["status"] == "ok"
    assert after["normalized_batch_entropy"] > before["normalized_batch_entropy"]
    assert after["same_batch_neighbor_fraction"] < before["same_batch_neighbor_fraction"]


def test_batch_correction_resolution_sweep_uses_corrected_representation(tmp_path, monkeypatch):
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.clustering as clustering_mod

    def _fake_leiden(adata, resolution=None, key_added=None, **kwargs):
        key = key_added or "leiden"
        labels = ["0", "0", "1", "1"] if float(resolution) < 0.9 else ["0", "1", "2", "3"]
        adata.obs[key] = pd.Categorical(labels)

    monkeypatch.setattr(clustering_mod.sc.tl, "leiden", _fake_leiden)
    adata = ad.AnnData(np.ones((4, 3), dtype=np.float32))
    adata.obs["leiden"] = pd.Categorical(["final", "final", "other", "other"])
    adata.uns["neighbors"] = {}
    adata.obsm["X_pca_harmony"] = np.ones((4, 2), dtype=np.float32)
    ctx = _ctx(
        tmp_path,
        adata,
        clustering=ClusteringConfig(leiden_resolution_sweep=(0.5, 1.0)),
    )

    BatchCorrectionModule._write_corrected_resolution_sweep(adata, ctx, "X_pca_harmony")

    table = pd.read_csv(tmp_path / "run" / "leiden_resolution_sweep.csv")
    assert table["resolution"].tolist() == [0.5, 1.0]
    assert ctx.metadata["leiden_resolution_sweep_representation"] == "X_pca_harmony"


def test_annotation_writes_epithelial_marker_qc(tmp_path, monkeypatch):
    from workflow.modular.modules.annotation import AnnotationModule
    import workflow.modular.modules.annotation as ann_mod

    monkeypatch.setattr(ann_mod.sc.pl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "close", lambda *a, **kw: None)
    monkeypatch.setattr(
        AnnotationModule,
        "_plot_composition",
        staticmethod(lambda adata, ctx: None),
    )
    monkeypatch.setattr(
        AnnotationModule,
        "_try_reference_mapping",
        staticmethod(lambda adata, ctx: None),
    )

    genes = ["EPCAM", "KRT8", "KRT18", "CD3D", "CD3E", "GENE1", "GENE2"]
    x = np.zeros((6, len(genes)), dtype=np.float32)
    x[:3, :3] = 4.0
    x[3:, 3:5] = 4.0
    adata = ad.AnnData(x)
    adata.var_names = genes
    adata.obs["leiden"] = pd.Categorical(["0", "0", "0", "1", "1", "1"])
    adata.obsm["X_umap"] = np.random.default_rng(0).random((6, 2)).astype(np.float32)
    ctx = _ctx(
        tmp_path,
        adata,
        markers={"Tumor epithelial": ["EPCAM", "KRT8", "KRT18"], "T cell": ["CD3D", "CD3E"]},
    )

    AnnotationModule().run(ctx)

    payload = json.loads((tmp_path / "run" / "epithelial_marker_qc.json").read_text())
    assert payload["present_markers"] == ["EPCAM", "KRT8", "KRT18"]
    assert payload["epithelial_cells"] == 3
    assert (tmp_path / "run" / "epithelial_marker_summary_by_leiden.csv").exists()
    assert (tmp_path / "run" / "epithelial_marker_summary_by_cell_type.csv").exists()


def test_cnv_depends_on_annotation_and_records_gate(tmp_path):
    from workflow.modular.module_catalog import MODULE_SPECS
    from workflow.modular.modules.cnv_inference import CNVInferenceModule

    assert MODULE_SPECS["cnv_inference"].depends_on == ("annotation",)
    assert CNVInferenceModule.requires_keys == {"obs": ["cell_type"]}

    adata = ad.AnnData(np.ones((4, 5), dtype=np.float32))
    adata.obs["cell_type"] = pd.Categorical(["Fibroblast", "Fibroblast", "Tumor", "Tumor"])
    ctx = _ctx(tmp_path, adata)

    CNVInferenceModule._write_annotation_gate(adata, ctx, "Fibroblast")

    gate = json.loads((tmp_path / "run" / "cnv_annotation_qc.json").read_text())
    assert gate["status"] == "present"
    assert gate["reference_group_cells"] == 2
    assert gate["reference_group_status"] == "insufficient_cells_global_mean"


def test_leiden_resolution_sweep_preserves_final_label(tmp_path, monkeypatch):
    from workflow.modular.modules.clustering import ClusteringModule
    import workflow.modular.modules.clustering as clustering_mod

    calls = []

    def _fake_leiden(adata, resolution=None, key_added=None, **kwargs):
        calls.append(resolution)
        key = key_added or "leiden"
        labels = ["0", "0", "1", "1"] if float(resolution) < 0.9 else ["0", "1", "2", "3"]
        adata.obs[key] = pd.Categorical(labels)

    monkeypatch.setattr(clustering_mod.sc.tl, "leiden", _fake_leiden)
    adata = ad.AnnData(np.ones((4, 3), dtype=np.float32))
    adata.obs["leiden"] = pd.Categorical(["final", "final", "other", "other"])
    adata.uns["neighbors"] = {}
    ctx = _ctx(
        tmp_path,
        adata,
        clustering=ClusteringConfig(leiden_resolution_sweep=(0.5, 1.0)),
    )

    ClusteringModule._write_resolution_sweep_audit(adata, ctx, ctx.cfg.clustering)

    assert adata.obs["leiden"].astype(str).tolist() == ["final", "final", "other", "other"]
    assert calls == [0.5, 1.0]
    table = pd.read_csv(tmp_path / "run" / "leiden_resolution_sweep.csv")
    assert table["resolution"].tolist() == [0.5, 1.0]
