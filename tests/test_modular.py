from pathlib import Path

import numpy as np
from anndata import AnnData

from workflow.modular.config import CellRangerConfig, PipelineConfig
from workflow.modular.pipeline import run_pipeline


def test_modular_pipeline_minimal(monkeypatch, tmp_path):
    import workflow.modular.modules.cellranger as cellranger_mod
    import workflow.modular.modules.qc as qc_mod

    def fake_read_10x_mtx(path, var_names="gene_symbols", cache=False):
        adata = AnnData(np.array([[1.0, 0.0], [2.0, 3.0]], dtype=float))
        adata.var_names = ["MT-CO1", "CD3D"]
        adata.obs_names = ["C1", "C2"]
        return adata

    def fake_qc_metrics(adata, qc_vars=None, inplace=True):
        adata.obs["n_genes_by_counts"] = np.array([500, 600])
        adata.obs["total_counts"] = np.array([1000.0, 2000.0])
        adata.obs["pct_counts_mt"] = np.array([1.0, 2.0])
        adata.obs["pct_counts_ribo"] = np.array([5.0, 6.0])
        adata.obs["pct_counts_hb"] = np.array([0.0, 0.0])

    monkeypatch.setattr(cellranger_mod.sc, "read_10x_mtx", fake_read_10x_mtx)
    monkeypatch.setattr(qc_mod.sc.pp, "calculate_qc_metrics", fake_qc_metrics)
    outs = tmp_path / "dataset" / "outs" / "filtered_feature_bc_matrix"
    outs.mkdir(parents=True, exist_ok=True)

    cfg = PipelineConfig(
        project="m1",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path / "dataset", outs_dir=outs),
        optional_modules=[],
    )
    manifest = run_pipeline(cfg)
    assert manifest.exists()
    assert manifest.name == "run_manifest.json"


def test_modular_cli_parse(monkeypatch):
    import workflow.modular.cli as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp/s",
        ],
    )
    args = mod.parse_args()
    assert args.project == "x"


def test_new_modules_in_dag():
    from workflow.modular.pipeline import MODULE_DEPENDENCIES

    # Verify all 21 modules are registered (including evolution)
    assert len(MODULE_DEPENDENCIES) == 21

    # Verify new modules exist with correct dependencies
    assert "immune_phenotyping" in MODULE_DEPENDENCIES
    assert "tumor_microenvironment" in MODULE_DEPENDENCIES
    assert "gene_signature_scoring" in MODULE_DEPENDENCIES
    assert "annotation" in MODULE_DEPENDENCIES["immune_phenotyping"]
    assert "annotation" in MODULE_DEPENDENCIES["tumor_microenvironment"]
    assert "clustering" in MODULE_DEPENDENCIES["gene_signature_scoring"]


def test_registry_includes_new_modules():
    from workflow.modular.pipeline import _build_registry

    registry = _build_registry()
    assert "immune_phenotyping" in registry
    assert "tumor_microenvironment" in registry
    assert "gene_signature_scoring" in registry
    assert registry["immune_phenotyping"].name == "immune_phenotyping"
    assert registry["tumor_microenvironment"].name == "tumor_microenvironment"
    assert registry["gene_signature_scoring"].name == "gene_signature_scoring"


def test_cli_scanorama_batch_method(monkeypatch):
    import workflow.modular.cli as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp/s",
            "--batch-method",
            "scanorama",
        ],
    )
    args = mod.parse_args()
    assert args.batch_method == "scanorama"


def test_cli_signature_json(monkeypatch):
    import workflow.modular.cli as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp/s",
            "--signature-json",
            "/tmp/sigs.json",
        ],
    )
    args = mod.parse_args()
    assert args.signature_json == "/tmp/sigs.json"


def test_dependency_resolution_auto_includes():
    from workflow.modular.pipeline import _resolve_execution_order

    # Requesting immune_phenotyping should auto-include annotation and clustering
    mandatory = ["cellranger", "qc", "doublet_detection"]
    optional = ["immune_phenotyping"]
    order = _resolve_execution_order(mandatory, optional)
    assert "annotation" in order
    assert "clustering" in order
    assert "immune_phenotyping" in order
    # Verify correct ordering
    assert order.index("clustering") < order.index("annotation")
    assert order.index("annotation") < order.index("immune_phenotyping")


def test_dependency_resolution_cycle_detection(monkeypatch):
    import workflow.modular.pipeline as pipe

    deps = dict(pipe.MODULE_DEPENDENCIES)
    deps["a"] = {"b"}
    deps["b"] = {"a"}
    monkeypatch.setattr(pipe, "MODULE_DEPENDENCIES", deps)

    try:
        pipe._resolve_execution_order([], ["a", "b"])
    except ValueError as exc:
        assert "Cyclic dependency detected" in str(exc)
        assert "a" in str(exc) and "b" in str(exc)
    else:
        raise AssertionError("Expected ValueError for cyclic dependencies")


def test_parallel_figure_pool_shutdown_on_exception(monkeypatch, tmp_path):
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext

    pools = []

    class DummyPool:
        def __init__(self, max_workers=2):
            self.max_workers = max_workers
            self.shutdown_called = False
            pools.append(self)

        def shutdown(self, wait=True):
            self.shutdown_called = True

    class OkModule:
        def run(self, ctx):
            return None

    class FailModule:
        def run(self, ctx):
            raise RuntimeError("boom")

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path / "dataset", outs_dir=tmp_path / "outs"),
        optional_modules=[],
        parallel_workers=2,
    )
    ctx = PipelineContext(cfg=cfg, run_dir=tmp_path / "run", figure_dir=tmp_path / "run", table_dir=tmp_path / "run")
    ctx.run_dir.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(pipe, "_prepare_output", lambda _cfg: ctx)
    monkeypatch.setattr(
        pipe,
        "_build_registry",
        lambda: {
            "cellranger": OkModule(),
            "qc": OkModule(),
            "doublet_detection": FailModule(),
        },
    )
    monkeypatch.setattr(pipe, "ThreadPoolExecutor", DummyPool)

    try:
        pipe.run_pipeline(cfg)
    except RuntimeError as exc:
        assert "boom" in str(exc)
    else:
        raise AssertionError("Expected RuntimeError from failing mandatory module")

    assert ctx._figure_pool is None
    assert pools and pools[0].shutdown_called is True


def test_parallel_memory_guard_disables_parallel(monkeypatch, tmp_path):
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext

    adata = AnnData(np.ones((4, 4), dtype=float))
    ctx = PipelineContext(
        cfg=PipelineConfig(
            project="p",
            output_dir=tmp_path / "out",
            cellranger=CellRangerConfig(sample_root=tmp_path / "dataset", outs_dir=tmp_path / "outs"),
        ),
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
        adata=adata,
    )
    monkeypatch.setattr(pipe, "_estimate_adata_copy_bytes", lambda _adata: pipe.PARALLEL_COPY_BUDGET_BYTES)
    assert pipe._should_parallelize_appending(ctx, ["annotation", "trajectory"]) is False
