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

    def fake_qc_metrics(adata, qc_vars=None, inplace=True, **kwargs):
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

    # Verify all modules are registered (3 mandatory + 21 optional = 24)
    assert len(MODULE_DEPENDENCIES) == 24

    # Verify new modules exist with correct dependencies
    assert "immune_phenotyping" in MODULE_DEPENDENCIES
    assert "tumor_microenvironment" in MODULE_DEPENDENCIES
    assert "gene_signature_scoring" in MODULE_DEPENDENCIES
    assert "annotation" in MODULE_DEPENDENCIES["immune_phenotyping"]
    assert "annotation" in MODULE_DEPENDENCIES["tumor_microenvironment"]
    assert "clustering" in MODULE_DEPENDENCIES["gene_signature_scoring"]


def test_rna_velocity_not_in_mutating_modules():
    from workflow.modular.pipeline import MUTATING_MODULES

    assert "rna_velocity" not in MUTATING_MODULES
    assert "batch_correction" in MUTATING_MODULES


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


def test_cli_de_options(monkeypatch):
    import workflow.modular.cli as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp/s",
            "--de-method",
            "t-test",
            "--de-n-genes",
            "123",
        ],
    )
    args = mod.parse_args()
    assert args.de_method == "t-test"
    assert args.de_n_genes == 123


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


def test_gpu_utils_returns_bool(monkeypatch):
    """gpu_available() returns True when cupy+rapids are importable, False otherwise."""
    import workflow.modular.modules._gpu_utils as gutils

    # Reset cache
    monkeypatch.setattr(gutils, "_gpu_ok", None)
    # When imports fail, should return False
    import builtins
    _real_import = builtins.__import__

    def _block_cupy(name, *args, **kwargs):
        if name == "cupy":
            raise ImportError("no cupy")
        return _real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _block_cupy)
    assert gutils.gpu_available() is False
    # Cached
    assert gutils._gpu_ok is False


def test_clustering_gpu_fallback(monkeypatch, tmp_path):
    """If GPU clustering raises, module falls back to CPU and records metadata."""
    from workflow.modular.modules.clustering import ClusteringModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", True)

    mod = ClusteringModule()
    call_log = []

    def mock_run_gpu(self, adata, cfg, ctx):
        call_log.append("gpu")
        raise RuntimeError("GPU OOM")

    def mock_run_cpu(self, adata, cfg, ctx):
        call_log.append("cpu")
        adata.obs["leiden"] = "0"

    monkeypatch.setattr(ClusteringModule, "_run_gpu", mock_run_gpu)
    monkeypatch.setattr(ClusteringModule, "_run_cpu", mock_run_cpu)
    monkeypatch.setattr(ClusteringModule, "_plot_umap_clusters", lambda self, a, c: None)

    from workflow.modular.context import PipelineContext

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
        adata=AnnData(np.ones((4, 4), dtype=float)),
    )
    mod.run(ctx)

    assert call_log == ["gpu", "cpu"]
    assert ctx.metadata["clustering_backend"] == "cpu"

    # Reset
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_clustering_gpu_fallback_does_not_reuse_mutated_gpu_adata(monkeypatch, tmp_path):
    """GPU failure should not leak partially mutated state into CPU fallback."""
    from workflow.modular.modules.clustering import ClusteringModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", True)
    mod = ClusteringModule()

    def mock_run_gpu(self, adata, cfg, ctx):
        adata.obs["gpu_only"] = "1"
        adata.X = np.zeros_like(adata.X)
        raise RuntimeError("GPU failure after mutation")

    def mock_run_cpu(self, adata, cfg, ctx):
        # CPU fallback must run on the original (unmutated) object.
        assert "gpu_only" not in adata.obs
        assert float(np.asarray(adata.X).sum()) > 0.0
        adata.obs["leiden"] = "0"

    monkeypatch.setattr(ClusteringModule, "_run_gpu", mock_run_gpu)
    monkeypatch.setattr(ClusteringModule, "_run_cpu", mock_run_cpu)
    monkeypatch.setattr(ClusteringModule, "_plot_umap_clusters", lambda self, a, c: None)

    from workflow.modular.context import PipelineContext

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
        adata=AnnData(np.ones((4, 4), dtype=float)),
    )
    mod.run(ctx)

    assert ctx.metadata["clustering_backend"] == "cpu"
    assert "gpu_only" not in ctx.adata.obs
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_de_backend_metadata(monkeypatch, tmp_path):
    """DE module records de_backend in metadata."""
    from workflow.modular.modules.differential_expression import DifferentialExpressionModule
    import workflow.modular.modules._gpu_utils as gutils

    # Force GPU off
    monkeypatch.setattr(gutils, "_gpu_ok", False)

    adata = AnnData(np.random.default_rng(42).random((20, 10)).astype(np.float32))
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.var_names = [f"G{i}" for i in range(10)]

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    run_dir = tmp_path / "run"
    fig_dir = run_dir / "differential_expression"
    tab_dir = run_dir / "differential_expression"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=fig_dir,
        table_dir=tab_dir,
        adata=adata,
    )

    mod = DifferentialExpressionModule()
    mod.run(ctx)
    assert ctx.metadata["de_backend"] == "cpu"
    assert ctx.metadata["de_significant_genes"] >= 0

    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_parallel_memory_guard_adjusts_worker_count(monkeypatch, tmp_path):
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
    monkeypatch.setattr(pipe, "_estimate_adata_copy_bytes", lambda _adata: 100)
    monkeypatch.setattr(pipe, "_get_available_memory_bytes", lambda: 350)
    monkeypatch.setattr(pipe, "MEMORY_RESERVE_BYTES", 100)
    workers = pipe._safe_parallel_worker_count(
        ctx,
        ["annotation", "trajectory", "immune_phenotyping"],
        requested_workers=8,
    )
    assert workers == 2
