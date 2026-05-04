import json
from pathlib import Path

import numpy as np
from anndata import AnnData

from workflow.modular.config import BatchConfig, CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext
from workflow.modular.pipeline import run_pipeline


class _CudaOutOfMemoryError(Exception):
    """Simulates cupy.cuda.memory.OutOfMemoryError for testing."""


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

    # Verify all modules are registered (3 mandatory + 26 optional = 29).
    # Phase B (v2.1) added the protein_adt ingest module as the first consumer
    # of the bundle extension API; the spatial lane (spatial_ingest +
    # spatial_neighborhoods) is the second; the multimodal_integration lane
    # (EXPERIMENTAL, off by default) is the third.
    assert len(MODULE_DEPENDENCIES) == 29

    # Verify new modules exist with correct dependencies
    assert "immune_phenotyping" in MODULE_DEPENDENCIES
    assert "tumor_microenvironment" in MODULE_DEPENDENCIES
    assert "gene_signature_scoring" in MODULE_DEPENDENCIES
    assert "paper_repro" in MODULE_DEPENDENCIES
    assert "protein_adt" in MODULE_DEPENDENCIES
    assert "spatial_ingest" in MODULE_DEPENDENCIES
    assert "spatial_neighborhoods" in MODULE_DEPENDENCIES
    assert "multimodal_integration" in MODULE_DEPENDENCIES
    assert "annotation" in MODULE_DEPENDENCIES["immune_phenotyping"]
    assert "annotation" in MODULE_DEPENDENCIES["tumor_microenvironment"]
    assert "clustering" in MODULE_DEPENDENCIES["gene_signature_scoring"]
    assert "clustering" in MODULE_DEPENDENCIES["paper_repro"]
    assert "qc" in MODULE_DEPENDENCIES["protein_adt"]
    assert "qc" in MODULE_DEPENDENCIES["spatial_ingest"]
    assert "spatial_ingest" in MODULE_DEPENDENCIES["spatial_neighborhoods"]
    assert "clustering" in MODULE_DEPENDENCIES["multimodal_integration"]


def test_rna_velocity_not_in_mutating_modules():
    from workflow.modular.pipeline import _discover_mutating, _build_registry

    mutating = _discover_mutating(_build_registry())
    assert "rna_velocity" not in mutating
    assert "batch_correction" in mutating


def test_registry_includes_new_modules():
    from workflow.modular.pipeline import _build_registry

    registry = _build_registry()
    assert "immune_phenotyping" in registry
    assert "tumor_microenvironment" in registry
    assert "gene_signature_scoring" in registry
    assert "paper_repro" in registry
    assert registry["immune_phenotyping"].name == "immune_phenotyping"
    assert registry["tumor_microenvironment"].name == "tumor_microenvironment"
    assert registry["gene_signature_scoring"].name == "gene_signature_scoring"
    assert registry["paper_repro"].name == "paper_repro"


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


def test_cli_scvi_batch_method(monkeypatch):
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
            "scvi",
        ],
    )
    args = mod.parse_args()
    assert args.batch_method == "scvi"
    assert args.scvi_max_epochs == 200
    assert args.scvi_n_latent == 30
    assert args.no_scvi_early_stopping is False


def test_cli_mnn_batch_method(monkeypatch):
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
            "mnn",
        ],
    )
    args = mod.parse_args()
    assert args.batch_method == "mnn"


def test_cli_fastmnn_batch_method(monkeypatch):
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
            "fastmnn",
        ],
    )
    args = mod.parse_args()
    assert args.batch_method == "fastmnn"


def test_cli_scvi_invalid_training_args(monkeypatch):
    import pytest
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
            "scvi",
            "--scvi-max-epochs",
            "0",
            "--scvi-n-latent",
            "-1",
        ],
    )
    args = mod.parse_args()
    with pytest.raises(SystemExit, match="--scvi-max-epochs"):
        mod._validate_args(args)


def test_cli_reference_mapping_args(monkeypatch):
    import workflow.modular.cli as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp/s",
            "--reference-adata",
            "/tmp/ref.h5ad",
            "--reference-label-key",
            "label",
            "--reference-k",
            "9",
            "--reference-min-confidence",
            "0.75",
            "--reference-override-mode",
            "all",
        ],
    )
    args = mod.parse_args()
    assert args.reference_adata == "/tmp/ref.h5ad"
    assert args.reference_label_key == "label"
    assert args.reference_k == 9
    assert abs(args.reference_min_confidence - 0.75) < 1e-8
    assert args.reference_override_mode == "all"


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


# ---------------------------------------------------------------------------
# GPU OOM Recovery Tests
# ---------------------------------------------------------------------------


def test_clustering_gpu_fallback_on_cuda_oom(monkeypatch, tmp_path):
    """Clustering catches CUDA-specific OOM (not just RuntimeError) and falls back to CPU."""
    from workflow.modular.modules.clustering import ClusteringModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", True)

    mod = ClusteringModule()
    call_log = []

    def mock_run_gpu(self, adata, cfg, ctx):
        call_log.append("gpu")
        raise _CudaOutOfMemoryError(
            "CUDA out of memory. Tried to allocate 2.00 GiB"
        )

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
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_clustering_gpu_copy_creates_separate_object(monkeypatch, tmp_path):
    """GPU path receives an adata.copy(), not the original object."""
    from workflow.modular.modules.clustering import ClusteringModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", True)

    mod = ClusteringModule()
    recorded_ids = []

    def mock_run_gpu(self, adata, cfg, ctx):
        recorded_ids.append(id(adata))
        adata.obs["leiden"] = "0"

    monkeypatch.setattr(ClusteringModule, "_run_gpu", mock_run_gpu)
    monkeypatch.setattr(ClusteringModule, "_plot_umap_clusters", lambda self, a, c: None)

    from workflow.modular.context import PipelineContext

    original_adata = AnnData(np.ones((4, 4), dtype=float))
    original_id = id(original_adata)

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
        adata=original_adata,
    )
    mod.run(ctx)

    assert len(recorded_ids) == 1
    assert recorded_ids[0] != original_id
    assert ctx.metadata["clustering_backend"] == "gpu"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_de_gpu_fallback_on_cuda_oom(monkeypatch, tmp_path):
    """DE module catches CUDA OOM from rapids and falls back to CPU scanpy."""
    import sys
    import types
    from workflow.modular.modules.differential_expression import DifferentialExpressionModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", True)

    def _raise_oom(*args, **kwargs):
        raise _CudaOutOfMemoryError("CUDA OOM during rank_genes_groups")

    mock_tl = types.SimpleNamespace(rank_genes_groups=_raise_oom)
    mock_rsc = types.SimpleNamespace(tl=mock_tl)
    monkeypatch.setitem(sys.modules, "rapids_singlecell", mock_rsc)

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
    assert (tab_dir / "marker_genes_all.csv").exists()
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_gpu_fallback_on_cuda_oom(monkeypatch, tmp_path):
    """Batch correction post-processing catches CUDA OOM and falls back to CPU."""
    import sys
    import types
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", True)

    # Mock rapids_singlecell to raise OOM on any call
    def _raise_oom(*args, **kwargs):
        raise _CudaOutOfMemoryError("CUDA OOM in neighbors")

    mock_pp = types.SimpleNamespace(neighbors=_raise_oom)
    mock_tl = types.SimpleNamespace(umap=_raise_oom, leiden=_raise_oom)
    mock_rsc = types.SimpleNamespace(pp=mock_pp, tl=mock_tl)
    monkeypatch.setitem(sys.modules, "rapids_singlecell", mock_rsc)

    # Mock harmony to no-op
    monkeypatch.setattr(
        BatchCorrectionModule, "_run_harmony",
        staticmethod(lambda adata, batch_key, ctx: None),
    )

    # Mock plotting
    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    mock_fig = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (mock_fig, [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(42)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = (["A"] * 10 + ["B"] * 10)
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"].copy()
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    run_dir = tmp_path / "run"
    fig_dir = run_dir / "batch_correction"
    tab_dir = run_dir / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=fig_dir,
        table_dir=tab_dir,
        adata=adata,
    )

    mod = BatchCorrectionModule()
    mod.run(ctx)

    assert ctx.metadata["batch_correction_status"] == "completed"
    assert "leiden" in ctx.adata.obs
    assert ctx.metadata["n_clusters_after_batch"] >= 1
    monkeypatch.setattr(gutils, "_gpu_ok", None)


# ---------------------------------------------------------------------------
# CPU-Only Environment Tests
# ---------------------------------------------------------------------------


def test_gpu_utils_cache_prevents_reimport(monkeypatch):
    """Cached _gpu_ok skips cupy/rapids import probes entirely."""
    import builtins
    import workflow.modular.modules._gpu_utils as gutils

    _real_import = builtins.__import__
    import_log = []

    def _tracking_import(name, *args, **kwargs):
        if name in ("cupy", "rapids_singlecell"):
            import_log.append(name)
        return _real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _tracking_import)

    # Cache set to True — should return True without importing
    monkeypatch.setattr(gutils, "_gpu_ok", True)
    assert gutils.gpu_available() is True
    assert import_log == []

    # Cache set to False — should return False without importing
    monkeypatch.setattr(gutils, "_gpu_ok", False)
    assert gutils.gpu_available() is False
    assert import_log == []

    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_clustering_cpu_only_no_gpu_imports(monkeypatch, tmp_path):
    """In CPU-only mode, clustering never attempts rapids_singlecell import."""
    import builtins
    from workflow.modular.modules.clustering import ClusteringModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)

    _real_import = builtins.__import__

    def _guard_import(name, *args, **kwargs):
        if name == "rapids_singlecell":
            raise AssertionError("rapids_singlecell should not be imported in CPU-only mode")
        return _real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _guard_import)

    mod = ClusteringModule()
    call_log = []

    def mock_run_cpu(self, adata, cfg, ctx):
        call_log.append("cpu")
        adata.obs["leiden"] = "0"

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

    assert call_log == ["cpu"]
    assert ctx.metadata["clustering_backend"] == "cpu"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_sequential_clustering_then_de_cpu_only(monkeypatch, tmp_path):
    """Integration: clustering then DE on CPU — no cross-module state leakage."""
    import builtins
    from workflow.modular.modules.clustering import ClusteringModule
    from workflow.modular.modules.differential_expression import DifferentialExpressionModule
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)

    _real_import = builtins.__import__

    def _guard_import(name, *args, **kwargs):
        if name == "rapids_singlecell":
            raise AssertionError("rapids_singlecell should not be imported in CPU-only mode")
        return _real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _guard_import)

    # --- Clustering ---
    def mock_run_cpu(self, adata, cfg, ctx):
        adata.obs["leiden"] = (np.arange(adata.n_obs) % 3).astype(str)

    monkeypatch.setattr(ClusteringModule, "_run_cpu", mock_run_cpu)
    monkeypatch.setattr(ClusteringModule, "_plot_umap_clusters", lambda self, a, c: None)

    from workflow.modular.context import PipelineContext

    adata = AnnData(np.random.default_rng(42).random((20, 10)).astype(np.float32))
    adata.var_names = [f"G{i}" for i in range(10)]

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    run_dir = tmp_path / "run"
    fig_dir = run_dir / "clustering"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=fig_dir,
        table_dir=run_dir,
        adata=adata,
    )

    ClusteringModule().run(ctx)
    assert ctx.metadata["clustering_backend"] == "cpu"

    # --- DE ---
    de_fig_dir = run_dir / "differential_expression"
    de_tab_dir = run_dir / "differential_expression"
    de_fig_dir.mkdir(parents=True, exist_ok=True)
    ctx.figure_dir = de_fig_dir
    ctx.table_dir = de_tab_dir

    DifferentialExpressionModule().run(ctx)
    assert ctx.metadata["de_backend"] == "cpu"

    # Cache not corrupted between modules
    assert gutils._gpu_ok is False
    assert "leiden" in ctx.adata.obs
    assert "rank_genes_groups" in ctx.adata.uns
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


# ---------------------------------------------------------------------------
# Fix 1: mutates_structure discovery
# ---------------------------------------------------------------------------


def test_discover_mutating_from_registry():
    """_discover_mutating picks up mutates_structure=True from module classes."""
    import workflow.modular.pipeline as pipe

    registry = pipe._build_registry()
    result = pipe._discover_mutating(registry)
    # batch_correction declares mutates_structure = True
    assert "batch_correction" in result
    # clustering does NOT
    assert "clustering" not in result


def test_mutating_modules_fallback():
    """Static fallback set is included even if class attribute is missing."""
    import workflow.modular.pipeline as pipe

    # Even with an empty registry, the static fallback should contain batch_correction
    result = pipe._discover_mutating({})
    assert "batch_correction" in result


# ---------------------------------------------------------------------------
# Fix 2: merge-back warnings
# ---------------------------------------------------------------------------


def test_warn_dropped_changes_logs_x_mutation(caplog):
    """_warn_dropped_changes detects X modifications in a parallel branch."""
    import logging
    import workflow.modular.pipeline as pipe

    adata_main = AnnData(np.ones((4, 4), dtype=float))
    adata_branch = AnnData(np.zeros((4, 4), dtype=float))  # X modified

    with caplog.at_level(logging.WARNING):
        pipe._warn_dropped_changes("fake_module", adata_branch, adata_main)

    assert any("modified adata.X" in msg for msg in caplog.messages)


def test_warn_dropped_changes_logs_new_layer(caplog):
    """_warn_dropped_changes detects new layers added in a parallel branch."""
    import logging
    import workflow.modular.pipeline as pipe

    adata_main = AnnData(np.ones((4, 4), dtype=float))
    adata_branch = AnnData(np.ones((4, 4), dtype=float))
    adata_branch.layers["spliced"] = np.ones((4, 4))

    with caplog.at_level(logging.WARNING):
        pipe._warn_dropped_changes("fake_module", adata_branch, adata_main)

    assert any("layer 'spliced'" in msg for msg in caplog.messages)


def test_warn_dropped_changes_silent_when_clean(caplog):
    """No warnings when branch didn't modify structure."""
    import logging
    import workflow.modular.pipeline as pipe

    adata_main = AnnData(np.ones((4, 4), dtype=float))
    adata_branch = AnnData(np.ones((4, 4), dtype=float))
    # Only add obs column (this is expected for appending modules)
    adata_branch.obs["new_col"] = "a"

    with caplog.at_level(logging.WARNING):
        pipe._warn_dropped_changes("fake_module", adata_branch, adata_main)

    assert not any("dropped" in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# Fix 3: batch_post_backend metadata
# ---------------------------------------------------------------------------


def test_batch_correction_records_post_backend(monkeypatch, tmp_path):
    """batch_correction records batch_post_backend in metadata."""
    import sys
    import types
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    # Force CPU path
    monkeypatch.setattr(gutils, "_gpu_ok", False)

    # Mock harmony to no-op
    monkeypatch.setattr(
        BatchCorrectionModule, "_run_harmony",
        staticmethod(lambda adata, batch_key, ctx: None),
    )

    # Mock plotting
    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(42)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"].copy()
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert ctx.metadata["batch_post_backend"] == "cpu"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_scvi_uses_scvi_embedding(monkeypatch, tmp_path):
    """scvi batch correction should use X_scvi for neighbors/UMAP/Leiden."""
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)
    monkeypatch.setattr(
        BatchCorrectionModule,
        "_run_scvi",
        staticmethod(lambda adata, batch_key, ctx: adata.obsm.__setitem__("X_scvi", adata.obsm["X_pca"].copy())),
    )

    captured = {}
    monkeypatch.setattr(
        bc_mod.sc.pp,
        "neighbors",
        lambda adata, use_rep=None, **kwargs: captured.setdefault("use_rep", use_rep),
    )
    monkeypatch.setattr(bc_mod.sc.tl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.sc.tl, "leiden", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(7)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="scvi"),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert captured["use_rep"] == "X_scvi"
    assert ctx.metadata["batch_method"] == "scvi"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_mnn_uses_mnn_embedding(monkeypatch, tmp_path):
    """mnn batch correction should use X_mnn for neighbors/UMAP/Leiden."""
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)
    monkeypatch.setattr(
        BatchCorrectionModule,
        "_run_mnn",
        staticmethod(lambda adata, batch_key, ctx: adata.obsm.__setitem__("X_mnn", adata.obsm["X_pca"].copy())),
    )

    captured = {}
    monkeypatch.setattr(
        bc_mod.sc.pp,
        "neighbors",
        lambda adata, use_rep=None, **kwargs: captured.setdefault("use_rep", use_rep),
    )
    monkeypatch.setattr(bc_mod.sc.tl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.sc.tl, "leiden", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(21)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="mnn"),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert captured["use_rep"] == "X_mnn"
    assert ctx.metadata["batch_method"] == "mnn"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_fastmnn_uses_fastmnn_embedding(monkeypatch, tmp_path):
    """fastmnn batch correction should use X_fastmnn for neighbors/UMAP/Leiden."""
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)
    monkeypatch.setattr(
        BatchCorrectionModule,
        "_run_fastmnn",
        staticmethod(lambda adata, batch_key, ctx: adata.obsm.__setitem__("X_fastmnn", adata.obsm["X_pca"].copy())),
    )

    captured = {}
    monkeypatch.setattr(
        bc_mod.sc.pp,
        "neighbors",
        lambda adata, use_rep=None, **kwargs: captured.setdefault("use_rep", use_rep),
    )
    monkeypatch.setattr(bc_mod.sc.tl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.sc.tl, "leiden", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(22)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="fastmnn"),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert captured["use_rep"] == "X_fastmnn"
    assert ctx.metadata["batch_method"] == "fastmnn"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_scvi_failure_skips_module(monkeypatch, tmp_path):
    """When scVI backend fails, batch_correction should be skipped (not failed)."""
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)
    monkeypatch.setattr(
        BatchCorrectionModule,
        "_run_scvi",
        staticmethod(lambda adata, batch_key, ctx: (_ for _ in ()).throw(ImportError("no scvi-tools"))),
    )

    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(11)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="scvi"),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert ctx.metadata["batch_correction_status"] == "skipped_scvi_unavailable_or_failed"
    assert "no scvi-tools" in ctx.metadata["batch_correction_skip_reason"]
    assert ctx.module_status[-1]["status"] == "skipped"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_mnn_failure_skips_module(monkeypatch, tmp_path):
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)
    monkeypatch.setattr(
        BatchCorrectionModule,
        "_run_mnn",
        staticmethod(lambda adata, batch_key, ctx: (_ for _ in ()).throw(ImportError("no mnnpy"))),
    )

    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(111)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="mnn"),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert ctx.metadata["batch_correction_status"] == "skipped_mnn_unavailable_or_failed"
    assert "no mnnpy" in ctx.metadata["batch_correction_skip_reason"]
    assert ctx.module_status[-1]["status"] == "skipped"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_batch_correction_fastmnn_failure_skips_module(monkeypatch, tmp_path):
    from unittest.mock import MagicMock
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    import workflow.modular.modules.batch_correction as bc_mod
    import workflow.modular.modules._gpu_utils as gutils

    monkeypatch.setattr(gutils, "_gpu_ok", False)
    monkeypatch.setattr(
        BatchCorrectionModule,
        "_run_fastmnn",
        staticmethod(lambda adata, batch_key, ctx: (_ for _ in ()).throw(ImportError("no mnnpy"))),
    )

    monkeypatch.setattr(bc_mod.sc.pl, "umap", lambda *a, **kw: None)
    mock_ax = MagicMock()
    monkeypatch.setattr(bc_mod.plt, "subplots", lambda *a, **kw: (MagicMock(), [mock_ax, mock_ax]))
    monkeypatch.setattr(bc_mod.plt, "tight_layout", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(bc_mod.plt, "close", lambda *a, **kw: None)

    rng = np.random.default_rng(222)
    adata = AnnData(rng.random((20, 10)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 10 + ["B"] * 10
    adata.obs["leiden"] = (np.arange(20) % 3).astype(str)
    adata.obsm["X_pca"] = rng.random((20, 10)).astype(np.float32)
    adata.obsm["X_umap"] = rng.random((20, 2)).astype(np.float32)

    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="fastmnn"),
    )
    fig_dir = tmp_path / "run" / "batch_correction"
    fig_dir.mkdir(parents=True, exist_ok=True)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=fig_dir,
        table_dir=fig_dir,
        adata=adata,
    )

    BatchCorrectionModule().run(ctx)
    assert ctx.metadata["batch_correction_status"] == "skipped_fastmnn_unavailable_or_failed"
    assert "no mnnpy" in ctx.metadata["batch_correction_skip_reason"]
    assert ctx.module_status[-1]["status"] == "skipped"
    monkeypatch.setattr(gutils, "_gpu_ok", None)


def test_run_scvi_reads_config_and_records_metadata(monkeypatch, tmp_path):
    """_run_scvi should consume config values and emit training metadata."""
    import sys
    import types
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    from workflow.modular.context import PipelineContext

    captured = {}

    class _DummySCVI:
        @staticmethod
        def setup_anndata(adata, batch_key=None, layer=None):
            captured["setup_batch_key"] = batch_key
            captured["setup_layer"] = layer

        def __init__(self, adata, n_latent=10):
            captured["n_latent"] = n_latent
            self._n_obs = adata.n_obs

        def train(self, **kwargs):
            captured["train_kwargs"] = kwargs

        def get_latent_representation(self):
            return np.ones((self._n_obs, 3), dtype=np.float32)

    fake_scvi = types.SimpleNamespace(model=types.SimpleNamespace(SCVI=_DummySCVI))
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    adata = AnnData(np.ones((8, 5), dtype=np.float32))
    adata.obs["sample"] = ["A"] * 4 + ["B"] * 4
    adata.layers["counts"] = np.ones((8, 5), dtype=np.float32)

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="scvi", scvi_max_epochs=77, scvi_n_latent=19, scvi_early_stopping=False),
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
        adata=adata,
    )

    BatchCorrectionModule._run_scvi(adata, "sample", ctx)

    assert captured["setup_batch_key"] == "sample"
    assert captured["setup_layer"] == "counts"
    assert captured["n_latent"] == 19
    assert captured["train_kwargs"]["max_epochs"] == 77
    assert captured["train_kwargs"]["early_stopping"] is False
    assert ctx.metadata["scvi_input_source"] == "layers.counts"
    assert ctx.metadata["scvi_train_config"]["n_latent"] == 19
    assert "X_scvi" in adata.obsm


def test_run_scvi_rejects_non_count_like_input(monkeypatch, tmp_path):
    """_run_scvi should guard against non-count normalized inputs."""
    import sys
    import types
    import pytest
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    from workflow.modular.context import PipelineContext

    class _NoopSCVI:
        @staticmethod
        def setup_anndata(adata, batch_key=None, layer=None):
            return None

        def __init__(self, adata, n_latent=10):
            self._n_obs = adata.n_obs

        def train(self, **kwargs):
            return None

        def get_latent_representation(self):
            return np.ones((self._n_obs, 3), dtype=np.float32)

    fake_scvi = types.SimpleNamespace(model=types.SimpleNamespace(SCVI=_NoopSCVI))
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    adata = AnnData(np.random.default_rng(0).random((8, 5)).astype(np.float32))
    adata.obs["sample"] = ["A"] * 4 + ["B"] * 4
    # No counts layer, X is non-integer normalized-like values.

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        batch=BatchConfig(method="scvi"),
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
        adata=adata,
    )

    with pytest.raises(ValueError, match="raw count-like input"):
        BatchCorrectionModule._run_scvi(adata, "sample", ctx)


def test_annotation_reference_mapping_overrides_labels(monkeypatch, tmp_path):
    from workflow.modular.modules.annotation import AnnotationModule
    import workflow.modular.modules.annotation as ann_mod
    from workflow.modular.context import PipelineContext

    # Avoid heavy plotting in test.
    monkeypatch.setattr(ann_mod.sc.pl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "close", lambda *a, **kw: None)
    monkeypatch.setattr(AnnotationModule, "_plot_composition", staticmethod(lambda adata, ctx: None))

    def _fake_score_genes(ad, genes, score_name, use_raw=False):
        ad.obs[score_name] = 0.0  # Force marker labels to collapse into first class.

    monkeypatch.setattr(ann_mod.sc.tl, "score_genes", _fake_score_genes)

    genes = [f"G{i}" for i in range(60)]
    qx = np.zeros((4, 60), dtype=np.float32)
    qx[0, :30] = 2.0
    qx[1, :30] = 2.0
    qx[2, 30:] = 2.0
    qx[3, 30:] = 2.0
    adata = AnnData(qx)
    adata.var_names = genes
    adata.obs["leiden"] = ["0", "0", "1", "1"]
    adata.obsm["X_umap"] = np.random.default_rng(123).random((4, 2)).astype(np.float32)

    rx = np.zeros((6, 60), dtype=np.float32)
    rx[:3, :30] = 2.0
    rx[3:, 30:] = 2.0
    ref = AnnData(rx)
    ref.var_names = genes
    ref.obs["cell_type"] = ["TypeA", "TypeA", "TypeA", "TypeB", "TypeB", "TypeB"]
    ref_path = tmp_path / "ref.h5ad"
    ref.write_h5ad(ref_path)

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        markers={"TypeA": ["G0"], "TypeB": ["G1"]},
        reference_adata=ref_path,
        reference_label_key="cell_type",
        reference_k=3,
        reference_min_confidence=0.6,
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

    AnnotationModule().run(ctx)
    assert ctx.metadata["reference_mapping_status"] == "completed"
    assert "reference_cell_type" in adata.obs.columns
    assert "reference_confidence" in adata.obs.columns
    assert int((adata.obs["cell_type"] == "TypeB").sum()) >= 2


def test_annotation_reference_mapping_missing_file_skips(monkeypatch, tmp_path):
    from workflow.modular.modules.annotation import AnnotationModule
    import workflow.modular.modules.annotation as ann_mod
    from workflow.modular.context import PipelineContext

    monkeypatch.setattr(ann_mod.sc.pl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "close", lambda *a, **kw: None)
    monkeypatch.setattr(AnnotationModule, "_plot_composition", staticmethod(lambda adata, ctx: None))

    def _fake_score_genes(ad, genes, score_name, use_raw=False):
        ad.obs[score_name] = 0.0

    monkeypatch.setattr(ann_mod.sc.tl, "score_genes", _fake_score_genes)

    genes = [f"G{i}" for i in range(60)]
    adata = AnnData(np.ones((4, 60), dtype=np.float32))
    adata.var_names = genes
    adata.obs["leiden"] = ["0", "0", "1", "1"]
    adata.obsm["X_umap"] = np.random.default_rng(321).random((4, 2)).astype(np.float32)

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        markers={"TypeA": ["G0"], "TypeB": ["G1"]},
        reference_adata=tmp_path / "not_exists.h5ad",
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

    AnnotationModule().run(ctx)
    assert ctx.metadata["reference_mapping_status"] == "skipped_missing_reference_file"


def test_annotation_reference_mapping_conservative_preserves_high_conf_marker(monkeypatch, tmp_path):
    from workflow.modular.modules.annotation import AnnotationModule
    import workflow.modular.modules.annotation as ann_mod
    from workflow.modular.context import PipelineContext

    monkeypatch.setattr(ann_mod.sc.pl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "close", lambda *a, **kw: None)
    monkeypatch.setattr(AnnotationModule, "_plot_composition", staticmethod(lambda adata, ctx: None))

    scores = {
        "score_TypeA": np.array([5.0, 5.0, 0.2, 0.2], dtype=float),
        "score_TypeB": np.array([0.1, 0.1, 4.5, 4.5], dtype=float),
    }

    def _fake_score_genes(ad, genes, score_name, use_raw=False):
        ad.obs[score_name] = scores[score_name]

    monkeypatch.setattr(ann_mod.sc.tl, "score_genes", _fake_score_genes)

    genes = [f"G{i}" for i in range(60)]
    qx = np.zeros((4, 60), dtype=np.float32)
    qx[0, :30] = 2.0
    qx[1, :30] = 2.0
    qx[2, 30:] = 2.0
    qx[3, 30:] = 2.0
    adata = AnnData(qx)
    adata.var_names = genes
    adata.obs["leiden"] = ["0", "0", "1", "1"]
    adata.obsm["X_umap"] = np.random.default_rng(555).random((4, 2)).astype(np.float32)

    rx = np.zeros((6, 60), dtype=np.float32)
    rx[:3, 30:] = 2.0  # reference votes TypeB for first two query cells
    rx[3:, :30] = 2.0
    ref = AnnData(rx)
    ref.var_names = genes
    ref.obs["cell_type"] = ["TypeA", "TypeA", "TypeA", "TypeB", "TypeB", "TypeB"]
    ref_path = tmp_path / "ref.h5ad"
    ref.write_h5ad(ref_path)

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        markers={"TypeA": ["G0"], "TypeB": ["G1"]},
        reference_adata=ref_path,
        reference_label_key="cell_type",
        reference_k=3,
        reference_min_confidence=0.6,
        reference_override_mode="conservative",
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

    AnnotationModule().run(ctx)
    # High-confidence marker labels should be preserved in conservative mode.
    assert adata.obs["cell_type"].iloc[0] == "TypeA"
    assert adata.obs["cell_type"].iloc[1] == "TypeA"
    assert ctx.metadata["reference_mapping_override_mode"] == "conservative"


def test_annotation_reference_mapping_all_overrides_high_conf_marker(monkeypatch, tmp_path):
    from workflow.modular.modules.annotation import AnnotationModule
    import workflow.modular.modules.annotation as ann_mod
    from workflow.modular.context import PipelineContext

    monkeypatch.setattr(ann_mod.sc.pl, "umap", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "savefig", lambda *a, **kw: None)
    monkeypatch.setattr(ann_mod.plt, "close", lambda *a, **kw: None)
    monkeypatch.setattr(AnnotationModule, "_plot_composition", staticmethod(lambda adata, ctx: None))

    scores = {
        "score_TypeA": np.array([5.0, 5.0, 0.2, 0.2], dtype=float),
        "score_TypeB": np.array([0.1, 0.1, 4.5, 4.5], dtype=float),
    }

    def _fake_score_genes(ad, genes, score_name, use_raw=False):
        ad.obs[score_name] = scores[score_name]

    monkeypatch.setattr(ann_mod.sc.tl, "score_genes", _fake_score_genes)

    genes = [f"G{i}" for i in range(60)]
    qx = np.zeros((4, 60), dtype=np.float32)
    qx[0, :30] = 2.0
    qx[1, :30] = 2.0
    qx[2, 30:] = 2.0
    qx[3, 30:] = 2.0
    adata = AnnData(qx)
    adata.var_names = genes
    adata.obs["leiden"] = ["0", "0", "1", "1"]
    adata.obsm["X_umap"] = np.random.default_rng(556).random((4, 2)).astype(np.float32)

    rx = np.zeros((6, 60), dtype=np.float32)
    rx[:3, 30:] = 2.0
    rx[3:, :30] = 2.0
    ref = AnnData(rx)
    ref.var_names = genes
    ref.obs["cell_type"] = ["TypeA", "TypeA", "TypeA", "TypeB", "TypeB", "TypeB"]
    ref_path = tmp_path / "ref.h5ad"
    ref.write_h5ad(ref_path)

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        markers={"TypeA": ["G0"], "TypeB": ["G1"]},
        reference_adata=ref_path,
        reference_label_key="cell_type",
        reference_k=3,
        reference_min_confidence=0.6,
        reference_override_mode="all",
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

    AnnotationModule().run(ctx)
    # In all mode, high-confidence marker labels can be overridden.
    assert adata.obs["cell_type"].iloc[0] == "TypeB"
    assert adata.obs["cell_type"].iloc[1] == "TypeB"
    assert ctx.metadata["reference_mapping_override_mode"] == "all"


# ---------------------------------------------------------------------------
# Fix 4: requires_keys contract validation
# ---------------------------------------------------------------------------


def test_check_requires_detects_missing_obs_key():
    """_check_requires returns missing obs keys."""
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext

    class FakeMod:
        name = "fake"
        requires_keys = {"obs": ["leiden", "cell_type"]}

    adata = AnnData(np.ones((4, 4), dtype=float))
    adata.obs["leiden"] = "0"
    # cell_type is missing

    cfg = PipelineConfig(
        project="p", output_dir=Path("/tmp/out"),
        cellranger=CellRangerConfig(sample_root=Path("/tmp"), outs_dir=Path("/tmp")),
    )
    ctx = PipelineContext(
        cfg=cfg, run_dir=Path("/tmp"), figure_dir=Path("/tmp"),
        table_dir=Path("/tmp"), adata=adata,
    )
    missing = pipe._check_requires(FakeMod(), ctx)
    assert missing == ["obs.cell_type"]


def test_check_requires_detects_missing_obsm_key():
    """_check_requires returns missing obsm keys."""
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext

    class FakeMod:
        name = "fake"
        requires_keys = {"obsm": ["X_umap"]}

    adata = AnnData(np.ones((4, 4), dtype=float))

    cfg = PipelineConfig(
        project="p", output_dir=Path("/tmp/out"),
        cellranger=CellRangerConfig(sample_root=Path("/tmp"), outs_dir=Path("/tmp")),
    )
    ctx = PipelineContext(
        cfg=cfg, run_dir=Path("/tmp"), figure_dir=Path("/tmp"),
        table_dir=Path("/tmp"), adata=adata,
    )
    missing = pipe._check_requires(FakeMod(), ctx)
    assert missing == ["obsm.X_umap"]


def test_check_requires_passes_when_satisfied():
    """_check_requires returns empty list when all keys present."""
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext

    class FakeMod:
        name = "fake"
        requires_keys = {"obs": ["leiden"], "obsm": ["X_umap"]}

    adata = AnnData(np.ones((4, 4), dtype=float))
    adata.obs["leiden"] = "0"
    adata.obsm["X_umap"] = np.ones((4, 2))

    cfg = PipelineConfig(
        project="p", output_dir=Path("/tmp/out"),
        cellranger=CellRangerConfig(sample_root=Path("/tmp"), outs_dir=Path("/tmp")),
    )
    ctx = PipelineContext(
        cfg=cfg, run_dir=Path("/tmp"), figure_dir=Path("/tmp"),
        table_dir=Path("/tmp"), adata=adata,
    )
    missing = pipe._check_requires(FakeMod(), ctx)
    assert missing == []


def test_run_module_skips_optional_on_missing_keys():
    """_run_module raises _SkipModule for optional modules with missing keys."""
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext
    import pytest

    class FakeMod:
        name = "fake"
        requires_keys = {"obs": ["nonexistent"]}
        def run(self, ctx):
            pass

    adata = AnnData(np.ones((4, 4), dtype=float))
    cfg = PipelineConfig(
        project="p", output_dir=Path("/tmp/out"),
        cellranger=CellRangerConfig(sample_root=Path("/tmp"), outs_dir=Path("/tmp")),
    )
    ctx = PipelineContext(
        cfg=cfg, run_dir=Path("/tmp"), figure_dir=Path("/tmp"),
        table_dir=Path("/tmp"), adata=adata,
    )
    with pytest.raises(pipe._SkipModule):
        pipe._run_module(FakeMod(), ctx, mandatory=False)


def test_run_module_raises_for_mandatory_on_missing_keys():
    """_run_module raises ValueError for mandatory modules with missing keys."""
    import workflow.modular.pipeline as pipe
    from workflow.modular.context import PipelineContext
    import pytest

    class FakeMod:
        name = "fake"
        requires_keys = {"obs": ["nonexistent"]}
        def run(self, ctx):
            pass

    adata = AnnData(np.ones((4, 4), dtype=float))
    cfg = PipelineConfig(
        project="p", output_dir=Path("/tmp/out"),
        cellranger=CellRangerConfig(sample_root=Path("/tmp"), outs_dir=Path("/tmp")),
    )
    ctx = PipelineContext(
        cfg=cfg, run_dir=Path("/tmp"), figure_dir=Path("/tmp"),
        table_dir=Path("/tmp"), adata=adata,
    )
    with pytest.raises(ValueError, match="missing required keys"):
        pipe._run_module(FakeMod(), ctx, mandatory=True)


def test_module_declarations_consistent_with_dependencies():
    """Modules with requires_keys have matching MODULE_DEPENDENCIES."""
    import workflow.modular.pipeline as pipe

    registry = pipe._build_registry()
    # Modules that require "leiden" should depend on clustering (directly or transitively)
    leiden_requirers = [
        name for name, mod in registry.items()
        if "leiden" in getattr(mod, "requires_keys", {}).get("obs", [])
    ]
    for name in leiden_requirers:
        deps = pipe.MODULE_DEPENDENCIES.get(name, set())
        # Must depend on clustering or on a module that depends on clustering
        assert deps, f"{name} requires leiden but has no dependencies"


def test_paper_repro_skips_without_spec_and_writes_template(tmp_path):
    from workflow.modular.modules.paper_repro import PaperReproModule

    adata = AnnData(np.ones((4, 4), dtype=np.float32))
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

    PaperReproModule().run(ctx)
    assert ctx.metadata["paper_repro_status"] == "skipped_no_spec"
    assert ctx.module_status[-1]["module"] == "paper_repro"
    assert ctx.module_status[-1]["status"] == "skipped"
    assert (run_dir / "paper_repro_spec.template.json").exists()


def test_paper_repro_reads_spec_and_writes_report(tmp_path):
    import matplotlib.pyplot as plt

    from workflow.modular.config import PaperReproConfig
    from workflow.modular.modules.paper_repro import PaperReproModule

    adata = AnnData(np.ones((4, 4), dtype=np.float32))
    run_dir = tmp_path / "run"
    spec_dir = tmp_path / "spec"
    run_dir.mkdir(parents=True, exist_ok=True)
    spec_dir.mkdir(parents=True, exist_ok=True)

    reference_img = spec_dir / "fig_ref.png"
    output_img = run_dir / "clustering" / "umap_leiden.png"
    output_img.parent.mkdir(parents=True, exist_ok=True)
    img = np.tile(np.linspace(0, 1, 32, dtype=np.float32), (32, 1))
    plt.imsave(reference_img, img, cmap="gray")
    plt.imsave(output_img, img, cmap="gray")

    spec_path = spec_dir / "paper_spec.json"
    spec_path.write_text(
        json.dumps(
            {
                "papers": [
                    {
                        "paper_id": "p1",
                        "title": "Demo Paper",
                        "source": {
                            "paper_path": "paper.pdf",
                            "repo_url": "https://github.com/example/demo",
                            "repo_commit": "abc123",
                            "license": "MIT",
                        },
                        "adaptation_targets": ["clustering"],
                        "reproductions": [
                            {
                                "figure_id": "fig2a",
                                "reference_path": "fig_ref.png",
                                "pipeline_output": "clustering/umap_leiden.png",
                                "metric": "mae",
                                "min_score": 0.99,
                            }
                        ],
                    }
                ]
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        paper_repro=PaperReproConfig(spec_json=spec_path, strict=False),
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )

    PaperReproModule().run(ctx)
    assert ctx.metadata["paper_repro_status"] == "completed"
    assert (run_dir / "paper_repro_registry.csv").exists()
    assert (run_dir / "paper_repro_figures.csv").exists()
    report_path = run_dir / "paper_repro_report.json"
    assert report_path.exists()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["papers_total"] == 1
    assert report["figure_checks_passed"] == 1
