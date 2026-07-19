"""WS1 #1 leiden-race regression (2026-05-24 governance remediation).

Covers the annotation/batch_correction leiden race fix:
  ordering-only ``runs_after`` (never auto-includes, never feeds tiering) +
  non-destructive ``leiden_corrected`` + annotation leiden resolver, so
  ``cell_type`` is a deterministic function of the FINAL leiden in both
  sequential and parallel execution modes.

These tests depend only on the self-contained #1 surface
(module_catalog.py, pipeline.py, annotation.py). The #4 ComBat and #5 doublet
regressions live in test_ws1_science_fixes.py.

Plan: .omc/plans/2026-05-24-pipeline-governance-remediation.md (WS1 #1).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import pytest

import workflow.modular.pipeline as pipe
from workflow.modular.pipeline import (
    _resolve_execution_order,
    _compute_tiers,
    MODULE_RUNS_AFTER,
)
from workflow.modular.module_catalog import MODULE_SPECS
from workflow.modular.modules.annotation import AnnotationModule


# ---------------------------------------------------------------------------
# #1 runs_after ordering-only semantics
# ---------------------------------------------------------------------------


def test_runs_after_catalog_hint_present():
    """annotation declares runs_after=batch_correction (ordering-only)."""
    assert MODULE_SPECS["annotation"].runs_after == ("batch_correction",)
    assert "annotation" in MODULE_RUNS_AFTER
    assert "batch_correction" in MODULE_RUNS_AFTER["annotation"]
    # depends_on must NOT include batch_correction (that would force-pull it).
    assert "batch_correction" not in MODULE_SPECS["annotation"].depends_on


def test_runs_after_orders_annotation_after_batch_correction():
    """When BOTH requested, annotation is sequenced after batch_correction."""
    order = _resolve_execution_order(
        ["cellranger", "qc", "ambient_correction", "doublet_detection"],
        ["clustering", "batch_correction", "annotation"],
    )
    assert order.index("batch_correction") < order.index("annotation"), order


def test_runs_after_does_not_auto_include_batch_correction():
    """annotation WITHOUT batch_correction must NOT pull batch_correction in.

    Locks the resolved iter-1 CRITICAL: runs_after never feeds inclusion
    auto-pull (that is depends_on's job alone).
    """
    order = _resolve_execution_order(
        ["cellranger", "qc", "ambient_correction", "doublet_detection"],
        ["clustering", "annotation"],
    )
    assert "annotation" in order
    assert "batch_correction" not in order, order


def test_runs_after_never_seeds_compute_tiers_indegree():
    """_compute_tiers keys off depends_on only, never runs_after.

    If runs_after leaked into tiering, annotation and batch_correction (both
    depend only on clustering) would be split into different tiers. They must
    share a tier so the mutating-first ordering inside _execute_tier handles
    sequencing in parallel mode.
    """
    order = _resolve_execution_order(
        ["cellranger", "qc", "ambient_correction", "doublet_detection"],
        ["clustering", "batch_correction", "annotation"],
    )
    completed = {
        "cellranger", "qc", "ambient_correction", "doublet_detection", "clustering",
    }
    tiers = _compute_tiers(order, completed)
    tier_of = {m: i for i, t in enumerate(tiers) for m in t}
    assert tier_of["annotation"] == tier_of["batch_correction"], tiers


def test_runs_after_combined_stall_drops_hint_and_breadcrumbs(monkeypatch):
    """A combined-graph stall DROPS runs_after edges (never raises) + breadcrumb.

    Force a stall by making batch_correction runs_after annotation while
    annotation runs_after batch_correction. The depends_on graph is acyclic, so
    the base order is still emitted; only the ordering hints are dropped.
    """
    forced = dict(MODULE_RUNS_AFTER)
    forced["annotation"] = {"batch_correction"}
    forced["batch_correction"] = {"annotation"}
    monkeypatch.setattr(pipe, "MODULE_RUNS_AFTER", forced)

    sink: list[dict] = []
    order = _resolve_execution_order(
        ["cellranger", "qc", "ambient_correction", "doublet_detection"],
        ["clustering", "batch_correction", "annotation"],
        sink,
    )
    # Did NOT raise; base depends_on order is intact.
    assert "annotation" in order and "batch_correction" in order
    # Breadcrumb recorded.
    assert sink, "expected a dropped-hint breadcrumb on combined-graph stall"
    assert sink[0]["reason"] == "combined_graph_stall"
    assert "annotation" in sink[0]["dropped_runs_after"]


def test_depends_on_cycle_still_raises(monkeypatch):
    """The hard cycle raise stays bound to the depends_on-only graph."""
    deps = dict(pipe.MODULE_DEPENDENCIES)
    deps["a"] = {"b"}
    deps["b"] = {"a"}
    monkeypatch.setattr(pipe, "MODULE_DEPENDENCIES", deps)
    with pytest.raises(ValueError, match="Cyclic dependency"):
        _resolve_execution_order([], ["a", "b"])


# ---------------------------------------------------------------------------
# #1 annotation leiden resolver
# ---------------------------------------------------------------------------


def _annotation_adata(corrected: bool) -> ad.AnnData:
    """5-cluster fixture; optionally carry a divergent leiden_corrected."""
    n_per = 60
    n_clusters = 5
    n_cells = n_per * n_clusters
    marker_genes = [
        "MS4A1", "CD79A", "CD79B", "CD19",
        "CD3D", "CD3E", "TRAC", "CD4",
        "TPSAB1", "TPSB2", "CPA3", "KIT",
        "LYZ", "CD68", "CST3", "FCGR3A",
        "GNLY", "NKG7", "KLRD1", "NCAM1",
    ]
    background = [f"GENE{i:04d}" for i in range(160)]
    genes = list(dict.fromkeys(marker_genes + background))
    gene_idx = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(7)
    X = rng.uniform(0.0, 0.1, size=(n_cells, len(genes))).astype(np.float32)
    cluster_markers = [
        ["MS4A1", "CD79A", "CD79B", "CD19"],
        ["CD3D", "CD3E", "TRAC", "CD4"],
        ["TPSAB1", "TPSB2", "CPA3", "KIT"],
        ["LYZ", "CD68", "CST3", "FCGR3A"],
        ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    ]
    leiden = []
    for c in range(n_clusters):
        s, e = c * n_per, (c + 1) * n_per
        for g in cluster_markers[c]:
            X[s:e, gene_idx[g]] = rng.uniform(3.0, 5.0, size=n_per).astype(np.float32)
        leiden.extend([str(c)] * n_per)
    obs = pd.DataFrame(
        {"leiden": pd.Categorical(leiden)},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    if corrected:
        # leiden_corrected scrambles the cluster -> cell mapping so a resolver
        # that wrongly reads `leiden` would bind cell_type differently.
        shifted = [str((int(c) + 1) % n_clusters) for c in leiden]
        adata.obs["leiden_corrected"] = pd.Categorical(shifted)
    return adata


class _Cfg:
    markers: dict = {}
    annotation_confidence_threshold: float = -999.0
    annotation_strategy: str = "cluster_voting"
    reference_adata = None


class _Ctx:
    def __init__(self, adata, tmp_path):
        self.adata = adata
        self.cfg = _Cfg()
        self.metadata: dict = {}
        self.table_dir = tmp_path
        self.figure_dir = tmp_path


def _run_annotation(adata, tmp_path):
    import unittest.mock as mock
    import scanpy as sc_real
    import workflow.modular.modules.annotation as ann_mod

    ctx = _Ctx(adata, tmp_path)
    _orig_viz = AnnotationModule.__dict__["_plot_composition"]
    _orig_ref = AnnotationModule.__dict__["_try_reference_mapping"]
    _orig_epi = AnnotationModule.__dict__["_write_epithelial_marker_qc"]
    try:
        AnnotationModule._plot_composition = staticmethod(lambda adata, ctx: None)
        AnnotationModule._try_reference_mapping = staticmethod(lambda adata, ctx: None)
        AnnotationModule._write_epithelial_marker_qc = classmethod(
            lambda cls, adata, ctx: None
        )
        with mock.patch.object(sc_real.pl, "umap", return_value=None), \
             mock.patch.object(ann_mod.plt, "savefig", return_value=None), \
             mock.patch.object(ann_mod.plt, "close", return_value=None):
            AnnotationModule().run(ctx)
    finally:
        AnnotationModule._plot_composition = _orig_viz
        AnnotationModule._try_reference_mapping = _orig_ref
        AnnotationModule._write_epithelial_marker_qc = _orig_epi
    return ctx


def test_annotation_resolver_prefers_leiden_corrected():
    adata = _annotation_adata(corrected=True)
    assert AnnotationModule._resolve_leiden_key(adata) == "leiden_corrected"


def test_annotation_resolver_falls_back_to_leiden():
    adata = _annotation_adata(corrected=False)
    assert AnnotationModule._resolve_leiden_key(adata) == "leiden"


def test_annotation_binds_celltype_to_leiden_corrected(tmp_path):
    """cell_type is a deterministic function of leiden_corrected when present."""
    adata = _annotation_adata(corrected=True)
    ctx = _run_annotation(adata, tmp_path)
    assert ctx.metadata["annotation_leiden_source"] == "leiden_corrected"
    # Each leiden_corrected cluster maps to exactly one cell_type.
    grouped = adata.obs.groupby("leiden_corrected", observed=True)["cell_type"].nunique()
    assert (grouped == 1).all(), grouped.to_dict()


def test_annotation_recomputes_celltype_on_resume(tmp_path):
    """Resume-lock: annotation always rebinds cell_type from the resolved key.

    A stale cell_type column carried in on resume must not short-circuit the
    recompute. Plant a wrong cell_type, run annotation, and assert it is
    overwritten by a clean leiden_corrected-derived binding.
    """
    adata = _annotation_adata(corrected=True)
    adata.obs["cell_type"] = "STALE"  # simulate resumed/carried-in state
    ctx = _run_annotation(adata, tmp_path)
    assert (adata.obs["cell_type"] != "STALE").all()
    grouped = adata.obs.groupby("leiden_corrected", observed=True)["cell_type"].nunique()
    assert (grouped == 1).all(), grouped.to_dict()


def test_annotation_celltype_leiden_mapping_consistent(tmp_path):
    """cell_type is a pure function of the resolved cluster column.

    This is the mode-independence invariant the leiden race fix guarantees:
    binding the same clusters yields the same cell_type mapping.
    """
    adata = _annotation_adata(corrected=False)
    ctx = _run_annotation(adata, tmp_path)
    mapping = (
        adata.obs.groupby("leiden", observed=True)["cell_type"]
        .agg(lambda s: s.unique().tolist())
    )
    assert all(len(v) == 1 for v in mapping), mapping.to_dict()


# ---------------------------------------------------------------------------
# #1 end-to-end: cell_type<->leiden mapping identical across worker counts
# ---------------------------------------------------------------------------


def _ws1_integration_registry():
    """Stub registry mirroring the real ordering contract for the leiden race.

    clustering -> leiden; batch_correction (mutating) overwrites leiden with a
    *corrected* labelling and writes leiden_corrected; annotation binds
    cell_type from the resolved leiden key. The pre- vs post-correction leiden
    deliberately diverge so a mode-dependent binding (annotation seeing
    pre-correction leiden) would produce a DIFFERENT mapping.
    """
    from workflow.modular.modules.annotation import AnnotationModule as _Ann

    class _Cellranger:
        name = "cellranger"

        def run(self, ctx):
            ctx.status("cellranger", True, "ok")

    class _QC:
        name = "qc"

        def run(self, ctx):
            ctx.status("qc", True, "ok")

    class _Ambient:
        name = "ambient_correction"

        def run(self, ctx):
            ctx.status("ambient_correction", True, "ok")

    class _Doublet:
        name = "doublet_detection"

        def run(self, ctx):
            ctx.status("doublet_detection", True, "ok")

    class _Clustering:
        name = "clustering"

        def run(self, ctx):
            n = ctx.adata.n_obs
            # Pre-correction clusters: 4 blocks.
            ctx.adata.obs["leiden"] = pd.Categorical(
                [str(i % 4) for i in range(n)]
            )
            ctx.status("clustering", True, "ok")

    class _BatchCorrection:
        name = "batch_correction"
        mutates_structure = True

        def run(self, ctx):
            n = ctx.adata.n_obs
            # Post-correction clusters DIFFER from pre-correction leiden.
            corrected = [str((i // (n // 4)) % 4) for i in range(n)]
            ctx.adata.obs["leiden"] = pd.Categorical(corrected)
            ctx.adata.obs["leiden_corrected"] = pd.Categorical(corrected)
            ctx.status("batch_correction", True, "ok")

    class _Annotation:
        name = "annotation"
        requires_keys = {"obs": ["leiden"]}

        def run(self, ctx):
            key = _Ann._resolve_leiden_key(ctx.adata)
            ctx.metadata["annotation_leiden_source"] = key
            labels = {
                "0": "TypeA", "1": "TypeB", "2": "TypeC", "3": "TypeD",
            }
            ctx.adata.obs["cell_type"] = (
                ctx.adata.obs[key].astype(str).map(labels).astype(str)
            )
            ctx.status("annotation", True, "ok")

    return {
        "cellranger": _Cellranger(),
        "qc": _QC(),
        "ambient_correction": _Ambient(),
        "doublet_detection": _Doublet(),
        "clustering": _Clustering(),
        "batch_correction": _BatchCorrection(),
        "annotation": _Annotation(),
    }


def _run_ws1_pipeline(workers, tmp_path, monkeypatch):
    from workflow.modular.config import PipelineConfig, CellRangerConfig
    from workflow.modular.context import PipelineContext

    cfg = PipelineConfig(
        project=f"ws1_{workers}",
        output_dir=tmp_path / f"out{workers}",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        optional_modules=["clustering", "batch_correction", "annotation"],
        parallel_workers=workers,
    )
    rng = np.random.default_rng(1)
    adata = ad.AnnData(rng.random((48, 6)).astype(np.float32))
    adata.obs_names = [f"cell_{i}" for i in range(48)]
    adata.obs["sample"] = (["A"] * 24) + (["B"] * 24)
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / f"run{workers}",
        figure_dir=tmp_path / f"run{workers}",
        table_dir=tmp_path / f"run{workers}",
        adata=adata,
    )
    ctx.run_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(pipe, "_prepare_output", lambda _cfg: ctx)
    monkeypatch.setattr(pipe, "_build_registry", _ws1_integration_registry)
    monkeypatch.setattr(
        pipe,
        "_save_manifest",
        lambda _ctx, _requested, _planned: _ctx.run_dir / "run_manifest.json",
    )
    pipe.run_pipeline(cfg)
    obs = ctx.adata.obs
    mapping = (
        obs.groupby(obs["leiden_corrected"].astype(str), observed=True)["cell_type"]
        .agg(lambda s: sorted(set(s)))
        .to_dict()
    )
    return mapping, ctx.metadata.get("annotation_leiden_source")


def test_celltype_leiden_mapping_identical_seq_vs_parallel(tmp_path, monkeypatch):
    """cell_type<->leiden_corrected mapping is identical at workers 1 and N."""
    map1, src1 = _run_ws1_pipeline(1, tmp_path / "a", monkeypatch)
    map2, src2 = _run_ws1_pipeline(2, tmp_path / "b", monkeypatch)
    # Each corrected cluster binds to exactly one cell_type.
    assert all(len(v) == 1 for v in map1.values()), map1
    assert all(len(v) == 1 for v in map2.values()), map2
    # Mapping identity across modes (the acceptance invariant).
    assert map1 == map2, (map1, map2)
    # annotation resolved against the corrected clusters in both modes.
    assert src1 == "leiden_corrected"
    assert src2 == "leiden_corrected"
