"""Regression tests for replicate-aware pseudobulk inference boundaries."""

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.modular.config import CellRangerConfig, PipelineConfig, PseudobulkConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule
from workflow.modular.cli import _resolve_optional_modules


def _ctx(tmp_path, *, sample_col="sample", sample_condition_pairs=None):
    pairs = sample_condition_pairs or [("S1", "A"), ("S2", "A"), ("S3", "B"), ("S4", "B")]
    samples = []
    conditions = []
    for sample, condition in pairs:
        samples.extend([sample] * 3)
        conditions.extend([condition] * 3)
    counts = np.tile(np.array([[10, 2]], dtype=np.int32), (len(samples), 1))
    adata = ad.AnnData(np.log1p(counts.astype(np.float32)))
    adata.layers["counts"] = sparse.csr_matrix(counts)
    adata.var_names = ["g1", "g2"]
    adata.obs["sample"] = samples
    adata.obs["cell_type"] = "Tumor"
    adata.obs["condition"] = conditions
    cfg = PipelineConfig(
        project="p",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        pseudobulk=PseudobulkConfig(
            sample_col=sample_col,
            contrast_col="condition",
            contrast_a="A",
            contrast_b="B",
        ),
    )
    return PipelineContext(
        cfg=cfg,
        run_dir=tmp_path,
        figure_dir=tmp_path,
        table_dir=tmp_path,
        adata=adata,
    )


def _patch_result(monkeypatch, *, test_used):
    def _fake(self, counts, meta, cond_col, ca, cb, min_samples_per_condition=2):
        return pd.DataFrame([{
            "gene": "g1",
            "log2FC": 1.0,
            "pvalue": 0.01,
            "padj": 0.02,
            "test_used": test_used,
        }])

    monkeypatch.setattr(PseudobulkDEModule, "_run_de", _fake)
    monkeypatch.setattr(PseudobulkDEModule, "_plot_volcano", lambda *args, **kwargs: None)
    monkeypatch.setattr(PseudobulkDEModule, "_plot_heatmap", lambda *args, **kwargs: None)


def test_pydeseq2_confirmatory_contract_is_claimable(monkeypatch, tmp_path):
    ctx = _ctx(tmp_path)
    _patch_result(monkeypatch, test_used="pydeseq2")

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_class"] == "replicate_aware_pseudobulk"
    assert ctx.metadata["pseudobulk_de_inference_status"] == "supported_confirmatory"
    assert ctx.metadata["pseudobulk_de_claimable"] is True
    results = pd.read_csv(tmp_path / "pseudobulk_de_results.csv")
    assert set(results["inference_status"]) == {"supported_confirmatory"}
    assert results["claimable"].all()
    assert set(results["n_biological_replicates_a"]) == {2}
    assert set(results["n_biological_replicates_b"]) == {2}
    counts = pd.read_csv(tmp_path / "pseudobulk_counts.csv")
    assert counts["claimable"].all()


def test_rank_backend_is_visibly_nonclaimable(monkeypatch, tmp_path):
    ctx = _ctx(tmp_path)
    _patch_result(monkeypatch, test_used="mann_whitney_u")

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "exploratory_nonclaimable_backend_fallback"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False
    results = pd.read_csv(tmp_path / "pseudobulk_de_results.csv")
    assert not results["claimable"].any()


def test_confirmatory_contract_requires_explicit_biological_sample_col(tmp_path):
    ctx = _ctx(tmp_path, sample_col=None)

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_status"] == (
        "skipped_missing_biological_sample_col"
    )
    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_missing_biological_sample_col"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_confirmatory_contract_rejects_missing_named_sample_column(tmp_path):
    ctx = _ctx(tmp_path, sample_col="donor_id")

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_missing_biological_sample_col"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_sample_mapped_to_multiple_conditions_is_rejected(tmp_path):
    ctx = _ctx(
        tmp_path,
        sample_condition_pairs=[("S1", "A"), ("S1", "B"), ("S2", "A"), ("S3", "B")],
    )

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_sample_condition_nonunique"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_minimum_unique_biological_replicates_is_enforced(tmp_path):
    ctx = _ctx(
        tmp_path,
        sample_condition_pairs=[("S1", "A"), ("S2", "B"), ("S3", "B")],
    )

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_insufficient_biological_replicates"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_explicit_condition_contrast_routes_to_pseudobulk_module():
    args = SimpleNamespace(
        optional_modules="clustering,differential_expression",
        select_integration=False,
        pseudobulk_contrast_json="",
        pseudobulk_contrast_a="A",
        pseudobulk_contrast_b="B",
    )

    modules = _resolve_optional_modules(args)

    assert modules[-1] == "pseudobulk_de"


def test_partial_condition_contrast_fails_loud():
    args = SimpleNamespace(
        optional_modules="clustering",
        select_integration=False,
        pseudobulk_contrast_json="",
        pseudobulk_contrast_a="A",
        pseudobulk_contrast_b="",
    )

    try:
        _resolve_optional_modules(args)
    except SystemExit as exc:
        assert "requires both" in str(exc)
    else:
        raise AssertionError("partial condition contrast was silently ignored")
