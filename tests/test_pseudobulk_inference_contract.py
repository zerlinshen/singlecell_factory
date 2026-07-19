"""Regression tests for replicate-aware pseudobulk inference boundaries."""

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
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
    assert all(np.issubdtype(dtype, np.number) for dtype in counts.dtypes)
    metadata = pd.read_csv(tmp_path / "pseudobulk_metadata.csv")
    assert len(metadata) == len(counts)
    assert metadata["claimable"].all()


def test_rank_backend_is_visibly_nonclaimable(monkeypatch, tmp_path):
    ctx = _ctx(tmp_path)
    _patch_result(monkeypatch, test_used="mann_whitney_u")

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "exploratory_nonclaimable_backend_fallback"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False
    assert ctx.metadata["pseudobulk_de_status"] == (
        "completed_nonclaimable_backend_fallback"
    )
    results = pd.read_csv(tmp_path / "pseudobulk_de_results.csv")
    assert not results["claimable"].any()
    metadata = pd.read_csv(tmp_path / "pseudobulk_metadata.csv")
    assert not metadata["claimable"].any()
    assert PseudobulkDEModule._plot_claim_label(results).startswith("NON-CLAIMABLE:")


def test_unsupported_pydeseq2_contrast_api_downgrades_to_rank_backend(
    monkeypatch, tmp_path
):
    ctx = _ctx(tmp_path)

    def _unsupported(*args, **kwargs):
        raise TypeError("contrast keyword unsupported")

    def _rank(*args, **kwargs):
        return pd.DataFrame(
            [{"gene": "g1", "log2FC": 1.0, "pvalue": 0.01, "padj": 0.02}]
        )

    monkeypatch.setattr(PseudobulkDEModule, "_de_pydeseq2", _unsupported)
    monkeypatch.setattr(PseudobulkDEModule, "_de_ranktest", _rank)
    monkeypatch.setattr(PseudobulkDEModule, "_plot_volcano", lambda *a, **k: None)
    monkeypatch.setattr(PseudobulkDEModule, "_plot_heatmap", lambda *a, **k: None)

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_claimable"] is False
    assert ctx.metadata["pseudobulk_de_status"] == (
        "completed_nonclaimable_backend_fallback"
    )


def test_confirmatory_contract_requires_explicit_biological_sample_col(tmp_path):
    ctx = _ctx(tmp_path, sample_col=None)

    with pytest.raises(ValueError, match="explicit biological"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_status"] == "failed_missing_biological_sample_col"
    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_missing_biological_sample_col"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_confirmatory_contract_rejects_missing_named_sample_column(tmp_path):
    ctx = _ctx(tmp_path, sample_col="donor_id")

    with pytest.raises(ValueError, match="explicit biological"):
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

    with pytest.raises(ValueError, match="exactly one condition"):
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

    with pytest.raises(ValueError, match="biological replicates"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_insufficient_biological_replicates"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_explicit_contrast_missing_counts_records_then_raises(tmp_path):
    ctx = _ctx(tmp_path)
    del ctx.adata.layers["counts"]

    with pytest.raises(ValueError, match="raw UMI counts"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_status"] == "failed_missing_counts_layer"
    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_missing_raw_counts"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_explicit_contrast_missing_grouping_column_records_then_raises(tmp_path):
    ctx = _ctx(tmp_path)
    del ctx.adata.obs["cell_type"]

    with pytest.raises(ValueError, match="No grouping column"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_status"] == "failed_missing_grouping_columns"
    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_missing_grouping_columns"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_explicit_contrast_invalid_label_records_then_raises(tmp_path):
    ctx = _ctx(tmp_path)
    ctx.cfg.pseudobulk.contrast_b = "missing"

    with pytest.raises(ValueError, match="must both exist"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_invalid_contrast_labels"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_null_biological_identity_records_then_raises(tmp_path):
    ctx = _ctx(tmp_path)
    ctx.adata.obs.loc[ctx.adata.obs.index[0], "sample"] = None

    with pytest.raises(ValueError, match="non-null and non-blank"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_invalid_biological_identifiers"
    )
    assert ctx.metadata["pseudobulk_de_claimable"] is False


def test_numeric_condition_labels_are_normalized(monkeypatch, tmp_path):
    ctx = _ctx(tmp_path)
    ctx.adata.obs["condition"] = ctx.adata.obs["condition"].map({"A": 0, "B": 1})
    ctx.cfg.pseudobulk.contrast_a = "0"
    ctx.cfg.pseudobulk.contrast_b = "1"
    _patch_result(monkeypatch, test_used="pydeseq2")

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_claimable"] is True


def test_result_replicate_counts_are_exact_for_each_cell_group(monkeypatch):
    pb_meta = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3", "S4", "S1", "S2", "S5", "S3", "S4", "S6"],
            "cell_type": ["Tumor"] * 4 + ["Stromal"] * 6,
            "condition": ["A", "A", "B", "B", "A", "A", "A", "B", "B", "B"],
        }
    )
    pb_counts = pd.DataFrame({"g1": np.arange(1, len(pb_meta) + 1)})
    _patch_result(monkeypatch, test_used="pydeseq2")

    results = PseudobulkDEModule()._run_confirmatory(
        pb_counts,
        pb_meta,
        "sample",
        "cell_type",
        "condition",
        [{"name": "A_vs_B", "contrast_a": "A", "contrast_b": "B"}],
        min_samples_per_condition=2,
    )

    counts_by_group = {
        row.group: (row.n_biological_replicates_a, row.n_biological_replicates_b)
        for row in results.itertuples()
    }
    assert counts_by_group == {"Tumor|A_vs_B": (2, 2), "Stromal|A_vs_B": (3, 3)}


def test_json_contrast_col_is_honored_once(monkeypatch, tmp_path):
    ctx = _ctx(tmp_path)
    contract = tmp_path / "contrasts.json"
    contract.write_text(
        '[{"name":"A_vs_B","contrast_col":"condition",'
        '"contrast_a":"A","contrast_b":"B"}]',
        encoding="utf-8",
    )
    ctx.cfg.pseudobulk.contrast_json = contract
    ctx.cfg.pseudobulk.contrast_a = None
    ctx.cfg.pseudobulk.contrast_b = None
    ctx.cfg.pseudobulk.contrast_col = None
    _patch_result(monkeypatch, test_used="pydeseq2")

    PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_contrast_contract"]["contrast_col"] == "condition"
    assert ctx.metadata["pseudobulk_de_claimable"] is True


def test_json_mixed_contrast_columns_record_then_raise(tmp_path):
    ctx = _ctx(tmp_path)
    ctx.adata.obs["disease"] = ctx.adata.obs["condition"]
    contract = tmp_path / "contrasts.json"
    contract.write_text(
        '[{"contrast_col":"condition","contrast_a":"A","contrast_b":"B"},'
        '{"contrast_col":"disease","contrast_a":"A","contrast_b":"B"}]',
        encoding="utf-8",
    )
    ctx.cfg.pseudobulk.contrast_json = contract
    ctx.cfg.pseudobulk.contrast_a = None
    ctx.cfg.pseudobulk.contrast_b = None
    ctx.cfg.pseudobulk.contrast_col = None

    with pytest.raises(ValueError, match="one consistent contrast_col"):
        PseudobulkDEModule().run(ctx)

    assert ctx.metadata["pseudobulk_de_inference_status"] == (
        "not_testable_mixed_contrast_columns"
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
