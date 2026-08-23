"""Fail-closed contracts for staged scATAC pseudobulk DA.

Only the two tests marked ``r_contract`` execute the DESeq2/edgeR adapter.
Everything else exercises Python-side axis, replicate, sparse-integer, and
provenance contracts without treating a planted fixture as real-data evidence.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import jsonschema
import numpy as np
import pandas as pd
import pytest
from scipy import io as spio
from scipy import sparse

from workflow.modular.cli import build_parser, resolve_requested_modules
from workflow.modular.config import CellRangerConfig, PipelineConfig, ScatacPseudobulkDAConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules import scatac_pseudobulk_da as da_module
from workflow.modular.modules.scatac_pseudobulk_da import ScatacContractError, ScatacPseudobulkDAModule

SCHEMA_PATH = Path(__file__).resolve().parents[1] / "contracts" / "scatac_pseudobulk_da.schema.json"
INT64_MAX = np.iinfo(np.int64).max
R_INT_MAX = np.iinfo(np.int32).max


def _schema() -> dict:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def _definition_schema(name: str) -> dict:
    schema = _schema()
    return {"definitions": schema["definitions"], **schema["definitions"][name]}


def _axis_certificate(peak_ids: list[str]) -> dict[str, object]:
    return {
        "schema_version": "atac_peak_axis/v1",
        "producer": "atac_ingest",
        "source_status": "test_fixture_bound_axis",
        "peak_id_col": "peak_id",
        "n_peaks": len(peak_ids),
        "ordered_peak_id_sha256": hashlib.sha256("\0".join(peak_ids).encode("utf-8")).hexdigest(),
    }


def _from_sample_counts(
    sample_counts: np.ndarray,
    *,
    groups: list[str] | None = None,
    conditions: list[str] | None = None,
    cells_per_sample: int = 3,
) -> ad.AnnData:
    """Split exact sample-level peak totals across cells without changing sums."""
    sample_counts = np.asarray(sample_counts, dtype=np.int64)
    n_samples, n_peaks = sample_counts.shape
    if conditions is None:
        assert n_samples % 2 == 0
        conditions = ["CTRL"] * (n_samples // 2) + ["STIM"] * (n_samples // 2)
    if groups is None:
        groups = ["GroupA"] * n_samples
    rows: list[np.ndarray] = []
    obs_rows: list[dict[str, str]] = []
    for sample_index, totals in enumerate(sample_counts):
        quotient, remainder = np.divmod(totals, cells_per_sample)
        for cell_index in range(cells_per_sample):
            rows.append(quotient + (remainder > cell_index).astype(np.int64))
            obs_rows.append(
                {
                    "sample": f"S{sample_index + 1}",
                    "condition": conditions[sample_index],
                    "cell_type": groups[sample_index],
                }
            )
    peak_ids = [f"peak_{index:03d}" for index in range(n_peaks)]
    adata = ad.AnnData(
        X=sparse.csr_matrix((len(rows), 3), dtype=np.float32),
        obs=pd.DataFrame(obs_rows, index=[f"cell_{index:03d}" for index in range(len(rows))]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    adata.obsm["atac_peaks"] = sparse.csr_matrix(np.vstack(rows), dtype=np.int64)
    adata.uns["atac_var"] = pd.DataFrame(
        {
            "peak_id": peak_ids,
            "name": [f"chr1:{index * 100}-{index * 100 + 50}" for index in range(n_peaks)],
        }
    )
    adata.uns["atac_peak_axis"] = _axis_certificate(peak_ids)
    return adata


def _balanced_adata(n_peaks: int = 24) -> ad.AnnData:
    base = np.arange(20, 20 + n_peaks, dtype=np.int64)
    sample_counts = np.vstack([base + offset for offset in (-2, 0, 2, -2, 0, 2)])
    return _from_sample_counts(sample_counts)


def _config(**overrides) -> ScatacPseudobulkDAConfig:
    base = {
        "sample_col": "sample",
        "group_col": "cell_type",
        "condition_col": "condition",
        "peak_id_col": "peak_id",
        "test_level": "STIM",
        "reference_level": "CTRL",
        "mode": "confirmatory_da",
        "min_samples_per_condition": 2,
        "min_total_count": 5,
        "fdr_threshold": 0.05,
        "abs_log2fc_threshold": 1.0,
    }
    base.update(overrides)
    return ScatacPseudobulkDAConfig(**base)


def _context(tmp_path: Path, adata: ad.AnnData, cfg: ScatacPseudobulkDAConfig) -> PipelineContext:
    run_dir = tmp_path / "run"
    for path in (run_dir, run_dir / "figures", run_dir / "tables"):
        path.mkdir(parents=True, exist_ok=True)
    pipeline_cfg = PipelineConfig(
        project="scatac-contract-test",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        scatac_pseudobulk_da=cfg,
    )
    return PipelineContext(
        cfg=pipeline_cfg,
        run_dir=run_dir,
        figure_dir=run_dir / "figures",
        table_dir=run_dir / "tables",
        adata=adata,
    )


def test_schema_v2_parses_and_rejects_legacy_shape():
    schema = _schema()
    assert schema["definitions"]["aggregation_manifest"]["properties"]["schema_version"]["enum"] == ["2.0.0"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({"schema_version": "1.0.0"}, schema["definitions"]["aggregation_manifest"])


def test_explicit_peak_column_and_axis_certificate_are_required():
    adata = _balanced_adata()
    module = ScatacPseudobulkDAModule()
    missing_column = _config(peak_id_col="does_not_exist")
    with pytest.raises(ScatacContractError, match="Index fallback is prohibited"):
        module.validate_scatac_contract(adata, missing_column)

    absent_certificate = _balanced_adata()
    del absent_certificate.uns["atac_peak_axis"]
    with pytest.raises(ScatacContractError, match="producer-emitted"):
        module.validate_scatac_contract(absent_certificate, _config())

    reordered = _balanced_adata()
    reordered.uns["atac_var"] = reordered.uns["atac_var"].iloc[::-1].reset_index(drop=True)
    with pytest.raises(ScatacContractError, match="Ordered peak-ID digest"):
        module.validate_scatac_contract(reordered, _config())


def test_peak_metadata_rejects_duplicate_null_and_length_mismatch():
    module = ScatacPseudobulkDAModule()
    duplicate = _balanced_adata()
    duplicate.uns["atac_var"].loc[1, "peak_id"] = duplicate.uns["atac_var"].loc[0, "peak_id"]
    with pytest.raises(ScatacContractError, match="strictly unique"):
        module.validate_scatac_contract(duplicate, _config())

    null = _balanced_adata()
    null.uns["atac_var"].loc[0, "peak_id"] = None
    with pytest.raises(ScatacContractError, match="null"):
        module.validate_scatac_contract(null, _config())

    mismatch = _balanced_adata()
    mismatch.uns["atac_var"] = mismatch.uns["atac_var"].iloc[:-1].copy()
    with pytest.raises(ScatacContractError, match="row count"):
        module.validate_scatac_contract(mismatch, _config())


@pytest.mark.parametrize("field", ["sample", "cell_type", "condition"])
@pytest.mark.parametrize("bad_value", [None, "", "   ", "NA", "<NA>"])
def test_design_labels_reject_null_blank_and_na_like_values(field, bad_value):
    adata = _balanced_adata()
    adata.obs.loc[adata.obs.index[0], field] = bad_value
    with pytest.raises(ScatacContractError, match="null|blank|NA-like"):
        ScatacPseudobulkDAModule().validate_scatac_contract(adata, _config())


def test_exact_two_level_design_and_one_condition_per_sample_are_fail_closed():
    extra = _balanced_adata()
    extra.obs.loc[extra.obs.index[0], "condition"] = "THIRD"
    with pytest.raises(ScatacContractError, match="exactly the configured two"):
        ScatacPseudobulkDAModule().validate_scatac_contract(extra, _config())

    inconsistent = _balanced_adata()
    same_sample = inconsistent.obs.index[inconsistent.obs["sample"] == "S1"]
    inconsistent.obs.loc[same_sample[0], "condition"] = "STIM"
    with pytest.raises(ScatacContractError, match="Sample-to-condition integrity"):
        ScatacPseudobulkDAModule().validate_scatac_contract(inconsistent, _config())


def test_group_allowlist_precedes_exact_contrast_check():
    adata = _balanced_adata()
    adata.obs.loc[adata.obs["sample"] == "S1", "cell_type"] = "Excluded"
    adata.obs.loc[adata.obs["sample"] == "S1", "condition"] = "THIRD"
    # The out-of-scope group cannot introduce a third level after allowlisting.
    info = ScatacPseudobulkDAModule().validate_scatac_contract(adata, _config(groups=("GroupA",)))
    assert info["groups"] == ["GroupA"]


@pytest.mark.parametrize("value", [0, 1])
def test_config_rejects_replication_floor_below_two(value):
    with pytest.raises(ValueError, match=">= 2"):
        _config(min_samples_per_condition=value)


def test_runtime_rejects_bypassed_replication_floor():
    cfg = SimpleNamespace(**_config().__dict__)
    cfg.min_samples_per_condition = 1
    with pytest.raises(ScatacContractError, match="must be >= 2"):
        ScatacPseudobulkDAModule().validate_scatac_contract(_balanced_adata(), cfg)


def test_malformed_counts_and_signed_int64_range_fail_before_aggregation():
    module = ScatacPseudobulkDAModule()
    negative = _balanced_adata()
    negative.obsm["atac_peaks"].data[0] = -1
    with pytest.raises(ScatacContractError, match="negative"):
        module.validate_scatac_contract(negative, _config())

    fractional = _balanced_adata()
    fractional.obsm["atac_peaks"].data = fractional.obsm["atac_peaks"].data.astype(float)
    fractional.obsm["atac_peaks"].data[0] = 1.5
    with pytest.raises(ScatacContractError, match="fractional"):
        module.validate_scatac_contract(fractional, _config())

    nonfinite = _balanced_adata()
    nonfinite.obsm["atac_peaks"].data = nonfinite.obsm["atac_peaks"].data.astype(float)
    nonfinite.obsm["atac_peaks"].data[0] = np.inf
    with pytest.raises(ScatacContractError, match="non-finite"):
        module.validate_scatac_contract(nonfinite, _config())

    dense = _balanced_adata()
    dense.obsm["atac_peaks"] = dense.obsm["atac_peaks"].toarray()
    with pytest.raises(ScatacContractError, match="Dense matrix rejected"):
        module.validate_scatac_contract(dense, _config())

    outside_int64 = _balanced_adata()
    outside_int64.obsm["atac_peaks"].data = outside_int64.obsm["atac_peaks"].data.astype(np.uint64)
    outside_int64.obsm["atac_peaks"].data[0] = np.uint64(INT64_MAX) + np.uint64(1)
    with pytest.raises(ScatacContractError, match="outside signed int64"):
        module.validate_scatac_contract(outside_int64, _config())


def test_sparse_aggregation_matches_trusted_integer_reference_and_manifest(tmp_path):
    sample_counts = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]], dtype=np.int64)
    adata = _from_sample_counts(sample_counts, cells_per_sample=2)
    cfg = _config(mode="aggregation_only")
    ctx = _context(tmp_path, adata, cfg)
    ScatacPseudobulkDAModule().run(ctx)
    out_dir = ctx.run_dir / "scatac_pseudobulk_da"
    with gzip.open(out_dir / "aggregation" / "peak_by_pseudobulk_sample_counts.mtx.gz", "rb") as handle:
        observed = spio.mmread(handle).toarray()
    np.testing.assert_array_equal(observed, sample_counts.T)
    manifest = json.loads((out_dir / "aggregation" / "aggregation_manifest.json").read_text())
    jsonschema.validate(manifest, _definition_schema("aggregation_manifest"))
    assert manifest["count_invariants"]["input_total_exact"] == str(int(sample_counts.sum()))
    assert manifest["backend_observation"]["resolved"] == "cpu"
    claim = json.loads((out_dir / "claimability_manifest.json").read_text())
    jsonschema.validate(claim, _definition_schema("claimability_manifest"))
    assert claim["technical_inference_status"] == "aggregation_completed"
    assert not (out_dir / "groups").exists()


def test_sparse_overflow_fails_with_a_machine_readable_manifest(tmp_path):
    # Every individual stored value is representable, but their exact selected
    # total exceeds int64 and must be rejected before sparse multiplication.
    counts = np.array([[INT64_MAX // 2 + 1], [INT64_MAX // 2 + 1], [1], [1]], dtype=np.int64)
    adata = _from_sample_counts(counts, cells_per_sample=1)
    ctx = _context(tmp_path, adata, _config(mode="aggregation_only"))
    with pytest.raises(ScatacContractError, match="overflow"):
        ScatacPseudobulkDAModule().run(ctx)
    failure = json.loads((ctx.run_dir / "scatac_pseudobulk_da" / "failure_manifest.json").read_text())
    jsonschema.validate(failure, _definition_schema("failure_manifest"))
    assert failure["status"] == "aggregation_failed"


def test_aggregation_never_densifies_source_even_with_a_selected_subset(tmp_path, monkeypatch):
    adata = _balanced_adata()
    adata.obs.loc[adata.obs["sample"].isin(["S1", "S4"]), "cell_type"] = "Excluded"
    cfg = _config(mode="aggregation_only", groups=("GroupA",))
    original_toarray = sparse.csr_matrix.toarray
    original_todense = sparse.csr_matrix.todense

    def forbidden_toarray(self, *args, **kwargs):
        raise AssertionError("cell-by-peak aggregation called toarray()")

    def forbidden_todense(self, *args, **kwargs):
        raise AssertionError("cell-by-peak aggregation called todense()")

    monkeypatch.setattr(sparse.csr_matrix, "toarray", forbidden_toarray)
    monkeypatch.setattr(sparse.csr_matrix, "todense", forbidden_todense)
    try:
        ScatacPseudobulkDAModule().run(_context(tmp_path, adata, cfg))
    finally:
        monkeypatch.setattr(sparse.csr_matrix, "toarray", original_toarray)
        monkeypatch.setattr(sparse.csr_matrix, "todense", original_todense)


def test_r_integer_boundary_is_checked_before_any_subprocess(tmp_path, monkeypatch):
    adata = _balanced_adata()
    cfg = _config()
    module = ScatacPseudobulkDAModule()
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    info = module.validate_scatac_contract(adata, cfg)
    aggregation = module.aggregate_sparse_counts(adata, cfg, info, out_dir)
    aggregation["aggregate"].data[0] = R_INT_MAX + 1
    monkeypatch.setattr(da_module, "_resolve_conda_executable", lambda: Path("/not-used/conda"))
    with pytest.raises(da_module.ScatacInferenceError, match="integer-range"):
        module.run_inference(cfg, info, aggregation, out_dir)


def test_cpu_only_implementation_has_no_gpu_import_or_route():
    source = Path(da_module.__file__).read_text(encoding="utf-8")
    assert "import cupy" not in source
    assert "import cudf" not in source
    assert "aggregation_backend != \"cpu\"" in source
    assert _config().aggregation_backend == "cpu"


def test_pipeline_fails_closed_before_an_unrelated_module_runs(tmp_path):
    import workflow.modular.pipeline as pipeline

    adata = _balanced_adata()
    del adata.uns["atac_peak_axis"]
    ctx = _context(tmp_path, adata, _config(mode="aggregation_only"))
    unrelated_called = []

    class UnrelatedModule:
        name = "unrelated"

        def run(self, _ctx):
            unrelated_called.append(True)

    with pytest.raises(ValueError, match="missing required keys: uns.atac_peak_axis"):
        pipeline._execute_tier(
            ["scatac_pseudobulk_da", "unrelated"],
            {
                "scatac_pseudobulk_da": ScatacPseudobulkDAModule(),
                "unrelated": UnrelatedModule(),
            },
            ctx,
            mandatory=set(),
            max_workers=2,
        )
    assert unrelated_called == []


def test_pipeline_executes_scatac_module_sequentially_without_anndata_copy(tmp_path, monkeypatch):
    import workflow.modular.pipeline as pipeline

    adata = _balanced_adata()
    ctx = _context(tmp_path, adata, _config(mode="aggregation_only"))

    def forbidden_copy(self, *args, **kwargs):
        raise AssertionError("force-sequential scATAC module must not copy AnnData")

    monkeypatch.setattr(type(adata), "copy", forbidden_copy)
    pipeline._execute_tier(
        ["scatac_pseudobulk_da"],
        {"scatac_pseudobulk_da": ScatacPseudobulkDAModule()},
        ctx,
        mandatory=set(),
        max_workers=2,
    )
    assert ctx.metadata["scatac_pseudobulk_da_status"] == "completed"


def test_cli_requires_complete_explicit_contract_and_forwards_cpu_only_defaults():
    parser = build_parser()
    partial = parser.parse_args(["--project", "p", "--sample-root", "/tmp", "--scatac-da-sample-col", "sample"])
    with pytest.raises(SystemExit, match="scATAC pseudobulk DA requires explicit"):
        resolve_requested_modules(partial)

    complete = parser.parse_args([
        "--project", "p", "--sample-root", "/tmp",
        "--scatac-da-sample-col", "sample", "--scatac-da-group-col", "cell_type",
        "--scatac-da-condition-col", "condition", "--scatac-da-peak-id-col", "peak_id",
        "--scatac-da-test-level", "STIM", "--scatac-da-reference-level", "CTRL",
    ])
    assert "scatac_pseudobulk_da" in resolve_requested_modules(complete)
    invalid_floor = parser.parse_args([
        "--project", "p", "--sample-root", "/tmp", "--optional-modules", "scatac_pseudobulk_da",
        "--scatac-da-sample-col", "sample", "--scatac-da-group-col", "cell_type",
        "--scatac-da-condition-col", "condition", "--scatac-da-peak-id-col", "peak_id",
        "--scatac-da-test-level", "STIM", "--scatac-da-reference-level", "CTRL",
        "--scatac-da-min-samples-per-condition", "1",
    ])
    with pytest.raises(SystemExit, match="at least 2"):
        resolve_requested_modules(invalid_floor)


def _actual_engine_counts(*, null: bool) -> tuple[np.ndarray, list[str], list[str]]:
    n_peaks = 80
    base = np.array([80 + (index % 11) * 7 for index in range(n_peaks)], dtype=np.int64)
    # The condition-matched rows are deliberately identical in the null case,
    # while peak-specific biological-sample variation makes this a real
    # DESeq2/edgeR fit rather than a constant-count toy that cannot estimate a
    # dispersion trend.
    rng = np.random.default_rng(20260823)
    reference = np.vstack([
        np.maximum(2, base + rng.integers(-35, 36, size=n_peaks, dtype=np.int64))
        for _ in range(3)
    ])
    stimulated = reference.copy()
    if not null:
        stimulated[:, :5] *= 8
        stimulated[:, 5:10] = np.maximum(2, stimulated[:, 5:10] // 8)
    return np.vstack([reference, stimulated]), [f"peak_{index:03d}" for index in range(n_peaks)], ["UP"] * 5 + ["DOWN"] * 5


@pytest.mark.r_contract
def test_actual_engine_recovers_implanted_peak_ids_and_directions(tmp_path):
    counts, implanted_ids, directions = _actual_engine_counts(null=False)
    adata = _from_sample_counts(counts, cells_per_sample=4)
    ctx = _context(tmp_path, adata, _config(min_total_count=10, abs_log2fc_threshold=1.0))
    ScatacPseudobulkDAModule().run(ctx)
    group_dir = next((ctx.run_dir / "scatac_pseudobulk_da" / "groups").glob("grp_*"))
    with gzip.open(group_dir / "da_results.tsv.gz", "rt", encoding="utf-8") as handle:
        result = pd.read_csv(handle, sep="\t").set_index("peak_id")
    for peak_id, expected_direction in zip(implanted_ids[:10], directions, strict=True):
        row = result.loc[peak_id]
        assert bool(row["deseq2_significant"]), peak_id
        assert bool(row["edgeR_significant"]), peak_id
        assert bool(row["cross_engine_supported"]), peak_id
        if expected_direction == "UP":
            assert row["log2FC"] > 0 and row["edgeR_log2FC"] > 0
        else:
            assert row["log2FC"] < 0 and row["edgeR_log2FC"] < 0
    manifest = json.loads((group_dir / "inference_manifest.json").read_text())
    jsonschema.validate(manifest, _definition_schema("inference_manifest"))
    assert manifest["argv"]["nul_delimited_sha256"]


@pytest.mark.r_contract
def test_actual_engine_balanced_null_has_zero_individual_engine_false_positives(tmp_path):
    counts, _, _ = _actual_engine_counts(null=True)
    adata = _from_sample_counts(counts, cells_per_sample=4)
    ctx = _context(tmp_path, adata, _config(min_total_count=10, abs_log2fc_threshold=1.0))
    ScatacPseudobulkDAModule().run(ctx)
    group_dir = next((ctx.run_dir / "scatac_pseudobulk_da" / "groups").glob("grp_*"))
    with gzip.open(group_dir / "da_results.tsv.gz", "rt", encoding="utf-8") as handle:
        result = pd.read_csv(handle, sep="\t")
    assert int(result["deseq2_significant"].sum()) == 0
    assert int(result["edgeR_significant"].sum()) == 0
    assert int(result["cross_engine_supported"].sum()) == 0
