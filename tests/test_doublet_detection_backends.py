from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.modular.config import CellRangerConfig, DoubletConfig, PipelineConfig, QCConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.doublet_detection import DoubletDetectionModule


def _adata(n_obs: int = 30, n_vars: int = 60) -> ad.AnnData:
    rng = np.random.default_rng(7)
    x = sparse.csr_matrix(rng.poisson(1.2, size=(n_obs, n_vars)).astype("float32"))
    obs = pd.DataFrame(index=[f"cell{i:03d}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(n_vars)])
    return ad.AnnData(X=x, obs=obs, var=var)


def _ctx(tmp_path: Path, doublet: DoubletConfig) -> PipelineContext:
    cfg = PipelineConfig(
        project="doublet-test",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        qc=QCConfig(min_genes=0, min_counts=0),
        doublet=doublet,
        random_state=11,
    )
    fig_dir = tmp_path / "figures"
    table_dir = tmp_path / "tables"
    fig_dir.mkdir()
    table_dir.mkdir()
    return PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=fig_dir, table_dir=table_dir)


def test_scdblfinder_backend_records_effective_metadata(monkeypatch, tmp_path):
    adata = _adata()
    calls = np.zeros(adata.n_obs, dtype=bool)
    calls[[1, 3, 5]] = True
    scores = np.linspace(0.0, 1.0, adata.n_obs, dtype=np.float32)

    def fake_scdbl(self, adata_arg, cfg, ctx):
        assert adata_arg is adata
        return scores, calls, None

    monkeypatch.setattr(DoubletDetectionModule, "_run_scdblfinder_via_r", fake_scdbl)
    ctx = _ctx(tmp_path, DoubletConfig(backend="scdblfinder", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_backend_requested"] == "scdblfinder"
    assert ctx.metadata["doublet_backend"] == "scdblfinder"
    assert ctx.metadata["doublet_method"] == "scdblfinder"
    assert ctx.metadata["doublets_detected"] == 3
    assert ctx.adata.obs["predicted_doublet"].sum() == 3


def test_scdbl_alias_maps_to_scdblfinder(monkeypatch, tmp_path):
    adata = _adata()
    calls = np.zeros(adata.n_obs, dtype=bool)
    scores = np.linspace(0.0, 1.0, adata.n_obs, dtype=np.float32)

    def fake_scdbl(self, adata_arg, cfg, ctx):
        return scores, calls, None

    monkeypatch.setattr(DoubletDetectionModule, "_run_scdblfinder_via_r", fake_scdbl)
    ctx = _ctx(tmp_path, DoubletConfig(backend="scdbl", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_backend_requested"] == "scdbl"
    assert ctx.metadata["doublet_backend"] == "scdblfinder"
    assert ctx.metadata["doublet_method"] == "scdblfinder"


def test_unknown_backend_records_fallback_reason_on_small_dataset(tmp_path):
    adata = _adata(n_obs=10, n_vars=20)
    ctx = _ctx(tmp_path, DoubletConfig(backend="bogus", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_backend_requested"] == "bogus"
    assert ctx.metadata["doublet_backend"] == "scrublet"
    assert ctx.metadata["doublet_backend_fallback_reason"] == "unknown_backend:bogus"
    assert ctx.metadata["doublet_method"] == "fallback_all_singlets"


def test_consensus_can_use_scrublet_scdblfinder_pair(monkeypatch, tmp_path):
    adata = _adata()
    scrub_calls = np.zeros(adata.n_obs, dtype=bool)
    scrub_calls[[0, 1]] = True
    scdbl_calls = np.zeros(adata.n_obs, dtype=bool)
    scdbl_calls[[1, 2, 3]] = True

    def fake_scrub(self, adata_arg, cfg, ctx):
        return np.linspace(0.0, 1.0, adata_arg.n_obs, dtype=np.float32), scrub_calls, 0.5

    def fake_scdbl(self, adata_arg, cfg, ctx):
        return np.linspace(1.0, 0.0, adata_arg.n_obs, dtype=np.float32), scdbl_calls, None

    monkeypatch.setattr(DoubletDetectionModule, "_run_scrublet_only", fake_scrub)
    monkeypatch.setattr(DoubletDetectionModule, "_run_scdblfinder_via_r", fake_scdbl)
    ctx = _ctx(
        tmp_path,
        DoubletConfig(
            backend="consensus",
            consensus_pair="scrublet_scdblfinder",
            consensus_logic="or",
            remove_doublets=False,
        ),
    )
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_method"] == "consensus_or_scrublet_scdblfinder"
    assert ctx.metadata["doublet_consensus_pair"] == "scrublet_scdblfinder"
    assert ctx.metadata["doublet_consensus_backends"] == "scrublet,scdblfinder"
    assert ctx.metadata["doublet_consensus_per_backend_rates"]["final"] == 4 / adata.n_obs
    assert "scrublet_score" in ctx.adata.obs
    assert "scdblfinder_score" in ctx.adata.obs
    assert ctx.adata.obs["predicted_doublet"].sum() == 4


def test_resolve_scdblfinder_driver_from_script_dir(monkeypatch, tmp_path):
    driver = tmp_path / "scdblfinder_run.R"
    driver.write_text("# dummy\n", encoding="utf-8")
    monkeypatch.setenv("R_MULTIOMICS_SCRIPT_DIR", str(tmp_path))

    assert DoubletDetectionModule._resolve_scdblfinder_driver() == driver


def test_null_fallback_gate_raises_without_opt_in(monkeypatch):
    """Failure-based all-singlets fallback must raise unless explicitly opted in.

    Mirrors the SC_ALLOW_WELCH_FALLBACK / Principle 9 / Plan F-3 contract used
    by differential_expression.py: a silent algorithmic compromise on a
    mandatory QC stage is not allowed.
    """
    monkeypatch.delenv(DoubletDetectionModule._NULL_FALLBACK_ENV, raising=False)
    primary = RuntimeError("rsc.pp.scrublet boom")
    secondary = RuntimeError("cpu scrublet also boom")
    try:
        DoubletDetectionModule._require_null_fallback_opt_in(primary, secondary)
    except RuntimeError as e:
        msg = str(e)
        assert DoubletDetectionModule._NULL_FALLBACK_ENV in msg
        assert "banned by default" in msg
        assert "rsc.pp.scrublet boom" in msg
        assert "cpu scrublet also boom" in msg
    else:
        raise AssertionError("expected RuntimeError when opt-in env is unset")


def test_null_fallback_gate_allows_with_opt_in(monkeypatch):
    """Explicit SC_ALLOW_DOUBLET_NULL_FALLBACK=1 lets the gate return cleanly."""
    monkeypatch.setenv(DoubletDetectionModule._NULL_FALLBACK_ENV, "1")
    # Must not raise.
    DoubletDetectionModule._require_null_fallback_opt_in(
        RuntimeError("primary boom"), RuntimeError("secondary boom")
    )


def test_tiny_dataset_path_records_method_actually_used(tmp_path):
    """The legitimate tiny-dataset degenerate path is exempt from the gate
    but must still record doublet_method_actually_used in metadata so the
    run manifest unambiguously identifies the engine that produced calls."""
    adata = _adata(n_obs=10, n_vars=20)
    ctx = _ctx(tmp_path, DoubletConfig(backend="scrublet", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_method"] == "fallback_all_singlets"
    assert ctx.metadata["doublet_method_actually_used"] == "fallback_all_singlets"


def test_sample_aware_expected_rate_contract_is_recorded_once(tmp_path):
    adata = _adata(n_obs=80)
    adata.obs["sample"] = ["small"] * 20 + ["large"] * 60
    cfg = DoubletConfig(
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)

    effective = DoubletDetectionModule._resolve_expected_rate_contract(
        ctx, cfg, adata
    )

    expected_small = DoubletDetectionModule._resolved_expected_rate(cfg, 20)
    expected_large = DoubletDetectionModule._resolved_expected_rate(cfg, 60)
    expected_effective = (20 * expected_small + 60 * expected_large) / 80
    assert effective == expected_effective
    assert ctx.metadata["doublet_expected_rate_requested"] == 0.06
    assert ctx.metadata["doublet_expected_rate_resolved_per_sample"] == {
        "small": round(expected_small, 6),
        "large": round(expected_large, 6),
    }
    assert ctx.metadata["doublet_expected_rate_sample_sizes"] == {
        "small": 20,
        "large": 60,
    }
    assert ctx.metadata["doublet_expected_rate_effective"] == round(
        expected_effective, 6
    )
    assert ctx.metadata["doublet_expected_rate_aggregation_rule"] == (
        "cell_weighted_mean_of_resolved_per_sample"
    )


def test_backends_record_per_sample_or_documented_aggregate_prior(tmp_path):
    adata = _adata(n_obs=80)
    adata.obs["sample"] = ["small"] * 20 + ["large"] * 60
    cfg = DoubletConfig(
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)
    effective = DoubletDetectionModule._resolve_expected_rate_contract(
        ctx, cfg, adata
    )

    per_sample = DoubletDetectionModule._record_backend_expected_rate(
        ctx, "scrublet", mode="per_sample"
    )
    for backend in ("doubletfinder", "scdblfinder", "consensus"):
        recorded = DoubletDetectionModule._record_backend_expected_rate(
            ctx, backend, mode="aggregate"
        )
        assert recorded == effective

    priors = ctx.metadata["doublet_expected_rate_by_backend"]
    assert set(priors) == {
        "scrublet",
        "doubletfinder",
        "scdblfinder",
        "consensus",
    }
    assert per_sample == ctx.metadata["doublet_expected_rate_resolved_per_sample"]
    assert priors["scrublet"] == {
        "requested": 0.06,
        "mode": "per_sample",
        "resolved_per_sample": per_sample,
        "sample_sizes": {"small": 20, "large": 60},
    }
    aggregate_records = [
        priors[name] for name in ("doubletfinder", "scdblfinder", "consensus")
    ]
    assert {entry["mode"] for entry in aggregate_records} == {"aggregate"}
    assert {entry["effective"] for entry in aggregate_records} == {
        round(effective, 6)
    }
    assert {entry["requested"] for entry in aggregate_records} == {0.06}
    assert {entry["aggregation_rule"] for entry in aggregate_records} == {
        "cell_weighted_mean_of_resolved_per_sample"
    }


def test_cpu_grouped_scrublet_consumes_exact_resolved_per_sample_rates(tmp_path):
    adata = _adata(n_obs=80)
    adata.obs["sample"] = ["small"] * 20 + ["large"] * 60
    cfg = DoubletConfig(
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)
    DoubletDetectionModule._resolve_expected_rate_contract(ctx, cfg, adata)
    seen = {}

    class FakeScrublet:
        def __init__(self, matrix, expected_doublet_rate, random_state):
            self.n_obs = matrix.shape[0]
            seen[self.n_obs] = expected_doublet_rate
            self.threshold_ = 0.5

        def scrub_doublets(self, **_kwargs):
            return np.zeros(self.n_obs), np.zeros(self.n_obs, dtype=bool)

    DoubletDetectionModule()._run_cpu_scrublet_grouped(
        adata, FakeScrublet, cfg, ctx
    )

    assert set(seen) == {20, 60}
    np.testing.assert_allclose(
        [seen[20], seen[60]],
        [
            DoubletDetectionModule._resolved_expected_rate(cfg, 20),
            DoubletDetectionModule._resolved_expected_rate(cfg, 60),
        ],
    )
    assert ctx.metadata["doublet_expected_rate_by_backend"]["scrublet"][
        "mode"
    ] == "per_sample"


def test_unavailable_r_backend_still_records_attempted_effective_prior(
    monkeypatch, tmp_path
):
    adata = _adata(n_obs=80)
    cfg = DoubletConfig(
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)
    effective = DoubletDetectionModule._resolve_expected_rate_contract(
        ctx, cfg, adata
    )
    monkeypatch.setattr("shutil.which", lambda _name: None)

    with np.testing.assert_raises_regex(RuntimeError, "requires conda"):
        DoubletDetectionModule()._run_doubletfinder_via_r(adata, cfg, ctx)
    with np.testing.assert_raises_regex(RuntimeError, "requires conda"):
        DoubletDetectionModule()._run_scdblfinder_via_r(adata, cfg, ctx)

    assert ctx.metadata["doublet_expected_rate_by_backend"] == {
        "doubletfinder": {
            "requested": 0.06,
            "mode": "aggregate",
            "effective": effective,
            "aggregation_rule": "cell_weighted_mean_of_resolved_per_sample",
        },
        "scdblfinder": {
            "requested": 0.06,
            "mode": "aggregate",
            "effective": effective,
            "aggregation_rule": "cell_weighted_mean_of_resolved_per_sample",
        },
    }


def test_rank_consensus_uses_effective_not_requested_rate(monkeypatch, tmp_path):
    adata = _adata(n_obs=80)
    adata.obs["sample"] = ["small"] * 20 + ["large"] * 60
    cfg = DoubletConfig(
        backend="consensus",
        consensus_pair="scrublet_scdblfinder",
        consensus_logic="rank",
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)
    ctx.adata = adata
    scores = np.linspace(0.0, 1.0, adata.n_obs, dtype=np.float32)
    calls = np.zeros(adata.n_obs, dtype=bool)

    monkeypatch.setattr(
        DoubletDetectionModule,
        "_run_scrublet_only",
        lambda self, adata_arg, cfg_arg, ctx_arg: (scores, calls, 0.5),
    )
    monkeypatch.setattr(
        DoubletDetectionModule,
        "_run_scdblfinder_via_r",
        lambda self, adata_arg, cfg_arg, ctx_arg: (scores**2, calls, None),
    )

    DoubletDetectionModule().run(ctx)

    effective = ctx.metadata["doublet_expected_rate_effective"]
    expected_calls = max(1, int(round(effective * adata.n_obs)))
    assert ctx.adata.obs["predicted_doublet"].sum() == expected_calls
    assert expected_calls < int(round(cfg.expected_doublet_rate * adata.n_obs))


def test_undercall_uses_effective_prior_from_contract(tmp_path):
    adata = _adata(n_obs=80)
    adata.obs["sample"] = ["small"] * 20 + ["large"] * 60
    cfg = DoubletConfig(
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)
    effective = DoubletDetectionModule._resolve_expected_rate_contract(
        ctx, cfg, adata
    )

    DoubletDetectionModule._check_undercall(
        call_rate=effective * 0.25,
        expected_rate=effective,
        ctx=ctx,
        backend="scrublet",
    )

    assert ctx.metadata["doublet_undercall_expected_rate"] == round(effective, 6)
    assert ctx.metadata["doublet_undercall_ratio"] == 4.0

    sample_rate = ctx.metadata["doublet_expected_rate_resolved_per_sample"]["small"]
    DoubletDetectionModule._check_undercall(
        call_rate=sample_rate * 0.25,
        expected_rate=sample_rate,
        ctx=ctx,
        backend="scrublet",
        scope="small",
    )
    assert ctx.metadata["doublet_undercall_expected_rate_per_sample"] == {
        "small": sample_rate
    }
    assert ctx.metadata["doublet_undercall_ratio_per_sample"] == {"small": 4.0}


def test_grouped_finalize_runs_undercall_diagnostics_at_sample_scope(tmp_path):
    adata = _adata(n_obs=80)
    adata.obs["sample"] = ["small"] * 20 + ["large"] * 60
    cfg = DoubletConfig(
        expected_doublet_rate=0.06,
        scale_expected_doublet_rate=True,
        remove_doublets=False,
    )
    ctx = _ctx(tmp_path, cfg)
    module = DoubletDetectionModule()
    module._resolve_expected_rate_contract(ctx, cfg, adata)
    module._record_backend_expected_rate(ctx, "scrublet", mode="per_sample")

    module._finalize_doublets(
        ctx,
        adata,
        cfg,
        "scrublet",
        11,
        np.zeros(adata.n_obs, dtype=np.float32),
        np.zeros(adata.n_obs, dtype=bool),
        None,
    )

    assert ctx.metadata["doublet_undercall_expected_rate_per_sample"] == (
        ctx.metadata["doublet_expected_rate_resolved_per_sample"]
    )
    assert set(ctx.metadata["doublet_undercall_warning_per_sample"]) == {
        "small",
        "large",
    }
    assert "doublet_undercall_expected_rate" not in ctx.metadata
