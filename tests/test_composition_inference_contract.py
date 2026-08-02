"""Replicate-aware inference contract for the composition module.

Each test names a scientific break demonstrated by the 2026-08-02 review:
sample identifiers were reused as model covariates, so six biological samples
produced sample-level coefficients instead of the requested condition contrast.
"""
from __future__ import annotations

import sys
import subprocess
from types import SimpleNamespace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from workflow.modular.config import CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.composition import CompositionModule


def _adata(
    *,
    samples_by_condition: dict[str, tuple[str, ...]],
    cells_per_type: int = 3,
) -> ad.AnnData:
    rows: list[dict[str, str]] = []
    for condition, samples in samples_by_condition.items():
        for sample in samples:
            for cell_type in ("A", "B"):
                rows.extend(
                    {
                        "sample": sample,
                        "condition": condition,
                        "cell_type": cell_type,
                    }
                    for _ in range(cells_per_type)
                )
    obs = pd.DataFrame(rows, index=[f"cell_{i}" for i in range(len(rows))])
    return ad.AnnData(
        X=np.ones((len(obs), 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["dummy"]),
    )


def _ctx(tmp_path, adata: ad.AnnData, **composition) -> PipelineContext:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = PipelineConfig(
        project="composition-contract",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "sample", outs_dir=tmp_path / "outs"
        ),
        optional_modules=[],
    )
    # The tests intentionally attach the future config shape without importing a
    # not-yet-existing class. RED therefore exercises today's runtime behavior
    # instead of failing during test collection.
    cfg.composition = SimpleNamespace(
        sample_col=composition.get("sample_col", "sample"),
        condition_col=composition.get("condition_col", "condition"),
        contrast_a=composition.get("contrast_a", "CTRL"),
        contrast_b=composition.get("contrast_b", "KO"),
        covariates=composition.get("covariates", ()),
        min_samples_per_condition=composition.get("min_samples_per_condition", 3),
    )
    return PipelineContext(
        cfg=cfg,
        run_dir=run_dir,
        figure_dir=run_dir,
        table_dir=run_dir,
        adata=adata,
    )


def _silence_plots(monkeypatch) -> None:
    monkeypatch.setattr(
        CompositionModule, "_plot_barplot", staticmethod(lambda prop_df, ctx: None)
    )
    monkeypatch.setattr(
        CompositionModule, "_plot_boxplot", staticmethod(lambda prop_df, ctx: None)
    )


class _FakeSccoda:
    def __init__(self) -> None:
        self.load_kwargs: dict = {}
        self.formula = ""
        self.available_covariates: set[str] = set()
        self.rng_key: int | None = None

    def load(self, adata, **kwargs):
        self.load_kwargs = kwargs
        self.available_covariates = set(kwargs.get("covariate_obs") or ())
        return {"adata": adata}

    def prepare(self, data, *, formula, reference_cell_type):
        if "condition" not in self.available_covariates:
            raise ValueError("condition missing from generated sample observations")
        self.formula = formula

    def run_nuts(self, data, *, num_warmup, num_samples, rng_key=0):
        self.rng_key = rng_key
        return None

    def credible_effects(self, data):
        coefficient = f"{self.formula}[T.KO]"
        index = pd.MultiIndex.from_product(
            [[coefficient], ["A", "B"]], names=["Covariate", "Cell Type"]
        )
        return pd.Series([True, False], index=index, name="Final Parameter")


def test_three_vs_three_models_condition_not_sample(tmp_path, monkeypatch):
    """Changing the formula back to the sample ID must fail this test."""
    adata = _adata(
        samples_by_condition={
            "CTRL": ("C1", "C2", "C3"),
            "KO": ("K1", "K2", "K3"),
        }
    )
    ctx = _ctx(tmp_path, adata)
    ctx.cfg.random_state = 314
    fake = _FakeSccoda()
    monkeypatch.setitem(
        sys.modules,
        "pertpy",
        SimpleNamespace(tl=SimpleNamespace(Sccoda=lambda: fake)),
    )
    _silence_plots(monkeypatch)

    CompositionModule().run(ctx)

    assert fake.load_kwargs["sample_identifier"] == "sample"
    assert fake.load_kwargs["covariate_obs"] == ["condition"]
    assert fake.formula == "C(condition)"
    assert fake.rng_key == 314
    result = pd.read_csv(ctx.table_dir / "composition_test_results.csv")
    assert result["covariate"].str.contains("condition", regex=False).all()
    assert not result["covariate"].str.contains("C1|C2|C3|K1|K2|K3").any()
    assert ctx.metadata["composition_sample_col"] == "sample"
    assert ctx.metadata["composition_condition_col"] == "condition"
    assert ctx.metadata["composition_claimable"] is True


def test_real_pertpy_load_prepare_receives_sample_covariates() -> None:
    """Exercise the installed pertpy 1.x API without running expensive MCMC."""
    python = Path("/home/zerlinshen/conda/envs/sc_gpu/bin/python")
    if not python.is_file():
        pytest.skip("sc_gpu Python is unavailable")
    probe = """
import anndata as ad
import numpy as np
import pandas as pd
import pertpy as pt

samples = ["C1", "C2", "C3", "K1", "K2", "K3"]
conditions = ["CTRL"] * 3 + ["KO"] * 3
obs = pd.DataFrame({
    "sample": np.repeat(samples, 4),
    "condition": np.repeat(conditions, 4),
    "cell_type": ["A", "A", "B", "B"] * 6,
})
adata = ad.AnnData(X=np.ones((len(obs), 1)), obs=obs)
sccoda = pt.tl.Sccoda()
data = sccoda.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_type",
    sample_identifier="sample",
    covariate_obs=["condition"],
)
sccoda.prepare(data, formula="C(condition)", reference_cell_type="automatic")
assert "condition" in data["coda"].obs.columns
"""
    completed = subprocess.run(
        [str(python), "-c", probe],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout


def test_no_explicit_condition_is_descriptive_and_never_calls_sccoda(
    tmp_path, monkeypatch
):
    """A multi-level sample column is not an experimental condition."""
    adata = _adata(
        samples_by_condition={"unspecified": ("S1", "S2", "S3", "S4")}
    )
    ctx = _ctx(
        tmp_path,
        adata,
        condition_col=None,
        contrast_a=None,
        contrast_b=None,
        min_samples_per_condition=2,
    )

    def fail_if_called(*args, **kwargs):
        raise AssertionError("scCODA must not run without an explicit condition contract")

    monkeypatch.setattr(
        CompositionModule, "_try_pertpy", staticmethod(fail_if_called)
    )
    _silence_plots(monkeypatch)

    CompositionModule().run(ctx)

    assert ctx.metadata["composition_engine"] == "descriptive_only"
    assert (
        ctx.metadata["composition_inference_status"]
        == "descriptive_only_no_condition_contract"
    )
    assert ctx.metadata["composition_claimable"] is False
    result = pd.read_csv(ctx.table_dir / "composition_test_results.csv")
    assert result.empty


def test_sample_mapped_to_two_conditions_fails_before_inference(tmp_path, monkeypatch):
    """One biological replicate cannot contribute to both contrast arms."""
    adata = _adata(
        samples_by_condition={
            "CTRL": ("shared", "C2", "C3"),
            "KO": ("shared", "K2", "K3"),
        }
    )
    ctx = _ctx(tmp_path, adata)
    _silence_plots(monkeypatch)

    with pytest.raises(ValueError, match="one condition"):
        CompositionModule().run(ctx)


def test_under_replicated_condition_fails_before_inference(tmp_path, monkeypatch):
    """Dropping the biological-replicate floor must fail this test."""
    adata = _adata(samples_by_condition={"CTRL": ("C1",), "KO": ("K1",)})
    ctx = _ctx(tmp_path, adata, min_samples_per_condition=3)
    _silence_plots(monkeypatch)

    with pytest.raises(ValueError, match="biological samples"):
        CompositionModule().run(ctx)


def test_sample_identifier_cannot_be_reintroduced_as_a_covariate(
    tmp_path, monkeypatch
):
    """A future covariate option must not recreate the original C(sample) bug."""
    adata = _adata(
        samples_by_condition={
            "CTRL": ("C1", "C2", "C3"),
            "KO": ("K1", "K2", "K3"),
        }
    )
    ctx = _ctx(tmp_path, adata, covariates=("sample",))
    _silence_plots(monkeypatch)

    with pytest.raises(ValueError, match="sample identifier.*covariate"):
        CompositionModule().run(ctx)


def test_lusc_shaped_sample_metadata_stays_descriptive_without_condition(
    tmp_path, monkeypatch
):
    """The retained LUSC object has 87 sample levels; that is not 87 conditions."""
    adata = _adata(
        samples_by_condition={
            "not_a_contrast": tuple(f"LUSC_{i:03d}" for i in range(87))
        },
        cells_per_type=1,
    )
    adata.obs["dataset"] = "lusc_squamous"
    ctx = _ctx(
        tmp_path,
        adata,
        condition_col=None,
        contrast_a=None,
        contrast_b=None,
        min_samples_per_condition=2,
    )

    def fail_if_called(*args, **kwargs):
        raise AssertionError("87 sample IDs must not become scCODA coefficients")

    monkeypatch.setattr(
        CompositionModule, "_try_pertpy", staticmethod(fail_if_called)
    )
    _silence_plots(monkeypatch)

    CompositionModule().run(ctx)

    counts = pd.read_csv(
        ctx.table_dir / "composition_counts.csv", index_col=0
    )
    assert len(counts) == 87
    assert ctx.metadata["composition_sample_col"] == "sample"
    assert ctx.metadata["composition_condition_col"] is None
    assert ctx.metadata["composition_claimable"] is False


def test_cli_keeps_sample_and_condition_arguments_separate(monkeypatch):
    """Removing either CLI field must not collapse the design back to batch_key."""
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp",
            "--composition-sample-col",
            "sample",
            "--composition-condition-col",
            "condition",
            "--composition-contrast-a",
            "CTRL",
            "--composition-contrast-b",
            "KO",
            "--composition-covariates",
            "donor_sex,chemistry",
        ],
    )
    args = cli_mod.parse_args()
    modules = cli_mod._resolve_optional_modules(args)

    assert args.composition_sample_col == "sample"
    assert args.composition_condition_col == "condition"
    assert args.composition_contrast_a == "CTRL"
    assert args.composition_contrast_b == "KO"
    assert args.composition_covariates == "donor_sex,chemistry"
    assert "composition" in modules


def test_cli_rejects_half_a_composition_contrast(monkeypatch):
    import workflow.modular.cli as cli_mod

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            "/tmp",
            "--composition-condition-col",
            "condition",
            "--composition-contrast-a",
            "CTRL",
        ],
    )
    args = cli_mod.parse_args()
    with pytest.raises(SystemExit, match="requires both"):
        cli_mod._resolve_optional_modules(args)


def test_cli_wires_composition_design_into_pipeline_config(tmp_path, monkeypatch):
    """Parser-only coverage cannot catch dropping the fields during cfg assembly."""
    from workflow.modular import cli as cli_mod
    from workflow.modular import manifest_writer

    captured: dict = {}
    sample_root = tmp_path / "sample"
    project_root = tmp_path / "project"
    sample_root.mkdir()

    monkeypatch.setattr(
        manifest_writer,
        "factory_git_state",
        lambda root: {"sha": "aaaaaaa", "dirty": False, "branch": "test"},
    )

    def fake_run_pipeline(cfg, ledger=None):
        captured["cfg"] = cfg
        return {"modules_run": ["composition"], "bundle_sha256": "abc123"}

    monkeypatch.setattr(cli_mod, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prog",
            "--project",
            "x",
            "--sample-root",
            str(sample_root),
            "--project-root",
            str(project_root),
            "--run-id",
            "2026-08-02T0000Z-aaaaaaa",
            "--allow-dirty",
            "--composition-sample-col",
            "sample",
            "--composition-condition-col",
            "condition",
            "--composition-contrast-a",
            "CTRL",
            "--composition-contrast-b",
            "KO",
            "--composition-covariates",
            "donor_sex,chemistry",
            "--composition-min-samples-per-condition",
            "3",
        ],
    )
    monkeypatch.delenv("SC_REQUIRE_PROJECT_ROOT", raising=False)

    cli_mod.main()

    design = captured["cfg"].composition
    assert design.sample_col == "sample"
    assert design.condition_col == "condition"
    assert (design.contrast_a, design.contrast_b) == ("CTRL", "KO")
    assert design.covariates == ("donor_sex", "chemistry")
    assert design.min_samples_per_condition == 3
    assert "composition" in captured["cfg"].optional_modules
