"""Tests for scripts/scfactory.py — the user-friendly CLI wrapper.

These tests exercise the auto-detection planner and the doctor health check
without invoking the underlying pipeline. We import scfactory as a module
(via importlib from a file path, since scripts/ is not a package) and call
``main()`` with a synthesized argv to keep tests fast.
"""

from __future__ import annotations

import importlib.util
import json
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCFACTORY_PATH = ROOT / "scripts" / "scfactory.py"


def _load_scfactory():
    spec = importlib.util.spec_from_file_location("scfactory", SCFACTORY_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["scfactory"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def scfactory():
    return _load_scfactory()


def _make_minimal_adata(tmp_path: Path, *, with_protein=False, with_spatial=False) -> Path:
    """Synthesize a tiny .h5ad with the requested obsm slots."""
    anndata = pytest.importorskip("anndata")
    n_obs, n_vars = 12, 8
    X = np.zeros((n_obs, n_vars), dtype=np.float32)
    adata = anndata.AnnData(X=X)
    if with_protein:
        adata.obsm["protein_counts"] = np.zeros((n_obs, 3), dtype=np.float32)
    if with_spatial:
        adata.obsm["spatial"] = np.zeros((n_obs, 2), dtype=np.float32)
    out = tmp_path / "tiny.h5ad"
    adata.write_h5ad(out)
    return out


def _run_scfactory_subprocess(*args: str) -> subprocess.CompletedProcess[str]:
    """Run the public adapter exactly as a user would from the repository root."""
    return subprocess.run(
        [sys.executable, str(SCFACTORY_PATH), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def _planned_cli_tokens(stdout: str) -> list[str]:
    line = next(
        line for line in stdout.splitlines()
        if line.startswith("scfactory: would execute: ")
    )
    return shlex.split(line.partition("scfactory: would execute: ")[2])


# ---------------------------------------------------------------------------
# run --dry-run modality detection
# ---------------------------------------------------------------------------

def test_run_detects_rna_only_h5ad(tmp_path, capsys, scfactory):
    h5ad = _make_minimal_adata(tmp_path)
    rc = scfactory.main(["run", str(h5ad), "--dry-run"])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "modality=rna_only" in out
    assert "protein_adt" not in out
    assert "spatial_ingest" not in out
    # baseline RNA modules present
    assert "clustering" in out
    assert "differential_expression" in out


def test_run_detects_protein_h5ad(tmp_path, capsys, scfactory):
    h5ad = _make_minimal_adata(tmp_path, with_protein=True)
    rc = scfactory.main(["run", str(h5ad), "--dry-run"])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "modality=cite_seq" in out
    assert "protein_adt" in out
    assert "spatial_ingest" not in out


def test_run_detects_spatial_h5ad(tmp_path, monkeypatch, capsys, scfactory):
    # Wave-2 preflight requires squidpy whenever spatial_neighborhoods is
    # planned; this host may lack it. Stub package check so this test only
    # covers modality detection + planning (not install state).
    monkeypatch.setattr(scfactory, "_missing_python_packages", lambda names: [])
    h5ad = _make_minimal_adata(tmp_path, with_spatial=True)
    rc = scfactory.main(["run", str(h5ad), "--dry-run"])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "modality=spatial" in out
    assert "spatial_ingest" in out
    assert "spatial_neighborhoods" in out


def test_run_detection_failure_exits_2(tmp_path, capsys, scfactory):
    nonexistent = tmp_path / "does_not_exist.h5ad"
    rc = scfactory.main(["run", str(nonexistent), "--dry-run"])
    err = capsys.readouterr().err
    assert rc == 2
    assert "input not found" in err or "could not auto-detect" in err


def test_subprocess_forwards_arbitrary_named_h5ad_to_canonical_cli(tmp_path):
    h5ad = _make_minimal_adata(tmp_path)
    arbitrary = h5ad.with_name("patient-42.expression-matrix.h5ad")
    h5ad.rename(arbitrary)

    proc = _run_scfactory_subprocess("run", str(arbitrary), "--dry-run")

    assert proc.returncode == 0, proc.stderr
    tokens = _planned_cli_tokens(proc.stdout)
    assert "--input-h5ad" in tokens
    assert tokens[tokens.index("--input-h5ad") + 1] == str(arbitrary.resolve())
    assert "--sample-root" not in tokens


def test_subprocess_forwards_governed_project_root_and_run_id(tmp_path):
    h5ad = _make_minimal_adata(tmp_path)
    project_root = tmp_path / "governed-project"
    run_id = "2026-08-02T1200Z-abcdef0"

    proc = _run_scfactory_subprocess(
        "run",
        str(h5ad),
        "--project-root",
        str(project_root),
        "--run-id",
        run_id,
        "--dry-run",
    )

    assert proc.returncode == 0, proc.stderr
    tokens = _planned_cli_tokens(proc.stdout)
    assert tokens[tokens.index("--project-root") + 1] == str(project_root.resolve())
    assert tokens[tokens.index("--run-id") + 1] == run_id
    assert "--output-dir" not in tokens


def test_subprocess_forwards_scientific_non_equivalence_acknowledgement(tmp_path):
    h5ad = _make_minimal_adata(tmp_path)

    proc = _run_scfactory_subprocess(
        "run",
        str(h5ad),
        "--scientific-profile",
        "paper-15pc",
        "--acknowledge-scientific-non-equivalence",
        "--dry-run",
    )

    assert proc.returncode == 0, proc.stderr
    tokens = _planned_cli_tokens(proc.stdout)
    assert tokens[tokens.index("--scientific-profile") + 1] == "paper-15pc"
    assert "--acknowledge-scientific-non-equivalence" in tokens


def test_governed_bundle_reuses_run_id_and_project_layout(
    tmp_path, capsys, monkeypatch, scfactory
):
    h5ad = _make_minimal_adata(tmp_path)
    project_root = tmp_path / "governed-project"
    run_id = "2026-08-02T1200Z-abcdef0"
    final_h5ad = (
        project_root
        / "runs"
        / run_id
        / "python"
        / "adapter-test_20260802_120000"
        / "final_adata.h5ad"
    )
    calls: list[list[str]] = []

    def _fake_run(cmd, **kwargs):
        calls.append(list(cmd))
        if cmd[1:3] == ["-m", "workflow.modular.cli"]:
            final_h5ad.parent.mkdir(parents=True)
            final_h5ad.write_bytes(b"pipeline-output-placeholder")

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr(scfactory.subprocess, "run", _fake_run)

    rc = scfactory.main(
        [
            "run",
            str(h5ad),
            "--project",
            "adapter-test",
            "--project-root",
            str(project_root),
            "--run-id",
            run_id,
            "--bundle",
        ]
    )

    assert rc == 0, capsys.readouterr().err
    assert len(calls) == 2
    export_cmd = calls[1]
    assert export_cmd[export_cmd.index("--input") + 1] == str(final_h5ad)
    assert export_cmd[export_cmd.index("--project-root") + 1] == str(
        project_root.resolve()
    )
    assert export_cmd[export_cmd.index("--run-id") + 1] == run_id
    assert "--output" not in export_cmd


def test_subprocess_rna_plan_uses_canonical_catalog_defaults(tmp_path):
    from workflow.modular.module_catalog import DEFAULT_OPTIONAL_MODULES

    h5ad = _make_minimal_adata(tmp_path)
    proc = _run_scfactory_subprocess("run", str(h5ad), "--dry-run")

    assert proc.returncode == 0, proc.stderr
    tokens = _planned_cli_tokens(proc.stdout)
    actual = tokens[tokens.index("--optional-modules") + 1].split(",")
    assert actual == list(DEFAULT_OPTIONAL_MODULES)


def test_subprocess_dry_run_forwards_complete_explicit_scatac_contract(tmp_path):
    h5ad = _make_minimal_adata(tmp_path)
    proc = _run_scfactory_subprocess(
        "run", str(h5ad),
        "--optional-modules", "scatac_pseudobulk_da",
        "--scatac-da-sample-col", "sample",
        "--scatac-da-group-col", "cell_type",
        "--scatac-da-condition-col", "condition",
        "--scatac-da-peak-id-col", "peak_id",
        "--scatac-da-test-level", "STIM",
        "--scatac-da-reference-level", "CTRL",
        "--scatac-da-min-samples-per-condition", "2",
        "--scatac-da-min-total-count", "10",
        "--scatac-da-fdr-threshold", "0.05",
        "--scatac-da-abs-log2fc-threshold", "1.0",
        "--scatac-da-r-conda-env", "r_multiomics",
        "--scatac-da-timeout", "1800",
        "--scatac-da-aggregation-backend", "cpu",
        "--dry-run",
    )

    assert proc.returncode == 0, proc.stderr
    tokens = _planned_cli_tokens(proc.stdout)
    assert tokens[tokens.index("--optional-modules") + 1] == "scatac_pseudobulk_da"
    expected = {
        "--scatac-da-sample-col": "sample",
        "--scatac-da-group-col": "cell_type",
        "--scatac-da-condition-col": "condition",
        "--scatac-da-peak-id-col": "peak_id",
        "--scatac-da-test-level": "STIM",
        "--scatac-da-reference-level": "CTRL",
        "--scatac-da-min-samples-per-condition": "2",
        "--scatac-da-min-total-count": "10",
        "--scatac-da-fdr-threshold": "0.05",
        "--scatac-da-abs-log2fc-threshold": "1.0",
        "--scatac-da-r-conda-env": "r_multiomics",
        "--scatac-da-timeout": "1800",
        "--scatac-da-aggregation-backend": "cpu",
    }
    for option, value in expected.items():
        assert tokens[tokens.index(option) + 1] == value


def test_public_dry_run_rejects_incomplete_scatac_contract(tmp_path):
    h5ad = _make_minimal_adata(tmp_path)
    proc = _run_scfactory_subprocess(
        "run", str(h5ad), "--scatac-da-sample-col", "sample", "--dry-run"
    )

    assert proc.returncode == 2
    assert "requires explicit --scatac-da-group-col" in proc.stderr
    assert "would execute" not in proc.stdout


def test_dry_run_rejects_unknown_user_module_before_printing_a_plan(tmp_path):
    h5ad = _make_minimal_adata(tmp_path)
    proc = _run_scfactory_subprocess(
        "run",
        str(h5ad),
        "--optional-modules",
        "clustering,definitely_not_a_module",
        "--dry-run",
    )

    assert proc.returncode == 2
    assert "unknown module name" in proc.stderr
    assert "definitely_not_a_module" in proc.stderr
    assert "scfactory: would execute:" not in proc.stdout


def test_canonical_cli_accepts_input_h5ad_argument(tmp_path, monkeypatch):
    from workflow.modular import cli

    h5ad = _make_minimal_adata(tmp_path)
    arbitrary = h5ad.with_name("arbitrary-name.h5ad")
    h5ad.rename(arbitrary)
    monkeypatch.setattr(
        sys,
        "argv",
        ["cli", "--project", "adapter-test", "--input-h5ad", str(arbitrary)],
    )

    args = cli.parse_args()

    assert args.input_h5ad == str(arbitrary)
    assert args.sample_root is None


def test_canonical_cli_plumbs_exact_h5ad_into_pipeline_config(tmp_path, monkeypatch):
    from workflow.modular import cli

    h5ad = _make_minimal_adata(tmp_path)
    arbitrary = h5ad.with_name("pipeline-config-source.h5ad")
    h5ad.rename(arbitrary)
    captured: dict[str, object] = {}

    def _fake_run_pipeline(cfg, ledger=None):
        captured["cfg"] = cfg
        return {"modules_run": [], "failed_modules": []}

    monkeypatch.setattr(cli, "run_pipeline", _fake_run_pipeline)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "cli",
            "--project",
            "adapter-test",
            "--input-h5ad",
            str(arbitrary),
            "--output-dir",
            str(tmp_path / "legacy-output"),
        ],
    )
    monkeypatch.delenv("SC_REQUIRE_PROJECT_ROOT", raising=False)
    monkeypatch.setenv(
        "SC_LEGACY_OUTPUT_ACCESS_LOG", str(tmp_path / "legacy-access.jsonl")
    )

    cli.main()

    cfg = captured["cfg"]
    assert cfg.cellranger.input_h5ad == arbitrary.resolve()
    assert cfg.cellranger.sample_root == arbitrary.parent


def test_cellranger_loads_the_exact_direct_h5ad(tmp_path):
    from workflow.modular.config import CellRangerConfig, PipelineConfig
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.cellranger import CellRangerModule

    h5ad = _make_minimal_adata(tmp_path)
    arbitrary = h5ad.with_name("not-prepared-input.h5ad")
    h5ad.rename(arbitrary)
    cfg = PipelineConfig(
        project="adapter-test",
        output_dir=tmp_path / "output",
        cellranger=CellRangerConfig(
            sample_root=arbitrary.parent,
            outs_dir=arbitrary.parent / "unused-outs",
            input_h5ad=arbitrary,
        ),
        optional_modules=[],
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
    )

    CellRangerModule().run(ctx)

    assert ctx.adata is not None
    assert ctx.adata.shape == (12, 8)
    assert ctx.metadata["prepared_input_source"] == str(arbitrary)
    assert ctx.metadata["prepared_input_h5ad"] == str(arbitrary)


def test_direct_normalized_h5ad_never_synthesizes_a_counts_layer(tmp_path):
    """A requested pseudobulk module must not turn arbitrary X into raw counts."""
    from workflow.modular.config import CellRangerConfig, PipelineConfig
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.cellranger import CellRangerModule

    h5ad = _make_minimal_adata(tmp_path)
    cfg = PipelineConfig(
        project="adapter-test",
        output_dir=tmp_path / "output",
        cellranger=CellRangerConfig(
            sample_root=h5ad.parent,
            outs_dir=h5ad.parent / "unused-outs",
            input_h5ad=h5ad,
        ),
        optional_modules=["pseudobulk_de"],
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
    )

    CellRangerModule().run(ctx)

    assert "counts" not in ctx.adata.layers
    assert ctx.metadata["counts_layer_preserved"] is False
    assert ctx.metadata["counts_provenance_status"] == "missing_explicit_counts_layer"


def test_direct_h5ad_rejects_fractional_declared_counts_at_ingest(tmp_path):
    """A layer called counts is not sufficient evidence when its values are normalized."""
    from workflow.modular.config import CellRangerConfig, PipelineConfig
    from workflow.modular.context import PipelineContext
    from workflow.modular.modules.cellranger import CellRangerModule

    anndata = pytest.importorskip("anndata")
    h5ad = tmp_path / "fractional-counts.h5ad"
    adata = anndata.AnnData(X=np.ones((12, 8), dtype=np.float32))
    adata.layers["counts"] = np.full((12, 8), 0.25, dtype=np.float32)
    adata.uns["counts_provenance"] = {
        "schema_version": "1.0",
        "matrix": "layers/counts",
        "semantic": "raw_umi_counts",
        "source": "author_supplied_raw_umi_layer",
    }
    adata.write_h5ad(h5ad)
    cfg = PipelineConfig(
        project="adapter-test",
        output_dir=tmp_path / "output",
        cellranger=CellRangerConfig(
            sample_root=h5ad.parent,
            outs_dir=h5ad.parent / "unused-outs",
            input_h5ad=h5ad,
        ),
        optional_modules=["pseudobulk_de"],
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
    )

    with pytest.raises(ValueError, match="nonnegative integer raw UMI counts"):
        CellRangerModule().run(ctx)


# ---------------------------------------------------------------------------
# doctor
# ---------------------------------------------------------------------------

def _stub_claim_critical_probe(monkeypatch, scfactory) -> None:
    """Stub out the suite-root claim-critical contract lookup for doctor
    tests that are not specifically exercising claim-critical behavior.

    The real check shells out to the `sc_gpu` and `r_multiomics` conda envs
    (see `test_doctor_claim_critical_*` below, which intentionally exercise
    that real subprocess path). Every OTHER `doctor` test just needs a fast,
    deterministic no-op here so the suite doesn't pay that real probe cost
    (observed ~12s/call) once per unrelated test.
    """

    def fake_load(*args, **kwargs):
        return None, "stubbed for test speed (see _stub_claim_critical_probe)"

    monkeypatch.setattr(scfactory, "_load_claim_critical_contract", fake_load)


def test_doctor_emits_pass_or_warn(monkeypatch, capsys, scfactory):
    _stub_claim_critical_probe(monkeypatch, scfactory)
    rc = scfactory.main(["doctor"])
    out = capsys.readouterr().out
    # Lines we expect regardless of pass/warn outcome.
    assert "python:" in out
    assert "rscript:" in out
    assert "bridge_symlinks:" in out
    assert "summary:" in out
    # Allow either 0 (no fails) or 1 (fail). In this dev env we expect 0.
    assert rc in (0, 1)


def test_doctor_json_shape(monkeypatch, capsys, scfactory):
    _stub_claim_critical_probe(monkeypatch, scfactory)
    rc = scfactory.main(["doctor", "--json"])
    out = capsys.readouterr().out
    payload = json.loads(out)
    for key in ("python", "conda_envs", "rscript", "bridge_symlinks",
                "python_deps", "r_packages", "last_run", "summary",
                "lock_identity_match", "claim_critical_deps", "readiness"):
        assert key in payload, f"missing top-level key: {key}"
    summary = payload["summary"]
    for k in ("pass", "warn", "fail"):
        assert k in summary
        assert isinstance(summary[k], int)
    for key in (
        "core_environment_ready",
        "selected_profile_ready",
        "optional_capabilities",
        "spatial_analytics",
        "claim_critical_ready",
    ):
        assert key in payload["readiness"], f"missing readiness dimension: {key}"
        assert "status" in payload["readiness"][key]
    assert rc in (0, 1)


def test_doctor_searches_governed_project_manifest(tmp_path, monkeypatch, capsys, scfactory):
    _stub_claim_critical_probe(monkeypatch, scfactory)
    project_root = tmp_path / "project"
    run_id = "2026-08-02T1200Z-abcdef0"
    manifest = project_root / "runs" / run_id / "manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(
        json.dumps(
            {
                "overall_status": "complete",
                "completed_modules": ["cellranger", "qc"],
            }
        ),
        encoding="utf-8",
    )

    rc = scfactory.main(
        ["doctor", "--json", "--project-root", str(project_root)]
    )
    payload = json.loads(capsys.readouterr().out)

    assert rc in (0, 1)
    assert payload["last_run"]["found"] is True
    assert payload["last_run"]["info"]["path"] == str(manifest)


def test_doctor_lock_identity_mismatch_is_non_pass(monkeypatch, capsys, scfactory):
    """A running interpreter that diverges from the declared environment.yml
    pin must surface as a non-pass, even though it clears the >=3.10 floor
    that `_check_python` alone would report as PASS."""
    _stub_claim_critical_probe(monkeypatch, scfactory)

    def fake_lock(*args, **kwargs):
        return {
            "path": "fake/environment.yml",
            "found": True,
            # Guaranteed to differ from whatever interpreter runs this test.
            "python_version": (2, 7),
            "declared_packages": {"scanpy", "anndata", "scipy"},
        }

    monkeypatch.setattr(scfactory, "_read_environment_lock", fake_lock)

    rc = scfactory.main(["doctor", "--json"])
    payload = json.loads(capsys.readouterr().out)

    lock = payload["lock_identity_match"]
    assert lock["status"] == "warn"
    assert lock["python"]["match"] is False
    assert "MISMATCH" in lock["message"]
    # A version mismatch is a warning, never an automatic suite-level fail.
    assert payload["summary"]["fail"] == 0
    assert rc in (0, 1)


def test_doctor_missing_required_python_dep_is_non_pass(monkeypatch, capsys, scfactory):
    """A dep declared in environment.yml (tier=required) that fails to import
    must surface as non-pass at both the per-dep level and the rolled-up
    `core_environment_ready` readiness dimension — not collapse into the
    flat "pass" that the old hardcoded status produced."""
    _stub_claim_critical_probe(monkeypatch, scfactory)
    import importlib as _importlib

    def fake_lock(*args, **kwargs):
        return {
            "path": "fake/environment.yml",
            "found": True,
            "python_version": (sys.version_info.major, sys.version_info.minor),
            # Mark mofapy2 as lock-required for this test, regardless of the
            # real environment.yml (which does not declare it).
            "declared_packages": {"scanpy", "anndata", "scipy", "mofapy2"},
        }

    monkeypatch.setattr(scfactory, "_read_environment_lock", fake_lock)

    real_import_module = _importlib.import_module

    def fake_import_module(name, *args, **kwargs):
        if name == "mofapy2":
            raise ImportError("simulated: mofapy2 not installed")
        return real_import_module(name, *args, **kwargs)

    monkeypatch.setattr(_importlib, "import_module", fake_import_module)

    rc = scfactory.main(["doctor", "--json"])
    payload = json.loads(capsys.readouterr().out)

    dep = payload["python_deps"]["deps"]["mofapy2"]
    assert dep["tier"] == "required"
    assert dep["present"] is False
    assert dep["status"] == "warn"
    assert "mofapy2" in payload["python_deps"]["missing_required"]
    assert payload["python_deps"]["status"] == "warn"
    assert payload["readiness"]["core_environment_ready"]["status"] == "warn"
    # Still never a hard suite-level fail — doctor stays non-aggressive.
    assert payload["summary"]["fail"] == 0
    assert rc in (0, 1)


def test_doctor_missing_optional_dep_does_not_cause_global_fail(monkeypatch, capsys, scfactory):
    """A dep NOT declared in environment.yml (tier=optional, e.g. squidpy)
    that fails to import must be individually visible as a warning but must
    never push the suite-level status to fail — optional capabilities are
    legitimately absent on many hosts."""
    _stub_claim_critical_probe(monkeypatch, scfactory)
    import importlib as _importlib

    real_import_module = _importlib.import_module

    def fake_import_module(name, *args, **kwargs):
        if name == "squidpy":
            raise ImportError("simulated: squidpy not installed")
        return real_import_module(name, *args, **kwargs)

    monkeypatch.setattr(_importlib, "import_module", fake_import_module)

    rc = scfactory.main(["doctor", "--json"])
    payload = json.loads(capsys.readouterr().out)

    dep = payload["python_deps"]["deps"]["squidpy"]
    assert dep["tier"] == "optional"
    assert dep["present"] is False
    assert dep["status"] == "warn"
    assert "squidpy" in payload["readiness"]["optional_capabilities"]["missing"]
    assert payload["readiness"]["optional_capabilities"]["status"] == "warn"
    # Wave-2: explicit spatial readiness surface.
    assert payload["readiness"]["spatial_analytics"]["present"] is False
    assert payload["readiness"]["spatial_analytics"]["status"] == "warn"
    assert "squidpy" in payload["readiness"]["spatial_analytics"]["message"]
    # core readiness is untouched by an optional-tier absence.
    assert payload["readiness"]["core_environment_ready"]["status"] == "pass"
    assert payload["summary"]["fail"] == 0
    assert rc == 0


def test_visium_recipe_preflight_fails_without_squidpy(
    tmp_path, monkeypatch, capsys, scfactory
):
    """Wave-2 W2.1: visium_neighborhoods must fail before pipeline start."""
    import importlib as _importlib

    real_import_module = _importlib.import_module

    def fake_import_module(name, *args, **kwargs):
        if name == "squidpy":
            raise ImportError("simulated: squidpy not installed")
        return real_import_module(name, *args, **kwargs)

    monkeypatch.setattr(_importlib, "import_module", fake_import_module)
    h5ad = _make_minimal_adata(tmp_path)
    rc = scfactory.main(
        [
            "run",
            str(h5ad),
            "--recipe",
            "visium_neighborhoods",
            "--dry-run",
        ]
    )
    err = capsys.readouterr().err
    assert rc == 2, err
    assert "squidpy" in err
    assert "missing required Python package" in err


def test_doctor_claim_critical_contract_unavailable_is_explicit(monkeypatch, capsys, scfactory):
    """When the suite-root claim-critical contract cannot be read (e.g. this
    repo cloned standalone, outside the suite), doctor must say so explicitly
    rather than silently reporting healthy."""

    def fake_load(*args, **kwargs):
        return None, "simulated: contract not found"

    monkeypatch.setattr(scfactory, "_load_claim_critical_contract", fake_load)

    rc = scfactory.main(["doctor", "--json"])
    payload = json.loads(capsys.readouterr().out)

    cc = payload["claim_critical_deps"]
    assert cc["contract_available"] is False
    assert cc["status"] == "warn"
    assert "simulated: contract not found" in cc["message"]
    assert payload["readiness"]["claim_critical_ready"]["contract_available"] is False
    # Contract unavailability is a "cannot verify" state, not a confirmed gap.
    assert payload["summary"]["fail"] == 0
    assert rc in (0, 1)


def test_doctor_missing_claim_critical_dep_is_fail_and_not_maskable(monkeypatch, capsys, scfactory):
    """A declared claim-critical dependency that fails to import in its
    declared env must escalate all the way to "fail" — the one tier in
    doctor allowed to do that — and that fail must survive to the overall
    exit code even while unrelated optional deps are also (harmlessly)
    missing on this host. This exercises the REAL cross-env subprocess probe
    (not a mocked import), using a genuinely nonexistent module name in the
    real `sc_gpu` env so the probe pathway itself is under test."""

    def fake_load(*args, **kwargs):
        return [{
            "name": "definitely_not_a_real_package_xyz_123",
            "import_name": "definitely_not_a_real_package_xyz_123",
            "env": "sc_gpu",
            "version": "0.0.0",
            "declared_in": "singlecell_factory/environment_gpu.yml",
            "guards": "singlecell_factory/workflow/modular/modules/pseudobulk_de.py",
            "claim_effect": "test fixture",
            "without_it": "simulated claim-path loss for test coverage",
        }], None

    monkeypatch.setattr(scfactory, "_load_claim_critical_contract", fake_load)

    rc = scfactory.main(["doctor", "--json"])
    payload = json.loads(capsys.readouterr().out)

    cc = payload["claim_critical_deps"]
    dep = cc["deps"]["definitely_not_a_real_package_xyz_123"]
    assert dep["verifiable"] is True  # sc_gpu env is genuinely present on this host
    assert dep["present"] is False
    assert dep["status"] == "fail"
    assert cc["status"] == "fail"
    assert payload["readiness"]["claim_critical_ready"]["status"] == "fail"
    # The fail must reach the top-level tally and exit code, unmasked by any
    # coexisting optional-tier warnings (this host has real optional gaps,
    # e.g. squidpy/fastparquet, which must stay warn-only and coexist with
    # this fail rather than diluting it).
    assert payload["summary"]["fail"] >= 1
    assert rc == 1


def test_doctor_claim_critical_env_absent_is_unverifiable_not_fail(monkeypatch, capsys, scfactory):
    """A claim-critical dependency declared in a conda env that does not
    exist on this host is honestly "unverifiable from here" (warn), not a
    confirmed failure — doctor cannot claim to have tested something it
    could not even launch an interpreter for."""

    def fake_load(*args, **kwargs):
        return [{
            "name": "some_claim_critical_thing",
            "import_name": "some_claim_critical_thing",
            "env": "definitely_not_a_real_conda_env_xyz",
            "version": "1.0.0",
            "declared_in": "singlecell_factory/environment_gpu.yml",
            "guards": "singlecell_factory/workflow/modular/modules/pseudobulk_de.py",
            "claim_effect": "test fixture",
            "without_it": "simulated",
        }], None

    monkeypatch.setattr(scfactory, "_load_claim_critical_contract", fake_load)

    rc = scfactory.main(["doctor", "--json"])
    payload = json.loads(capsys.readouterr().out)

    dep = payload["claim_critical_deps"]["deps"]["some_claim_critical_thing"]
    assert dep["verifiable"] is False
    assert dep["status"] == "warn"
    assert payload["claim_critical_deps"]["status"] == "warn"
    assert payload["summary"]["fail"] == 0
    assert rc in (0, 1)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

# A tiny but valid 1x1 transparent PNG (67 bytes), decoded from a constant.
_TINY_PNG_B64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNkYAAAAAYAAjCB"
    "0C8AAAAASUVORK5CYII="
)


def _write_tiny_png(path: Path) -> None:
    import base64 as _b64
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_b64.b64decode(_TINY_PNG_B64))


def _make_run_dir(
    tmp_path: Path,
    *,
    project: str = "DEMO_RUN",
    claim_guard: str | None = "MUST_NOT_TRUST_AS_CANONICAL",
    modules: tuple[str, ...] = ("clustering", "differential_expression"),
    add_pdf: bool = False,
    add_oversize_png: bool = False,
) -> Path:
    """Synthesize a run-like directory layout for report tests."""
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    manifest: dict = {
        "project": project,
        "generated_at": "2026-01-01T00:00:00Z",
        "module_status": [
            {"module": m, "status": "ok", "message": "completed"} for m in modules
        ],
    }
    if claim_guard is not None:
        manifest["expression"] = {"claim_guard": claim_guard}
    (run_dir / "run_manifest.json").write_text(json.dumps(manifest))

    # module_status.csv
    csv_lines = ["module,status,message"]
    for m in modules:
        csv_lines.append(f"{m},ok,completed")
    (run_dir / "module_status.csv").write_text("\n".join(csv_lines) + "\n")

    # PNGs per module
    for m in modules:
        _write_tiny_png(run_dir / m / f"{m}_fig1.png")
        _write_tiny_png(run_dir / m / f"{m}_fig2.png")

    if add_pdf:
        # Use the first module's dir
        first = modules[0]
        (run_dir / first / "summary.pdf").write_bytes(b"%PDF-1.4\n%fake\n")

    if add_oversize_png:
        # Just over 10 MB to trip the per-image cap.
        big = run_dir / modules[0] / "huge.png"
        big.parent.mkdir(parents=True, exist_ok=True)
        with big.open("wb") as fh:
            fh.write(b"\x00" * (10 * 1024 * 1024 + 16))

    return run_dir


def test_report_writes_html(tmp_path, capsys, scfactory):
    run_dir = _make_run_dir(tmp_path)
    rc = scfactory.main(["report", str(run_dir)])
    out = capsys.readouterr().out
    assert rc == 0, capsys.readouterr().err
    out_path_line = out.strip().splitlines()[-1]
    out_path = Path(out_path_line)
    assert out_path.is_file(), f"report html not written: {out_path}"

    text = out_path.read_text()
    assert "DEMO_RUN" in text  # project name
    assert "clustering" in text
    assert "differential_expression" in text
    assert "MUST_NOT_TRUST_AS_CANONICAL" in text  # claim_guard string
    assert "data:image/png;base64," in text  # images embedded
    # claim_guard appears at top AND bottom (two banners).
    assert text.count('class="claim-guard-banner"') >= 2


def test_report_exit_2_when_manifest_missing(tmp_path, capsys, scfactory):
    empty = tmp_path / "empty_run"
    empty.mkdir()
    rc = scfactory.main(["report", str(empty)])
    err = capsys.readouterr().err
    assert rc == 2
    assert "run_manifest.json" in err


def test_report_links_pdfs_not_embedded(tmp_path, capsys, scfactory):
    run_dir = _make_run_dir(tmp_path, add_pdf=True)
    rc = scfactory.main(["report", str(run_dir)])
    out = capsys.readouterr().out
    assert rc == 0
    out_path = Path(out.strip().splitlines()[-1])
    text = out_path.read_text()
    # PDF must be referenced via <a href> with relative path, not embedded.
    assert "summary.pdf" in text
    assert 'href="clustering/summary.pdf"' in text
    assert "data:application/pdf;base64," not in text


def test_report_size_ceiling(tmp_path, capsys, scfactory):
    run_dir = _make_run_dir(
        tmp_path, modules=("clustering",), add_oversize_png=True
    )
    rc = scfactory.main(["report", str(run_dir)])
    out = capsys.readouterr().out
    assert rc == 0
    out_path = Path(out.strip().splitlines()[-1])
    text = out_path.read_text()
    # The oversized PNG must be linked, not embedded.
    assert 'href="clustering/huge.png"' in text
    # And a notice must mention it.
    assert "huge.png" in text
    assert "10 MB" in text or "exceeds per-image" in text
    # The two small PNGs are still embedded (under cap).
    assert "data:image/png;base64," in text


# ---------------------------------------------------------------------------
# recipes
# ---------------------------------------------------------------------------

# All starter recipes shipped under recipes/. Tests assert presence by name.
_STARTER_RECIPES = (
    "nc2024_paper",
    "quick_explore",
    "cite_seq_full",
    "visium_neighborhoods",
)


def test_list_recipes_prints_available(capsys, scfactory):
    rc = scfactory.main(["run", "--list-recipes"])
    out = capsys.readouterr().out
    assert rc == 0, out
    for name in _STARTER_RECIPES:
        assert name in out, f"recipe {name!r} missing from --list-recipes output"


def test_recipe_loads_and_applies_dry_run(tmp_path, capsys, scfactory):
    h5ad = _make_minimal_adata(tmp_path)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "quick_explore", "--dry-run",
    ])
    out = capsys.readouterr().out
    assert rc == 0, out
    # Recipe identity is surfaced.
    assert "recipe=quick_explore" in out or "recipe: quick_explore" in out
    # Recipe's optional modules show up in the planned list.
    assert "clustering" in out
    assert "differential_expression" in out
    # quick_explore omits annotation/trajectory; the dry-run module list
    # should reflect that (the line "optional modules: [...]" is the planned
    # list — assert annotation is not in that specific line).
    plan_lines = [
        ln for ln in out.splitlines() if "optional modules:" in ln
    ]
    assert plan_lines, "no 'optional modules:' line in dry-run output"
    plan_line = plan_lines[0]
    assert "annotation" not in plan_line
    assert "trajectory" not in plan_line
    # scale_preset reflected.
    assert "scale_preset: standard" in out


def test_recipe_with_explicit_modules_user_wins(tmp_path, capsys, scfactory):
    h5ad = _make_minimal_adata(tmp_path)
    rc = scfactory.main([
        "run", str(h5ad),
        "--recipe", "nc2024_paper",
        "--optional-modules", "clustering",
        "--dry-run",
    ])
    out = capsys.readouterr().out
    assert rc == 0, out
    plan_lines = [
        ln for ln in out.splitlines() if "optional modules:" in ln
    ]
    assert plan_lines
    plan_line = plan_lines[0]
    # Only `clustering` should be in the planned module list — assert that
    # the recipe-only modules (paper_repro, annotation, trajectory) are
    # absent from that planned list line.
    assert "clustering" in plan_line
    assert "paper_repro" not in plan_line
    assert "annotation" not in plan_line
    assert "trajectory" not in plan_line
    # Plan source attribution should also note user-provided override.
    assert "user-provided --optional-modules" in out


def test_recipe_env_applied_to_subprocess(tmp_path, capsys, monkeypatch, scfactory):
    h5ad = _make_minimal_adata(tmp_path)
    captured: dict = {}

    def _fake_run(cmd, **kwargs):
        # Record the env passed to the subprocess invocation.
        captured["cmd"] = cmd
        captured["env"] = kwargs.get("env")

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr(scfactory.subprocess, "run", _fake_run)

    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "nc2024_paper",
        # NOTE: not --dry-run, so subprocess.run is actually invoked.
    ])
    assert rc == 0, capsys.readouterr().out
    assert captured.get("env") is not None
    assert captured["env"].get("SC_CLUSTERING_ENGINE") == "sparse_exact"
    # Recipe env override must NOT mutate the parent process env.
    import os as _os
    assert _os.environ.get("SC_CLUSTERING_ENGINE") != "sparse_exact" or \
        _os.environ.get("SC_CLUSTERING_ENGINE") is None
    # And the underlying CLI command should pass --scale-mode massive
    # (from the recipe's scale_preset).
    cmd = captured["cmd"]
    assert "--scale-mode" in cmd
    assert cmd[cmd.index("--scale-mode") + 1] == "massive"


def test_recipe_unknown_name_exits_2(tmp_path, capsys, scfactory):
    h5ad = _make_minimal_adata(tmp_path)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "does_not_exist", "--dry-run",
    ])
    err = capsys.readouterr().err
    assert rc == 2
    assert "unknown recipe" in err
    # Helpful listing of available recipes.
    for name in _STARTER_RECIPES:
        assert name in err, f"available-recipes listing should mention {name}"


def test_recipe_yaml_missing_pyyaml_clean_error(tmp_path, capsys, monkeypatch, scfactory):
    # Force the lazy yaml import inside scfactory to fail with ImportError.
    real_import = __builtins__["__import__"] if isinstance(__builtins__, dict) \
        else __builtins__.__import__

    def _bad_import(name, *args, **kwargs):
        if name == "yaml":
            raise ImportError("no module named 'yaml' (simulated)")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", _bad_import)

    h5ad = _make_minimal_adata(tmp_path)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "quick_explore", "--dry-run",
    ])
    err = capsys.readouterr().err
    assert rc != 0
    # Error message must include the install hint.
    assert "pip install pyyaml" in err.lower()
