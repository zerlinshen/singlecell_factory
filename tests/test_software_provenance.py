"""Tests for workflow.modular.software_provenance and its manifest wiring.

Regression target: before this recorder existed, a produced manifest named the
git SHA of the factory but not a single scientific package version, so a run
could not be reproduced from its own record. Verified on a real production
manifest (round9-singlecell-comparison, 2026-05-22): zero of scanpy, anndata,
numpy, scipy, scikit-learn, leidenalg, igraph, harmonypy, scrublet,
rapids-singlecell, pandas, matplotlib appeared anywhere in the file.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import warnings
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as dist_version
from pathlib import Path

import pytest

from workflow.modular.config import CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext
from workflow.modular.manifest_writer import write_manifest
from workflow.modular.pipeline import _save_manifest
from workflow.modular.software_provenance import (
    CORE_PACKAGE_NAMES,
    OPTIONAL_PACKAGE_NAMES,
    RECORDED_PACKAGE_NAMES,
    collect_software_versions,
    resolve_package_version,
)

ROOT = Path(__file__).resolve().parents[1]

# The packages whose version can move a published number. Kept as a literal list
# rather than importing CORE_PACKAGE_NAMES so that silently dropping a package
# from the recorder fails a test instead of weakening it.
EXPECTED_CORE = {
    "scanpy",
    "anndata",
    "numpy",
    "scipy",
    "scikit-learn",
    "pandas",
    "matplotlib",
    "leidenalg",
    "igraph",
}
EXPECTED_OPTIONAL = {"rapids-singlecell", "harmonypy", "scrublet"}


# ---------------------------------------------------------------------------
# collect_software_versions
# ---------------------------------------------------------------------------

def test_block_has_expected_shape():
    block = collect_software_versions()
    for key in (
        "python",
        "python_implementation",
        "platform",
        "executable",
        "packages",
        "core_packages",
        "unresolved",
        "core_complete",
    ):
        assert key in block, f"software_versions missing {key!r}"
    assert re.match(r"^\d+\.\d+\.\d+", block["python"]), block["python"]
    assert block["executable"] == sys.executable


def test_core_and_optional_packages_are_all_recorded():
    """Every package in the contract has a key, present or absent."""
    packages = collect_software_versions()["packages"]
    assert EXPECTED_CORE <= set(packages), EXPECTED_CORE - set(packages)
    assert EXPECTED_OPTIONAL <= set(packages), EXPECTED_OPTIONAL - set(packages)
    assert set(CORE_PACKAGE_NAMES) == EXPECTED_CORE
    assert set(OPTIONAL_PACKAGE_NAMES) == EXPECTED_OPTIONAL
    assert set(RECORDED_PACKAGE_NAMES) == EXPECTED_CORE | EXPECTED_OPTIONAL


def test_core_packages_resolve_to_real_versions():
    """The block is populated, not a dict of nulls.

    Any interpreter able to import this test has the core stack installed, so
    an unresolved core package means the resolver — not the environment — is
    broken.
    """
    block = collect_software_versions()
    packages = block["packages"]
    for name in EXPECTED_CORE:
        assert packages[name], f"core package {name} resolved to {packages[name]!r}"
        assert re.match(r"^\d+", str(packages[name])), packages[name]
    assert block["core_complete"] is True
    assert not (set(block["unresolved"]) & EXPECTED_CORE)


@pytest.mark.parametrize("name", sorted(EXPECTED_CORE))
def test_recorded_version_matches_installed_distribution(name):
    """Recorded value equals what importlib.metadata reports for the same env."""
    packages = collect_software_versions()["packages"]
    candidates = [name, name.replace("-", "_")]
    if name == "scikit-learn":
        candidates.append("sklearn")
    if name == "igraph":
        candidates.append("python-igraph")
    truth = None
    for candidate in candidates:
        try:
            truth = dist_version(candidate)
            break
        except PackageNotFoundError:
            continue
    assert truth is not None, f"{name} not installed; cannot validate recorder"
    assert packages[name] == truth


def test_absent_package_is_recorded_as_none_not_omitted():
    """A missing backend must be provable from the manifest, not merely absent."""
    assert resolve_package_version(("definitely-not-installed-xyz",)) is None
    assert resolve_package_version("definitely-not-installed-xyz", "no_such_module") is None


def test_resolver_falls_back_to_already_imported_module(monkeypatch):
    """Packages installed without metadata still resolve, without new imports."""
    fake = type(sys)("fake_prov_mod")
    fake.__version__ = "9.9.9"
    monkeypatch.setitem(sys.modules, "fake_prov_mod", fake)
    assert resolve_package_version(("no-such-dist",), "fake_prov_mod") == "9.9.9"


def test_resolver_never_raises_on_broken_metadata(monkeypatch):
    """Bookkeeping must not be able to abort a completed run."""
    def boom(_name):
        raise RuntimeError("corrupt dist-info")

    monkeypatch.setattr(
        "workflow.modular.software_provenance._distribution_version", boom
    )
    block = collect_software_versions()  # must not raise
    assert set(block["packages"]) == set(RECORDED_PACKAGE_NAMES)

    # With metadata dead, only packages already imported by an earlier test can
    # still resolve via the __version__ fallback, so the expectation is stated
    # per package rather than as a blanket "everything is unresolved" — the
    # latter would pass alone and fail inside the full suite.
    import_names = {"scikit-learn": "sklearn", "rapids-singlecell": "rapids_singlecell"}
    for name, value in block["packages"].items():
        module = sys.modules.get(import_names.get(name, name))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # anndata.__version__ is deprecated
            fallback = getattr(module, "__version__", None) if module is not None else None
        if fallback:
            assert value == str(fallback)
        else:
            assert value is None, f"{name} resolved to {value!r} with metadata broken"
            assert name in block["unresolved"]


def test_block_is_json_serialisable():
    json.dumps(collect_software_versions())


def test_recorder_does_not_import_the_packages_it_reports():
    """Metadata lookup only — importing scanpy/cupy here would cost seconds and
    could initialise a CUDA context in a process that never asked for one."""
    code = (
        "import sys, json;"
        "from workflow.modular.software_provenance import collect_software_versions;"
        "b = collect_software_versions();"
        "print(json.dumps({'scanpy': b['packages']['scanpy'],"
        " 'imported': [m for m in ('scanpy','anndata','sklearn','matplotlib')"
        " if m in sys.modules]}))"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=120
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout.strip().splitlines()[-1])
    assert payload["scanpy"], "scanpy version not resolved in a clean interpreter"
    assert payload["imported"] == [], f"recorder imported {payload['imported']}"


# ---------------------------------------------------------------------------
# Manifest wiring — these fail without the software_versions block
# ---------------------------------------------------------------------------

def test_project_root_manifest_carries_software_versions(tmp_path):
    out = write_manifest(
        tmp_path / "run",
        project_id="prov",
        run_id="2026-07-29T0000Z-0000000",
        modules_run=["qc"],
    )
    block = json.loads(out.read_text(encoding="utf-8"))["software_versions"]
    assert EXPECTED_CORE <= set(block["packages"])
    assert block["packages"]["scanpy"]
    assert block["packages"]["anndata"]
    assert block["packages"]["numpy"]


def test_run_manifest_carries_software_versions(tmp_path):
    """run_manifest.json is written for every run, including runs without
    --project-root, so it is the manifest most likely to be the only record."""
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    cfg = PipelineConfig(
        project="prov",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(
            sample_root=tmp_path / "input",
            outs_dir=tmp_path / "input" / "outs" / "filtered_feature_bc_matrix",
        ),
        optional_modules=[],
    )
    ctx = PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)
    ctx.status("qc", True, "ok")

    payload = json.loads(_save_manifest(ctx).read_text(encoding="utf-8"))
    block = payload["software_versions"]
    assert EXPECTED_CORE <= set(block["packages"])
    assert block["packages"]["scanpy"]
    assert block["core_complete"] is True


# ---------------------------------------------------------------------------
# Environment-file floors
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("env_file", ["environment.yml", "environment_gpu.yml"])
def test_declared_floors_are_met_by_this_interpreter(env_file):
    """Guard against declaring a floor the working environments do not satisfy.

    A wrong pin breaks the solve, which is worse than no pin, so every `>=`
    written into an environment file must hold for the stack these tests run on.
    """
    conda_to_dist = {"python-igraph": ("igraph", "python-igraph")}
    text = (ROOT / env_file).read_text(encoding="utf-8")
    floors = re.findall(r"^\s*-\s*([A-Za-z0-9_.-]+)>=([0-9][0-9A-Za-z.]*)\s*$", text, re.M)
    assert floors, f"{env_file} declares no >= floors; update this test if intentional"

    for pkg, floor in floors:
        installed = resolve_package_version(conda_to_dist.get(pkg, (pkg,)))
        if installed is None:
            continue  # optional backend absent from this env — not this test's business
        got = tuple(int(p) for p in re.findall(r"\d+", installed)[:3])
        want = tuple(int(p) for p in re.findall(r"\d+", floor)[:3])
        got += (0,) * (len(want) - len(got))
        assert got >= want, f"{env_file}: {pkg}>={floor} declared but {installed} installed"
