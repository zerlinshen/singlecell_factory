"""scfactory wiring for named analysis profiles and batch declarations.

Kept separate from tests/test_scfactory.py so the two surfaces can evolve
independently. Covers the launcher-side half of the 2026-08-03 multi-batch fix:
the warning must reach the operator on the DRY-RUN path (previewing a plan is
the last moment the plan can still be changed), and a contradicted declaration
must stop before any runnable command is printed.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
RECIPES_DIR = REPO_ROOT / "recipes"


@pytest.fixture(scope="module")
def scfactory():
    spec = importlib.util.spec_from_file_location(
        "scfactory_batch_profiles_under_test", REPO_ROOT / "scripts" / "scfactory.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _h5ad(tmp_path: Path, name: str, n_samples: int) -> Path:
    n = 40
    obs = pd.DataFrame(
        {"sample": pd.Categorical([f"S{i % n_samples}" for i in range(n)])},
        index=[f"c{i}" for i in range(n)],
    )
    path = tmp_path / name
    ad.AnnData(
        X=np.zeros((n, 3), dtype="float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    ).write_h5ad(path)
    return path


# ---------------------------------------------------------------------------
# recipes resolve through the catalog
# ---------------------------------------------------------------------------


def test_profile_recipes_are_shipped_and_listed(capsys, scfactory):
    assert (RECIPES_DIR / "single_batch.yaml").is_file()
    assert (RECIPES_DIR / "multi_batch_harmony.yaml").is_file()
    assert scfactory.main(["run", "--list-recipes"]) == 0
    out = capsys.readouterr().out
    assert "single_batch" in out
    assert "multi_batch_harmony" in out


def test_profile_recipes_carry_no_copied_module_list(scfactory):
    """The catalog owns the module names; the YAML must not restate them."""
    yaml = scfactory._import_yaml()
    for name in ("single_batch", "multi_batch_harmony"):
        data = yaml.safe_load((RECIPES_DIR / f"{name}.yaml").read_text())
        assert data["profile"] == name
        assert "optional_modules" not in data


def test_recipe_may_not_set_both_profile_and_modules(tmp_path, capsys, scfactory):
    recipe = tmp_path / "recipes" / "both.yaml"
    recipe.parent.mkdir(parents=True)
    recipe.write_text(
        "name: both\nprofile: single_batch\noptional_modules:\n  - clustering\n"
    )
    original = scfactory.RECIPES_DIR
    scfactory.RECIPES_DIR = recipe.parent
    try:
        with pytest.raises(SystemExit) as excinfo:
            scfactory._load_recipe("both")
    finally:
        scfactory.RECIPES_DIR = original
    assert excinfo.value.code == 2
    assert "not both" in capsys.readouterr().err


def test_recipe_naming_an_unknown_profile_is_rejected(tmp_path, capsys, scfactory):
    recipe = tmp_path / "recipes" / "bogus.yaml"
    recipe.parent.mkdir(parents=True)
    recipe.write_text("name: bogus\nprofile: no_such_profile\n")
    original = scfactory.RECIPES_DIR
    scfactory.RECIPES_DIR = recipe.parent
    try:
        with pytest.raises(SystemExit) as excinfo:
            scfactory._load_recipe("bogus")
    finally:
        scfactory.RECIPES_DIR = original
    assert excinfo.value.code == 2
    assert "unknown analysis profile" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# dry-run behaviour on multi-batch input
# ---------------------------------------------------------------------------


def test_dry_run_warns_on_multi_batch_input_with_default_plan(
    tmp_path, capsys, scfactory
):
    """The exact defect: previewing the default plan said nothing about batches."""
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    assert scfactory.main(["run", str(h5ad), "--dry-run"]) == 0
    captured = capsys.readouterr()
    assert "BATCH_RISK:" in captured.err
    assert "5 batches" in captured.err
    assert "claim=exploratory" in captured.out


def test_dry_run_is_silent_on_single_batch_input(tmp_path, capsys, scfactory):
    h5ad = _h5ad(tmp_path, "single.h5ad", n_samples=1)
    assert scfactory.main(["run", str(h5ad), "--dry-run"]) == 0
    captured = capsys.readouterr()
    assert "BATCH_RISK:" not in captured.err
    assert "claim=claimable" in captured.out


def test_multi_batch_harmony_recipe_plans_integration_and_clears_the_claim(
    tmp_path, capsys, scfactory
):
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "multi_batch_harmony", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 0, captured.err
    assert "BATCH_RISK:" not in captured.err
    assert "batch_correction" in captured.out
    assert "claim=claimable" in captured.out
    assert "profile=multi_batch_harmony" in captured.out
    # The declaration is forwarded to the canonical CLI, not just printed.
    assert "--batch-strategy integrate" in captured.out


def test_single_batch_recipe_on_multi_batch_input_refuses_to_plan(
    tmp_path, capsys, scfactory
):
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "single_batch", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 2
    assert "single-batch" in captured.err
    assert "no dry-run or real run was started" in captured.err
    assert "would execute" not in captured.out, (
        "a contradicted declaration must not print a runnable command"
    )


def test_single_batch_recipe_accepts_genuinely_single_batch_input(
    tmp_path, capsys, scfactory
):
    h5ad = _h5ad(tmp_path, "single.h5ad", n_samples=1)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "single_batch", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 0, captured.err
    assert "--batch-strategy single-batch" in captured.out


def test_accept_uncorrected_declares_but_still_warns(tmp_path, capsys, scfactory):
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    rc = scfactory.main([
        "run", str(h5ad), "--batch-strategy", "accept-uncorrected", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 0, captured.err
    assert "BATCH_RISK:" in captured.err
    assert "claim=exploratory" in captured.out
    assert "--batch-strategy accept-uncorrected" in captured.out


# ---------------------------------------------------------------------------
# precedence
# ---------------------------------------------------------------------------


def test_explicit_batch_strategy_wins_over_the_recipe(tmp_path, capsys, scfactory):
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "single_batch",
        "--batch-strategy", "accept-uncorrected", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 0, captured.err
    assert "accept-uncorrected (user-provided --batch-strategy)" in captured.out


def test_explicit_optional_modules_still_wins_over_a_recipe(
    tmp_path, capsys, scfactory
):
    """Documented precedence is unchanged: --optional-modules > --recipe."""
    h5ad = _h5ad(tmp_path, "single.h5ad", n_samples=1)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "quick_explore",
        "--optional-modules", "clustering", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 0, captured.err
    assert "user-provided --optional-modules" in captured.out
    assert "optional modules: ['clustering']" in captured.out


def test_module_override_of_a_profile_recipe_needs_its_own_declaration(
    tmp_path, capsys, scfactory
):
    """Modules and the declaration are overridden independently.

    --optional-modules wins for modules, but it cannot silently satisfy the
    profile's inherited 'integrate' declaration; releasing that takes an
    explicit --batch-strategy, which is then what the manifest records.
    """
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "multi_batch_harmony",
        "--optional-modules", "clustering",
        "--batch-strategy", "accept-uncorrected", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 0, captured.err
    assert "optional modules: ['clustering']" in captured.out
    assert "claim=exploratory" in captured.out
    assert "BATCH_RISK:" in captured.err


def test_overriding_modules_away_from_integration_fails_the_declaration(
    tmp_path, capsys, scfactory
):
    """Dropping batch_correction while the recipe declares integrate must fail."""
    h5ad = _h5ad(tmp_path, "multi.h5ad", n_samples=5)
    rc = scfactory.main([
        "run", str(h5ad), "--recipe", "multi_batch_harmony",
        "--optional-modules", "clustering,annotation", "--dry-run",
    ])
    captured = capsys.readouterr()
    assert rc == 2
    assert "batch_correction" in captured.err
    assert "would execute" not in captured.out
