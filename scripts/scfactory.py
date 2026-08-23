#!/usr/bin/env python3
"""scfactory — thin user-friendly CLI wrapper over the modular pipeline.

This is a UX layer ONLY: it dispatches to the canonical
``python -m workflow.modular.cli`` invocation and to
``scripts/export_singlecell_r_bundle.py``. It never alters their behavior.

Subcommands:
  run     auto-detect modality, build optional-modules list, run pipeline
  doctor  read-only environment health check
  report  walk a run dir and render a single self-contained HTML report

Hard contracts:
  * The dependency-light canonical module catalog is the only source for
    auto-selected modules. anndata is imported lazily for ``run`` modality
    detection on .h5ad inputs.
  * No state mutation in ``doctor`` or ``report``.
  * Every planned argument is forwarded visibly to the canonical CLI.
"""

from __future__ import annotations

import argparse
import base64
import csv
import datetime as _dt
import html
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPORT_SCRIPT = REPO_ROOT / "scripts" / "export_singlecell_r_bundle.py"
RECIPES_DIR = REPO_ROOT / "recipes"
DEFAULT_PROJECTS_ROOT = Path.home() / "projects"
SCFACTORY_VERSION = "0.4.0"

# Report constants (HARD limits per spec).
REPORT_PNG_EMBED_MAX_BYTES = 10 * 1024 * 1024  # 10 MB per PNG
REPORT_TOTAL_EMBED_MAX_BYTES = 50 * 1024 * 1024  # 50 MB total embedded
REPORT_WALK_MAX_DEPTH = 3

# Import the dependency-light canonical catalog from the repository even when
# this file is invoked as an absolute script path (where sys.path[0] is
# ``scripts/``). Import failures are fatal: a stale local fallback would make
# this adapter advertise a plan different from the canonical CLI.
_repo_path = str(REPO_ROOT)
_added_repo_path = _repo_path not in sys.path
if _added_repo_path:
    sys.path.insert(0, _repo_path)
try:
    from workflow.modular.batch_risk import (
        BATCH_STRATEGY_AUTO,
        BATCH_STRATEGY_CHOICES,
        BatchStrategyConflict,
        announce_batch_risk,
        detect_batch_structure,
        resolve_batch_risk,
        resolve_obs_source,
    )
    from workflow.modular.legacy_output import LEGACY_DEFAULT_OUTPUT_DIR
    from workflow.modular.module_catalog import MODULE_SPECS
    from workflow.modular.module_catalog import analysis_profile
    from workflow.modular.module_catalog import analysis_profile_names
    from workflow.modular.module_catalog import optional_modules_for_modality
finally:
    if _added_repo_path:
        sys.path.remove(_repo_path)


# ---------------------------------------------------------------------------
# Modality detection
# ---------------------------------------------------------------------------

def _detect_h5ad(path: Path) -> tuple[str, list[str]]:
    """Return (modality, evidence) for an .h5ad input.

    Modalities: "rna_only", "cite_seq", "spatial", "multimodal".
    """
    try:
        import anndata as ad  # type: ignore
    except ImportError as exc:  # pragma: no cover - tested env always has it
        raise SystemExit(
            "anndata is required for .h5ad detection. "
            "Run inside the singlecell_factory python env."
        ) from exc

    adata = ad.read_h5ad(path, backed="r")
    try:
        obsm_keys = list(adata.obsm.keys())
    finally:
        # Close backing file when possible.
        try:
            adata.file.close()  # type: ignore[attr-defined]
        except Exception:
            pass

    has_protein = any(k in obsm_keys for k in ("protein_counts", "protein_clr"))
    has_spatial = "spatial" in obsm_keys
    evidence = [f"obsm keys: {obsm_keys}"]

    if has_protein and has_spatial:
        return "multimodal", evidence
    if has_protein:
        return "cite_seq", evidence
    if has_spatial:
        return "spatial", evidence
    return "rna_only", evidence


def _detect_sample_root(path: Path) -> tuple[str, list[str]]:
    """Return (modality, evidence) for a sample-root directory."""
    evidence: list[str] = []

    # Spaceranger layout probe.
    spatial_marker = path / "spatial" / "tissue_positions_list.csv"
    spatial_marker_alt = path / "outs" / "spatial" / "tissue_positions_list.csv"
    if spatial_marker.exists() or spatial_marker_alt.exists():
        hit = spatial_marker if spatial_marker.exists() else spatial_marker_alt
        evidence.append(f"spaceranger tissue_positions_list at {hit}")
        return "spatial", evidence

    # CITE-seq probe: cellranger filtered_feature_bc_matrix.h5 with an
    # 'Antibody Capture' modality string. Probe both at sample-root and
    # within outs/.
    candidates = [
        path / "filtered_feature_bc_matrix.h5",
        path / "outs" / "filtered_feature_bc_matrix.h5",
    ]
    h5_path: Path | None = None
    for cand in candidates:
        if cand.exists():
            h5_path = cand
            break

    if h5_path is not None:
        evidence.append(f"cellranger h5: {h5_path}")
        try:
            import h5py  # type: ignore

            with h5py.File(h5_path, "r") as fh:
                feat = fh.get("matrix/features/feature_type")
                if feat is not None:
                    types = {x.decode() if isinstance(x, bytes) else str(x)
                             for x in feat[:]}
                    evidence.append(f"feature_types={sorted(types)}")
                    if "Antibody Capture" in types:
                        return "cite_seq", evidence
        except ImportError:
            evidence.append("h5py unavailable; cannot probe feature_type")
        except Exception as exc:
            evidence.append(f"h5 probe failed: {exc}")

        return "rna_only", evidence

    # No matrix file found.
    evidence.append("no filtered_feature_bc_matrix.h5 or spaceranger layout")
    return "unknown", evidence


def detect_modality(input_path: Path) -> tuple[str, list[str]]:
    """Dispatch on input type. Returns ("unknown", evidence) if not detectable."""
    if not input_path.exists():
        return "unknown", [f"path does not exist: {input_path}"]
    if input_path.is_file() and input_path.suffix.lower() == ".h5ad":
        return _detect_h5ad(input_path)
    if input_path.is_dir():
        return _detect_sample_root(input_path)
    return "unknown", [f"unsupported input type: {input_path}"]


def plan_optional_modules(modality: str) -> list[str]:
    """Map modality -> optional module list."""
    return list(optional_modules_for_modality(modality))


# ---------------------------------------------------------------------------
# Recipes (preset bundles of optional_modules + env + bundle config)
# ---------------------------------------------------------------------------

VALID_MODALITY_HINTS = ("rna_only", "cite_seq", "spatial", "multimodal")
VALID_SCALE_PRESETS = ("standard", "large", "massive")


def _import_yaml():
    """Lazy-import PyYAML; raise a clear error pointing to install command.

    PyYAML is the ONLY recipe-related optional dependency. It is gated behind
    `--recipe` / `--list-recipes` so the default no-recipe path needs nothing.
    """
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "scfactory: recipes require PyYAML, but it is not installed.\n"
            "  install with: pip install pyyaml\n"
            "  (PyYAML is only needed when using --recipe / --list-recipes;"
            " the default scfactory run path does not need it.)"
        ) from exc
    return yaml


def _list_recipe_files() -> list[Path]:
    """Return all *.yaml/*.yml files in the recipes directory, sorted."""
    if not RECIPES_DIR.is_dir():
        return []
    files = list(RECIPES_DIR.glob("*.yaml")) + list(RECIPES_DIR.glob("*.yml"))
    return sorted(files)


def _recipe_path_for_name(name: str) -> Path | None:
    """Locate a recipe file by stem name. Returns None if not found."""
    for ext in (".yaml", ".yml"):
        cand = RECIPES_DIR / f"{name}{ext}"
        if cand.is_file():
            return cand
    return None


def _load_recipe(name: str) -> dict[str, Any]:
    """Load and validate a recipe by name. Raises SystemExit on errors."""
    yaml = _import_yaml()
    path = _recipe_path_for_name(name)
    if path is None:
        available = sorted(p.stem for p in _list_recipe_files())
        listing = ", ".join(available) if available else "(none found)"
        print(
            f"scfactory: unknown recipe '{name}'.\n"
            f"  available recipes: {listing}\n"
            f"  list with: scfactory run --list-recipes",
            file=sys.stderr,
        )
        raise SystemExit(2)

    try:
        with path.open() as fh:
            data = yaml.safe_load(fh)
    except Exception as exc:
        print(
            f"scfactory: failed to parse recipe {path}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(2) from exc

    if not isinstance(data, dict):
        print(
            f"scfactory: recipe {path} did not parse to a mapping "
            f"(got {type(data).__name__}).",
            file=sys.stderr,
        )
        raise SystemExit(2)

    _validate_recipe(data, path, name)
    return data


def _known_module_names() -> set[str]:
    """Return canonical catalog names used to validate recipe modules."""

    return set(MODULE_SPECS)


def _validate_recipe(data: dict[str, Any], path: Path, expected_name: str) -> None:
    """Validate recipe schema. Raises SystemExit(2) with a clear error on failure."""
    where = f"recipe {path}"

    # name must match filename stem
    name = data.get("name")
    if not isinstance(name, str) or not name:
        print(f"scfactory: {where}: missing required 'name' string", file=sys.stderr)
        raise SystemExit(2)
    if name != expected_name:
        print(
            f"scfactory: {where}: 'name' field ({name!r}) does not match "
            f"filename stem ({expected_name!r}).",
            file=sys.stderr,
        )
        raise SystemExit(2)

    # description optional but must be str if present
    if "description" in data and not isinstance(data["description"], str):
        print(f"scfactory: {where}: 'description' must be a string", file=sys.stderr)
        raise SystemExit(2)

    # modality_hint optional, must be in VALID_MODALITY_HINTS if present
    mh = data.get("modality_hint")
    if mh is not None and mh not in VALID_MODALITY_HINTS:
        print(
            f"scfactory: {where}: 'modality_hint' must be one of "
            f"{list(VALID_MODALITY_HINTS)}, got {mh!r}",
            file=sys.stderr,
        )
        raise SystemExit(2)

    # Exactly one module source: a named catalog profile, or an inline list.
    # `profile` is preferred for designs the catalog owns — it resolves through
    # module_catalog.ANALYSIS_PROFILES instead of copying module names into
    # YAML, so the catalog stays the single source of truth.
    profile_name = data.get("profile")
    mods = data.get("optional_modules")
    if profile_name is not None and mods is not None:
        print(
            f"scfactory: {where}: set either 'profile' or 'optional_modules', "
            f"not both.",
            file=sys.stderr,
        )
        raise SystemExit(2)
    if profile_name is not None:
        if not isinstance(profile_name, str):
            print(
                f"scfactory: {where}: 'profile' must be a string", file=sys.stderr
            )
            raise SystemExit(2)
        try:
            analysis_profile(profile_name)
        except ValueError as exc:
            print(f"scfactory: {where}: {exc}", file=sys.stderr)
            raise SystemExit(2) from exc
        mods = []
    if not isinstance(mods, list) or not all(isinstance(m, str) for m in mods):
        print(
            f"scfactory: {where}: 'optional_modules' is required (or a "
            f"'profile' naming one of: {', '.join(analysis_profile_names())}) "
            f"and must be a list of strings.",
            file=sys.stderr,
        )
        raise SystemExit(2)
    known = _known_module_names()
    unknown = [m for m in mods if m not in known]
    if unknown:
        print(
            f"scfactory: {where}: unknown module name(s) {unknown!r}.\n"
            f"  recipe references modules that are not in "
            f"workflow.modular.module_catalog.MODULE_SPECS.\n"
            f"  fix the offending recipe at {path}.",
            file=sys.stderr,
        )
        raise SystemExit(2)

    # scale_preset optional, must be in VALID_SCALE_PRESETS if present
    sp = data.get("scale_preset")
    if sp is not None and sp not in VALID_SCALE_PRESETS:
        print(
            f"scfactory: {where}: 'scale_preset' must be one of "
            f"{list(VALID_SCALE_PRESETS)}, got {sp!r}",
            file=sys.stderr,
        )
        raise SystemExit(2)

    # batch_strategy optional; must be a declared strategy when present. A
    # recipe that names a profile inherits the profile's strategy instead.
    bs = data.get("batch_strategy")
    if bs is not None and bs not in BATCH_STRATEGY_CHOICES:
        print(
            f"scfactory: {where}: 'batch_strategy' must be one of "
            f"{list(BATCH_STRATEGY_CHOICES)}, got {bs!r}",
            file=sys.stderr,
        )
        raise SystemExit(2)

    # env optional, dict[str, str]
    env = data.get("env")
    if env is not None:
        if not isinstance(env, dict) or not all(
            isinstance(k, str) and isinstance(v, (str, int, float, bool))
            for k, v in env.items()
        ):
            print(
                f"scfactory: {where}: 'env' must be a flat mapping of "
                f"string keys to scalar values.",
                file=sys.stderr,
            )
            raise SystemExit(2)

    # bundle optional, dict
    bundle = data.get("bundle")
    if bundle is not None and not isinstance(bundle, dict):
        print(f"scfactory: {where}: 'bundle' must be a mapping", file=sys.stderr)
        raise SystemExit(2)

    # required_python_packages optional: list of importable module names that
    # must be present before the recipe may run (Wave-2 spatial preflight).
    rpp = data.get("required_python_packages")
    if rpp is not None:
        if not isinstance(rpp, list) or not all(isinstance(x, str) and x for x in rpp):
            print(
                f"scfactory: {where}: 'required_python_packages' must be a "
                f"list of non-empty strings.",
                file=sys.stderr,
            )
            raise SystemExit(2)


# Modules that require squidpy for their primary analytics path.
_MODULES_REQUIRING_SQUIDPY = frozenset({"spatial_neighborhoods"})


def _missing_python_packages(names: list[str]) -> list[str]:
    """Return package names that fail importlib.import_module."""
    import importlib

    missing: list[str] = []
    for name in names:
        try:
            importlib.import_module(name)
        except Exception:
            missing.append(name)
    return missing


def _recipe_and_module_required_packages(
    recipe: dict[str, Any] | None, planned_modules: list[str]
) -> list[str]:
    """Union of recipe-declared required packages and module-implied deps."""
    required: list[str] = []
    if recipe is not None:
        for name in recipe.get("required_python_packages") or []:
            if isinstance(name, str) and name and name not in required:
                required.append(name)
    if any(m in _MODULES_REQUIRING_SQUIDPY for m in planned_modules):
        if "squidpy" not in required:
            required.append("squidpy")
    return required


def cmd_list_recipes() -> int:
    """Print available recipes (one per line: 'name: description'). Exit 0."""
    files = _list_recipe_files()
    if not files:
        print(f"scfactory: no recipes found in {RECIPES_DIR}")
        return 0
    yaml = _import_yaml()
    for path in files:
        try:
            with path.open() as fh:
                data = yaml.safe_load(fh) or {}
        except Exception as exc:
            print(f"{path.stem}: <unparseable: {exc}>")
            continue
        name = data.get("name", path.stem)
        desc = data.get("description", "")
        print(f"{name}: {desc}")
    return 0


# ---------------------------------------------------------------------------
# `run` subcommand
# ---------------------------------------------------------------------------

def cmd_run(args: argparse.Namespace) -> int:
    # `--list-recipes` is a self-contained mode that ignores `input`.
    if getattr(args, "list_recipes", False):
        return cmd_list_recipes()

    # Load recipe early so that "unknown recipe" errors fire before any other
    # validation (matches user mental model: misspelled --recipe fails fast,
    # regardless of whether `input` was provided).
    recipe: dict[str, Any] | None = None
    if args.recipe:
        recipe = _load_recipe(args.recipe)

    if not args.input:
        print(
            "scfactory: 'input' is required (.h5ad file or sample-root dir).\n"
            "  to list available recipes: scfactory run --list-recipes",
            file=sys.stderr,
        )
        return 2

    for option, attribute in (
        ("--atac-peak-matrix-path", "atac_peak_matrix_path"),
        ("--atac-peaks-bed-path", "atac_peaks_bed_path"),
    ):
        raw_path = getattr(args, attribute)
        if not raw_path:
            continue
        resolved_path = Path(raw_path).expanduser().resolve()
        if not resolved_path.is_file():
            print(f"scfactory: {option} file not found: {resolved_path}", file=sys.stderr)
            return 2
        setattr(args, attribute, str(resolved_path))
    if args.atac_n_components is not None and args.atac_n_components < 2:
        print("scfactory: --atac-n-components must be at least 2.", file=sys.stderr)
        return 2

    input_path = Path(args.input).resolve()

    if not input_path.exists():
        print(
            f"scfactory: input not found: {input_path}\n"
            "  pass an existing .h5ad file or sample-root dir, or run "
            "'scfactory doctor' for environment checks.",
            file=sys.stderr,
        )
        return 2

    # Recipe was already loaded above; its modality_hint can compensate when
    # auto-detection is "unknown" (e.g. a directory with neither cellranger H5
    # nor spaceranger layout, but the user has named the workflow via --recipe).
    modality, evidence = detect_modality(input_path)
    if modality == "unknown":
        if recipe and recipe.get("modality_hint"):
            modality = str(recipe["modality_hint"])
            evidence.append(f"recipe modality_hint={modality}")
        elif not args.optional_modules:
            print(
                f"scfactory: could not auto-detect modality from {input_path}\n"
                "  evidence: " + "; ".join(evidence) + "\n"
                "  pass --optional-modules explicitly, --recipe NAME, or run "
                "scfactory doctor.",
                file=sys.stderr,
            )
            return 2

    # Precedence (documented in --help):
    #   1. --optional-modules (escape hatch) overrides everything.
    #   2. --recipe provides the optional_modules list.
    #   3. Auto-detect (existing logic) fills the gap.
    # A recipe naming a catalog profile resolves its modules AND its batch
    # declaration from module_catalog, never from a copy in the YAML.
    recipe_profile = (
        analysis_profile(str(recipe["profile"]))
        if recipe is not None and recipe.get("profile")
        else None
    )
    if args.optional_modules:
        planned = [m.strip() for m in args.optional_modules.split(",") if m.strip()]
        plan_source = "user-provided --optional-modules"
    elif recipe_profile is not None:
        planned = list(recipe_profile.optional_modules)
        plan_source = f"recipe={recipe['name']} profile={recipe_profile.name}"
    elif recipe is not None:
        planned = [str(m) for m in recipe.get("optional_modules", [])]
        plan_source = f"recipe={recipe['name']}"
    else:
        planned = plan_optional_modules(modality)
        plan_source = f"auto for modality={modality}"

    # scATAC pseudobulk DA is an explicit contract, never an inference from a
    # generic ATAC modality.  Keep the public dry-run plan identical to the
    # canonical modular CLI selection and fail before starting either route.
    scatac_flag_values = (
        args.scatac_da_sample_col,
        args.scatac_da_group_col,
        args.scatac_da_condition_col,
        args.scatac_da_peak_id_col,
        args.scatac_da_test_level,
        args.scatac_da_reference_level,
        args.scatac_da_groups,
        args.scatac_da_mode,
        args.scatac_da_min_samples_per_condition,
        args.scatac_da_min_total_count,
        args.scatac_da_fdr_threshold,
        args.scatac_da_abs_log2fc_threshold,
        args.scatac_da_r_conda_env,
        args.scatac_da_timeout,
        args.scatac_da_aggregation_backend,
    )
    scatac_requested = any(value is not None for value in scatac_flag_values)
    scatac_selected = scatac_requested or "scatac_pseudobulk_da" in planned
    if scatac_requested and "scatac_pseudobulk_da" not in planned:
        planned.append("scatac_pseudobulk_da")
        plan_source += " + explicit scatac DA contract"
    if scatac_selected:
        mode = args.scatac_da_mode or "confirmatory_da"
        missing = [
            name for name, value in (
                ("--scatac-da-sample-col", args.scatac_da_sample_col),
                ("--scatac-da-group-col", args.scatac_da_group_col),
                ("--scatac-da-peak-id-col", args.scatac_da_peak_id_col),
            ) if not value
        ]
        if mode == "confirmatory_da":
            missing.extend(
                name for name, value in (
                    ("--scatac-da-condition-col", args.scatac_da_condition_col),
                    ("--scatac-da-test-level", args.scatac_da_test_level),
                    ("--scatac-da-reference-level", args.scatac_da_reference_level),
                ) if not value
            )
        if missing:
            print(
                "scfactory: scATAC pseudobulk DA requires explicit " + ", ".join(missing) + ".\n"
                "  no dry-run or real run was started.",
                file=sys.stderr,
            )
            return 2
        if args.scatac_da_min_samples_per_condition is not None and args.scatac_da_min_samples_per_condition < 2:
            print("scfactory: --scatac-da-min-samples-per-condition must be at least 2.", file=sys.stderr)
            return 2
        if args.scatac_da_timeout is not None and args.scatac_da_timeout <= 0:
            print("scfactory: --scatac-da-timeout must be positive.", file=sys.stderr)
            return 2
        if args.scatac_da_r_conda_env not in (None, "r_multiomics"):
            print("scfactory: --scatac-da-r-conda-env must be r_multiomics.", file=sys.stderr)
            return 2
        if args.scatac_da_aggregation_backend not in (None, "cpu"):
            print("scfactory: --scatac-da-aggregation-backend must be cpu.", file=sys.stderr)
            return 2

    # Wave-2 W2.1: fail early when recipe/modules require packages that are
    # not installed (e.g. squidpy for spatial_neighborhoods / visium recipe).
    # Applies to dry-run too so operators discover the gap before a long run.
    required_pkgs = _recipe_and_module_required_packages(recipe, planned)
    missing_pkgs = _missing_python_packages(required_pkgs)
    if missing_pkgs:
        hint_lines = [
            f"scfactory: missing required Python package(s) for this run: "
            f"{missing_pkgs}",
            f"  planned modules: {planned}",
        ]
        if recipe is not None:
            hint_lines.append(f"  recipe: {recipe.get('name')}")
        if "squidpy" in missing_pkgs:
            hint_lines.append(
                "  spatial_neighborhoods / visium_neighborhoods need squidpy "
                "(Palla et al. 2022). Install into the active env, e.g.:"
            )
            hint_lines.append("    pip install 'squidpy>=1.2'   # or conda-forge squidpy")
            hint_lines.append(
                "  then re-run `scfactory doctor --json` and check "
                "readiness.spatial_analytics."
            )
        else:
            hint_lines.append(
                "  install the packages into the active Python env, then "
                "re-run `scfactory doctor --json`."
            )
        print("\n".join(hint_lines), file=sys.stderr)
        return 2

    # Batch declaration precedence mirrors the module precedence above:
    # explicit flag > recipe (profile or literal) > undeclared.
    if args.batch_strategy and args.batch_strategy != BATCH_STRATEGY_AUTO:
        batch_strategy = args.batch_strategy
        batch_strategy_source = "user-provided --batch-strategy"
    elif recipe_profile is not None:
        batch_strategy = recipe_profile.batch_strategy
        batch_strategy_source = f"profile={recipe_profile.name}"
    elif recipe is not None and recipe.get("batch_strategy"):
        batch_strategy = str(recipe["batch_strategy"])
        batch_strategy_source = f"recipe={recipe['name']}"
    else:
        batch_strategy = BATCH_STRATEGY_AUTO
        batch_strategy_source = "undeclared"

    unknown_modules = sorted(set(planned) - _known_module_names())
    if unknown_modules:
        print(
            f"scfactory: unknown module name(s) {unknown_modules!r} in {plan_source}.\n"
            "  choose names from workflow.modular.module_catalog.MODULE_SPECS; "
            "no dry-run or real run was started.",
            file=sys.stderr,
        )
        return 2

    # Plan-time batch accounting. Runs on the dry-run path too: previewing a
    # plan that will silently produce batch-driven clusters is exactly the
    # moment the operator can still change it. The canonical CLI re-resolves
    # this itself (cheap: obs-only) and owns the manifest record.
    detection = detect_batch_structure(
        resolve_obs_source(
            input_h5ad=input_path if input_path.is_file() else None,
            sample_root=None if input_path.is_file() else input_path,
        ),
        preferred_key="sample",
    )
    try:
        batch_risk = resolve_batch_risk(
            detection=detection,
            declared_strategy=batch_strategy,
            planned_modules=planned,
        )
    except BatchStrategyConflict as exc:
        print(
            f"scfactory: {exc}\n"
            f"  batch strategy came from: {batch_strategy_source}\n"
            "  no dry-run or real run was started.",
            file=sys.stderr,
        )
        return 2
    announce_batch_risk(batch_risk)

    project = args.project or (
        f"scfactory_{recipe['name']}" if recipe else f"scfactory_{modality}"
    )
    project_root = (
        Path(args.project_root).expanduser().resolve()
        if args.project_root
        else None
    )
    if args.run_id and project_root is None:
        print(
            "scfactory: --run-id requires --project-root because legacy output "
            "directories do not use governed run identifiers.",
            file=sys.stderr,
        )
        return 2
    output_dir = str(Path(args.out).expanduser().resolve()) if args.out else str(
        LEGACY_DEFAULT_OUTPUT_DIR
    )

    cli_cmd = [
        sys.executable, "-m", "workflow.modular.cli",
        "--project", project,
    ]
    if input_path.is_file():
        cli_cmd += ["--input-h5ad", str(input_path)]
    else:
        cli_cmd += ["--sample-root", str(input_path)]

    if project_root is not None:
        cli_cmd += ["--project-root", str(project_root)]
        if args.run_id:
            cli_cmd += ["--run-id", args.run_id]
    else:
        cli_cmd += ["--output-dir", output_dir]

    cli_cmd += ["--optional-modules", ",".join(planned)]

    if args.atac_peak_matrix_path:
        cli_cmd += ["--atac-peak-matrix-path", args.atac_peak_matrix_path]
    if args.atac_peaks_bed_path:
        cli_cmd += ["--atac-peaks-bed-path", args.atac_peaks_bed_path]
    if args.atac_n_components is not None:
        cli_cmd += ["--atac-n-components", str(args.atac_n_components)]

    if scatac_selected:
        cli_cmd += [
            "--scatac-da-sample-col", args.scatac_da_sample_col,
            "--scatac-da-group-col", args.scatac_da_group_col,
            "--scatac-da-peak-id-col", args.scatac_da_peak_id_col,
            "--scatac-da-mode", args.scatac_da_mode or "confirmatory_da",
            "--scatac-da-min-samples-per-condition", str(args.scatac_da_min_samples_per_condition or 2),
            "--scatac-da-min-total-count", str(args.scatac_da_min_total_count or 10),
            "--scatac-da-fdr-threshold", str(args.scatac_da_fdr_threshold or 0.05),
            "--scatac-da-abs-log2fc-threshold", str(args.scatac_da_abs_log2fc_threshold or 1.0),
            "--scatac-da-r-conda-env", args.scatac_da_r_conda_env or "r_multiomics",
            "--scatac-da-timeout", str(args.scatac_da_timeout or 1800),
            "--scatac-da-aggregation-backend", args.scatac_da_aggregation_backend or "cpu",
        ]
        if args.scatac_da_condition_col:
            cli_cmd += ["--scatac-da-condition-col", args.scatac_da_condition_col]
        if args.scatac_da_test_level:
            cli_cmd += ["--scatac-da-test-level", args.scatac_da_test_level]
        if args.scatac_da_reference_level:
            cli_cmd += ["--scatac-da-reference-level", args.scatac_da_reference_level]
        if args.scatac_da_groups:
            cli_cmd += ["--scatac-da-groups", args.scatac_da_groups]

    if batch_strategy != BATCH_STRATEGY_AUTO:
        cli_cmd += ["--batch-strategy", batch_strategy]

    if args.scientific_profile:
        cli_cmd += ["--scientific-profile", args.scientific_profile]
    if args.acknowledge_scientific_non_equivalence:
        cli_cmd.append("--acknowledge-scientific-non-equivalence")
    if getattr(args, "allow_dirty", False):
        cli_cmd.append("--allow-dirty")

    # Recipe scale_preset -> --scale-mode passthrough (only if not overridden
    # by a user --optional-modules escape, which is purely about modules).
    scale_preset: str | None = None
    if recipe is not None and recipe.get("scale_preset"):
        scale_preset = str(recipe["scale_preset"])
        cli_cmd += ["--scale-mode", scale_preset]

    # Resolve bundle config: recipe sets defaults, --bundle CLI flag forces on.
    bundle_cfg: dict[str, Any] = {}
    if recipe is not None and isinstance(recipe.get("bundle"), dict):
        bundle_cfg = dict(recipe["bundle"])
    bundle_enabled = bool(bundle_cfg.get("enabled", False)) or bool(args.bundle)
    if bundle_enabled and project_root is not None and not args.run_id:
        print(
            "scfactory: governed bundle export requires an explicit --run-id so "
            "the pipeline and bundle are written to the same run.",
            file=sys.stderr,
        )
        return 2
    if project_root is not None and args.bundle_out:
        print(
            "scfactory: --bundle-out is a legacy-layout option and cannot be "
            "combined with --project-root; governed bundles live under the run.",
            file=sys.stderr,
        )
        return 2

    # Resolve env overrides from recipe (subprocess env, NOT parent env).
    recipe_env: dict[str, str] = {}
    if recipe is not None and isinstance(recipe.get("env"), dict):
        recipe_env = {str(k): str(v) for k, v in recipe["env"].items()}

    print(
        f"scfactory: detected modality={modality} ({plan_source}).\n"
        f"  evidence: {'; '.join(evidence)}\n"
        f"  optional modules: {planned}\n"
        + (f"  recipe: {recipe['name']} ({recipe.get('description', '')})\n"
           if recipe else "")
        + (f"  scale_preset: {scale_preset}\n" if scale_preset else "")
        + (f"  batch strategy: {batch_strategy} ({batch_strategy_source}); "
           f"detected_batches={batch_risk['detected_batches']} "
           f"key={batch_risk['key'] or '-'}; "
           f"claim={batch_risk['clustering_claim_status']}\n")
        + (f"  env overrides: {recipe_env}\n" if recipe_env else "")
        + (f"  bundle: enabled={bundle_enabled} cfg={bundle_cfg}\n"
           if recipe or args.bundle else "")
        + "  to override: --optional-modules a,b,c    "
        "to preview only: --dry-run"
    )
    print(f"scfactory: would execute: {shlex.join(cli_cmd)}")

    if args.dry_run:
        return 0

    # Build subprocess env: parent env + recipe env (recipe wins on conflict).
    sub_env = dict(os.environ)
    if recipe_env:
        sub_env.update(recipe_env)

    proc = subprocess.run(cli_cmd, cwd=str(REPO_ROOT), env=sub_env)
    if proc.returncode != 0:
        return 1

    # Bundle export.
    if bundle_enabled:
        if project_root is not None:
            governed_python_dir = (
                project_root / "runs" / args.run_id / "python"
            )
            candidates = list(governed_python_dir.glob("*/final_adata.h5ad"))
            bundle_out = governed_python_dir / "bundle"
        else:
            legacy_root = Path(output_dir)
            candidates = [
                path
                for path in (
                    legacy_root / project / "final_adata.h5ad",
                    *legacy_root.glob(f"{project}_*/final_adata.h5ad"),
                )
                if path.is_file()
            ]
            bundle_out = (
                Path(args.bundle_out).resolve()
                if args.bundle_out
                else legacy_root / project / "r_bundle"
            )

        if not candidates:
            print(
                "scfactory: bundle requested but no final_adata.h5ad was produced "
                "in the selected run; bundle export cannot proceed.",
                file=sys.stderr,
            )
            return 1
        candidate_h5ad = max(candidates, key=lambda path: path.stat().st_mtime)
        bundle_cmd: list[str] = [
            sys.executable, str(EXPORT_SCRIPT),
            "--input", str(candidate_h5ad),
        ]
        if project_root is not None:
            bundle_cmd += [
                "--project-root", str(project_root),
                "--run-id", args.run_id,
            ]
        else:
            bundle_cmd += ["--output", str(bundle_out)]
        # Bundle protein/spatial/multimodal flags: recipe wins, modality is fallback.
        include_protein = bundle_cfg.get(
            "include_protein", modality in ("cite_seq", "multimodal")
        )
        include_spatial = bundle_cfg.get(
            "include_spatial", modality in ("spatial", "multimodal")
        )
        include_multimodal_obsm = bundle_cfg.get("include_multimodal_obsm", False)
        if include_protein:
            bundle_cmd.append("--include-protein")
        if include_spatial:
            bundle_cmd.append("--include-spatial")
        if include_multimodal_obsm:
            bundle_cmd.append("--include-multimodal-obsm")
        print(f"scfactory: bundling -> {bundle_out}")
        print(f"scfactory: would execute: {shlex.join(bundle_cmd)}")
        b_proc = subprocess.run(bundle_cmd, cwd=str(REPO_ROOT), env=sub_env)
        if b_proc.returncode != 0:
            return 1

    return 0


# ---------------------------------------------------------------------------
# `doctor` subcommand
# ---------------------------------------------------------------------------

def _resolve_rscript() -> tuple[str | None, str]:
    """Mirror tests/conftest.py rscript resolution. Returns (path, source)."""
    env_override = os.environ.get("RSCRIPT_BIN")
    if env_override:
        if os.path.isfile(env_override) and os.access(env_override, os.X_OK):
            return env_override, "RSCRIPT_BIN env"
        return None, f"RSCRIPT_BIN set but not executable: {env_override}"
    candidates = [
        Path.home() / "conda/envs/r_multiomics_arrow/bin/Rscript",
        Path.home() / "conda/envs/r_multiomics/bin/Rscript",
    ]
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return str(c), f"conda env: {c.parent.parent.name}"
    which = shutil.which("Rscript")
    if which:
        return which, "PATH"
    return None, "not found"


ENVIRONMENT_LOCK_PATH = REPO_ROOT / "environment.yml"
# The R-side lock lives in the sibling factory repo (this repo has no R
# dependency file of its own). README.md's "installed and validated" list and
# PROTOCOL.md's Tier-3 table both point at that renv-managed environment, so
# it — not a hardcoded guess — is the source of truth for which R packages
# are actually declared.
R_LOCK_PATH = REPO_ROOT.parent / "r_multiomics_factory" / "renv.lock"


def _read_environment_lock(path: Path = ENVIRONMENT_LOCK_PATH) -> dict[str, Any]:
    """Best-effort parse of environment.yml: declared python pin + package names.

    Deliberately NOT a full YAML parse — doctor stays dependency-light and
    must not gain a hard PyYAML requirement just to run. Only extracts what
    doctor needs: the flat list of `dependencies:` entries (conda-level and
    nested `pip:` entries alike) and the `python=X.Y` pin. Read-only; never
    raises. A missing/unparseable lock degrades to ``found=False`` so callers
    can surface an honest "could not verify" warning instead of a false pass.
    """
    result: dict[str, Any] = {
        "path": str(path),
        "found": False,
        "python_version": None,  # (major, minor) or None
        "declared_packages": set(),
    }
    if not path.is_file():
        return result
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return result
    result["found"] = True
    in_deps = False
    for raw_line in lines:
        stripped = raw_line.strip()
        if not in_deps:
            if stripped == "dependencies:":
                in_deps = True
            continue
        if not stripped.startswith("-"):
            continue
        entry = stripped[1:].strip().split("#", 1)[0].strip().rstrip(":")
        if not entry:
            continue
        if entry.startswith("python"):
            m = re.search(r"python\s*=\s*(\d+)\.(\d+)", entry)
            if m:
                result["python_version"] = (int(m.group(1)), int(m.group(2)))
            continue
        name = re.split(r"[=<>! ]", entry, maxsplit=1)[0].strip().lower()
        if name:
            result["declared_packages"].add(name)
    return result


def _read_r_lock(path: Path = R_LOCK_PATH) -> dict[str, Any]:
    """Best-effort parse of the sibling r_multiomics_factory renv.lock.

    Read-only; absence (e.g. the sibling repo is not checked out on this
    host) degrades to ``found=False`` rather than raising — doctor must never
    crash the host it is diagnosing.
    """
    result: dict[str, Any] = {
        "path": str(path),
        "found": False,
        "r_version": None,
        "declared_packages": set(),
    }
    if not path.is_file():
        return result
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return result
    result["found"] = True
    r_info = data.get("R")
    if isinstance(r_info, dict):
        ver = r_info.get("Version")
        if isinstance(ver, str):
            result["r_version"] = ver
    packages = data.get("Packages")
    if isinstance(packages, dict):
        result["declared_packages"] = {str(k) for k in packages}
    return result


def _check_python() -> dict[str, Any]:
    v = sys.version_info
    status = "pass" if (v.major, v.minor) >= (3, 10) else "warn"
    return {
        "status": status,
        "version": f"{v.major}.{v.minor}.{v.micro}",
        "message": (
            f"Python {v.major}.{v.minor}.{v.micro}"
            + ("" if status == "pass" else " (recommend >= 3.10)")
        ),
    }


def _check_conda_envs() -> dict[str, Any]:
    envs = {}
    home = Path.home()
    for name in ("r_multiomics_arrow", "r_multiomics"):
        rscript = home / "conda" / "envs" / name / "bin" / "Rscript"
        envs[name] = {"present": rscript.is_file(), "path": str(rscript)}
    any_present = any(e["present"] for e in envs.values())
    return {
        "status": "pass" if any_present else "warn",
        "envs": envs,
        "message": (
            "conda envs: "
            + ", ".join(f"{n}={'ok' if e['present'] else 'missing'}"
                        for n, e in envs.items())
        ),
    }


def _check_rscript() -> dict[str, Any]:
    path, source = _resolve_rscript()
    if path is None:
        return {"status": "warn", "path": None, "source": source,
                "message": f"Rscript: not available ({source})"}
    return {"status": "pass", "path": path, "source": source,
            "message": f"Rscript: {path} ({source})"}


def _check_bridges() -> dict[str, Any]:
    bridge_dir = REPO_ROOT / "bridges" / "local_r_pipeline_macbook"
    results: dict[str, Any] = {}
    overall = "pass"
    for name, expected_target_substr in (
        ("R", "r_multiomics_factory/R"),
        ("R_bundle", "r_multiomics_factory/R_bundle"),
    ):
        link = bridge_dir / name
        info: dict[str, Any] = {"path": str(link)}
        if not link.exists() and not link.is_symlink():
            info["status"] = "fail"
            info["message"] = f"missing: {link}"
            overall = "fail"
        elif not link.is_symlink():
            info["status"] = "fail"
            info["message"] = f"{link} is not a symlink (must be)"
            overall = "fail"
        else:
            target = os.readlink(link)
            resolved = link.resolve()
            info["target"] = target
            info["resolved"] = str(resolved)
            if expected_target_substr in target or expected_target_substr in str(resolved):
                info["status"] = "pass"
                info["message"] = f"{name} -> {target}"
            else:
                info["status"] = "fail"
                info["message"] = (
                    f"{name} symlink target unexpected: {target} "
                    f"(expected substring '{expected_target_substr}')"
                )
                overall = "fail"
        results[name] = info
    return {"status": overall, "links": results,
            "message": "; ".join(v["message"] for v in results.values())}


def _check_python_deps() -> dict[str, Any]:
    """Probe optional/scientific python deps and classify each pass/warn.

    Each dep is tiered "required" or "optional" by looking it up in
    environment.yml's own declared dependency list (see
    ``_read_environment_lock``) rather than a hardcoded guess: scanpy/anndata/
    scipy are pinned there, so a host missing them is materially broken;
    pyarrow/fastparquet/squidpy/muon/mofapy2 are NOT declared there — they
    back genuinely optional, module-gated capabilities (bundle parquet
    export, spatial_neighborhoods, multimodal_integration MOFA) that clean-
    skip at runtime when absent. A missing dep — required or optional — is
    always individually visible and never silently collapsed into "pass";
    this section's own status is "warn" whenever anything is missing so a
    wrong environment cannot summarize as all-green. It still never escalates
    to "fail": optional capabilities are legitimately absent on many hosts,
    and the fail-worthy distinction (required vs optional) is surfaced
    separately in `readiness.core_environment_ready` /
    `readiness.optional_capabilities`.
    """
    deps = ("scanpy", "anndata", "pyarrow", "fastparquet", "scipy",
            "squidpy", "muon", "mofapy2")
    lock = _read_environment_lock()
    declared = lock["declared_packages"]
    # Listed in environment.yml for install discoverability but recipe-gated
    # (Wave-2 spatial preflight; Wave-4 env pin). Always capability-optional
    # for doctor core rollup even when declared.
    capability_optional = frozenset({"squidpy"})

    out: dict[str, Any] = {}
    missing_required: list[str] = []
    missing_optional: list[str] = []
    import importlib
    import importlib.metadata as _imd
    for name in deps:
        if name.lower() in capability_optional:
            tier = "optional"
        else:
            tier = "required" if name.lower() in declared else "optional"
        try:
            importlib.import_module(name)
            present = True
        except Exception:
            present = False
        ver: str | None = None
        if present:
            try:
                ver = _imd.version(name)
            except Exception:
                ver = "unknown"
        else:
            (missing_required if tier == "required" else missing_optional).append(name)
        out[name] = {
            "present": present,
            "version": ver,
            "tier": tier,
            "status": "pass" if present else "warn",
        }

    any_missing = bool(missing_required or missing_optional)
    return {
        "status": "warn" if any_missing else "pass",
        "deps": out,
        "missing_required": missing_required,
        "missing_optional": missing_optional,
        "lock_source": lock["path"] if lock["found"] else None,
        "message": "python deps: "
                   + ", ".join(
                       f"{n}=ok" if v["present"] else f"{n}=MISSING({v['tier']})"
                       for n, v in out.items()
                   ),
    }


def _check_r_packages(rscript_path: str | None) -> dict[str, Any]:
    """Probe R bridge/bundle packages and classify each pass/warn.

    Tiering comes from the sibling repo's renv.lock (see ``_read_r_lock``):
    packages it declares (Seurat, arrow, jsonlite, Matrix, zellkonverter,
    harmony) are required for the R-bundle/bridge lane; SeuratDisk is
    deliberately absent from that lock — README.md documents it as
    intentionally not installed (conda-forge SeuratDisk currently conflicts
    with R 4.5 / zellkonverter via old spatstat requirements; zellkonverter/
    hdf5r cover the same h5ad-bridge need) — so it is tiered "optional", not
    a defect. When renv.lock cannot be read at all (e.g. the sibling repo
    isn't checked out on this host), tiering degrades to "unknown" rather
    than guessing.
    """
    pkgs = ("Seurat", "arrow", "jsonlite", "Matrix", "zellkonverter",
            "SeuratDisk", "harmony")
    r_lock = _read_r_lock()
    declared = r_lock["declared_packages"]

    def _tier(pkg: str) -> str:
        if not r_lock["found"]:
            return "unknown"
        return "required" if pkg in declared else "optional"

    lock_source = r_lock["path"] if r_lock["found"] else None

    if rscript_path is None:
        packages = {
            p: {"present": False, "version": None, "tier": _tier(p), "status": "warn"}
            for p in pkgs
        }
        return {
            "status": "warn",
            "packages": packages,
            "missing_required": [p for p, v in packages.items() if v["tier"] == "required"],
            "missing_optional": [p for p, v in packages.items() if v["tier"] != "required"],
            "lock_source": lock_source,
            "message": "R: not available; skipping package probe",
        }

    pkg_vec = ",".join(f"'{p}'" for p in pkgs)
    code = (
        f"for (p in c({pkg_vec})) "
        "cat(sprintf('%s|%s\\n', p, "
        "if (requireNamespace(p, quietly=TRUE)) "
        "as.character(packageVersion(p)) else 'MISSING'))"
    )
    try:
        proc = subprocess.run(
            [rscript_path, "-e", code],
            capture_output=True, text=True, timeout=5,
        )
    except subprocess.TimeoutExpired:
        return {"status": "warn", "packages": {},
                "missing_required": [], "missing_optional": [],
                "lock_source": lock_source,
                "message": "R package probe timed out (>5s)"}
    except Exception as exc:
        return {"status": "warn", "packages": {},
                "missing_required": [], "missing_optional": [],
                "lock_source": lock_source,
                "message": f"R package probe failed: {exc}"}

    parsed: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if "|" in line:
            name, ver = line.split("|", 1)
            parsed[name.strip()] = ver.strip()

    packages: dict[str, Any] = {}
    missing_required: list[str] = []
    missing_optional: list[str] = []
    for p in pkgs:
        ver = parsed.get(p, "MISSING")
        present = ver != "MISSING"
        tier = _tier(p)
        if not present:
            (missing_required if tier == "required" else missing_optional).append(p)
        packages[p] = {
            "present": present,
            "version": ver if present else None,
            "tier": tier,
            "status": "pass" if present else "warn",
        }

    any_missing = bool(missing_required or missing_optional)
    return {
        "status": "warn" if any_missing else "pass",
        "packages": packages,
        "missing_required": missing_required,
        "missing_optional": missing_optional,
        "lock_source": lock_source,
        "message": "R packages: "
                   + ", ".join(
                       f"{n}={v['version']}" if v["present"] else f"{n}=MISSING({v['tier']})"
                       for n, v in packages.items()
                   ),
    }


def _check_last_run(project_roots: list[Path] | None = None) -> dict[str, Any]:
    """Inspect governed run manifests plus the deprecated local layout.

    Explicit project roots keep the search bounded to paths named by the user.
    Without one, doctor checks each immediate project under ``~/projects``.
    Producer-native ``run_manifest.json`` files remain a fallback, while the
    cross-factory run-root ``manifest.json`` is the governed source of truth.
    """

    candidates: set[Path] = set()
    if project_roots:
        governed_roots = project_roots
        for root in governed_roots:
            candidates.update(root.glob("runs/*/manifest.json"))
            candidates.update(root.glob("runs/*/python/**/run_manifest.json"))
    elif DEFAULT_PROJECTS_ROOT.is_dir():
        candidates.update(DEFAULT_PROJECTS_ROOT.glob("*/runs/*/manifest.json"))
        candidates.update(
            DEFAULT_PROJECTS_ROOT.glob("*/runs/*/python/**/run_manifest.json")
        )

    for sub in ("results", "output"):
        d = REPO_ROOT / sub
        if d.is_dir():
            candidates.update(d.glob("*/run_manifest.json"))
    if not candidates:
        return {
            "status": "warn",
            "found": False,
            "message": (
                "no governed manifest.json found under project runs and no "
                "legacy run_manifest.json found in results/ or output/"
            ),
        }
    newest = max(candidates, key=lambda p: p.stat().st_mtime)
    info: dict[str, Any] = {"path": str(newest)}
    try:
        data = json.loads(newest.read_text())
        modules = (
            data.get("completed_modules")
            or data.get("modules_run")
            or data.get("modules")
            or data.get("module_status")
            or []
        )
        info["module_count"] = len(modules) if hasattr(modules, "__len__") else 0
        info["status_field"] = data.get(
            "overall_status", data.get("status", "unknown")
        )
        info["layout"] = (
            "governed" if newest.name == "manifest.json" else "producer_or_legacy"
        )
        return {"status": "pass", "found": True, "info": info,
                "message": (f"last run: {newest.parent.name} "
                            f"(modules={info['module_count']}, "
                            f"status={info['status_field']})")}
    except Exception as exc:
        return {"status": "warn", "found": True, "info": info,
                "message": f"could not parse {newest}: {exc}"}


# The claim-critical dependency contract lives at the SUITE root (one level
# above every factory repo, including this one), not inside any single
# factory. It records dependencies whose presence is what separates a
# scientific result the suite will let you CLAIM from one it marks
# exploratory (see the contract file's own rationale: the 2026-08-02 audit
# found pydeseq2/pertpy/multiHiCcompare installed but declared in NO
# committed spec, and pertpy was additionally installed-but-unimportable for
# an unknown period, silently downgrading every scCODA composition analysis
# to a compositionally-invalid fallback). This is a SEPARATE, higher-severity
# tier above the environment.yml/renv.lock-derived "required" tier used by
# `_check_python_deps` / `_check_r_packages`.
CLAIM_CRITICAL_CONTRACT_PATH = REPO_ROOT.parent / "contracts" / "claim_critical_dependencies.yaml"
CONDA_ENVS_ROOT = Path.home() / "conda" / "envs"
# sc_gpu-hosted claim-critical deps (pydeseq2, pertpy) transitively import
# jax/numpyro/torch; a cold-cache import can legitimately take well past the
# 5s budget used for the lightweight R package probe elsewhere in this file.
CLAIM_CRITICAL_PROBE_TIMEOUT = 180


def _load_claim_critical_contract(
    path: Path = CLAIM_CRITICAL_CONTRACT_PATH,
) -> tuple[list[dict[str, Any]] | None, str | None]:
    """Best-effort load of the suite-root claim-critical dependency contract.

    Returns ``(dependencies, error)``. ``dependencies`` is None when the
    contract cannot be read at all — this repo cloned standalone outside the
    suite, PyYAML unavailable, or the file malformed — so callers can surface
    an explicit "contract unavailable" state instead of silently reporting
    healthy. This does not replace the dedicated suite gate at
    ``scripts/validate_claim_critical_deps.py`` (which enforces the contract's
    own structure — C1-C4: schema completeness, declared_in actually
    mentioning the dep, guards file existing); doctor only asks the
    environment-readiness question (that gate's C5), escalated to failure
    severity per the 2026-08-02 pertpy incident (see module docstring above).
    """
    if not path.is_file():
        return None, (
            f"contract not found at {path} "
            "(suite root not checked out alongside this repo?)"
        )
    try:
        import yaml  # type: ignore
    except ImportError:
        return None, "PyYAML not installed; cannot parse the claim-critical dependency contract"
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return None, f"could not parse {path}: {exc}"
    if not isinstance(data, dict):
        return None, f"{path} did not parse to a mapping"
    deps = data.get("dependencies")
    if not isinstance(deps, list) or not deps:
        return None, f"{path} declares no dependencies"
    return deps, None


def _probe_python_env_imports(env_prefix: Path, import_names: list[str]) -> dict[str, str | None]:
    """Actually IMPORT each name in one subprocess under ``env_prefix``'s python.

    Deliberately a real ``import <name>``, not merely
    ``importlib.metadata.version()``: the 2026-08-02 pertpy incident was
    exactly a package whose metadata was present but whose import raised
    (jax removed ``xla_pmap_p``, numpyro 0.20.1 still referenced it) — a
    metadata-only probe would have kept reporting healthy throughout that
    incident. Batches all names for one env into a single subprocess launch
    (mirrors ``_check_r_packages`` batching multiple R packages into one
    Rscript call). Returns ``{import_name: version-or-None}``; None means
    "not importable" (or the probe subprocess itself could not run) — this
    function does not distinguish the two, since either way the dependency
    cannot be relied on from this host right now.
    """
    python = env_prefix / "bin" / "python"
    if not python.is_file():
        return {n: None for n in import_names}
    safe_names = [n for n in import_names if n.isidentifier()]
    lines = ["import json", "out = {}"]
    for n in safe_names:
        lines += [
            "try:",
            f"    import {n} as _mod",
            "    try:",
            "        import importlib.metadata as _md",
            f"        out[{n!r}] = _md.version({n!r})",
            "    except Exception:",
            f"        out[{n!r}] = getattr(_mod, '__version__', 'unknown')",
            "except Exception:",
            f"    out[{n!r}] = None",
        ]
    lines.append("print(json.dumps(out))")
    try:
        proc = subprocess.run(
            [str(python), "-c", "\n".join(lines)],
            capture_output=True, text=True, timeout=CLAIM_CRITICAL_PROBE_TIMEOUT,
        )
    except (subprocess.TimeoutExpired, OSError):
        return {n: None for n in import_names}
    out_lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    try:
        parsed = json.loads(out_lines[-1]) if out_lines else {}
    except Exception:
        parsed = {}
    return {n: parsed.get(n) for n in import_names}


def _probe_r_env_package(env_prefix: Path, r_package: str) -> str | None:
    """``requireNamespace()``-based probe for one R package under ``env_prefix``.

    Like ``_check_r_packages``, ``requireNamespace`` actually attempts to
    load the package's namespace rather than only consulting installed
    package metadata.
    """
    rscript = env_prefix / "bin" / "Rscript"
    if not rscript.is_file():
        return None
    code = (
        f'if (requireNamespace("{r_package}", quietly=TRUE)) '
        f'cat(as.character(packageVersion("{r_package}"))) else cat("ABSENT")'
    )
    try:
        proc = subprocess.run(
            [str(rscript), "-e", code],
            capture_output=True, text=True, timeout=CLAIM_CRITICAL_PROBE_TIMEOUT,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None
    value = proc.stdout.strip()
    return None if (not value or value == "ABSENT") else value


def _check_claim_critical_deps(envs_root: Path = CONDA_ENVS_ROOT) -> dict[str, Any]:
    """Environment-readiness check for suite-declared claim-critical deps.

    A missing/unimportable claim-critical dependency escalates this check's
    status all the way to "fail" — the one place in doctor where a missing
    dep is allowed to do that. Every other dependency tier in this file tops
    out at "warn" because optional/required-but-absent capabilities are
    still an honest, loud degradation the pipeline handles gracefully; a
    claim-critical gap is different in kind — it silently converts a
    confirmatory scientific result into an unclaimable one (or, in the
    pertpy case, produced compositionally-invalid statistics
    indistinguishable from the valid path for an unknown period). Doctor
    should be loud here even though it stays quiet elsewhere.

    When the target conda env is not present on this host at all, the
    per-dependency status is "warn", not "fail": that is an honest
    "unverifiable from here" rather than a confirmed gap (mirrors
    ``validate_claim_critical_deps.py``'s ``env_unavailable`` info finding).
    """
    deps, load_error = _load_claim_critical_contract()
    if deps is None:
        return {
            "status": "warn",
            "contract_available": False,
            "contract_path": str(CLAIM_CRITICAL_CONTRACT_PATH),
            "deps": {},
            "message": f"claim-critical dependency contract unavailable: {load_error}",
        }

    by_env_python: dict[str, list[dict[str, Any]]] = {}
    r_deps: list[dict[str, Any]] = []
    for dep in deps:
        if not isinstance(dep, dict) or not dep.get("name") or not dep.get("env"):
            continue
        if dep.get("r_package"):
            r_deps.append(dep)
        else:
            by_env_python.setdefault(dep["env"], []).append(dep)

    results: dict[str, Any] = {}

    for env_name, env_deps in by_env_python.items():
        env_prefix = envs_root / env_name
        env_present = env_prefix.is_dir()
        import_names = [d.get("import_name", d["name"]) for d in env_deps]
        versions = _probe_python_env_imports(env_prefix, import_names) if env_present else {}
        for dep in env_deps:
            iname = dep.get("import_name", dep["name"])
            declared_version = str(dep.get("version") or "") or None
            if not env_present:
                results[dep["name"]] = {
                    "env": env_name, "verifiable": False, "present": None,
                    "version": None, "declared_version": declared_version,
                    "status": "warn",
                    "message": (
                        f"env '{env_name}' not present on this host; claim-critical "
                        f"status of '{dep['name']}' cannot be verified from here"
                    ),
                }
                continue
            installed_version = versions.get(iname)
            present = installed_version is not None
            results[dep["name"]] = {
                "env": env_name, "verifiable": True, "present": present,
                "version": installed_version, "declared_version": declared_version,
                "status": "pass" if present else "fail",
                "message": (
                    f"{dep['name']} {installed_version} importable in env '{env_name}'"
                    if present else
                    f"{dep['name']} NOT importable in env '{env_name}' — "
                    f"{dep.get('without_it', '').strip()[:160]}"
                ),
            }

    for dep in r_deps:
        env_name = dep["env"]
        env_prefix = envs_root / env_name
        env_present = env_prefix.is_dir()
        declared_version = str(dep.get("version") or "") or None
        if not env_present:
            results[dep["name"]] = {
                "env": env_name, "verifiable": False, "present": None,
                "version": None, "declared_version": declared_version,
                "status": "warn",
                "message": (
                    f"env '{env_name}' not present on this host; claim-critical "
                    f"status of '{dep['name']}' cannot be verified from here"
                ),
            }
            continue
        installed_version = _probe_r_env_package(env_prefix, dep["name"])
        present = installed_version is not None
        results[dep["name"]] = {
            "env": env_name, "verifiable": True, "present": present,
            "version": installed_version, "declared_version": declared_version,
            "status": "pass" if present else "fail",
            "message": (
                f"{dep['name']} {installed_version} loadable in env '{env_name}'"
                if present else
                f"{dep['name']} NOT loadable in env '{env_name}' — "
                f"{dep.get('without_it', '').strip()[:160]}"
            ),
        }

    if any(v["status"] == "fail" for v in results.values()):
        overall = "fail"
    elif any(v["status"] == "warn" for v in results.values()):
        overall = "warn"
    else:
        overall = "pass"

    return {
        "status": overall,
        "contract_available": True,
        "contract_path": str(CLAIM_CRITICAL_CONTRACT_PATH),
        "deps": results,
        "message": "claim-critical deps: " + ", ".join(
            f"{n}=ok({v['version']})" if v["present"]
            else f"{n}={'UNVERIFIABLE' if not v['verifiable'] else 'MISSING'}"
            for n, v in results.items()
        ),
    }


def _check_lock_identity(rscript_path: str | None) -> dict[str, Any]:
    """Compare the RUNNING interpreter/R against the repo's declared locks.

    Distinct from ``_check_python``'s ``>= 3.10`` floor check: that check
    passes for any modern python. This check compares against the ACTUAL
    pin recorded in environment.yml (``python=3.11``) and, when the sibling
    r_multiomics_factory renv.lock is reachable, the R version it records. A
    version-mismatched interpreter can pass the floor check while running
    code the lock was never validated against — that mismatch is invisible
    unless something explicitly diffs against the declared lock, which is
    what this check exists to do.
    """
    py_lock = _read_environment_lock()
    running_py = (sys.version_info.major, sys.version_info.minor)
    py_declared = py_lock["python_version"]
    py_match: bool | None = (running_py == py_declared) if py_declared else None

    r_lock = _read_r_lock()
    r_declared = r_lock["r_version"]
    r_running: str | None = None
    r_match: bool | None = None
    if rscript_path is not None and r_declared is not None:
        try:
            proc = subprocess.run(
                [rscript_path, "-e", "cat(as.character(getRversion()))"],
                capture_output=True, text=True, timeout=5,
            )
            r_running = proc.stdout.strip() or None
        except Exception:
            r_running = None
        if r_running:
            # Compare major.minor only; R patch releases are routinely mixed
            # across hosts without breaking package ABI.
            r_match = tuple(r_running.split(".")[:2]) == tuple(r_declared.split(".")[:2])

    unverifiable = py_declared is None or (rscript_path is not None and r_lock["found"] and r_declared is not None and r_running is None)
    mismatched = (py_match is False) or (r_match is False)
    status = "warn" if (mismatched or unverifiable) else "pass"

    parts: list[str] = []
    if py_declared is None:
        parts.append(
            f"python lock unverifiable ({py_lock['path']} "
            + ("not found" if not py_lock["found"] else "no python= pin found")
            + ")"
        )
    else:
        parts.append(
            f"python: running {running_py[0]}.{running_py[1]} vs locked "
            f"{py_declared[0]}.{py_declared[1]} "
            f"({'match' if py_match else 'MISMATCH'})"
        )
    if r_declared is None:
        parts.append(
            f"R lock unverifiable ({r_lock['path']} "
            + ("not found" if not r_lock["found"] else "no R version recorded")
            + ")"
        )
    elif r_running is None:
        parts.append(f"R: locked {r_declared}, running version could not be probed")
    else:
        parts.append(
            f"R: running {r_running} vs locked {r_declared} "
            f"({'match' if r_match else 'MISMATCH'})"
        )

    return {
        "status": status,
        "python": {
            "running": f"{running_py[0]}.{running_py[1]}",
            "declared": (f"{py_declared[0]}.{py_declared[1]}" if py_declared else None),
            "lock_path": py_lock["path"],
            "match": py_match,
        },
        "r": {
            "running": r_running,
            "declared": r_declared,
            "lock_path": r_lock["path"],
            "match": r_match,
        },
        "message": "; ".join(parts),
    }


def _rollup_readiness(report: dict[str, Any]) -> dict[str, Any]:
    """Roll the atomic checks up into legible, decision-relevant questions.

    These are DERIVED views for human/automation legibility — they
    deliberately do NOT feed the top-level pass/warn/fail tally (their
    inputs already do, via ``python``, ``python_deps``, ``rscript``,
    ``r_packages``, ``lock_identity_match``); counting them again would
    double-count the same underlying signal.
    """
    python_deps = report["python_deps"]["deps"]
    r_packages = report["r_packages"]["packages"]

    core_missing = [
        n for n, v in python_deps.items() if v["tier"] == "required" and not v["present"]
    ]
    core_ok = report["python"]["status"] == "pass" and not core_missing
    core_reasons: list[str] = []
    if report["python"]["status"] != "pass":
        core_reasons.append(report["python"]["message"])
    if core_missing:
        core_reasons.append(f"missing required deps (per environment.yml): {core_missing}")
    core = {
        "status": "pass" if core_ok else "warn",
        "message": (
            "core python interpreter + required deps (scanpy/anndata/scipy "
            "per environment.yml) are usable"
            if core_ok else
            "core environment is degraded: " + "; ".join(core_reasons)
        ),
    }

    bundle_missing_r = [
        n for n, v in r_packages.items() if v["tier"] == "required" and not v["present"]
    ]
    have_parquet_writer = (
        python_deps.get("pyarrow", {}).get("present", False)
        or python_deps.get("fastparquet", {}).get("present", False)
    )
    rscript_ok = report["rscript"]["status"] == "pass"
    profile_ok = rscript_ok and not bundle_missing_r and have_parquet_writer
    profile_reasons: list[str] = []
    if not rscript_ok:
        profile_reasons.append("Rscript not resolved")
    if bundle_missing_r:
        profile_reasons.append(f"missing required R packages (per renv.lock): {bundle_missing_r}")
    if not have_parquet_writer:
        profile_reasons.append("no parquet writer (pyarrow or fastparquet) available")
    _profile_heuristic_note = (
        " (heuristic: this is NOT tied to whatever --recipe/profile you are "
        "about to run — doctor has no notion of an active profile today and "
        "reports readiness of the R-bundle/bridge lane, the most commonly "
        "exercised non-default lane in this repo, as a stand-in)"
    )
    profile = {
        "status": "pass" if profile_ok else "warn",
        "rscript_ok": rscript_ok,
        "missing_r_packages": bundle_missing_r,
        "parquet_writer_available": have_parquet_writer,
        "message": (
            "R-bundle export lane (--bundle / bridge_ready modules) is ready"
            + _profile_heuristic_note
            if profile_ok else
            "R-bundle export lane is degraded: " + "; ".join(profile_reasons)
            + _profile_heuristic_note
        ),
    }

    optional_missing = [
        n for n, v in python_deps.items() if v["tier"] == "optional" and not v["present"]
    ] + [
        n for n, v in r_packages.items() if v["tier"] != "required" and not v["present"]
    ]
    optional = {
        "status": "warn" if optional_missing else "pass",
        "missing": optional_missing,
        "message": (
            f"optional capabilities not installed (legitimately absent unless "
            f"you need them): {optional_missing}"
            if optional_missing else
            "all probed optional capabilities are present"
        ),
    }

    # Wave-2: explicit spatial analytics readiness (squidpy).
    squidpy_present = bool(python_deps.get("squidpy", {}).get("present"))
    spatial = {
        "status": "pass" if squidpy_present else "warn",
        "package": "squidpy",
        "present": squidpy_present,
        "required_by": [
            "spatial_neighborhoods",
            "recipe:visium_neighborhoods",
        ],
        "message": (
            "squidpy present — spatial_neighborhoods / visium_neighborhoods ready"
            if squidpy_present else
            "squidpy MISSING — spatial_neighborhoods and recipe "
            "visium_neighborhoods will fail preflight (install squidpy, then "
            "re-check). Other scRNA modules are unaffected."
        ),
    }

    claim_critical = report["claim_critical_deps"]
    if not claim_critical.get("contract_available", False):
        cc_message = claim_critical["message"]
    elif claim_critical["status"] == "pass":
        cc_message = (
            "all suite-declared claim-critical dependencies import cleanly "
            "in their declared envs"
        )
    else:
        cc_deps = claim_critical["deps"]
        failing = [n for n, v in cc_deps.items() if v["status"] == "fail"]
        unverifiable = [n for n, v in cc_deps.items() if v["status"] == "warn"]
        parts = []
        if failing:
            parts.append(f"NOT importable, claim path broken: {failing}")
        if unverifiable:
            parts.append(f"unverifiable from this host (env absent): {unverifiable}")
        cc_message = "claim-critical dependency gap — " + "; ".join(parts)
    claim_critical_ready = {
        "status": claim_critical["status"],
        "contract_available": claim_critical.get("contract_available", False),
        "message": cc_message,
    }

    return {
        "core_environment_ready": core,
        "selected_profile_ready": profile,
        "optional_capabilities": optional,
        "spatial_analytics": spatial,
        "claim_critical_ready": claim_critical_ready,
    }


def cmd_doctor(args: argparse.Namespace) -> int:
    rscript = _check_rscript()
    project_roots = [
        Path(value).expanduser().resolve() for value in args.project_root
    ]
    report = {
        "python": _check_python(),
        "conda_envs": _check_conda_envs(),
        "rscript": rscript,
        "bridge_symlinks": _check_bridges(),
        "python_deps": _check_python_deps(),
        "r_packages": _check_r_packages(rscript.get("path")),
        "last_run": _check_last_run(project_roots or None),
        "lock_identity_match": _check_lock_identity(rscript.get("path")),
        "claim_critical_deps": _check_claim_critical_deps(),
    }

    # `counts` intentionally tallies only the atomic checks above — the
    # readiness rollup added below is a derived view over the same signals
    # and is deliberately excluded so nothing is double-counted.
    counts = {"pass": 0, "warn": 0, "fail": 0}
    for v in report.values():
        s = v.get("status", "warn")
        counts[s] = counts.get(s, 0) + 1
    report["summary"] = counts
    report["readiness"] = _rollup_readiness(report)

    if args.json:
        print(json.dumps(report, indent=2, default=str))
    else:
        print("scfactory doctor — read-only environment health check")
        print("=" * 56)
        for key in ("python", "conda_envs", "rscript", "bridge_symlinks",
                    "python_deps", "r_packages", "last_run", "lock_identity_match",
                    "claim_critical_deps"):
            v = report[key]
            tag = v["status"].upper()
            print(f"[{tag:4}] {key}: {v['message']}")
        print("-" * 56)
        print("readiness dimensions:")
        for key, v in report["readiness"].items():
            tag = v["status"].upper()
            print(f"  [{tag:4}] {key}: {v['message']}")
        print("-" * 56)
        print(f"summary: pass={counts['pass']} warn={counts['warn']} "
              f"fail={counts['fail']}")
        if counts["fail"]:
            print("fix hints:")
            for k, v in report.items():
                if isinstance(v, dict) and v.get("status") == "fail":
                    print(f"  - {k}: {v['message']}")

    return 1 if counts["fail"] else 0


# ---------------------------------------------------------------------------
# `report` subcommand
# ---------------------------------------------------------------------------

def _find_claim_guard(manifest: dict[str, Any], run_dir: Path) -> str | None:
    """Locate a claim_guard string in the manifest or sibling bundle.

    Search order:
      1. ``manifest["expression"]["claim_guard"]`` (canonical bundle layout)
      2. ``manifest["claim_guard"]`` (top-level shorthand)
      3. ``manifest["bundle"]["claim_guard"]``
      4. Any sibling ``r_bundle/manifest.json`` under ``run_dir``
         (looks for ``expression.claim_guard``).
    """
    expr = manifest.get("expression")
    if isinstance(expr, dict):
        cg = expr.get("claim_guard")
        if isinstance(cg, str) and cg:
            return cg
    cg = manifest.get("claim_guard")
    if isinstance(cg, str) and cg:
        return cg
    bundle = manifest.get("bundle")
    if isinstance(bundle, dict):
        cg = bundle.get("claim_guard")
        if isinstance(cg, str) and cg:
            return cg
    # Fallback: probe a colocated bundle manifest.
    bundle_manifest = run_dir / "r_bundle" / "manifest.json"
    if bundle_manifest.is_file():
        try:
            data = json.loads(bundle_manifest.read_text())
            expr = data.get("expression")
            if isinstance(expr, dict):
                cg = expr.get("claim_guard")
                if isinstance(cg, str) and cg:
                    return cg
        except Exception:
            pass
    return None


def _read_module_status_csv(path: Path) -> list[dict[str, str]]:
    """Parse module_status.csv. Returns [] if missing or unreadable."""
    if not path.is_file():
        return []
    rows: list[dict[str, str]] = []
    try:
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            for r in reader:
                rows.append({k: (v or "") for k, v in r.items()})
    except Exception:
        return []
    return rows


def _walk_module_assets(
    run_dir: Path, max_depth: int = REPORT_WALK_MAX_DEPTH
) -> dict[str, dict[str, list[Path]]]:
    """Group PNGs and PDFs under each top-level subdirectory of run_dir.

    Returns: {module_name: {"png": [Path, ...], "pdf": [Path, ...]}}.
    Walks each top-level subdir up to ``max_depth`` levels deep.
    """
    grouped: dict[str, dict[str, list[Path]]] = {}
    if not run_dir.is_dir():
        return grouped
    for child in sorted(run_dir.iterdir()):
        if not child.is_dir():
            continue
        pngs: list[Path] = []
        pdfs: list[Path] = []
        # Bounded walk: depth relative to `child`.
        for dirpath, _dirs, files in os.walk(child):
            rel_depth = len(Path(dirpath).relative_to(child).parts)
            if rel_depth > max_depth:
                # Stop descending: clear dirs to prune walk.
                _dirs[:] = []
                continue
            for fname in files:
                lower = fname.lower()
                if lower.endswith(".png"):
                    pngs.append(Path(dirpath) / fname)
                elif lower.endswith(".pdf"):
                    pdfs.append(Path(dirpath) / fname)
        if pngs or pdfs:
            grouped[child.name] = {"png": sorted(pngs), "pdf": sorted(pdfs)}
    return grouped


def _plan_png_embedding(
    grouped: dict[str, dict[str, list[Path]]],
    *,
    per_image_max: int = REPORT_PNG_EMBED_MAX_BYTES,
    total_embed_max: int = REPORT_TOTAL_EMBED_MAX_BYTES,
) -> tuple[dict[Path, str], list[str]]:
    """Decide which PNGs to embed vs link.

    Returns:
      decisions: {png_path: "embed" | "link-oversize" | "link-total-cap"}
      notices:   list of human-readable notice strings to surface in HTML.
    """
    decisions: dict[Path, str] = {}
    notices: list[str] = []

    # First pass: per-image ceiling. Collect candidates with sizes.
    candidates: list[tuple[Path, int]] = []
    for mod, assets in grouped.items():
        for p in assets["png"]:
            try:
                size = p.stat().st_size
            except OSError:
                decisions[p] = "link-oversize"
                notices.append(f"could not stat {p}; linking instead of embedding")
                continue
            if size > per_image_max:
                decisions[p] = "link-oversize"
                notices.append(
                    f"{mod}/{p.name} is {size/1_048_576:.1f} MB "
                    f"(> {per_image_max/1_048_576:.0f} MB); linked, not embedded"
                )
            else:
                candidates.append((p, size))

    # Second pass: total embedded cap. base64 inflates ~4/3.
    total_b64 = 0
    for p, size in candidates:
        total_b64 += (size + 2) // 3 * 4

    if total_b64 > total_embed_max:
        # Switch ALL candidates to links.
        notices.insert(
            0,
            f"total embedded payload would be ~{total_b64/1_048_576:.1f} MB "
            f"(> {total_embed_max/1_048_576:.0f} MB cap); "
            f"ALL images switched to relative links",
        )
        for p, _ in candidates:
            decisions[p] = "link-total-cap"
    else:
        for p, _ in candidates:
            decisions[p] = "embed"

    return decisions, notices


def _claim_guard_banner_html(claim_guard: str | None, *, position: str) -> str:
    """Render the claim_guard banner. Always rendered, even when missing."""
    pos_label = html.escape(position)
    if claim_guard:
        msg = html.escape(claim_guard)
        return (
            f'<div class="claim-guard-banner" data-pos="{pos_label}">'
            f'<div class="cg-title">CLAIM GUARD</div>'
            f'<div class="cg-body">{msg}</div>'
            f"</div>"
        )
    return (
        f'<div class="claim-guard-banner cg-missing" data-pos="{pos_label}">'
        f'<div class="cg-title">CLAIM GUARD</div>'
        f'<div class="cg-body">'
        f"no claim_guard recorded for this run "
        f"(no R bundle exported or manifest lacks claim_guard field)"
        f"</div>"
        f"</div>"
    )


_REPORT_CSS = """
body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif;
       margin: 0; padding: 0; color: #222; background: #fafafa; }
header, footer { padding: 16px 24px; }
header { background: #1f2937; color: #f9fafb; }
header h1 { margin: 0 0 4px 0; font-size: 22px; }
header .meta { font-size: 12px; opacity: 0.85; }
main { padding: 16px 24px; max-width: 1100px; margin: 0 auto; }
.claim-guard-banner {
    background: #fff1f0; border: 3px solid #c0392b; color: #7b1d12;
    padding: 14px 18px; margin: 12px 24px; border-radius: 6px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
}
.claim-guard-banner.cg-missing { background: #fff7e6; border-color: #d97706;
    color: #8a4b00; }
.cg-title { font-weight: 700; font-size: 13px; letter-spacing: 0.08em;
    text-transform: uppercase; margin-bottom: 4px; }
.cg-body { font-size: 14px; white-space: pre-wrap; font-family:
    SFMono-Regular, Consolas, monospace; }
table.summary { border-collapse: collapse; width: 100%; margin: 12px 0; }
table.summary th, table.summary td { padding: 6px 10px; border: 1px solid #ddd;
    font-size: 13px; text-align: left; }
table.summary th { background: #f3f4f6; }
.status-ok { color: #166534; font-weight: 600; }
.status-failed, .status-fail, .status-error { color: #991b1b; font-weight: 600; }
.status-skipped, .status-skip { color: #78716c; }
section.module { background: #fff; border: 1px solid #e5e7eb; border-radius: 6px;
    padding: 14px 16px; margin: 14px 0; }
section.module h2 { margin: 0 0 8px 0; font-size: 16px; }
.figure { margin: 10px 0; }
.figure img { max-width: 100%; height: auto; border: 1px solid #e5e7eb;
    border-radius: 4px; background: #fff; }
.figure .caption { font-size: 12px; color: #4b5563; margin-top: 4px; }
.notice { background: #fef9c3; border-left: 4px solid #ca8a04; padding: 8px 12px;
    margin: 12px 0; font-size: 13px; }
.pdf-list ul { padding-left: 20px; }
footer { color: #6b7280; font-size: 12px; text-align: center; }
"""


def _render_module_section(
    module: str,
    assets: dict[str, list[Path]],
    decisions: dict[Path, str],
    run_dir: Path,
    status_lookup: dict[str, dict[str, str]],
) -> str:
    """Render an HTML section for a single module."""
    parts: list[str] = []
    safe_mod = html.escape(module)
    status_row = status_lookup.get(module, {})
    status = html.escape(status_row.get("status", "—"))
    message = html.escape(status_row.get("message", ""))
    parts.append(f'<section class="module" id="mod-{safe_mod}">')
    parts.append(f"<h2>{safe_mod}</h2>")
    parts.append(
        f'<div class="meta">status: '
        f'<span class="status-{status}">{status}</span>'
        + (f' <span class="msg">— {message}</span>' if message else "")
        + "</div>"
    )

    for png_path in assets["png"]:
        decision = decisions.get(png_path, "embed")
        try:
            rel = png_path.relative_to(run_dir)
        except ValueError:
            rel = png_path
        rel_href = html.escape(str(rel))
        caption = html.escape(str(rel))
        if decision == "embed":
            try:
                data = png_path.read_bytes()
                b64 = base64.b64encode(data).decode("ascii")
                parts.append(
                    '<div class="figure">'
                    f'<img alt="{caption}" '
                    f'src="data:image/png;base64,{b64}">'
                    f'<div class="caption">{caption}</div>'
                    "</div>"
                )
            except OSError:
                parts.append(
                    '<div class="figure">'
                    f'<a href="{rel_href}">{caption}</a>'
                    f'<div class="caption">read failed; linked instead</div>'
                    "</div>"
                )
        else:
            note = (
                "exceeds per-image 10 MB ceiling"
                if decision == "link-oversize"
                else "switched to link due to total embed cap"
            )
            parts.append(
                '<div class="figure">'
                f'<a href="{rel_href}">{caption}</a>'
                f'<div class="caption">{html.escape(note)}</div>'
                "</div>"
            )

    if assets["pdf"]:
        parts.append('<div class="pdf-list"><strong>PDFs:</strong><ul>')
        for pdf_path in assets["pdf"]:
            try:
                rel = pdf_path.relative_to(run_dir)
            except ValueError:
                rel = pdf_path
            rel_href = html.escape(str(rel))
            label = html.escape(str(rel))
            parts.append(f'<li><a href="{rel_href}">{label}</a></li>')
        parts.append("</ul></div>")

    parts.append("</section>")
    return "\n".join(parts)


def _render_report_html(
    *,
    title: str,
    project: str,
    generated_at: str,
    manifest_modules: list[dict[str, str]],
    csv_modules: list[dict[str, str]],
    grouped: dict[str, dict[str, list[Path]]],
    decisions: dict[Path, str],
    notices: list[str],
    claim_guard: str | None,
    run_dir: Path,
) -> str:
    """Compose the full HTML report as a string."""
    # Merge module status from manifest + CSV (CSV wins on conflict).
    status_lookup: dict[str, dict[str, str]] = {}
    for row in manifest_modules:
        name = row.get("module") or row.get("name")
        if not name:
            continue
        status_lookup[name] = {
            "status": str(row.get("status", "")),
            "message": str(row.get("message", "")),
        }
    for row in csv_modules:
        name = row.get("module") or row.get("name")
        if not name:
            continue
        status_lookup[name] = {
            "status": str(row.get("status", "")),
            "message": str(row.get("message", "")),
        }

    # Summary table.
    summary_rows: list[str] = []
    all_module_names = sorted(set(status_lookup.keys()) | set(grouped.keys()))
    for name in all_module_names:
        status = html.escape(status_lookup.get(name, {}).get("status", "—"))
        n_fig = len(grouped.get(name, {"png": []})["png"]) + \
                len(grouped.get(name, {"pdf": []}).get("pdf", []))
        summary_rows.append(
            f"<tr><td>{html.escape(name)}</td>"
            f'<td><span class="status-{status}">{status}</span></td>'
            f"<td>{n_fig}</td></tr>"
        )

    # Per-module sections, in the same sorted order.
    module_sections = "\n".join(
        _render_module_section(name, grouped[name], decisions, run_dir, status_lookup)
        for name in all_module_names
        if name in grouped
    )

    # Notices block.
    notices_html = ""
    if notices:
        items = "".join(f"<li>{html.escape(n)}</li>" for n in notices)
        notices_html = f'<div class="notice"><strong>Notices:</strong><ul>{items}</ul></div>'

    now_iso = _dt.datetime.now().isoformat(timespec="seconds")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(title)}</title>
<style>{_REPORT_CSS}</style>
</head>
<body>
<header>
<h1>{html.escape(title)}</h1>
<div class="meta">project: {html.escape(project)} &middot; generated_at: {html.escape(generated_at)}</div>
</header>
{_claim_guard_banner_html(claim_guard, position="top")}
<main>
{notices_html}
<h2>Module summary</h2>
<table class="summary">
<thead><tr><th>module</th><th>status</th><th>n_figures</th></tr></thead>
<tbody>
{''.join(summary_rows) if summary_rows else '<tr><td colspan="3">no modules found</td></tr>'}
</tbody>
</table>
{module_sections}
</main>
{_claim_guard_banner_html(claim_guard, position="bottom")}
<footer>generated by scfactory v{html.escape(SCFACTORY_VERSION)} at {html.escape(now_iso)}</footer>
</body>
</html>
"""


def cmd_report(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).resolve()
    manifest_path = run_dir / "run_manifest.json"

    if not manifest_path.is_file():
        print(
            f"scfactory: run_manifest.json not found at {manifest_path}\n"
            "  pass a run directory that contains run_manifest.json "
            "(e.g. results/<project>/).",
            file=sys.stderr,
        )
        return 2

    try:
        manifest = json.loads(manifest_path.read_text())
    except Exception as exc:
        print(
            f"scfactory: could not parse {manifest_path}: {exc}",
            file=sys.stderr,
        )
        return 2

    project = str(manifest.get("project") or run_dir.name)
    generated_at = str(manifest.get("generated_at") or "")
    title = args.title or project

    manifest_modules_raw = (
        manifest.get("module_status")
        or manifest.get("modules")
        or []
    )
    manifest_modules: list[dict[str, str]] = []
    if isinstance(manifest_modules_raw, list):
        for row in manifest_modules_raw:
            if isinstance(row, dict):
                manifest_modules.append({k: str(v) for k, v in row.items()})

    csv_modules = _read_module_status_csv(run_dir / "module_status.csv")
    claim_guard = _find_claim_guard(manifest, run_dir)
    grouped = _walk_module_assets(run_dir)
    decisions, notices = _plan_png_embedding(grouped)

    html_doc = _render_report_html(
        title=title,
        project=project,
        generated_at=generated_at,
        manifest_modules=manifest_modules,
        csv_modules=csv_modules,
        grouped=grouped,
        decisions=decisions,
        notices=notices,
        claim_guard=claim_guard,
        run_dir=run_dir,
    )

    out_path = (
        Path(args.out).resolve() if args.out else (run_dir / "report.html")
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html_doc, encoding="utf-8")
    print(str(out_path))
    return 0


# ---------------------------------------------------------------------------
# Argparse wiring
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scfactory",
        description="User-friendly wrapper around the modular pipeline.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser(
        "run",
        help="Auto-detect modality and dispatch to workflow.modular.cli",
        description=(
            "Run the modular pipeline on an .h5ad or sample-root.\n\n"
            "Module precedence (highest wins):\n"
            "  1. --optional-modules <list>   (escape hatch; overrides all)\n"
            "  2. --recipe <name>             (preset modules + env + bundle)\n"
            "  3. auto-detected modality      (fills any remaining gap)\n\n"
            "Batch declaration follows the same precedence:\n"
            "  1. --batch-strategy <name>\n"
            "  2. --recipe (its named profile, or its batch_strategy field)\n"
            "  3. undeclared (multi-batch input is warned about and the "
            "clustering claim is downgraded to exploratory)\n\n"
            "--bundle and recipe.bundle.enabled are OR-ed: passing --bundle "
            "always forces bundle on even if the recipe disables it."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_run.add_argument("input", nargs="?", default=None,
                       help=".h5ad file or sample-root directory "
                            "(omit only with --list-recipes)")
    p_run.add_argument("--project", default=None,
                       help="Run name (default: scfactory_<modality>)")
    output_group = p_run.add_mutually_exclusive_group()
    output_group.add_argument(
        "--project-root",
        default=None,
        metavar="PATH",
        help=(
            "Governed project root; forwards outputs to "
            "<project-root>/runs/<run-id>/python/"
        ),
    )
    output_group.add_argument(
        "--out",
        default=None,
        help=(
            "Deprecated compatibility output directory (default: <repo>/results); "
            "real launches are access-logged by the canonical CLI"
        ),
    )
    p_run.add_argument(
        "--run-id",
        default=None,
        metavar="STR",
        help="Governed run identifier; requires --project-root",
    )
    p_run.add_argument("--bundle", action="store_true",
                       help="After run, export an R bundle via "
                            "scripts/export_singlecell_r_bundle.py "
                            "(forces bundle on even if recipe disables it)")
    p_run.add_argument("--bundle-out", default=None,
                       help="Bundle output path (default: <out>/<project>/r_bundle)")
    p_run.add_argument("--optional-modules", default=None,
                       help="Override auto-detected optional module list "
                            "(comma-separated; passthrough escape hatch — "
                            "wins over --recipe)")
    p_run.add_argument(
        "--atac-peak-matrix-path",
        default=None,
        help="Sparse peak-count matrix consumed by the canonical atac_ingest producer",
    )
    p_run.add_argument(
        "--atac-peaks-bed-path",
        default=None,
        help="Ordered peaks BED consumed by the canonical atac_ingest producer",
    )
    p_run.add_argument(
        "--atac-n-components",
        type=int,
        default=None,
        help="TF-IDF/LSI components for atac_ingest (canonical default: 30)",
    )
    p_run.add_argument("--scatac-da-sample-col", default=None,
                       help="Explicit biological-sample obs column for scATAC pseudobulk DA")
    p_run.add_argument("--scatac-da-group-col", default=None,
                       help="Explicit cell-group obs column for scATAC pseudobulk DA")
    p_run.add_argument("--scatac-da-condition-col", default=None,
                       help="Explicit two-level condition obs column for confirmatory scATAC DA")
    p_run.add_argument("--scatac-da-peak-id-col", default=None,
                       help="Explicit ordered peak-ID column in atac_var")
    p_run.add_argument("--scatac-da-test-level", default=None,
                       help="Test/numerator condition level for confirmatory scATAC DA")
    p_run.add_argument("--scatac-da-reference-level", default=None,
                       help="Reference/denominator condition level for confirmatory scATAC DA")
    p_run.add_argument("--scatac-da-mode", choices=["confirmatory_da", "aggregation_only"], default=None,
                       help="scATAC route: confirmatory_da or aggregation_only")
    p_run.add_argument("--scatac-da-min-samples-per-condition", type=int, default=None,
                       help="Minimum biological samples per condition (must be at least 2)")
    p_run.add_argument("--scatac-da-min-total-count", type=int, default=None,
                       help="Minimum aggregate peak count passed to the R engine")
    p_run.add_argument("--scatac-da-fdr-threshold", type=float, default=None,
                       help="FDR threshold for scATAC DA")
    p_run.add_argument("--scatac-da-abs-log2fc-threshold", type=float, default=None,
                       help="Absolute log2 fold-change threshold for scATAC DA")
    p_run.add_argument("--scatac-da-groups", default=None,
                       help="Optional comma-separated group allowlist for scATAC DA")
    p_run.add_argument("--scatac-da-r-conda-env", default=None,
                       help="Declared R environment for scATAC DA (r_multiomics)")
    p_run.add_argument("--scatac-da-timeout", type=int, default=None,
                       help="Per-group R inference timeout in seconds")
    p_run.add_argument("--scatac-da-aggregation-backend", choices=["cpu"], default=None,
                       help="Exact scATAC aggregation backend (CPU only)")
    p_run.add_argument("--recipe", default=None,
                       help=f"Apply a preset recipe from {RECIPES_DIR.name}/<name>.yaml "
                            "(modules + scale-mode + env + bundle config)")
    p_run.add_argument(
        "--batch-strategy",
        default=BATCH_STRATEGY_AUTO,
        choices=list(BATCH_STRATEGY_CHOICES),
        help=(
            "Declare the batch design of the input; forwarded to the canonical "
            "CLI. 'auto' is the absence of a declaration and leaves multi-batch "
            "clustering marked exploratory. Wins over a recipe's declaration."
        ),
    )
    p_run.add_argument("--list-recipes", action="store_true",
                       help="List available recipes (one per line) and exit 0")
    p_run.add_argument(
        "--scientific-profile",
        default=None,
        metavar="NAME",
        help="Forward a named scientific profile to the canonical CLI",
    )
    p_run.add_argument(
        "--acknowledge-scientific-non-equivalence",
        action="store_true",
        help=(
            "Forward explicit acknowledgement for non-canonical scientific "
            "settings"
        ),
    )
    p_run.add_argument(
        "--allow-dirty",
        action="store_true",
        help=(
            "Forward to the modular CLI: allow a dirty factory git tree and "
            "record the diff in the run manifest (required when uncommitted "
            "factory changes are present)"
        ),
    )
    p_run.add_argument("--dry-run", action="store_true",
                       help="Print what would run, do not execute")
    p_run.set_defaults(func=cmd_run)

    p_doc = sub.add_parser(
        "doctor",
        help="Read-only environment health check",
    )
    p_doc.add_argument("--json", action="store_true",
                       help="Emit machine-readable JSON")
    p_doc.add_argument(
        "--project-root",
        action="append",
        default=[],
        metavar="PATH",
        help=(
            "Governed project root to inspect for runs (repeatable; default: "
            "all immediate projects under ~/projects)"
        ),
    )
    p_doc.set_defaults(func=cmd_doctor)

    p_rep = sub.add_parser(
        "report",
        help="Render a single self-contained HTML report from a run directory",
    )
    p_rep.add_argument("run_dir",
                       help="Path to a run directory containing run_manifest.json")
    p_rep.add_argument("--out", default=None,
                       help="Output HTML path (default: <run_dir>/report.html)")
    p_rep.add_argument("--title", default=None,
                       help="Report title (default: project name from manifest)")
    p_rep.set_defaults(func=cmd_report)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except SystemExit as exc:
        # Recipe validation / PyYAML import errors raise SystemExit with an
        # int code and a stderr message. Convert that into a normal return so
        # programmatic callers (and tests) can assert on the exit code without
        # wrapping with pytest.raises.
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        # Non-int message: emit it then return 1.
        print(str(code), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
