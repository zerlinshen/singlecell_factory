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
    from workflow.modular.module_catalog import MODULE_SPECS
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

    # optional_modules required, must be list[str], all known
    mods = data.get("optional_modules")
    if not isinstance(mods, list) or not all(isinstance(m, str) for m in mods):
        print(
            f"scfactory: {where}: 'optional_modules' is required and must be "
            f"a list of strings.",
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
    if args.optional_modules:
        planned = [m.strip() for m in args.optional_modules.split(",") if m.strip()]
        plan_source = "user-provided --optional-modules"
    elif recipe is not None:
        planned = [str(m) for m in recipe.get("optional_modules", [])]
        plan_source = f"recipe={recipe['name']}"
    else:
        planned = plan_optional_modules(modality)
        plan_source = f"auto for modality={modality}"

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
        REPO_ROOT / "results"
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

    if args.scientific_profile:
        cli_cmd += ["--scientific-profile", args.scientific_profile]
    if args.acknowledge_scientific_non_equivalence:
        cli_cmd.append("--acknowledge-scientific-non-equivalence")

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
    deps = ("scanpy", "anndata", "pyarrow", "fastparquet", "scipy",
            "squidpy", "muon", "mofapy2")
    out: dict[str, Any] = {}
    import importlib
    import importlib.metadata as _imd
    for name in deps:
        try:
            importlib.import_module(name)
            try:
                ver = _imd.version(name)
            except Exception:
                ver = "unknown"
            out[name] = {"present": True, "version": ver}
        except Exception:
            out[name] = {"present": False, "version": None}
    return {
        "status": "pass",  # informational only
        "deps": out,
        "message": "python deps: "
                   + ", ".join(f"{n}={'ok' if v['present'] else 'missing'}"
                               for n, v in out.items()),
    }


def _check_r_packages(rscript_path: str | None) -> dict[str, Any]:
    if rscript_path is None:
        return {"status": "warn", "packages": {},
                "message": "R: not available; skipping package probe"}
    pkgs = ("Seurat", "arrow", "jsonlite", "Matrix", "zellkonverter",
            "SeuratDisk", "harmony")
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
                "message": "R package probe timed out (>5s)"}
    except Exception as exc:
        return {"status": "warn", "packages": {},
                "message": f"R package probe failed: {exc}"}

    parsed: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if "|" in line:
            name, ver = line.split("|", 1)
            parsed[name.strip()] = ver.strip()
    return {
        "status": "pass",
        "packages": parsed,
        "message": "R packages: "
                   + ", ".join(f"{n}={v}" for n, v in parsed.items()),
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
    }

    counts = {"pass": 0, "warn": 0, "fail": 0}
    for v in report.values():
        s = v.get("status", "warn")
        counts[s] = counts.get(s, 0) + 1
    report["summary"] = counts

    if args.json:
        print(json.dumps(report, indent=2, default=str))
    else:
        print("scfactory doctor — read-only environment health check")
        print("=" * 56)
        for key in ("python", "conda_envs", "rscript", "bridge_symlinks",
                    "python_deps", "r_packages", "last_run"):
            v = report[key]
            tag = v["status"].upper()
            print(f"[{tag:4}] {key}: {v['message']}")
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
        help="Deprecated factory-local output directory (default: <repo>/results)",
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
    p_run.add_argument("--recipe", default=None,
                       help=f"Apply a preset recipe from {RECIPES_DIR.name}/<name>.yaml "
                            "(modules + scale-mode + env + bundle config)")
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
