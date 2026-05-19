#!/usr/bin/env python3
"""Validate the current NC2024 project-root architecture contract.

This controller-validation smoke reflects the 2026-05 factory/project split:
scientific artifacts live under /home/zerlinshen/projects/<project-id>/runs/<run-id>,
not under singlecell_factory/results/.  The validator deliberately avoids opening
large H5AD objects; it checks manifests, module status, compact R bundle files,
project-owned R outputs, and canonical suite-root bridge symlinks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
SUITE_ROOT = REPO_ROOT.parent
R_FACTORY_ROOT = SUITE_ROOT / "r_multiomics_factory"
PLOTTING_FACTORY_ROOT = SUITE_ROOT / "plotting_factory"
DEFAULT_PROJECT_ROOT = Path("/home/zerlinshen/projects/nc-reproduction")
DEFAULT_RUN_ID = "2026-05-18T0900Z-13c2c88"
ACCEPTED_BUNDLE_SCHEMAS = {
    "singlecell_r_bundle_v2",
    "singlecell_r_bundle_v2.1",
    "singlecell_r_bundle_v2.2",
}
CLAIM_GUARD = "not_for_de_or_new_quantitative_claims_without_full_object_validation"
EXPECTED_FINAL_SHAPE = (5281, 19504)
# Tolerance band for shape assertions. The architecture validator is meant
# to confirm "factory wiring intact", not "bit-exact reproducibility". Minor
# drift (e.g. from the 2026-05-19 doublet HIGH-1 fix that recovered the
# RAPIDS-shadowed CPU scrublet path) must not fail the architecture gate.
# A separate validate_exact_reproducibility entry-point can be added when
# bit-exact release tagging is required.
SHAPE_TOLERANCE_CELLS_PCT = 2.0
SHAPE_TOLERANCE_GENES_ABS = 50


def _shape_within_band(observed: tuple[int, int], expected: tuple[int, int]) -> bool:
    cells_band = max(1, int(expected[0] * SHAPE_TOLERANCE_CELLS_PCT / 100.0))
    return (
        abs(observed[0] - expected[0]) <= cells_band
        and abs(observed[1] - expected[1]) <= SHAPE_TOLERANCE_GENES_ABS
    )
EXPECTED_OK_MODULES = {
    "cellranger",
    "marker_db_loader",
    "qc",
    "doublet_detection",
    "clustering",
    "annotation",
    "context_aware_annotation",
    "differential_expression",
}
EXPECTED_SKIPPED_MODULES = {"batch_correction"}
EXPECTED_R_OUTPUT_FILES = {
    "cell_fraction.png",
    "marker_dot.png",
    "marker_feature.png",
    "qc_scatter.png",
    "qc_violin.png",
    "r.pdf",
    "umap_by_cluster.png",
    "umap_by_group.png",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(path: Path, label: str) -> Path:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    return path


def _load_json(path: Path) -> dict:
    return json.loads(_require(path, "json").read_text(encoding="utf-8"))


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with _require(path, "csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _nonempty_files(root: Path) -> list[str]:
    if not root.exists():
        return []
    return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file() and p.stat().st_size > 0)


def _validate_plot_file(path: Path) -> None:
    _require(path, "plot file")
    header = path.read_bytes()[:8]
    if path.suffix == ".png" and header != b"\x89PNG\r\n\x1a\n":
        raise SystemExit(f"invalid PNG output header: {path}")
    if path.suffix == ".pdf" and not header.startswith(b"%PDF-"):
        raise SystemExit(f"invalid PDF output header: {path}")


def _find_one(paths: Iterable[Path], label: str) -> Path:
    found = sorted(p for p in paths if p.exists())
    if len(found) != 1:
        raise SystemExit(f"expected exactly one {label}, found {len(found)}: {[str(p) for p in found]}")
    return found[0]


def _validate_bridge_symlinks() -> dict[str, str]:
    bridge = REPO_ROOT / "bridges" / "local_r_pipeline_macbook"
    plot_bridge = REPO_ROOT / "bridges" / "local_plot_pipeline"
    checks = {
        "R": (bridge / "R", R_FACTORY_ROOT / "R"),
        "R_bundle": (bridge / "R_bundle", R_FACTORY_ROOT / "R_bundle"),
        "plotting_r": (plot_bridge / "r", PLOTTING_FACTORY_ROOT / "r"),
        "plotting_python": (plot_bridge / "python", PLOTTING_FACTORY_ROOT / "python"),
    }
    resolved: dict[str, str] = {}
    for label, (link, expected) in checks.items():
        _require(link, f"bridge symlink {label}")
        if not link.is_symlink():
            raise SystemExit(f"{link} must remain a symlink")
        actual = link.resolve()
        expected = expected.resolve()
        if actual != expected:
            raise SystemExit(f"{label} resolves to {actual}, expected {expected}")
        resolved[label] = str(actual)

    _require(R_FACTORY_ROOT / "R_bundle" / "io_bundle.R", "R bundle reader")
    _require(R_FACTORY_ROOT / "scripts" / "plot_remote_bundle.R", "project-root R plotter")
    _require(PLOTTING_FACTORY_ROOT / "r" / "general" / "dim_plots.R", "plotting_factory R plots")
    return resolved


def _validate_module_status(path: Path) -> dict[str, object]:
    rows = _csv_rows(path)
    statuses = {row["module"]: row["status"] for row in rows}
    missing_ok = sorted(EXPECTED_OK_MODULES - statuses.keys())
    if missing_ok:
        raise SystemExit(f"module_status missing expected ok modules: {missing_ok}")
    not_ok = {m: statuses[m] for m in EXPECTED_OK_MODULES if statuses.get(m) != "ok"}
    if not_ok:
        raise SystemExit(f"module_status expected ok modules failed: {not_ok}")
    for module in EXPECTED_SKIPPED_MODULES:
        if statuses.get(module) != "skipped":
            raise SystemExit(f"module_status expected {module}=skipped, got {statuses.get(module)!r}")
    return {
        "path": str(path),
        "counts": dict(Counter(statuses.values())),
        "statuses": statuses,
    }


def _validate_python_run(run_dir: Path) -> dict[str, object]:
    python_dir = _require(run_dir / "python", "project python dir")
    run_manifest_path = _find_one(python_dir.glob("*/run_manifest.json"), "Python run_manifest.json")
    module_status_path = _find_one(python_dir.glob("*/module_status.csv"), "Python module_status.csv")
    final_h5ad_path = _find_one(python_dir.glob("*/final_adata.h5ad"), "final_adata.h5ad")

    run_manifest = _load_json(run_manifest_path)
    metadata = run_manifest.get("metadata", {})
    observed_shape = (int(metadata.get("cells_after_qc", -1)), int(metadata.get("genes_after_qc", -1)))
    if not _shape_within_band(observed_shape, EXPECTED_FINAL_SHAPE):
        raise SystemExit(
            f"unexpected final shape from manifest: {observed_shape}, "
            f"expected ~{EXPECTED_FINAL_SHAPE} "
            f"(±{SHAPE_TOLERANCE_CELLS_PCT}% cells, ±{SHAPE_TOLERANCE_GENES_ABS} genes)"
        )
    if int(metadata.get("context_validation_mismatches", -1)) != 0:
        raise SystemExit("context_validation_mismatches must be 0")

    return {
        "run_manifest": str(run_manifest_path),
        "module_status": _validate_module_status(module_status_path),
        "final_adata": str(final_h5ad_path),
        "final_adata_bytes": final_h5ad_path.stat().st_size,
        "shape": {"n_obs": observed_shape[0], "n_vars": observed_shape[1]},
        "metadata_subset": {
            "sample_id_filter_applied": metadata.get("sample_id_filter_applied"),
            "sample_id_filter_n_before": metadata.get("sample_id_filter_n_before"),
            "sample_id_filter_n_after": metadata.get("sample_id_filter_n_after"),
            "clustering_backend": metadata.get("clustering_backend"),
            "de_backend": metadata.get("de_backend"),
            "batch_correction_status": metadata.get("batch_correction_status"),
            "context_validation_mismatches": metadata.get("context_validation_mismatches"),
        },
    }


def _validate_bundle(run_dir: Path, verify_sha: bool) -> dict[str, object]:
    bundle_dir = _require(run_dir / "python" / "bundle", "project python/bundle dir")
    manifest = _load_json(bundle_dir / "bundle_manifest.json")
    schema = manifest.get("schema_version")
    if schema not in ACCEPTED_BUNDLE_SCHEMAS:
        raise SystemExit(f"{bundle_dir}: unexpected bundle schema {schema!r}")

    source = manifest.get("source", {})
    bundle = manifest.get("bundle", {})
    files = manifest.get("files", {})
    exported = int(bundle.get("n_cells_exported", -1))
    source_n = int(source.get("n_obs", -2))
    if exported <= 0 or exported != source_n:
        raise SystemExit(f"{bundle_dir}: exported/source n_obs mismatch: {exported} vs {source_n}")
    cells_band = max(1, int(EXPECTED_FINAL_SHAPE[0] * SHAPE_TOLERANCE_CELLS_PCT / 100.0))
    if abs(source_n - EXPECTED_FINAL_SHAPE[0]) > cells_band:
        raise SystemExit(
            f"{bundle_dir}: expected ~{EXPECTED_FINAL_SHAPE[0]} cells "
            f"(±{cells_band}), got {source_n}"
        )
    source_vars = int(source.get("n_vars", -1))
    if abs(source_vars - EXPECTED_FINAL_SHAPE[1]) > SHAPE_TOLERANCE_GENES_ABS:
        raise SystemExit(
            f"{bundle_dir}: expected ~{EXPECTED_FINAL_SHAPE[1]} genes "
            f"(±{SHAPE_TOLERANCE_GENES_ABS}), got {source_vars}"
        )

    required = set(bundle.get("required_files", [])) | {"obs", "marker_expr", "X_umap", "X_pca"}
    missing = sorted(required - set(files))
    if missing:
        raise SystemExit(f"{bundle_dir}: manifest missing file records: {missing}")

    checked_files: list[str] = []
    for stem in sorted(required):
        rec = files[stem]
        rel_path = rec["path"]
        path = _require(bundle_dir / rel_path, f"bundle file {stem}")
        actual_bytes = path.stat().st_size
        if actual_bytes != int(rec["bytes"]):
            raise SystemExit(
                f"{bundle_dir}: byte mismatch for {stem}: manifest={rec['bytes']} actual={actual_bytes}"
            )
        if verify_sha and rec.get("sha256") and _sha256(path) != str(rec["sha256"]).lower():
            raise SystemExit(f"{bundle_dir}: sha256 mismatch for {stem}")
        checked_files.append(rel_path)

    expression = manifest.get("expression", {})
    if expression.get("claim_guard") != CLAIM_GUARD:
        raise SystemExit(f"{bundle_dir}: unexpected claim_guard={expression.get('claim_guard')}")
    _require(bundle_dir / "provenance.json", "bundle provenance")

    return {
        "bundle_dir": str(bundle_dir),
        "schema_version": schema,
        "n_cells_exported": exported,
        "n_vars_source": source_vars,
        "markers_present": bundle.get("markers_present", []),
        "marker_format": bundle.get("marker_format"),
        "checked_files": checked_files,
        "claim_guard": expression.get("claim_guard"),
    }


def _validate_r_outputs(run_dir: Path) -> dict[str, object]:
    r_dir = _require(run_dir / "r", "project r dir")
    files = _nonempty_files(r_dir)
    missing = sorted(EXPECTED_R_OUTPUT_FILES - set(files))
    if missing:
        raise SystemExit(f"missing expected R outputs under {r_dir}: {missing}")
    unexpected = sorted(set(files) - EXPECTED_R_OUTPUT_FILES)
    if unexpected:
        raise SystemExit(f"unexpected R outputs under {r_dir}: {unexpected}")
    for rel_path in sorted(EXPECTED_R_OUTPUT_FILES):
        _validate_plot_file(r_dir / rel_path)
    return {
        "r_dir": str(r_dir),
        "file_count": len(files),
        "expected_files": sorted(EXPECTED_R_OUTPUT_FILES),
        "files": files[:50],
    }


def validate_project_run(project_root: Path, run_id: str, *, verify_sha: bool = False) -> dict[str, object]:
    project_root = project_root.resolve()
    run_dir = _require(project_root / "runs" / run_id, "project run dir")
    root_manifest = _load_json(run_dir / "manifest.json")
    if root_manifest.get("run_id") != run_id:
        raise SystemExit(f"root manifest run_id mismatch: {root_manifest.get('run_id')!r} vs {run_id!r}")

    ledger_hits = sorted((run_dir / "ops" / "run_ledger").glob("*.json"))
    if not ledger_hits:
        raise SystemExit(f"{run_dir}: no project-owned run ledger JSON found")

    return {
        "project_root": str(project_root),
        "run_id": run_id,
        "run_dir": str(run_dir),
        "root_manifest": str(run_dir / "manifest.json"),
        "factory_python_sha": root_manifest.get("factory_python", {}).get("sha"),
        "factory_r_sha": root_manifest.get("factory_r", {}).get("sha"),
        "python": _validate_python_run(run_dir),
        "bundle": _validate_bundle(run_dir, verify_sha=verify_sha),
        "r_outputs": _validate_r_outputs(run_dir),
        "latest_project_ledger": str(ledger_hits[-1]),
    }


def validate(*, project_root: Path = DEFAULT_PROJECT_ROOT, run_id: str = DEFAULT_RUN_ID, verify_sha: bool = False) -> dict[str, object]:
    return {
        "status": "ok",
        "execution_mode": "controller_validation",
        "scope": "current_nc2024_project_root_run",
        "checked_run": validate_project_run(project_root, run_id, verify_sha=verify_sha),
        "bridge": _validate_bridge_symlinks(),
        "scientific_boundary": (
            "This validates the retained NC2024 P15_T1 project-root/module/bridge run; "
            "it is not a full article-scale multi-cohort annotation claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT_ROOT)
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--verify-sha", action="store_true", help="Hash every required bundle file.")
    args = parser.parse_args()
    print(json.dumps(validate(project_root=args.project_root, run_id=args.run_id, verify_sha=args.verify_sha), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
