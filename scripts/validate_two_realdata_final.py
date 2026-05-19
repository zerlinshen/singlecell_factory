#!/usr/bin/env python3
"""Validate final two-real-dataset evidence for the three-factory architecture."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# 2026-05-19 CRITICAL fix: run-id and output-dir arguments come from the
# command line and are used to construct filesystem paths. Validate them
# against a strict allowlist + containment check before use so a typo or a
# pathological argument cannot escape the project tree.
_SAFE_RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9._:T-]{1,128}$")


def _safe_run_id(value: str) -> str:
    if not _SAFE_RUN_ID_PATTERN.match(value):
        raise argparse.ArgumentTypeError(
            f"run id {value!r} is not in the allowed character class "
            f"[A-Za-z0-9._:T-] (max 128 chars)"
        )
    return value


def _safe_output_dir(value: str) -> Path:
    path = Path(value).resolve()
    suite_root = Path(__file__).resolve().parents[2]
    project_roots = (
        Path("/home/zerlinshen/projects").resolve(),
        suite_root,
        Path("/tmp").resolve(),
    )
    if not any(_is_relative_to(path, root) for root in project_roots):
        raise argparse.ArgumentTypeError(
            f"output dir {path} must live under projects/, suite root, or /tmp"
        )
    return path


def _is_relative_to(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False

import validate_nc2024_architecture_contract as nc_validator


DEFAULT_CELL_PROJECT_ROOT = Path("/home/zerlinshen/projects/wave5-trevino")
DEFAULT_CELL_PIPELINE_RUN_ID = "2026-05-17T2004Z-13c2c88"
DEFAULT_CELL_EVIDENCE_RUN_ID = "20260517T1436Z-13c2c88"
EXPECTED_CELL_SHAPE = (55653, 25519)
# Tolerance band — see validate_nc2024_architecture_contract.SHAPE_TOLERANCE_*
SHAPE_TOLERANCE_CELLS_PCT = 2.0
SHAPE_TOLERANCE_GENES_ABS = 50


def _shape_within_band(observed: tuple[int, int], expected: tuple[int, int]) -> bool:
    cells_band = max(1, int(expected[0] * SHAPE_TOLERANCE_CELLS_PCT / 100.0))
    return (
        abs(observed[0] - expected[0]) <= cells_band
        and abs(observed[1] - expected[1]) <= SHAPE_TOLERANCE_GENES_ABS
    )


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing json: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing csv: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _nonempty_files(root: Path) -> list[str]:
    if not root.exists():
        return []
    return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file() and p.stat().st_size > 0)


def _validate_plot_file(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"missing plot file: {path}")
    header = path.read_bytes()[:8]
    if path.suffix == ".png" and header != b"\x89PNG\r\n\x1a\n":
        raise SystemExit(f"invalid PNG output header: {path}")
    if path.suffix == ".pdf" and not header.startswith(b"%PDF-"):
        raise SystemExit(f"invalid PDF output header: {path}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _untracked_inventory(repo: Path) -> dict[str, object]:
    # 2026-05-19 MEDIUM fix: cap per-file hashing at 50 MB. Untracked
    # multi-GB AnnData / Zarr blobs were being SHA-hashed unconditionally,
    # turning a few seconds of git status into minutes of disk I/O.
    sha_size_cap = 50 * 1024 * 1024
    raw = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=all"], cwd=repo, text=True)
    entries: list[dict[str, object]] = []
    for line in raw.splitlines():
        if not line.startswith("?? "):
            continue
        rel = line[3:]
        path = repo / rel
        if not path.is_file():
            continue
        size = path.stat().st_size
        entry: dict[str, object] = {"path": rel, "bytes": size}
        if size <= sha_size_cap:
            entry["sha256"] = _sha256(path)
        else:
            entry["sha256"] = None
            entry["sha256_skipped_reason"] = f"file_above_{sha_size_cap}_byte_cap"
        entries.append(entry)
    payload = "\n".join(f"{entry['path']}\t{entry['bytes']}\t{entry['sha256']}" for entry in entries)
    return {
        "count": len(entries),
        "aggregate_sha256": _sha256_bytes(payload.encode("utf-8")),
        "files": entries,
    }


def _git_state(repo: Path) -> dict[str, object]:
    def run(*args: str) -> str:
        return subprocess.check_output(args, cwd=repo, text=True).strip()

    status = run("git", "status", "--short")
    tracked_diff = subprocess.check_output(["git", "diff", "--binary"], cwd=repo)
    return {
        "path": str(repo),
        "branch": run("git", "branch", "--show-current"),
        "head": run("git", "rev-parse", "--short", "HEAD"),
        "dirty": bool(status),
        "status_short": status.splitlines(),
        "tracked_diff_sha256": _sha256_bytes(tracked_diff),
        "untracked_inventory": _untracked_inventory(repo),
    }


def _validate_bundle(run_dir: Path, expected_shape: tuple[int, int], *, verify_sha: bool = False) -> dict[str, object]:
    bundle_dir = run_dir / "python" / "bundle"
    manifest = _load_json(bundle_dir / "bundle_manifest.json")
    schema = manifest.get("schema_version")
    if schema not in nc_validator.ACCEPTED_BUNDLE_SCHEMAS:
        raise SystemExit(f"{bundle_dir}: unexpected bundle schema {schema!r}")
    source_n = int(manifest.get("source", {}).get("n_obs", -1))
    source_vars = int(manifest.get("source", {}).get("n_vars", -1))
    exported = int(manifest.get("bundle", {}).get("n_cells_exported", -2))
    expected_cells, expected_genes = expected_shape
    if source_n != expected_cells or exported != expected_cells:
        raise SystemExit(f"{bundle_dir}: expected {expected_cells} cells, got source={source_n} exported={exported}")
    if source_vars != expected_genes:
        raise SystemExit(f"{bundle_dir}: expected {expected_genes} genes, got source={source_vars}")
    if manifest.get("expression", {}).get("claim_guard") != nc_validator.CLAIM_GUARD:
        raise SystemExit(f"{bundle_dir}: invalid claim_guard")
    files = manifest.get("files", {})
    required = set(manifest.get("bundle", {}).get("required_files", [])) | {"obs", "marker_expr", "X_umap", "X_pca"}
    missing_records = sorted(required - set(files))
    if missing_records:
        raise SystemExit(f"{bundle_dir}: manifest missing file records: {missing_records}")
    invalid = []
    for stem in sorted(required):
        rec = files[stem]
        path = bundle_dir / rec.get("path", "")
        path_ok = path.exists() and path.stat().st_size == int(rec.get("bytes", -1))
        if not path_ok:
            invalid.append(stem)
        if path_ok and verify_sha and rec.get("sha256") and _sha256(path) != str(rec["sha256"]).lower():
            raise SystemExit(f"{bundle_dir}: sha256 mismatch for {stem}")
    if invalid:
        raise SystemExit(f"{bundle_dir}: missing/truncated bundle file records: {invalid}")
    _load_json(bundle_dir / "provenance.json")
    return {
        "bundle_dir": str(bundle_dir),
        "schema_version": schema,
        "n_cells_exported": exported,
        "n_vars_source": source_vars,
        "markers_present": manifest.get("bundle", {}).get("markers_present", []),
        "checked_file_count": len(required),
        "provenance": str(bundle_dir / "provenance.json"),
    }


def _validate_r_outputs(run_dir: Path) -> dict[str, object]:
    r_dir = run_dir / "r"
    files = _nonempty_files(r_dir)
    missing = sorted(nc_validator.EXPECTED_R_OUTPUT_FILES - set(files))
    if missing:
        raise SystemExit(f"missing expected R outputs under {r_dir}: {missing}")
    unexpected = sorted(set(files) - nc_validator.EXPECTED_R_OUTPUT_FILES)
    if unexpected:
        raise SystemExit(f"unexpected R outputs under {r_dir}: {unexpected}")
    for rel_path in sorted(nc_validator.EXPECTED_R_OUTPUT_FILES):
        _validate_plot_file(r_dir / rel_path)
    return {
        "r_dir": str(r_dir),
        "file_count": len(files),
        "expected_files": sorted(nc_validator.EXPECTED_R_OUTPUT_FILES),
        "files": files[:50],
    }


def _find_cell_pipeline_output(run_dir: Path) -> Path:
    hits = sorted((run_dir / "python").glob("*/run_manifest.json"))
    if len(hits) != 1:
        raise SystemExit(f"expected exactly one Cell pipeline run_manifest, found {len(hits)}")
    return hits[0].parent


def _parse_shape(value: Any) -> tuple[int, int] | None:
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return int(value[0]), int(value[1])
    if isinstance(value, str):
        normalized = value.lower().replace("×", "x")
        if "x" in normalized:
            left, right = normalized.split("x", 1)
            return int(left.strip()), int(right.strip())
    return None


def _cell_pipeline_shape(manifest: dict[str, object]) -> tuple[int, int]:
    metadata = manifest.get("metadata", {})
    if not isinstance(metadata, dict):
        raise SystemExit("Cell/Trevino pipeline run_manifest metadata must be an object")
    for key in ("final_adata_shape", "final_shape"):
        shape = _parse_shape(metadata.get(key))
        if shape is not None:
            return shape
    cells = metadata.get("n_obs", metadata.get("cells_after_doublet_removal", metadata.get("cells_after_qc")))
    genes = metadata.get("n_vars", metadata.get("genes_after_qc"))
    if cells is None or genes is None:
        raise SystemExit("Cell/Trevino pipeline run_manifest missing final shape metadata")
    return int(cells), int(genes)


def validate_cell_trevino(
    project_root: Path,
    pipeline_run_id: str,
    evidence_run_id: str,
    *,
    verify_sha: bool = False,
) -> dict[str, object]:
    pipeline_run = project_root / "runs" / pipeline_run_id
    evidence_run = project_root / "runs" / evidence_run_id
    output = _find_cell_pipeline_output(pipeline_run)
    manifest = _load_json(output / "run_manifest.json")
    pipeline_shape = _cell_pipeline_shape(manifest)
    if not _shape_within_band(pipeline_shape, EXPECTED_CELL_SHAPE):
        raise SystemExit(
            f"Cell/Trevino linked pipeline shape outside band: {pipeline_shape} "
            f"vs ~{EXPECTED_CELL_SHAPE} "
            f"(±{SHAPE_TOLERANCE_CELLS_PCT}% cells, ±{SHAPE_TOLERANCE_GENES_ABS} genes)"
        )
    module_rows = _csv_rows(output / "module_status.csv")
    statuses = {row["module"]: row["status"] for row in module_rows}
    not_ok = {module: status for module, status in statuses.items() if status != "ok"}
    if not_ok:
        raise SystemExit(f"Cell/Trevino pipeline modules not ok: {not_ok}")

    summary = _load_json(evidence_run / "python" / "pipeline_viability" / "pipeline_viability_summary.json")
    adata = summary.get("adata", {})
    observed_shape = (int(adata.get("n_obs", -1)), int(adata.get("n_vars", -1)))
    if not _shape_within_band(observed_shape, EXPECTED_CELL_SHAPE):
        raise SystemExit(
            f"Cell/Trevino shape outside band: {observed_shape} vs ~{EXPECTED_CELL_SHAPE} "
            f"(±{SHAPE_TOLERANCE_CELLS_PCT}% cells, ±{SHAPE_TOLERANCE_GENES_ABS} genes)"
        )
    if not summary.get("all_requested_modules_ok"):
        raise SystemExit("Cell/Trevino pipeline_viability_summary all_requested_modules_ok is false")

    quality = _load_json(evidence_run / "python" / "figure_reproduction_evidence" / "all_figures_final_quality_gate_review.json")
    if quality.get("decision") != "conditional":
        raise SystemExit(f"unexpected Cell/Trevino quality decision: {quality.get('decision')}")
    if quality.get("unexpected_false_checks") not in ([], None):
        raise SystemExit(f"unexpected Cell/Trevino false checks: {quality.get('unexpected_false_checks')}")
    if quality.get("missing_referenced_paths") not in ([], None):
        raise SystemExit(f"Cell/Trevino missing referenced paths: {quality.get('missing_referenced_paths')}")
    if quality.get("hash_mismatches") not in ([], None):
        raise SystemExit(f"Cell/Trevino hash mismatches: {quality.get('hash_mismatches')}")
    direct_panel_count = quality.get("direct_evidence_panel_count")
    if direct_panel_count is None:
        direct_panel_count = quality.get("evidence_panel_count_direct")
    if direct_panel_count is None:
        raise SystemExit("Cell/Trevino quality gate missing direct evidence panel count")

    final_h5ad = output / "final_adata.h5ad"
    if not final_h5ad.exists():
        raise SystemExit(f"missing Cell/Trevino final_adata: {final_h5ad}")

    return {
        "project_root": str(project_root),
        "pipeline_run_id": pipeline_run_id,
        "evidence_run_id": evidence_run_id,
        "pipeline_output_dir": str(output),
        "run_manifest": str(output / "run_manifest.json"),
        "module_status": {"path": str(output / "module_status.csv"), "counts": dict(Counter(statuses.values())), "statuses": statuses},
        "final_adata": str(final_h5ad),
        "shape": {"n_obs": observed_shape[0], "n_vars": observed_shape[1]},
        "bundle": _validate_bundle(pipeline_run, EXPECTED_CELL_SHAPE, verify_sha=verify_sha),
        "r_outputs": _validate_r_outputs(pipeline_run),
        "quality_gate": {
            "path": str(evidence_run / "python" / "figure_reproduction_evidence" / "all_figures_final_quality_gate_review.json"),
            "decision": quality.get("decision"),
            "evidence_json_count": quality.get("evidence_json_count"),
            "direct_evidence_panel_count": direct_panel_count,
            "expected_false_checks": quality.get("expected_false_checks", []),
            "unexpected_false_checks": quality.get("unexpected_false_checks", []),
        },
        "scientific_boundary": (
            "Conditional public-resource reproduction from GEO processed matrices, author public RDS/tables, "
            "and Cell supplement tables; not FASTQ/fragments/BPNet exact parity."
        ),
    }


def _write_report(report: dict, output_dir: Path) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "two_real_dataset_final_validation.json"
    md_path = output_dir / "REPORT.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(
        "# Two Real-Dataset Final Factory Validation\n\n"
        "Control-plane validation report; scientific outputs remain under project roots.\n\n"
        f"- Generated UTC: `{report['generated_at_utc']}`\n"
        f"- Verdict: `{report['verdict']}`\n"
        "- Datasets: NC2024 NSCLC P15_T1 and Cell/Trevino public-resource reproduction.\n\n"
        "## NC2024\n\n"
        f"- Run: `{report['datasets']['nc2024']['checked_run']['run_dir']}`\n"
        f"- Shape: `{report['datasets']['nc2024']['checked_run']['python']['shape']['n_obs']} x {report['datasets']['nc2024']['checked_run']['python']['shape']['n_vars']}`\n"
        f"- Bundle: `{report['datasets']['nc2024']['checked_run']['bundle']['bundle_dir']}`\n"
        f"- R outputs: `{report['datasets']['nc2024']['checked_run']['r_outputs']['r_dir']}` ({report['datasets']['nc2024']['checked_run']['r_outputs']['file_count']} files)\n"
        "- Boundary: project/module architecture validation, not full article-scale cohort annotation.\n\n"
        "## Cell/Trevino\n\n"
        f"- Pipeline run: `{report['datasets']['cell_trevino']['pipeline_run_id']}`\n"
        f"- Evidence run: `{report['datasets']['cell_trevino']['evidence_run_id']}`\n"
        f"- Shape: `{report['datasets']['cell_trevino']['shape']['n_obs']} x {report['datasets']['cell_trevino']['shape']['n_vars']}`\n"
        f"- Bundle: `{report['datasets']['cell_trevino']['bundle']['bundle_dir']}`\n"
        f"- R outputs: `{report['datasets']['cell_trevino']['r_outputs']['r_dir']}` ({report['datasets']['cell_trevino']['r_outputs']['file_count']} files)\n"
        f"- Quality gate: `{report['datasets']['cell_trevino']['quality_gate']['decision']}` with expected resource gaps documented.\n\n"
        "## Three-Factory Result\n\n"
        "- `singlecell_factory`: upstream run evidence, bundle export, governance validator.\n"
        "- `r_multiomics_factory`: R bundle reading and plotting entrypoint.\n"
        "- `plotting_factory`: source modules used by R plotting and visual smoke tests.\n\n"
        "## Ownership Governance\n\n"
        "- Rule: `owner-by-primary-output`.\n"
        "- `singlecell_factory`: default global control plane for Python-heavy upstream and cross-factory validation.\n"
        "- `r_multiomics_factory`: local owner for R-heavy/spatial primary scientific truth when R produces primary objects and biological interpretation.\n"
        "- `plotting_factory`: presentation-only plotting surface; it must not own biological conclusions.\n"
        "- Canonical policy: `docs/OWNER_BY_PRIMARY_OUTPUT_GOVERNANCE.md`.\n\n"
        "## Hardening Gates\n\n"
        "- Clean-room minimal real-data gate: `cleanroom_minimal_realdata.json`.\n"
        "- Figure parity gate: `figure_parity_gate.json` (conditional until curated paper/process-data references are registered).\n"
        "- CI governance gate: `ci_governance_gate.json`.\n\n"
        "## Remaining Boundaries\n\n"
        "- NC2024 retained P15_T1 is not full multi-cohort paper-scale truth.\n"
        "- Cell/Trevino remains conditional public-resource reproduction, not raw FASTQ/fragments/BPNet parity.\n",
        encoding="utf-8",
    )
    return {"json": str(json_path), "markdown": str(md_path)}


def validate(args: argparse.Namespace) -> dict[str, object]:
    nc = nc_validator.validate(project_root=args.nc_project_root, run_id=args.nc_run_id, verify_sha=args.verify_sha)
    cell = validate_cell_trevino(
        args.cell_project_root,
        args.cell_pipeline_run_id,
        args.cell_evidence_run_id,
        verify_sha=args.verify_sha,
    )
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "PASS_TWO_REAL_DATASET_FACTORY_BRIDGE_VALIDATED",
        "repo_states": {
            "singlecell_factory": _git_state(nc_validator.REPO_ROOT),
            "r_multiomics_factory": _git_state(nc_validator.R_FACTORY_ROOT),
            "plotting_factory": _git_state(nc_validator.PLOTTING_FACTORY_ROOT),
        },
        "datasets": {"nc2024": nc, "cell_trevino": cell},
    }
    if args.output_dir:
        report["report_paths"] = _write_report(report, args.output_dir)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nc-project-root", type=Path, default=nc_validator.DEFAULT_PROJECT_ROOT)
    parser.add_argument("--nc-run-id", type=_safe_run_id, default=nc_validator.DEFAULT_RUN_ID)
    parser.add_argument("--cell-project-root", type=Path, default=DEFAULT_CELL_PROJECT_ROOT)
    parser.add_argument("--cell-pipeline-run-id", type=_safe_run_id, default=DEFAULT_CELL_PIPELINE_RUN_ID)
    parser.add_argument("--cell-evidence-run-id", type=_safe_run_id, default=DEFAULT_CELL_EVIDENCE_RUN_ID)
    parser.add_argument("--verify-sha", action="store_true")
    parser.add_argument("--output-dir", type=_safe_output_dir, default=None)
    args = parser.parse_args()
    print(json.dumps(validate(args), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
