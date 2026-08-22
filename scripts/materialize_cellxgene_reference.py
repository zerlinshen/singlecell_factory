#!/usr/bin/env python3
"""Materialize and verify a version-pinned CELLxGENE Census reference.

The official CELLxGENE Census / TileDB-SOMA client is intentionally confined
to the Python 3.12 ``sc_census_io`` lane.  Every output belongs to a
project-owned run directory and receives a hash-linked, independently
verifiable receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import re
import shlex
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from workflow.modular.manifest_writer import factory_git_state, write_manifest
from workflow.modular.project_paths import resolve_run_dir

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("materialize_census_reference")

_DATE_BUILD_RE = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
_REQUIRED_OBS_COLUMNS = (
    "soma_joinid",
    "dataset_id",
    "donor_id",
    "cell_type",
    "cell_type_ontology_term_id",
    "assay",
    "tissue",
    "tissue_general",
    "development_stage",
    "disease",
    "is_primary_data",
)


def validate_census_build_date(build_date: str) -> None:
    """Reject unresolved Census aliases before a network request is possible."""
    value = str(build_date).strip()
    if not value:
        raise ValueError("Census build date cannot be empty. Must be YYYY-MM-DD (e.g. 2025-11-08).")
    if value.casefold() in {"latest", "stable", "default"}:
        raise ValueError(
            f"Dynamic alias {build_date!r} is forbidden. "
            "An explicit YYYY-MM-DD build date is required (e.g. 2025-11-08)."
        )
    if not _DATE_BUILD_RE.fullmatch(value):
        raise ValueError(f"Census build date must match YYYY-MM-DD format; got {build_date!r}.")


def compute_sha256(file_path: Path) -> str:
    """Return the SHA-256 of a regular file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(file_path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(file_path: Path, payload: dict[str, Any]) -> None:
    """Atomically publish a JSON artifact in its destination directory."""
    destination = Path(file_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, destination)


def atomic_write_text(file_path: Path, text: str) -> None:
    destination = Path(file_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, destination)


def installed_version(distribution: str, fallback: str = "unknown") -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return fallback


def required_observation_columns(label_key: str, group_key: str) -> list[str]:
    """Return base Census columns plus user-selected label/group columns once."""
    selected = [str(label_key).strip(), str(group_key).strip()]
    if not all(selected):
        raise ValueError("label_key and group_key must both be non-empty Census observation columns.")
    return list(dict.fromkeys([*_REQUIRED_OBS_COLUMNS, *selected]))


def _distribution_counts(frame: pd.DataFrame, column: str) -> dict[str, int]:
    if column not in frame.columns:
        return {}
    values = frame[column].astype(str).value_counts().sort_index()
    return {str(key): int(value) for key, value in values.items()}


def _lock_hashes() -> dict[str, str]:
    lock_dir = _REPO_ROOT / "environments"
    lock_paths = {
        "conda_base": lock_dir / "sc_census_io.lock.conda",
        "pip_layer": lock_dir / "sc_census_io.lock.txt",
        "environment_spec": _REPO_ROOT / "environment_census_io.yml",
    }
    missing = [name for name, file_path in lock_paths.items() if not file_path.is_file()]
    if missing:
        raise RuntimeError(f"Required Census environment lock files are missing: {', '.join(missing)}")
    return {name: compute_sha256(file_path) for name, file_path in lock_paths.items()}


def _artifact_record(file_path: Path) -> dict[str, Any]:
    return {
        "path": str(file_path),
        "sha256": compute_sha256(file_path),
        "size_bytes": int(file_path.stat().st_size),
    }


def _sanitize_command(argv: list[str]) -> str:
    """Persist operator command shape without retaining likely secret values."""
    secret_markers = ("token", "secret", "password", "credential", "key")
    sanitized: list[str] = []
    redact_next = False
    for token in argv:
        if redact_next:
            sanitized.append("<redacted>")
            redact_next = False
            continue
        lowered = token.casefold()
        if lowered.startswith("--") and any(marker in lowered for marker in secret_markers):
            if "=" in token:
                sanitized.append(token.split("=", 1)[0] + "=<redacted>")
            else:
                sanitized.append(token)
                redact_next = True
        else:
            sanitized.append(token)
    return shlex.join(sanitized)


def _failure_receipt(
    *,
    status: str,
    census_build: str,
    error: BaseException,
    git_state: dict[str, Any],
    run_id: str,
) -> dict[str, Any]:
    return {
        "receipt_schema_version": "2.0",
        "status": status,
        "requested_census_build": census_build,
        "run_id": run_id,
        "error_type": error.__class__.__name__,
        "error_message": str(error),
        "factory_git_state": git_state,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def select_stratified_soma_joinids(
    obs_df: pd.DataFrame,
    label_key: str,
    max_cells: int = 12000,
    max_cells_per_label: int = 1500,
    seed: int = 42,
) -> tuple[np.ndarray, dict[str, int]]:
    """Select a deterministic label-stratified, coordinate-sorted Census sample."""
    if "soma_joinid" not in obs_df.columns:
        raise ValueError("obs_df missing required 'soma_joinid' column.")
    if label_key not in obs_df.columns:
        raise ValueError(f"obs_df missing requested label_key {label_key!r}.")

    labels = obs_df[label_key]
    if labels.isna().any():
        raise ValueError(f"obs_df label_key {label_key!r} contains missing values.")
    label_strings = labels.astype(str)
    invalid_tokens = {"", "unknown", "unassigned", "na", "n/a", "nan", "none", "null"}
    invalid_mask = label_strings.str.strip().str.casefold().isin(invalid_tokens)
    if invalid_mask.any():
        raise ValueError(
            f"obs_df label_key {label_key!r} contains {int(invalid_mask.sum())} invalid labels."
        )
    if max_cells < 1 or max_cells_per_label < 1:
        raise ValueError("max_cells and max_cells_per_label must be >= 1.")

    rng = np.random.default_rng(seed)
    selected_joinids: list[int] = []
    label_counts: dict[str, int] = {}
    for label in sorted(label_strings.unique().tolist()):
        joinids = np.sort(obs_df.loc[label_strings == label, "soma_joinid"].to_numpy(dtype=np.int64))
        rng.shuffle(joinids)
        chosen = joinids[: min(len(joinids), max_cells_per_label)]
        selected_joinids.extend(chosen.tolist())
        label_counts[str(label)] = int(len(chosen))

    selected_arr = np.asarray(selected_joinids, dtype=np.int64)
    if len(selected_arr) > max_cells:
        rng.shuffle(selected_arr)
        selected_arr = np.sort(selected_arr[:max_cells])
        sampled = obs_df.loc[obs_df["soma_joinid"].isin(selected_arr), label_key].astype(str)
        label_counts = {str(key): int(value) for key, value in sampled.value_counts().sort_index().items()}
    else:
        selected_arr = np.sort(selected_arr)
    return selected_arr, dict(sorted(label_counts.items()))


def _verify_materialization_artifacts(receipt: dict[str, Any]) -> dict[str, Any]:
    """Verify hashes, sizes, sparse H5AD structure, and exact coordinate order."""
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ValueError("Receipt is missing an artifacts object.")

    h5ad_path = Path(artifacts.get("reference_h5ad_path", ""))
    joinids_path = Path(artifacts.get("selected_joinids_path", ""))
    if not h5ad_path.is_file() or not joinids_path.is_file():
        raise FileNotFoundError("Receipt reference H5AD or selected soma_joinids artifact is missing.")

    expected_artifacts = (
        (h5ad_path, "reference_h5ad_sha256", "reference_h5ad_size_bytes"),
        (joinids_path, "selected_joinids_sha256", "selected_joinids_size_bytes"),
    )
    for file_path, hash_key, size_key in expected_artifacts:
        if artifacts.get(hash_key) != compute_sha256(file_path):
            raise ValueError(f"Artifact SHA-256 mismatch for {file_path.name}.")
        if int(artifacts.get(size_key, -1)) != int(file_path.stat().st_size):
            raise ValueError(f"Artifact size mismatch for {file_path.name}.")

    reopened = ad.read_h5ad(h5ad_path)
    if not sp.issparse(reopened.X):
        raise ValueError("Materialized H5AD matrix is not sparse.")
    if reopened.X.getformat() != "csr":
        raise ValueError(f"Materialized H5AD matrix must be CSR, found {reopened.X.getformat()!r}.")
    expected_shape = tuple(int(value) for value in receipt.get("shape", []))
    if expected_shape != tuple(reopened.shape):
        raise ValueError(f"H5AD shape mismatch: receipt={expected_shape}, reopened={tuple(reopened.shape)}.")
    if int(receipt.get("n_cells_selected", -1)) != int(reopened.n_obs):
        raise ValueError("Receipt selected-cell count does not match reopened H5AD.")
    if "soma_joinid" not in reopened.obs.columns:
        raise ValueError("Reopened H5AD is missing soma_joinid.")
    expected_joinids = np.loadtxt(joinids_path, dtype=np.int64, ndmin=1)
    observed_joinids = reopened.obs["soma_joinid"].to_numpy(dtype=np.int64)
    if not np.array_equal(observed_joinids, expected_joinids):
        raise ValueError("Reopened H5AD soma_joinid order does not exactly match the coordinate receipt.")
    return {
        "verified_shape": [int(reopened.n_obs), int(reopened.n_vars)],
        "matrix_storage": "csr",
        "row_order_matches_selected_soma_joinids": True,
    }


def verify_materialization_receipt(receipt_path: Path) -> dict[str, Any]:
    """Verify a completed receipt and its manifest-bound receipt checksum.

    A receipt cannot authenticate itself.  The canonical run manifest therefore
    stores its SHA-256 after artifact verification, so later receipt edits are
    detected rather than silently treated as a new source of truth.
    """
    receipt_path = Path(receipt_path)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") != "completed":
        raise ValueError(f"Receipt status must be 'completed', found {receipt.get('status')!r}.")
    artifact_verification = _verify_materialization_artifacts(receipt)

    manifest_path = receipt_path.parents[2] / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Canonical run manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_binding = manifest.get("extra", {}).get("census_reference", {})
    expected_receipt_sha = manifest_binding.get("receipt_sha256")
    if not expected_receipt_sha:
        raise ValueError("Run manifest does not bind census_reference receipt_sha256.")
    if expected_receipt_sha != compute_sha256(receipt_path):
        raise ValueError("Receipt SHA-256 does not match the canonical manifest binding.")
    if manifest_binding.get("receipt_path") != str(receipt_path):
        raise ValueError("Run manifest receipt path does not match the verified receipt path.")

    return {
        "status": "verified",
        "receipt_path": str(receipt_path),
        "manifest_path": str(manifest_path),
        "receipt_sha256": expected_receipt_sha,
        **artifact_verification,
    }


def _write_run_manifest(
    *,
    project_root: Path,
    run_dir: Path,
    git_state: dict[str, Any],
    receipt_path: Path,
    receipt: dict[str, Any],
) -> Path:
    receipt_sha256 = compute_sha256(receipt_path)
    return write_manifest(
        run_dir,
        project_id=project_root.name,
        run_id=run_dir.name,
        modules_run=["census_reference_materialization"],
        produced_on=platform.node(),
        factory_python_path=_REPO_ROOT,
        factory_python_state=git_state,
        extra={
            "census_reference": {
                "receipt_path": str(receipt_path),
                "receipt_sha256": receipt_sha256,
                "reference_h5ad_sha256": receipt["artifacts"]["reference_h5ad_sha256"],
                "selected_joinids_sha256": receipt["artifacts"]["selected_joinids_sha256"],
                "claim_class": "technical_reference_materialization_not_biological_validation",
            }
        },
    )


def materialize_reference(
    project_root: Path,
    census_build: str,
    organism: str = "Homo sapiens",
    obs_filter: str = "is_primary_data == True and tissue_general == 'central nervous system'",
    label_key: str = "cell_type",
    group_key: str = "donor_id",
    max_cells: int = 12000,
    max_cells_per_label: int = 1500,
    seed: int = 42,
    run_id: Optional[str] = None,
    census_client: Any = None,
    sanitized_command: Optional[str] = None,
) -> dict[str, Any]:
    """Materialize an official Census reference and publish a verified receipt."""
    validate_census_build_date(census_build)
    if max_cells < 1 or max_cells_per_label < 1:
        raise ValueError(
            f"max_cells and max_cells_per_label must be >= 1; got {max_cells} and {max_cells_per_label}."
        )
    obs_columns = required_observation_columns(label_key, group_key)
    project_root = Path(project_root)
    git_state = factory_git_state(_REPO_ROOT)
    run_dir = resolve_run_dir(project_root, run_id=run_id, factory_sha=git_state["sha"] or None)
    reference_dir = run_dir / "python" / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = reference_dir / "census_reference_receipt.json"
    log_path = logs_dir / "materialize_cellxgene_reference.log"
    command_path = logs_dir / "materialize_cellxgene_reference.command.txt"
    if sanitized_command:
        atomic_write_text(command_path, sanitized_command + "\n")

    h5ad_path = reference_dir / "census_reference.h5ad"
    temporary_h5ad_path = reference_dir / ".tmp_census_reference.h5ad"
    joinids_path = reference_dir / "selected_soma_joinids.txt"
    started_at = datetime.now(timezone.utc)
    started_clock = time.perf_counter()
    atomic_write_text(log_path, f"started_at_utc={started_at.isoformat()}\nrequested_build={census_build}\n")

    try:
        if census_client is None:
            try:
                import cellxgene_census
                import tiledbsoma  # noqa: F401 - import validates the pinned direct stack
            except ImportError as exc:
                raise RuntimeError(
                    f"CELLxGENE Census dependencies not installed in current environment: {exc}"
                ) from exc
            client = cellxgene_census
        else:
            client = census_client

        version_description = client.get_census_version_description(census_build)
        resolved_build = (
            version_description.get("release_build")
            or version_description.get("release_date")
            or version_description.get("build_date")
            or version_description.get("census_build")
        )
        if resolved_build != census_build:
            raise RuntimeError(
                f"Census version description resolved {resolved_build!r} for requested build {census_build!r}."
            )

        logger.info("Opening Census %s and reading configured obs columns.", census_build)
        with client.open_soma(census_version=census_build) as census:
            obs_df = client.get_obs(
                census,
                organism=organism,
                value_filter=obs_filter,
                column_names=obs_columns,
            )
            if obs_df.empty:
                raise ValueError(f"No cells returned from Census query matching filter: {obs_filter}")
            missing_columns = [column for column in obs_columns if column not in obs_df.columns]
            if missing_columns:
                raise ValueError(f"Census obs response is missing requested columns: {missing_columns}")

            selected_joinids, label_counts = select_stratified_soma_joinids(
                obs_df,
                label_key=label_key,
                max_cells=max_cells,
                max_cells_per_label=max_cells_per_label,
                seed=seed,
            )
            adata = client.get_anndata(
                census,
                organism=organism,
                X_name="raw",
                obs_coords=selected_joinids,
                obs_column_names=obs_columns,
            )

        if "soma_joinid" not in adata.obs.columns:
            raise ValueError("Materialized AnnData is missing required obs column 'soma_joinid'.")
        materialized_joinids = adata.obs["soma_joinid"].to_numpy(dtype=np.int64)
        if len(np.unique(materialized_joinids)) != len(materialized_joinids):
            raise ValueError("Materialized AnnData contains duplicate soma_joinid values.")
        if set(materialized_joinids.tolist()) != set(selected_joinids.tolist()):
            raise ValueError("Materialized AnnData coordinates do not match the selected soma_joinids.")
        row_by_joinid = {int(joinid): position for position, joinid in enumerate(materialized_joinids)}
        adata = adata[[row_by_joinid[int(joinid)] for joinid in selected_joinids]].copy()
        adata.X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

        if temporary_h5ad_path.exists():
            temporary_h5ad_path.unlink()
        adata.write_h5ad(temporary_h5ad_path)
        os.replace(temporary_h5ad_path, h5ad_path)
        temporary_joinids_path = joinids_path.with_name(f".{joinids_path.name}.tmp")
        np.savetxt(temporary_joinids_path, selected_joinids, fmt="%d")
        os.replace(temporary_joinids_path, joinids_path)

        h5ad_artifact = _artifact_record(h5ad_path)
        joinids_artifact = _artifact_record(joinids_path)
        lock_hashes = _lock_hashes()
        completed_at = datetime.now(timezone.utc)
        receipt = {
            "receipt_schema_version": "2.0",
            "status": "completed",
            "run_id": run_dir.name,
            "requested_census_build": census_build,
            "resolved_census_build": str(resolved_build),
            "census_version_description": version_description,
            "source_uri": version_description.get("soma", {}).get("uri") or "https://cellxgene.cziscience.com",
            "organism": organism,
            "query": {
                "obs_filter": obs_filter,
                "obs_columns": obs_columns,
                "x_name": "raw",
                "max_cells": int(max_cells),
                "max_cells_per_label": int(max_cells_per_label),
            },
            "label_key": label_key,
            "group_key": group_key,
            "seed": int(seed),
            "n_cells_selected": int(len(selected_joinids)),
            "n_genes": int(adata.n_vars),
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "matrix_storage": "csr",
            "label_distribution": label_counts,
            "dataset_counts": _distribution_counts(adata.obs, "dataset_id"),
            "group_counts": _distribution_counts(adata.obs, group_key),
            "schema": {
                "anndata_version": installed_version("anndata"),
                "var_names_name": adata.var_names.name or "var_names",
                "gene_namespace": (
                    "gene_id" if "gene_id" in adata.var.columns else "var_names"
                ),
                "var_columns": sorted(str(column) for column in adata.var.columns),
            },
            "artifacts": {
                "reference_h5ad_path": h5ad_artifact["path"],
                "reference_h5ad_sha256": h5ad_artifact["sha256"],
                "reference_h5ad_size_bytes": h5ad_artifact["size_bytes"],
                "selected_joinids_path": joinids_artifact["path"],
                "selected_joinids_sha256": joinids_artifact["sha256"],
                "selected_joinids_size_bytes": joinids_artifact["size_bytes"],
            },
            "software_provenance": {
                "python_version": sys.version,
                "python_executable": sys.executable,
                "platform": platform.platform(),
                "cellxgene_census_version": installed_version(
                    "cellxgene-census", "injected_client" if census_client is not None else "unknown"
                ),
                "tiledbsoma_version": installed_version(
                    "tiledbsoma", "injected_client" if census_client is not None else "unknown"
                ),
                "anndata_version": installed_version("anndata"),
                "numpy_version": np.__version__,
                "pandas_version": pd.__version__,
                "scipy_version": scipy.__version__,
                "materializer_script_sha256": compute_sha256(Path(__file__).resolve()),
            },
            "environment_locks": lock_hashes,
            "factory_git_state": git_state,
            "license_and_attribution": {
                "code_license": "MIT",
                "data_license": "CC-BY 4.0",
                "attribution": "CELLxGENE Census (Chan Zuckerberg Initiative)",
                "citation": "https://doi.org/10.1101/2023.05.08.539536",
            },
            "claim_class": "technical_reference_materialization_not_biological_validation",
            "timestamps": {
                "started_at_utc": started_at.isoformat(),
                "completed_at_utc": completed_at.isoformat(),
                "duration_seconds": round(time.perf_counter() - started_clock, 3),
            },
        }

        # Verify the exact in-memory receipt claims before exposing a completed
        # receipt, then bind the immutable serialized receipt through manifest.
        _verify_materialization_artifacts(receipt)
        atomic_write_json(receipt_path, receipt)
        manifest_path = _write_run_manifest(
            project_root=project_root,
            run_dir=run_dir,
            git_state=git_state,
            receipt_path=receipt_path,
            receipt=receipt,
        )
        verification = verify_materialization_receipt(receipt_path)
        atomic_write_text(
            log_path,
            (
                f"started_at_utc={started_at.isoformat()}\n"
                f"completed_at_utc={completed_at.isoformat()}\n"
                f"receipt={receipt_path}\nmanifest={manifest_path}\n"
                f"verification={json.dumps(verification, sort_keys=True)}\n"
            ),
        )
        return receipt
    except Exception as exc:
        failure = _failure_receipt(
            status="failed_materialization",
            census_build=census_build,
            error=exc,
            git_state=git_state,
            run_id=run_dir.name,
        )
        atomic_write_json(receipt_path, failure)
        atomic_write_text(
            log_path,
            f"started_at_utc={started_at.isoformat()}\nstatus=failed\nerror={exc.__class__.__name__}: {exc}\n",
        )
        raise


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Materialize or verify a deterministic CELLxGENE Census reference")
    parser.add_argument("--project-root", required=True, help="Project root containing runs/")
    parser.add_argument("--run-id", default=None, help="Run ID for materialization or --verify-only")
    parser.add_argument("--census-build", default=None, help="Explicit Census build YYYY-MM-DD")
    parser.add_argument("--organism", default="Homo sapiens")
    parser.add_argument(
        "--obs-filter",
        default="is_primary_data == True and tissue_general == 'central nervous system'",
    )
    parser.add_argument("--label-key", default="cell_type")
    parser.add_argument("--group-key", default="donor_id")
    parser.add_argument("--max-cells", type=int, default=12000)
    parser.add_argument("--max-cells-per-label", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--verify-only", action="store_true", help="Verify a completed receipt without network access")
    parser.add_argument("--receipt-path", default="", help="Optional explicit receipt path for --verify-only")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    project_root = Path(args.project_root)
    try:
        if args.verify_only:
            if args.receipt_path:
                receipt_path = Path(args.receipt_path)
            elif args.run_id:
                receipt_path = project_root / "runs" / args.run_id / "python" / "reference" / "census_reference_receipt.json"
            else:
                raise ValueError("--verify-only requires --receipt-path or --run-id.")
            print(json.dumps(verify_materialization_receipt(receipt_path), indent=2, sort_keys=True))
            return 0
        if not args.census_build:
            raise ValueError("--census-build is required unless --verify-only is used.")
        materialize_reference(
            project_root=project_root,
            census_build=args.census_build,
            organism=args.organism,
            obs_filter=args.obs_filter,
            label_key=args.label_key,
            group_key=args.group_key,
            max_cells=args.max_cells,
            max_cells_per_label=args.max_cells_per_label,
            seed=args.seed,
            run_id=args.run_id,
            sanitized_command=_sanitize_command([sys.executable, str(Path(__file__).resolve()), *(argv or sys.argv[1:])]),
        )
        return 0
    except Exception as exc:
        logger.error("Census materialization/verification failed: %s", exc, exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
