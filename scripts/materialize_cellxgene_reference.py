#!/usr/bin/env python3
"""Materialize a version-pinned CELLxGENE Census reference AnnData with provenance receipt.

Operates in an isolated Census I/O environment (Python 3.12) using official
`cellxgene-census` and `tiledbsoma` APIs.
Outputs are written strictly under the designated project run directory:
    <project-root>/runs/<run-id>/python/reference/
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import logging
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Add repo root to path for project_paths helper
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    from workflow.modular.manifest_writer import factory_git_state
    from workflow.modular.project_paths import generate_run_id, resolve_run_dir, validate_run_id
except ImportError:
    # Fallback if invoked stand-alone
    def generate_run_id(sha: Optional[str] = None) -> str:
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%MZ")
        return f"{ts}-{(sha or '0000000')[:7]}"

    def validate_run_id(rid: str) -> bool:
        return bool(re.match(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{4}Z-[0-9a-f]{7}$", rid))

    def resolve_run_dir(project_root: Path, run_id: Optional[str] = None, factory_sha: Optional[str] = None) -> Path:
        rid = run_id or generate_run_id(factory_sha)
        rd = Path(project_root) / "runs" / rid
        for d in [rd / "python" / "reference", rd / "logs"]:
            d.mkdir(parents=True, exist_ok=True)
        return rd

    def factory_git_state(path: Path) -> dict[str, Any]:
        return {"sha": "unknown", "dirty": False, "diff_sha256": None}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("materialize_census_reference")

_DATE_BUILD_RE = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")


def validate_census_build_date(build_date: str) -> None:
    """Ensure explicit YYYY-MM-DD build date; reject 'latest', 'stable', or empty aliases."""
    b = str(build_date).strip()
    if not b:
        raise ValueError("Census build date cannot be empty. Must be YYYY-MM-DD (e.g. 2025-11-08).")
    if b.lower() in {"latest", "stable", "default"}:
        raise ValueError(
            f"Dynamic alias {build_date!r} is forbidden. "
            "An explicit YYYY-MM-DD build date is required (e.g. 2025-11-08)."
        )
    if not _DATE_BUILD_RE.match(b):
        raise ValueError(f"Census build date must match YYYY-MM-DD format; got {build_date!r}.")


def compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write a JSON receipt atomically in the destination directory."""
    tmp_path = path.with_name(f".{path.name}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp_path, path)


def installed_version(distribution: str, fallback: str = "unknown") -> str:
    """Return installed distribution version without requiring module __version__."""
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return fallback


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Materialize a deterministic CELLxGENE Census reference H5AD"
    )
    parser.add_argument("--project-root", required=True, help="Root path of the project")
    parser.add_argument("--run-id", default=None, help="Optional run ID (YYYY-MM-DDTHHMMZ-<short-sha>)")
    parser.add_argument(
        "--census-build",
        required=True,
        help="Explicit Census build date (YYYY-MM-DD, e.g. 2025-11-08)",
    )
    parser.add_argument("--organism", default="Homo sapiens", help="Organism name (default: Homo sapiens)")
    parser.add_argument(
        "--obs-filter",
        default="is_primary_data == True and tissue_general == 'central nervous system'",
        help="TileDB-SOMA value filter for obs query",
    )
    parser.add_argument("--label-key", default="cell_type", help="Column in obs for cell type labels")
    parser.add_argument("--group-key", default="donor_id", help="Column in obs for donor/sample groups")
    parser.add_argument("--max-cells", type=int, default=12000, help="Maximum total cells to extract")
    parser.add_argument(
        "--max-cells-per-label",
        type=int,
        default=1500,
        help="Maximum cells per label class for balanced stratification",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for deterministic sampling")
    parser.add_argument("--allow-dirty", action="store_true", help="Allow dirty working tree for run ID creation")

    return parser.parse_args(argv)


def select_stratified_soma_joinids(
    obs_df: pd.DataFrame,
    label_key: str,
    max_cells: int = 12000,
    max_cells_per_label: int = 1500,
    seed: int = 42,
) -> tuple[np.ndarray, dict[str, int]]:
    """Deterministically select soma_joinids stratified by label."""
    if "soma_joinid" not in obs_df.columns:
        raise ValueError("obs_df missing required 'soma_joinid' column.")
    if label_key not in obs_df.columns:
        raise ValueError(f"obs_df missing requested label_key {label_key!r}.")

    rng = np.random.default_rng(seed)
    selected_joinids: list[int] = []
    label_counts: dict[str, int] = {}

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

    # Avoid pandas categorical get_group edge cases: filter the normalized
    # string labels directly and sort the closed observed vocabulary.
    for lbl in sorted(label_strings.unique().tolist()):
        group_df = obs_df.loc[label_strings == lbl]
        # Sort by soma_joinid before shuffling for stable permutation
        joinids = np.sort(group_df["soma_joinid"].to_numpy())
        rng.shuffle(joinids)

        take = min(len(joinids), max_cells_per_label)
        chosen = joinids[:take]
        selected_joinids.extend(chosen.tolist())
        label_counts[str(lbl)] = int(take)

    selected_arr = np.array(selected_joinids, dtype=np.int64)
    # If total exceeds max_cells, downsample deterministically
    if len(selected_arr) > max_cells:
        rng.shuffle(selected_arr)
        selected_arr = np.sort(selected_arr[:max_cells])
        # Recompute label counts
        sub_obs = obs_df[obs_df["soma_joinid"].isin(selected_arr)]
        observed_labels = sub_obs[label_key].astype(str)
        label_counts = {str(k): int(v) for k, v in observed_labels.value_counts().items()}
    else:
        selected_arr = np.sort(selected_arr)

    return selected_arr, label_counts


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
) -> dict[str, Any]:
    """Materialize reference H5AD and provenance receipt."""
    validate_census_build_date(census_build)
    if max_cells < 1 or max_cells_per_label < 1:
        raise ValueError(
            f"max_cells and max_cells_per_label must be >= 1; got {max_cells} and {max_cells_per_label}."
        )

    t_start = datetime.now(timezone.utc)
    t0 = time.perf_counter()

    git_state = factory_git_state(_REPO_ROOT)
    run_dir = resolve_run_dir(project_root, run_id=run_id, factory_sha=git_state["sha"] or None)
    effective_run_id = run_dir.name

    ref_dir = run_dir / "python" / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    h5ad_out = ref_dir / "census_reference.h5ad"
    tmp_h5ad_out = ref_dir / ".tmp_census_reference.h5ad"
    joinids_out = ref_dir / "selected_soma_joinids.txt"
    receipt_out = ref_dir / "census_reference_receipt.json"

    # Injected fake client or official cellxgene_census module
    if census_client is None:
        try:
            import cellxgene_census
            import tiledbsoma
            import anndata as ad
        except ImportError as exc:
            err_receipt = {
                "receipt_schema_version": "1.0",
                "status": "failed_missing_dependencies",
                "requested_census_build": census_build,
                "error": str(exc),
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            }
            atomic_write_json(receipt_out, err_receipt)
            raise RuntimeError(f"CELLxGENE Census dependencies not installed in current environment: {exc}") from exc
        client = cellxgene_census
    else:
        client = census_client
        import anndata as ad

    logger.info("Verifying Census version description for build %s...", census_build)
    try:
        ver_desc = client.get_census_version_description(census_build)
    except Exception as exc:
        err_receipt = {
            "receipt_schema_version": "1.0",
            "status": "failed_version_query",
            "requested_census_build": census_build,
            "error_type": exc.__class__.__name__,
            "error_message": str(exc),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        atomic_write_json(receipt_out, err_receipt)
        raise RuntimeError(f"Failed to query Census version description for {census_build}: {exc}") from exc

    resolved_build = (
        ver_desc.get("release_build")
        or ver_desc.get("release_date")
        or ver_desc.get("build_date")
        or ver_desc.get("census_build")
    )
    if resolved_build != census_build:
        err_receipt = {
            "receipt_schema_version": "1.0",
            "status": "failed_build_mismatch",
            "requested_census_build": census_build,
            "resolved_census_build": resolved_build,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        atomic_write_json(receipt_out, err_receipt)
        raise RuntimeError(
            f"Census version description resolved {resolved_build!r} for requested build {census_build!r}."
        )

    logger.info("Opening Census %s and querying obs with filter: %s", census_build, obs_filter)
    try:
        with client.open_soma(census_version=census_build) as census:
            obs_df = client.get_obs(
                census,
                organism=organism,
                value_filter=obs_filter,
                column_names=[
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
                ],
            )
            if obs_df.empty:
                raise ValueError(f"No cells returned from Census query matching filter: {obs_filter}")

            logger.info("Retrieved %d matching obs rows. Selecting stratified cells (cap=%d)...", len(obs_df), max_cells)
            selected_joinids, label_counts = select_stratified_soma_joinids(
                obs_df=obs_df,
                label_key=label_key,
                max_cells=max_cells,
                max_cells_per_label=max_cells_per_label,
                seed=seed,
            )

            logger.info("Materializing AnnData for %d selected cells (raw counts)...", len(selected_joinids))
            adata = client.get_anndata(
                census,
                organism=organism,
                X_name="raw",
                obs_coords=selected_joinids,
                obs_column_names=list(obs_df.columns),
            )

            if "soma_joinid" not in adata.obs.columns:
                raise ValueError("Materialized AnnData is missing required obs column 'soma_joinid'.")
            materialized_joinids = adata.obs["soma_joinid"].to_numpy(dtype=np.int64)
            if len(np.unique(materialized_joinids)) != len(materialized_joinids):
                raise ValueError("Materialized AnnData contains duplicate soma_joinid values.")
            if set(materialized_joinids.tolist()) != set(selected_joinids.tolist()):
                raise ValueError(
                    "Materialized AnnData soma_joinid set does not exactly match the selected coordinate receipt."
                )
            position_by_joinid = {int(joinid): pos for pos, joinid in enumerate(materialized_joinids)}
            adata = adata[[position_by_joinid[int(joinid)] for joinid in selected_joinids]].copy()
    except Exception as exc:
        err_receipt = {
            "receipt_schema_version": "1.0",
            "status": "failed_data_retrieval",
            "requested_census_build": census_build,
            "obs_filter": obs_filter,
            "error_type": exc.__class__.__name__,
            "error_message": str(exc),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        atomic_write_json(receipt_out, err_receipt)
        raise

    # Ensure sparse CSR matrix
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    else:
        adata.X = adata.X.tocsr()

    logger.info("Writing materialized reference to temporary file %s...", tmp_h5ad_out)
    if tmp_h5ad_out.exists():
        tmp_h5ad_out.unlink()
    adata.write_h5ad(tmp_h5ad_out)

    # Verification: re-open and check integrity
    logger.info("Verifying serialized H5AD...")
    test_ad = ad.read_h5ad(tmp_h5ad_out)
    if test_ad.n_obs != len(selected_joinids):
        raise RuntimeError(f"Obs count mismatch: {test_ad.n_obs} vs {len(selected_joinids)}")
    if test_ad.n_vars < 1:
        raise RuntimeError("No genes found in serialized H5AD")
    del test_ad

    # Atomically rename to final path
    os.replace(tmp_h5ad_out, h5ad_out)

    # Write selected joinids list
    tmp_joinids_out = joinids_out.with_name(f".{joinids_out.name}.tmp")
    np.savetxt(tmp_joinids_out, selected_joinids, fmt="%d")
    os.replace(tmp_joinids_out, joinids_out)

    h5ad_sha = compute_sha256(h5ad_out)
    joinids_sha = compute_sha256(joinids_out)

    t_end = datetime.now(timezone.utc)
    duration_sec = round(time.perf_counter() - t0, 2)

    receipt = {
        "receipt_schema_version": "1.0",
        "status": "completed",
        "requested_census_build": census_build,
        "resolved_census_build": str(resolved_build),
        "census_version_description": ver_desc,
        "source_uri": ver_desc.get("soma", {}).get("uri") or "https://cellxgene.cziscience.com",
        "organism": organism,
        "obs_filter": obs_filter,
        "label_key": label_key,
        "group_key": group_key,
        "seed": seed,
        "n_cells_selected": int(len(selected_joinids)),
        "n_genes": int(adata.n_vars),
        "is_sparse": bool(sp.issparse(adata.X)),
        "label_distribution": label_counts,
        "unique_datasets_count": int(adata.obs["dataset_id"].nunique()) if "dataset_id" in adata.obs else 0,
        "unique_donors_count": int(adata.obs[group_key].nunique()) if group_key in adata.obs else 0,
        "artifacts": {
            "reference_h5ad_path": str(h5ad_out),
            "reference_h5ad_sha256": h5ad_sha,
            "reference_h5ad_size_bytes": os.path.getsize(h5ad_out),
            "selected_joinids_path": str(joinids_out),
            "selected_joinids_sha256": joinids_sha,
        },
        "software_provenance": {
            "python_version": sys.version,
            "python_executable": sys.executable,
            "cellxgene_census_version": installed_version(
                "cellxgene-census", "injected_client" if census_client is not None else "unknown"
            ),
            "tiledbsoma_version": installed_version(
                "tiledbsoma", "injected_client" if census_client is not None else "unknown"
            ),
            "anndata_version": installed_version("anndata"),
            "factory_git_sha": git_state["sha"],
            "factory_git_dirty": git_state["dirty"],
            "materializer_script_sha256": compute_sha256(Path(__file__).resolve()),
        },
        "license_and_attribution": {
            "code_license": "MIT",
            "data_license": "CC-BY 4.0",
            "attribution": "CELLxGENE Census (Chan Zuckerberg Initiative)",
            "citation": "https://doi.org/10.1101/2023.05.08.539536",
        },
        "timestamps": {
            "started_at_utc": t_start.isoformat(),
            "completed_at_utc": t_end.isoformat(),
            "duration_seconds": duration_sec,
        },
    }

    atomic_write_json(receipt_out, receipt)
    logger.info("Reference materialization complete. Receipt written to %s", receipt_out)
    return receipt


def main() -> int:
    args = parse_args()
    try:
        materialize_reference(
            project_root=Path(args.project_root),
            census_build=args.census_build,
            organism=args.organism,
            obs_filter=args.obs_filter,
            label_key=args.label_key,
            group_key=args.group_key,
            max_cells=args.max_cells,
            max_cells_per_label=args.max_cells_per_label,
            seed=args.seed,
            run_id=args.run_id,
        )
        return 0
    except Exception as exc:
        logger.error("Materialization failed: %s", exc, exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
