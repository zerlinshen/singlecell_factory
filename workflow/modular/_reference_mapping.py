"""Reference mapping, gene alignment, OOD calibration, and deterministic KNN transfer.

Provides:
- One-to-one gene alignment across query and reference AnnData objects.
- Sparse L2 normalization.
- Reference-only OOD distance calibration (quantile or fixed mode).
- Deterministic KNN voting with stable tie-breaking and closed status vocabulary.
- Explicit CPU (scikit-learn) and GPU (cuML RAPIDS) execution.
- Analytical parity and performance metric calculations.
"""

from __future__ import annotations

import logging
import importlib
import importlib.metadata
import json
import hashlib
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import f1_score
from sklearn.neighbors import NearestNeighbors as SklearnNearestNeighbors
from sklearn.preprocessing import normalize

logger = logging.getLogger(__name__)

# Closed vocabulary for reference assignment status
ASSIGNMENT_STATUS_ACCEPTED = "accepted"
ASSIGNMENT_STATUS_REJECTED_LOW_CONF = "rejected_low_confidence"
ASSIGNMENT_STATUS_REJECTED_OOD_DIST = "rejected_ood_distance"
ASSIGNMENT_STATUS_REJECTED_BOTH = "rejected_low_confidence_and_ood_distance"

def resolve_reference_device(
    requested_device: str,
    *,
    gpu_mode: str = "auto",
) -> tuple[str, dict[str, Any]]:
    """Resolve an explicit CPU/GPU request without certificate-driven routing.

    The offline validation registry is evidence-only.  It must never select a
    production backend, so this function intentionally accepts only ``cpu`` or
    ``gpu`` and performs the same fail-closed GPU preflight for CLI and
    programmatic callers.
    """
    requested = str(requested_device).strip().lower()
    if requested not in {"cpu", "gpu"}:
        raise ValueError(f"reference device must be cpu or gpu; got {requested_device!r}")
    resolved_gpu_mode = str(gpu_mode or "auto").strip().lower()
    receipt: dict[str, Any] = {
        "requested_device": requested,
        "gpu_mode": resolved_gpu_mode,
        "decision_basis": "explicit_operator_selection",
        "validation_registry_role": "evidence_only_not_runtime_routing",
    }

    if requested == "cpu":
        receipt.update({"effective_device": "cpu", "decision": "explicit_cpu"})
        return "cpu", receipt

    if resolved_gpu_mode == "off":
        raise ValueError("reference_device='gpu' cannot be used when gpu_mode='off'.")
    try:
        cp = importlib.import_module("cupy")
        cuml = importlib.import_module("cuml")
        importlib.import_module("cuml.neighbors")
    except Exception as exc:
        raise RuntimeError(
            "Explicit GPU reference mapping requires importable CuPy and cuML. "
            "No CPU fallback is allowed."
        ) from exc

    try:
        device_count = int(cp.cuda.runtime.getDeviceCount())
        if device_count < 1:
            raise RuntimeError("CuPy reports zero CUDA devices.")
        free_bytes, total_bytes = cp.cuda.runtime.memGetInfo()
        properties = cp.cuda.runtime.getDeviceProperties(0)
        raw_name = properties.get("name", b"unknown")
        device_name = (
            raw_name.decode("utf-8", errors="replace")
            if isinstance(raw_name, bytes)
            else str(raw_name)
        )
    except Exception as exc:
        raise RuntimeError(
            "Explicit GPU reference mapping requires an available CUDA device. "
            "No CPU fallback is allowed."
        ) from exc

    receipt.update({
        "effective_device": "gpu",
        "decision": "explicit_gpu_preflight_passed",
        "cuml_version": str(getattr(cuml, "__version__", "unknown")),
        "cuda_device_count": device_count,
        "cuda_device": device_name,
        "cuda_free_mb": round(float(free_bytes) / (1024.0 ** 2), 2),
        "cuda_total_mb": round(float(total_bytes) / (1024.0 ** 2), 2),
        "cuda_preflight_verified": True,
    })
    return "gpu", receipt

VALID_ASSIGNMENT_STATUSES = {
    ASSIGNMENT_STATUS_ACCEPTED,
    ASSIGNMENT_STATUS_REJECTED_LOW_CONF,
    ASSIGNMENT_STATUS_REJECTED_OOD_DIST,
    ASSIGNMENT_STATUS_REJECTED_BOTH,
}


def _stable_json_sha256(value: Any) -> str:
    """Hash JSON-safe calibration identities with a canonical representation."""
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def normalize_l2(matrix: Union[sp.spmatrix, np.ndarray]) -> Union[sp.csr_matrix, np.ndarray]:
    """Perform row-wise L2 normalization preserving sparsity."""
    if sp.issparse(matrix):
        csr = matrix.tocsr().astype(np.float32)
        return normalize(csr, norm="l2", axis=1, copy=False)
    arr = np.asarray(matrix, dtype=np.float32)
    return normalize(arr, norm="l2", axis=1, copy=False)


def align_reference_genes(
    query_adata: Any,
    ref_adata: Any,
    min_shared: int = 50,
) -> tuple[list[str], list[str], str, int, float, float]:
    """Perform unambiguous 1-to-1 gene alignment between query and reference.

    Routes:
    1. Exact unique `var_names` intersection.
    2. Query stable ID (`gene_id` or `ensembl_id`) vs Reference ID (`feature_id`, `gene_id`, `ensembl_id`).
    3. Canonical symbol match.

    Returns:
        tuple: (
            query_genes (list of str in query.var_names order or matched),
            ref_genes (list of str in ref.var_names aligned 1:1),
            alignment_route (str),
            n_shared (int),
            query_fraction (float),
            ref_fraction (float),
        )
    """
    q_vars = pd.Index(query_adata.var_names.astype(str))
    r_vars = pd.Index(ref_adata.var_names.astype(str))

    if not q_vars.is_unique or not r_vars.is_unique:
        raise ValueError(
            "Reference mapping requires unique query and reference var_names; "
            "make feature identifiers unique before requesting mapping."
        )

    # Route 1: Exact var_names
    shared_direct = q_vars.intersection(r_vars)
    if len(shared_direct) >= min_shared:
        # Preserve query var_names order for determinism
        aligned = [g for g in q_vars if g in shared_direct]
        n_shared = len(aligned)
        q_frac = float(n_shared / max(len(q_vars), 1))
        r_frac = float(n_shared / max(len(r_vars), 1))
        return aligned, aligned, "exact_var_names", n_shared, q_frac, r_frac

    # Route 2: Stable feature/gene IDs in var tables
    q_id_cols = [c for c in ["gene_id", "ensembl_id", "feature_id"] if c in query_adata.var.columns]
    r_id_cols = [c for c in ["feature_id", "gene_id", "ensembl_id"] if c in ref_adata.var.columns]

    if q_id_cols and r_id_cols:
        q_col = q_id_cols[0]
        r_col = r_id_cols[0]
        q_map = query_adata.var[q_col].dropna().astype(str)
        r_map = ref_adata.var[r_col].dropna().astype(str)

        if not q_map.is_unique or not r_map.is_unique:
            raise ValueError(
                f"Stable-ID alignment is ambiguous because {q_col!r} or {r_col!r} contains duplicate values; "
                "clean the identifier columns explicitly instead of falling back to symbols."
            )

        q_rev = pd.Series(q_map.index.values, index=q_map.values)
        r_rev = pd.Series(r_map.index.values, index=r_map.values)
        common_ids = q_rev.index.intersection(r_rev.index)
        if len(common_ids) >= min_shared:
            q_matched = q_rev.loc[common_ids].tolist()
            r_matched = r_rev.loc[common_ids].tolist()
            n_shared = len(common_ids)
            q_frac = float(n_shared / max(len(q_vars), 1))
            r_frac = float(n_shared / max(len(r_vars), 1))
            return q_matched, r_matched, f"stable_ids_{q_col}_{r_col}", n_shared, q_frac, r_frac

    # Route 3: Uppercase symbol match fallback
    q_upper = pd.Series(q_vars.values, index=q_vars.str.upper())
    r_upper = pd.Series(r_vars.values, index=r_vars.str.upper())
    if q_upper.index.is_unique and r_upper.index.is_unique:
        common_upper = q_upper.index.intersection(r_upper.index)
        if len(common_upper) >= min_shared:
            q_matched = q_upper.loc[common_upper].tolist()
            r_matched = r_upper.loc[common_upper].tolist()
            n_shared = len(common_upper)
            q_frac = float(n_shared / max(len(q_vars), 1))
            r_frac = float(n_shared / max(len(r_vars), 1))
            return q_matched, r_matched, "symbol_upper", n_shared, q_frac, r_frac

    raise ValueError(
        f"Insufficient unambiguous shared genes between query ({len(q_vars)}) and reference ({len(r_vars)}). "
        f"Found {len(shared_direct)} direct matches (minimum required: {min_shared})."
    )


def calibrate_reference_ood_threshold(
    ref_adata: Any,
    ref_genes: list[str],
    mode: str = "reference_quantile",
    quantile: float = 0.95,
    fixed_threshold: Optional[float] = None,
    group_key: Optional[str] = None,
    cal_fraction: float = 0.2,
    seed: int = 42,
    k: int = 15,
) -> tuple[float, dict[str, Any]]:
    """Derive reference-only OOD cosine distance threshold.

    Guarantees:
    - Never uses query cells or held-out query labels.
    - Reference-quantile mode requires a whole-group split to prevent intra-group leakage.
    - If fixed mode is chosen, validates and records explicit threshold.
    """
    mode = str(mode).strip().lower()
    if int(k) < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if not (0.0 < float(cal_fraction) < 1.0):
        raise ValueError(f"cal_fraction must be in (0.0, 1.0), got {cal_fraction}")
    if not ref_genes:
        raise ValueError("Reference OOD calibration requires at least one aligned feature.")
    if mode == "fixed":
        if fixed_threshold is None:
            raise ValueError("OOD mode 'fixed' requires an explicit float fixed_threshold.")
        thresh = float(fixed_threshold)
        if not (0.0 <= thresh <= 2.0):
            raise ValueError(f"fixed_threshold must be in [0.0, 2.0], got {thresh}")
        receipt = {
            "ood_calibration_mode": "fixed",
            "distance_threshold": thresh,
            "distance_threshold_display": round(thresh, 6),
            "quantile": None,
            "group_key": None,
            "cal_fraction": None,
            "seed": seed,
            "n_reference_total": int(ref_adata.n_obs),
            "evidence": "operator_fixed_contract",
        }
        return thresh, receipt

    if mode != "reference_quantile":
        raise ValueError(f"Unsupported OOD calibration mode: {mode}. Must be 'reference_quantile' or 'fixed'.")

    if not (0.0 < float(quantile) < 1.0):
        raise ValueError(f"Quantile must be in (0.0, 1.0), got {quantile}")

    if not group_key:
        raise ValueError(
            "Reference-quantile calibration requires an explicit whole-group obs key "
            "(for example sample or donor_id); random cell partitioning is forbidden."
        )

    n_cells = ref_adata.n_obs
    if n_cells < 20:
        raise ValueError(
            f"Reference quantile calibration requires at least 20 cells; got {n_cells}. "
            "Use a larger reference or an explicit fixed threshold contract."
        )

    rng = np.random.default_rng(seed)
    fit_indices: np.ndarray
    cal_indices: np.ndarray

    if group_key and group_key not in ref_adata.obs.columns:
        raise ValueError(f"Requested reference calibration group key {group_key!r} is absent from reference .obs.")

    groups = np.array(sorted(ref_adata.obs[group_key].astype(str).unique().tolist()))
    if len(groups) >= 2:
        rng.shuffle(groups)
        n_cal_groups = max(1, int(round(len(groups) * cal_fraction)))
        cal_groups = sorted(str(value) for value in groups[:n_cal_groups])
        fit_groups = sorted(str(value) for value in groups[n_cal_groups:])
        if not fit_groups:
            raise ValueError(
                "Reference calibration split left no fit groups; use a smaller calibration fraction."
            )
        cal_group_set = set(cal_groups)
        fit_group_set = set(fit_groups)
        if (cal_group_set & fit_group_set) or ((cal_group_set | fit_group_set) != set(groups.astype(str))):
            raise RuntimeError("Reference calibration group split is not a disjoint exhaustive partition.")
        group_series = ref_adata.obs[group_key].astype(str)
        fit_mask = group_series.isin(fit_groups)
        cal_mask = group_series.isin(cal_groups)
        fit_indices = np.where(fit_mask)[0]
        cal_indices = np.where(cal_mask)[0]
        split_method = f"whole_group_split_on_{group_key}"
    else:
        raise ValueError(
            f"Requested whole-group calibration on {group_key!r} requires at least two groups; got {len(groups)}."
        )

    if len(fit_indices) < 5 or len(cal_indices) < 5:
        raise ValueError(
            f"Reference calibration partition is too small (fit={len(fit_indices)}, cal={len(cal_indices)}); "
            "use a larger reference or an explicit fixed threshold contract."
        )
    if int(k) > len(fit_indices):
        raise ValueError(
            f"Requested k={k} exceeds the reference calibration fit partition ({len(fit_indices)} cells); "
            "k is never reduced silently."
        )

    # Extract submatrices for calibration
    ref_sub = ref_adata[:, ref_genes]
    fit_x = normalize_l2(ref_sub.X[fit_indices])
    cal_x = normalize_l2(ref_sub.X[cal_indices])

    nn = SklearnNearestNeighbors(n_neighbors=int(k), metric="cosine", algorithm="brute")
    nn.fit(fit_x)
    cal_distances, _ = nn.kneighbors(cal_x)
    mean_cal_dists = np.mean(cal_distances, axis=1)

    derived_threshold = float(np.quantile(mean_cal_dists, float(quantile)))

    group_counts = {
        str(group): int(count)
        for group, count in group_series.value_counts().sort_index().items()
    }
    fit_group_counts = {group: group_counts[group] for group in fit_groups}
    cal_group_counts = {group: group_counts[group] for group in cal_groups}

    receipt = {
        "ood_calibration_mode": "reference_quantile",
        # Mapping consumes this full-precision value.  The display value is
        # explicitly separate so a rendered receipt cannot change behavior.
        "distance_threshold": derived_threshold,
        "distance_threshold_display": round(derived_threshold, 6),
        "quantile": float(quantile),
        "split_method": split_method,
        "group_key": group_key,
        "cal_fraction": float(cal_fraction),
        "seed": seed,
        "n_reference_total": int(n_cells),
        "n_fit_cells": int(len(fit_indices)),
        "n_cal_cells": int(len(cal_indices)),
        "fit_groups": fit_groups,
        "calibration_groups": cal_groups,
        "fit_groups_sha256": _stable_json_sha256(fit_groups),
        "calibration_groups_sha256": _stable_json_sha256(cal_groups),
        "fit_group_counts": fit_group_counts,
        "calibration_group_counts": cal_group_counts,
        "all_group_counts": group_counts,
        "n_aligned_features": int(len(ref_genes)),
        "cal_distance_min": float(np.min(mean_cal_dists)),
        "cal_distance_median": float(np.median(mean_cal_dists)),
        "cal_distance_mean": float(np.mean(mean_cal_dists)),
        "cal_distance_max": float(np.max(mean_cal_dists)),
        "cal_distance_std": float(np.std(mean_cal_dists)),
    }
    return derived_threshold, receipt


def map_knn_reference(
    query_adata: Any,
    ref_adata: Any,
    label_key: str,
    k: int = 15,
    min_confidence: float = 0.6,
    distance_threshold: float = 0.5,
    device: str = "cpu",
    seed: int = 42,
    min_shared_genes: int = 50,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Execute deterministic KNN reference mapping with strict OOD status classification.

    Supports:
    - CPU: sklearn NearestNeighbors(metric="cosine", algorithm="brute")
    - GPU: cuML NearestNeighbors(metric="cosine", algorithm="brute")

    Guarantees:
    - GPU failure raises RuntimeError; no silent CPU fallback.
    - Returns full per-cell diagnostic predictions and status classifications.
    - Tie-breaking is deterministic (by vote count desc, label asc).
    """
    device = str(device).strip().lower()
    if device not in {"cpu", "gpu"}:
        raise ValueError(f"device must be 'cpu' or 'gpu', got {device!r}")

    if label_key not in ref_adata.obs.columns:
        raise ValueError(f"Requested reference label_key {label_key!r} not found in reference .obs columns.")

    if query_adata.n_obs < 1 or ref_adata.n_obs < 1:
        raise ValueError(
            f"Reference mapping requires non-empty query and reference objects; "
            f"got query={query_adata.n_obs}, reference={ref_adata.n_obs}."
        )

    raw_labels = ref_adata.obs[label_key]
    if raw_labels.isna().any():
        raise ValueError(f"Reference label column {label_key!r} contains missing labels.")
    ref_labels_series = raw_labels.astype(str)
    invalid_tokens = {"", "unknown", "unassigned", "na", "n/a", "nan", "none", "null"}
    invalid_label_mask = ref_labels_series.str.strip().str.casefold().isin(invalid_tokens)
    if invalid_label_mask.any():
        raise ValueError(
            f"Reference label column {label_key!r} contains {int(invalid_label_mask.sum())} invalid labels; "
            "clean the reference explicitly before mapping."
        )
    unique_labels = ref_labels_series.unique()
    valid_labels = [lbl for lbl in unique_labels if lbl.strip().casefold() not in invalid_tokens]
    if len(valid_labels) < 2:
        raise ValueError(f"Reference label column {label_key!r} must contain at least 2 distinct valid labels; got {valid_labels}")

    q_genes, r_genes, route, n_shared, q_frac, r_frac = align_reference_genes(
        query_adata, ref_adata, min_shared=min_shared_genes
    )

    query_sub = query_adata[:, q_genes]
    ref_sub = ref_adata[:, r_genes]

    t0_prep = time.perf_counter()
    query_x = normalize_l2(query_sub.X)
    ref_x = normalize_l2(ref_sub.X)
    t_prep = time.perf_counter() - t0_prep

    effective_k = int(k)
    if effective_k < 1:
        raise ValueError(f"k must be >= 1, got {k}.")
    if effective_k > int(ref_adata.n_obs):
        raise ValueError(
            f"Requested k={effective_k} exceeds reference size ({ref_adata.n_obs} cells); "
            "k is never reduced silently."
        )
    min_conf = float(min_confidence)
    dist_thresh = float(distance_threshold)
    if not (0.0 <= min_conf <= 1.0):
        raise ValueError(f"min_confidence must be in [0.0, 1.0], got {min_confidence}")
    if not (0.0 <= dist_thresh <= 2.0):
        raise ValueError(f"distance_threshold must be in [0.0, 2.0], got {distance_threshold}")

    ref_labels = ref_labels_series.to_numpy()
    timing_metrics: dict[str, float] = {"prep_sec": t_prep}
    device_info: dict[str, Any] = {
        "device": device,
        "backend_residency": "cuda" if device == "gpu" else "host_cpu",
        "cuda_residency_verified": False,
        "numpy_version": importlib.metadata.version("numpy"),
        "scipy_version": importlib.metadata.version("scipy"),
    }

    if device == "gpu":
        # Ensure CUDA_PATH is set for CuPy JIT compilation if targets directory exists
        for candidate_dir in [
            Path(sys.prefix) / "targets" / "x86_64-linux",
            Path(os.environ.get("CONDA_PREFIX", "")) / "targets" / "x86_64-linux",
        ]:
            if candidate_dir.exists():
                os.environ.setdefault("CUDA_PATH", str(candidate_dir))
                break

        try:
            import cupy as cp
            import cupyx.scipy.sparse as cp_sp
            import cuml.neighbors
            from cuml.neighbors import NearestNeighbors as CuMLNearestNeighbors
        except ImportError as exc:
            raise RuntimeError(
                f"GPU reference mapping requested but cuML/CuPy cannot be imported: {exc}. "
                "No silent CPU fallback allowed."
            ) from exc

        try:
            free_before, total_vram = cp.cuda.runtime.memGetInfo()
            used_before = total_vram - free_before
            # Transfer to GPU device
            t0_h2d = time.perf_counter()
            if sp.issparse(ref_x):
                ref_x_dev = cp_sp.csr_matrix(ref_x)
                query_x_dev = cp_sp.csr_matrix(query_x)
            else:
                ref_x_dev = cp.asarray(ref_x, dtype=cp.float32)
                query_x_dev = cp.asarray(query_x, dtype=cp.float32)
            cp.cuda.Stream.null.synchronize()
            t_h2d = time.perf_counter() - t0_h2d

            # Fit and Query on GPU
            t0_fit_query = time.perf_counter()
            gpu_nn = CuMLNearestNeighbors(n_neighbors=effective_k, metric="cosine", algorithm="brute", output_type="cupy")
            gpu_nn.fit(ref_x_dev)
            dists_dev, indices_dev = gpu_nn.kneighbors(query_x_dev)
            cp.cuda.Stream.null.synchronize()
            t_fit_query = time.perf_counter() - t0_fit_query

            # Transfer back to Host CPU
            t0_d2h = time.perf_counter()
            distances = cp.asnumpy(dists_dev)
            indices = cp.asnumpy(indices_dev)
            t_d2h = time.perf_counter() - t0_d2h
            free_after, _ = cp.cuda.runtime.memGetInfo()
            used_after = total_vram - free_after

            timing_metrics.update({
                "h2d_sec": t_h2d,
                "fit_query_sec": t_fit_query,
                "d2h_sec": t_d2h,
                "total_gpu_sec": t_h2d + t_fit_query + t_d2h,
            })
            device_info.update({
                "cuml_version": cuml.__version__,
                "backend_class": f"{gpu_nn.__class__.__module__}.{gpu_nn.__class__.__qualname__}",
                "cuda_device": cp.cuda.runtime.getDeviceProperties(0)["name"].decode("utf-8", errors="ignore"),
                "vram_total_mb": round(float(total_vram) / (1024.0 ** 2), 2),
                "vram_used_before_mb": round(float(used_before) / (1024.0 ** 2), 2),
                "vram_used_after_query_mb": round(float(used_after) / (1024.0 ** 2), 2),
                "vram_delta_mb": round(float(used_after - used_before) / (1024.0 ** 2), 2),
                "cuda_residency_verified": True,
            })
        except Exception as exc:
            raise RuntimeError(f"cuML GPU KNN execution failed: {exc}. No silent CPU fallback allowed.") from exc

    else:
        # CPU sklearn execution
        t0_fit = time.perf_counter()
        cpu_nn = SklearnNearestNeighbors(n_neighbors=effective_k, metric="cosine", algorithm="brute", n_jobs=-1)
        cpu_nn.fit(ref_x)
        distances, indices = cpu_nn.kneighbors(query_x)
        t_fit_query = time.perf_counter() - t0_fit
        timing_metrics["fit_query_sec"] = t_fit_query
        timing_metrics["total_cpu_sec"] = t_fit_query
        device_info.update({
            "backend_class": f"{cpu_nn.__class__.__module__}.{cpu_nn.__class__.__qualname__}",
            "n_jobs": -1,
            "scikit_learn_version": importlib.metadata.version("scikit-learn"),
        })

    # Derive predictions, confidences, mean distances, and status
    pred_labels: list[str] = []
    confidences: list[float] = []
    mean_distances: list[float] = []
    statuses: list[str] = []
    ref_cell_types: list[str] = []
    ood_flags: list[bool] = []

    for i in range(len(indices)):
        row_indices = indices[i]
        row_dists = distances[i]
        neighbor_labels = ref_labels[row_indices]

        # Voting with deterministic tie-breaking (by count desc, then label asc)
        vals, counts = np.unique(neighbor_labels, return_counts=True)
        sort_keys = sorted(zip(vals, counts), key=lambda item: (-item[1], str(item[0])))
        top_label, top_count = sort_keys[0]

        conf = float(top_count) / float(effective_k)
        mean_dist = float(np.mean(row_dists))

        is_conf = conf >= min_conf
        is_in_dist = mean_dist <= dist_thresh
        is_ood = mean_dist > dist_thresh

        if is_conf and is_in_dist:
            status = ASSIGNMENT_STATUS_ACCEPTED
            final_type = str(top_label)
        elif (not is_conf) and is_in_dist:
            status = ASSIGNMENT_STATUS_REJECTED_LOW_CONF
            final_type = "Unknown"
        elif is_conf and (not is_in_dist):
            status = ASSIGNMENT_STATUS_REJECTED_OOD_DIST
            final_type = "Unknown"
        else:
            status = ASSIGNMENT_STATUS_REJECTED_BOTH
            final_type = "Unknown"

        pred_labels.append(str(top_label))
        confidences.append(conf)
        mean_distances.append(mean_dist)
        statuses.append(status)
        ref_cell_types.append(final_type)
        ood_flags.append(bool(is_ood))

    df = pd.DataFrame(
        {
            "reference_predicted_label": pred_labels,
            "reference_confidence": confidences,
            "reference_distance": mean_distances,
            "reference_assignment_status": statuses,
            "reference_cell_type": ref_cell_types,
            "reference_ood": ood_flags,
        },
        index=query_adata.obs_names,
    )

    n_total = len(df)
    n_accepted = int((df["reference_assignment_status"] == ASSIGNMENT_STATUS_ACCEPTED).sum())
    n_rej_conf = int((df["reference_assignment_status"] == ASSIGNMENT_STATUS_REJECTED_LOW_CONF).sum())
    n_rej_dist = int((df["reference_assignment_status"] == ASSIGNMENT_STATUS_REJECTED_OOD_DIST).sum())
    n_rej_both = int((df["reference_assignment_status"] == ASSIGNMENT_STATUS_REJECTED_BOTH).sum())
    n_total_rej = n_total - n_accepted

    summary_metadata = {
        "status": "completed",
        "alignment_route": route,
        "n_shared_genes": n_shared,
        "query_gene_fraction": float(q_frac),
        "ref_gene_fraction": float(r_frac),
        "requested_k": int(k),
        "effective_k": int(effective_k),
        "min_confidence": min_conf,
        "distance_threshold": dist_thresh,
        "distance_threshold_display": round(dist_thresh, 6),
        "device_info": device_info,
        "timings": timing_metrics,
        "total_query_cells": n_total,
        "accepted_cells": n_accepted,
        "accepted_pct": round(100.0 * n_accepted / max(n_total, 1), 2),
        "rejected_cells": n_total_rej,
        "rejected_pct": round(100.0 * n_total_rej / max(n_total, 1), 2),
        "rejected_low_confidence": n_rej_conf,
        "rejected_ood_distance": n_rej_dist,
        "rejected_low_confidence_and_ood_distance": n_rej_both,
    }

    return df, summary_metadata


def compute_parity_metrics(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    ground_truth: Optional[pd.Series] = None,
    ood_label: Optional[str] = None,
) -> dict[str, Any]:
    """Calculate rigorous analytical parity metrics between two mapping results (e.g. CPU vs GPU)."""
    if len(df_a) != len(df_b):
        raise ValueError(f"Length mismatch: {len(df_a)} vs {len(df_b)}")
    if not df_a.index.equals(df_b.index):
        raise ValueError("Index mismatch between mapping DataFrames")

    n_cells = len(df_a)

    # 1. Prediction Agreement
    pred_match = (df_a["reference_predicted_label"] == df_b["reference_predicted_label"]).sum()
    pred_agreement = float(pred_match / max(n_cells, 1))

    # 2. Status Agreement
    status_match = (df_a["reference_assignment_status"] == df_b["reference_assignment_status"]).sum()
    status_agreement = float(status_match / max(n_cells, 1))

    # 3. Reference Cell Type Agreement
    cell_type_match = (df_a["reference_cell_type"] == df_b["reference_cell_type"]).sum()
    cell_type_agreement = float(cell_type_match / max(n_cells, 1))

    # 4. Rejection Rates
    rej_a = float((df_a["reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED).mean())
    rej_b = float((df_b["reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED).mean())
    delta_rejection_rate = abs(rej_a - rej_b)

    # 5. Distance comparison
    dists_a = df_a["reference_distance"].to_numpy(dtype=float)
    dists_b = df_b["reference_distance"].to_numpy(dtype=float)
    dist_max_abs_diff = float(np.max(np.abs(dists_a - dists_b)))
    dist_mean_abs_diff = float(np.mean(np.abs(dists_a - dists_b)))
    dist_close = bool(np.allclose(dists_a, dists_b, rtol=1e-4, atol=1e-4))

    metrics: dict[str, Any] = {
        "n_cells": n_cells,
        "candidate_label_agreement": round(pred_agreement, 6),
        "assignment_status_agreement": round(status_agreement, 6),
        "reference_cell_type_agreement": round(cell_type_agreement, 6),
        "rejection_rate_a": round(rej_a, 6),
        "rejection_rate_b": round(rej_b, 6),
        "delta_rejection_rate": round(delta_rejection_rate, 6),
        "distance_max_abs_diff": round(dist_max_abs_diff, 8),
        "distance_mean_abs_diff": round(dist_mean_abs_diff, 8),
        "distance_allclose_1e4": dist_close,
    }

    # Biological / Falsifiability metrics if ground_truth is provided
    if ground_truth is not None:
        gt = ground_truth.loc[df_a.index].astype(str)

        # Held-out OOD metrics
        if ood_label is not None:
            ood_mask = gt == ood_label
            n_ood = int(ood_mask.sum())
            if n_ood > 0:
                rej_ood_a = float((df_a.loc[ood_mask, "reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED).mean())
                rej_ood_b = float((df_b.loc[ood_mask, "reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED).mean())
                metrics.update({
                    "n_held_out_ood_cells": n_ood,
                    "ood_rejection_recall_a": round(rej_ood_a, 6),
                    "ood_rejection_recall_b": round(rej_ood_b, 6),
                    "delta_ood_rejection_recall": round(abs(rej_ood_a - rej_ood_b), 6),
                })

            known_mask = ~ood_mask
            n_known = int(known_mask.sum())
            if n_known > 0:
                rej_known_a = float((df_a.loc[known_mask, "reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED).mean())
                rej_known_b = float((df_b.loc[known_mask, "reference_assignment_status"] != ASSIGNMENT_STATUS_ACCEPTED).mean())
                metrics.update({
                    "n_known_cells": n_known,
                    "known_rejection_rate_a": round(rej_known_a, 6),
                    "known_rejection_rate_b": round(rej_known_b, 6),
                    "known_acceptance_coverage_a": round(1.0 - rej_known_a, 6),
                    "known_acceptance_coverage_b": round(1.0 - rej_known_b, 6),
                })

                # Macro-F1 on accepted known cells vs ground truth
                for label, df, key in [("a", df_a, "macro_f1_known_a"), ("b", df_b, "macro_f1_known_b")]:
                    acc_mask = known_mask & (df["reference_assignment_status"] == ASSIGNMENT_STATUS_ACCEPTED)
                    if acc_mask.sum() > 0:
                        f1 = float(f1_score(gt[acc_mask], df.loc[acc_mask, "reference_cell_type"], average="macro", zero_division=0))
                    else:
                        f1 = 0.0
                    metrics[key] = round(f1, 6)

                metrics["delta_known_macro_f1"] = round(abs(metrics["macro_f1_known_a"] - metrics["macro_f1_known_b"]), 6)

    # Predeclared Parity Gates check
    gates_passed = bool(
        metrics["candidate_label_agreement"] >= 0.995
        and metrics["assignment_status_agreement"] >= 0.995
        and metrics["delta_rejection_rate"] <= 0.005
        and metrics["distance_allclose_1e4"]
    )
    if "delta_known_macro_f1" in metrics:
        gates_passed = gates_passed and (metrics["delta_known_macro_f1"] <= 0.01)
    if "delta_ood_rejection_recall" in metrics:
        gates_passed = gates_passed and (metrics["delta_ood_rejection_recall"] <= 0.02)

    metrics["parity_gates_passed"] = gates_passed
    metrics["parity_verdict"] = "PASS" if gates_passed else "FAIL"

    return metrics
