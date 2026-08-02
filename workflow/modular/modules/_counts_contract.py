"""Fail-closed raw-count semantics shared by ingestion and pseudobulk DE."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse


COUNTS_PROVENANCE_KEY = "counts_provenance"
COUNTS_SCHEMA_VERSION = "1.0"
COUNTS_MATRIX = "layers/counts"
COUNTS_SEMANTIC = "raw_umi_counts"


class CountsContractError(ValueError):
    """The declared raw-count matrix is missing provenance or is not count data."""


def _materialize(matrix: Any) -> Any:
    if hasattr(matrix, "compute"):
        return matrix.compute()
    if hasattr(matrix, "to_memory"):
        return matrix.to_memory()
    return matrix


def validate_nonnegative_integer_counts(
    matrix: Any,
    *,
    expected_shape: tuple[int, int] | None = None,
) -> dict[str, Any]:
    """Validate count semantics without rounding away evidence of normalization."""
    matrix = _materialize(matrix)
    shape = tuple(int(value) for value in matrix.shape)
    if expected_shape is not None and shape != tuple(expected_shape):
        raise CountsContractError(
            f"counts shape {shape} does not match AnnData shape {tuple(expected_shape)}"
        )

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    values = np.asarray(values)
    if values.size:
        if not np.isfinite(values).all():
            raise CountsContractError(
                "counts must be finite nonnegative integer raw UMI counts"
            )
        if bool((values < 0).any()):
            raise CountsContractError(
                "counts must be finite nonnegative integer raw UMI counts"
            )
        if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6):
            raise CountsContractError(
                "counts must be finite nonnegative integer raw UMI counts"
            )

    return {
        "shape": list(shape),
        "stored_values_checked": int(values.size),
        "finite": True,
        "nonnegative": True,
        "integer_like": True,
    }


def validate_counts_layer_contract(adata: Any) -> dict[str, Any]:
    """Require explicit semantics plus a numerically valid ``layers['counts']``."""
    if "counts" not in adata.layers:
        raise CountsContractError(
            "missing provenance-qualified adata.layers['counts'] raw UMI matrix"
        )

    provenance = adata.uns.get(COUNTS_PROVENANCE_KEY)
    if not isinstance(provenance, dict):
        raise CountsContractError(
            "confirmatory pseudobulk requires a provenance-qualified counts layer"
        )
    expected = {
        "schema_version": COUNTS_SCHEMA_VERSION,
        "matrix": COUNTS_MATRIX,
        "semantic": COUNTS_SEMANTIC,
    }
    mismatches = {
        field: {"expected": value, "observed": provenance.get(field)}
        for field, value in expected.items()
        if provenance.get(field) != value
    }
    source = provenance.get("source")
    if not isinstance(source, str) or not source.strip():
        mismatches["source"] = {"expected": "non-empty string", "observed": source}
    if mismatches:
        raise CountsContractError(
            "confirmatory pseudobulk requires a provenance-qualified counts layer; "
            f"invalid fields: {mismatches}"
        )

    validation = validate_nonnegative_integer_counts(
        adata.layers["counts"], expected_shape=adata.shape
    )
    return {
        "provenance": dict(provenance),
        "validation": validation,
    }


def stamp_factory_counts_provenance(adata: Any, *, source: str) -> None:
    """Stamp counts whose raw semantics are established by a factory-owned reader."""
    adata.uns[COUNTS_PROVENANCE_KEY] = {
        "schema_version": COUNTS_SCHEMA_VERSION,
        "matrix": COUNTS_MATRIX,
        "semantic": COUNTS_SEMANTIC,
        "source": source,
    }
