"""Global densify policy — single decision point for all .toarray() / .todense() in modules.

Every caller must call plan_densify() and act on the returned DensifyDecision before
materialising a sparse matrix. Direct .toarray() in workflow/modular/modules/ without a
# densify-allowed comment is banned by tests/test_densify_audit.py.
"""
from __future__ import annotations

import os
from enum import Enum

import numpy as np
import scipy.sparse as sp


class DensifyDecision(str, Enum):
    GO = "go"       # planned_bytes < soft cap; safe to toarray()
    CHUNK = "chunk" # soft_cap <= planned < hard_cap; caller must chunk
    ABORT = "abort" # planned >= hard_cap; raise MemoryGuardError


# Caps in bytes — overridable via env for testing / large-RAM machines
_DEFAULT_SOFT_CAP = 4 * 1024**3   # 4 GiB
_DEFAULT_HARD_CAP = 12 * 1024**3  # 12 GiB


def _soft_cap() -> int:
    return int(os.environ.get("SC_DENSIFY_SOFT_CAP_BYTES", _DEFAULT_SOFT_CAP))


def _hard_cap() -> int:
    return int(os.environ.get("SC_DENSIFY_HARD_CAP_BYTES", _DEFAULT_HARD_CAP))


def _free_bytes() -> int:
    """Best-effort free physical memory; falls back to a conservative constant."""
    try:
        import psutil
        return psutil.virtual_memory().available
    except Exception:
        pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except Exception:
        pass
    return 8 * 1024**3  # conservative fallback: 8 GiB


def planned_bytes(shape: tuple[int, int], dtype) -> int:
    """Estimate bytes required to materialise a dense array of given shape/dtype."""
    return int(np.prod(shape)) * np.dtype(dtype).itemsize


def plan_densify(
    shape: tuple[int, int],
    dtype,
    *,
    reason: str,
    hint_chunk_rows: int | None = None,
) -> DensifyDecision:
    """Return a DensifyDecision based on planned bytes vs memory caps.

    This function NEVER densifies; it only decides. Callers own the execution.

    Args:
        shape: (n_rows, n_cols) of the target dense array.
        dtype: numpy dtype of the target array.
        reason: human-readable explanation of why densification is needed.
        hint_chunk_rows: ignored; reserved for future adaptive chunking.

    Returns:
        DensifyDecision.GO / CHUNK / ABORT
    """
    nbytes = planned_bytes(shape, dtype)
    soft = _soft_cap()
    hard = _hard_cap()

    if nbytes >= hard:
        return DensifyDecision.ABORT
    if nbytes >= soft:
        return DensifyDecision.CHUNK
    # Also reject if planned would consume > 80 % of current free RAM
    free = _free_bytes()
    if nbytes > 0.80 * free:
        if nbytes >= hard:
            return DensifyDecision.ABORT
        return DensifyDecision.CHUNK
    return DensifyDecision.GO
