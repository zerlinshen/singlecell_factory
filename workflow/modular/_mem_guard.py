# Memory guard — dormant until Phase 3 imports this module.
from __future__ import annotations

from workflow.modular.pipeline import _get_available_memory_bytes, _estimate_adata_copy_bytes

_ABORT_FRACTION = 0.90
_CHUNK_FRACTION = 0.50


class MemoryGuard:
    """Decision gate for large memory allocations.

    Decisions: "go" | "chunk" | "backed" | "abort"
    """

    def __init__(self, ctx=None, label: str = ""):
        self._ctx = ctx
        self._label = label

    @staticmethod
    def _quick_check(planned_bytes: int) -> str:
        available = _get_available_memory_bytes()
        if available is None:
            return "go"
        if planned_bytes >= available * _ABORT_FRACTION:
            return "abort"
        if planned_bytes >= available * _CHUNK_FRACTION:
            return "chunk"
        return "go"

    def check(self, label: str = "") -> str:
        if self._ctx is None or self._ctx.adata is None:
            return "go"
        est = _estimate_adata_copy_bytes(self._ctx.adata)
        return self._quick_check(est)

    def estimate_and_decide(self, planned_bytes: int) -> str:
        available = _get_available_memory_bytes()
        if available is None:
            return "backed"
        if planned_bytes >= available * _ABORT_FRACTION:
            return "abort"
        if planned_bytes >= available * _CHUNK_FRACTION:
            return "chunk"
        return "go"
