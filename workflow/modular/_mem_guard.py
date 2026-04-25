# Memory guard — dormant until Phase 3 imports this module.
from __future__ import annotations

from workflow.modular.pipeline import _get_available_memory_bytes, _estimate_adata_copy_bytes

_ABORT_FRACTION = 0.90
_CHUNK_FRACTION = 0.50


class MemoryGuardError(Exception):
    """Raised by MemoryGuard.check() when memory pressure exceeds thresholds."""


class MemoryGuard:
    """Decision gate for large memory allocations.

    Decisions: "go" | "chunk" | "backed" | "abort"

    Used as a context manager: records entry/exit in ctx.metadata['mem_traces'].
    check() raises MemoryGuardError when decision is "abort".
    """

    def __init__(self, ctx=None, label: str = ""):
        self._ctx = ctx
        self._label = label

    def __enter__(self) -> "MemoryGuard":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        return False

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
            decision = "go"
        else:
            est = _estimate_adata_copy_bytes(self._ctx.adata)
            decision = self._quick_check(est)

        if self._ctx is not None:
            trace = {"module": self._label, "label": label, "decision": decision}
            self._ctx.metadata.setdefault("mem_traces", []).append(trace)

        if decision == "abort":
            raise MemoryGuardError(
                f"MemoryGuard abort: module={self._label!r} label={label!r} decision={decision!r}"
            )
        return decision

    def estimate_and_decide(self, planned_bytes: int) -> str:
        available = _get_available_memory_bytes()
        if available is None:
            return "backed"
        if planned_bytes >= available * _ABORT_FRACTION:
            return "abort"
        if planned_bytes >= available * _CHUNK_FRACTION:
            return "chunk"
        return "go"
