from __future__ import annotations

import os
import threading

from workflow.modular.pipeline import _get_available_memory_bytes, _estimate_adata_copy_bytes

# Legacy fractions kept for existing estimate_and_decide / _quick_check callers.
_ABORT_FRACTION = 0.90
_CHUNK_FRACTION = 0.50

# Env-overridable thresholds for check_budget (fraction of total RAM).
# hard_cap: watchdog triggers abort above this RSS fraction.
# soft_cap: check_budget returns "chunk" above this fraction.
_HARD_CAP = float(os.environ.get("SC_MEM_HARD_CAP", "0.85"))
_SOFT_CAP = float(os.environ.get("SC_MEM_SOFT_CAP", "0.65"))


def _get_total_memory_bytes() -> int | None:
    """Return total physical RAM in bytes."""
    try:
        import psutil  # type: ignore
        return int(psutil.virtual_memory().total)
    except Exception:
        pass
    # /proc/self/status fallback (Linux, no psutil)
    try:
        with open("/proc/meminfo", "r") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return kb * 1024
    except Exception:
        pass
    return None


def _get_rss_bytes() -> int | None:
    """Return resident set size of current process in bytes."""
    try:
        import psutil  # type: ignore
        return int(psutil.Process().memory_info().rss)
    except Exception:
        pass
    # /proc/self/status fallback
    try:
        with open("/proc/self/status", "r") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    kb = int(line.split()[1])
                    return kb * 1024
    except Exception:
        pass
    return None


class MemoryGuardError(Exception):
    """Raised by MemoryGuard.check() when memory pressure exceeds thresholds."""


class MemoryAbortError(Exception):
    """Raised cooperatively by module loops when watchdog sets abort signal."""


class MemoryGuard:
    """Decision gate for large memory allocations.

    Decisions: "go" | "chunk" | "backed" | "abort"

    check_budget(planned_bytes, label) — pre-flight budget check.
      Returns "go" | "chunk" | "abort". Never raises.

    check() — legacy entry-point used by existing modules.
      Raises MemoryGuardError on abort.

    Used as a context manager: records entry/exit in ctx.metadata['mem_traces'].
    """

    # Module-level abort event shared between watchdog and cooperative loop checks.
    _global_abort: threading.Event = threading.Event()

    def __init__(self, ctx=None, label: str = ""):
        self._ctx = ctx
        self._label = label

    def __enter__(self) -> "MemoryGuard":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        return False

    # ------------------------------------------------------------------
    # Pre-flight budget check (new API — returns, never raises)
    # ------------------------------------------------------------------

    @staticmethod
    def check_budget(planned_bytes: int, label: str = "") -> str:
        """Pre-flight budget check.

        Compares planned_bytes against current available memory using
        hard_cap / soft_cap thresholds (env-configurable).

        Returns: "go" | "chunk" | "abort"
        Does NOT raise.
        """
        available = _get_available_memory_bytes()
        if available is None:
            return "go"
        total = _get_total_memory_bytes()
        if total is None:
            # Fall back to available-only heuristic
            if planned_bytes >= available * _ABORT_FRACTION:
                return "abort"
            if planned_bytes >= available * _CHUNK_FRACTION:
                return "chunk"
            return "go"

        # Compute projected RSS fraction after allocation
        rss = _get_rss_bytes() or 0
        projected = rss + planned_bytes
        projected_fraction = projected / total

        if projected_fraction >= _HARD_CAP:
            return "abort"
        if projected_fraction >= _SOFT_CAP:
            return "chunk"
        return "go"

    # ------------------------------------------------------------------
    # Watchdog abort signal helpers
    # ------------------------------------------------------------------

    @classmethod
    def clear_abort(cls) -> None:
        cls._global_abort.clear()

    @classmethod
    def request_abort(cls) -> None:
        cls._global_abort.set()

    @classmethod
    def abort_requested(cls) -> bool:
        return cls._global_abort.is_set()

    # ------------------------------------------------------------------
    # Legacy API
    # ------------------------------------------------------------------

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
