from __future__ import annotations

import threading
import traceback

from workflow.modular.pipeline import _get_available_memory_bytes
from workflow.modular._mem_guard import MemoryGuard, _get_rss_bytes, _get_total_memory_bytes, _HARD_CAP


def start(ctx, period: float = 5.0, hard: int | None = None) -> threading.Thread:
    """Start a background thread monitoring memory usage.

    When RSS exceeds hard_cap fraction of total RAM (or `hard` bytes if supplied),
    sets MemoryGuard._global_abort so that cooperative module loops can exit cleanly
    and dumps stack traces into ctx.metadata["oom_traces"].
    """
    stop_event = threading.Event()

    def _watch():
        import logging
        log = logging.getLogger(__name__)
        while not stop_event.wait(timeout=period):
            available = _get_available_memory_bytes()
            if available is None:
                continue

            # Hard-bytes threshold (legacy)
            hard_triggered = hard is not None and available < hard

            # RSS fraction threshold (new)
            rss_triggered = False
            rss = _get_rss_bytes()
            total = _get_total_memory_bytes()
            if rss is not None and total is not None and total > 0:
                rss_triggered = (rss / total) >= _HARD_CAP

            if hard_triggered or rss_triggered:
                frames = traceback.format_stack()
                traces = ctx.metadata.setdefault("oom_traces", [])
                traces.append({"available_bytes": available, "rss_bytes": rss, "stack": frames})
                log.warning(
                    "MemoryWatchdog: available=%d bytes rss_fraction=%.1f%% — requesting cooperative abort",
                    available,
                    (rss / total * 100) if (rss and total) else 0.0,
                )
                MemoryGuard.request_abort()

    t = threading.Thread(target=_watch, daemon=True, name="mem-watchdog")
    t._stop_event = stop_event  # type: ignore[attr-defined]
    t.start()
    return t


def stop(t: threading.Thread) -> None:
    """Stop a watchdog thread started by start()."""
    stop_event = getattr(t, "_stop_event", None)
    if stop_event is not None:
        stop_event.set()
    t.join(timeout=10)
