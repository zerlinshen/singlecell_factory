# Memory watchdog — dormant until Phase 3 imports this module.
from __future__ import annotations

import threading
import traceback

from workflow.modular.pipeline import _get_available_memory_bytes


def start(ctx, period: float = 5.0, hard: int | None = None) -> threading.Thread:
    """Start a background thread monitoring memory usage.

    When available memory drops below `hard` bytes, dumps stack traces
    into ctx.metadata["oom_traces"] and logs a warning.
    """
    stop_event = threading.Event()

    def _watch():
        import logging
        log = logging.getLogger(__name__)
        while not stop_event.wait(timeout=period):
            available = _get_available_memory_bytes()
            if available is None:
                continue
            if hard is not None and available < hard:
                frames = traceback.format_stack()
                traces = ctx.metadata.setdefault("oom_traces", [])
                traces.append({"available_bytes": available, "stack": frames})
                log.warning("MemoryWatchdog: available=%d bytes < hard=%d", available, hard)

    t = threading.Thread(target=_watch, daemon=True, name="mem-watchdog")
    t._stop_event = stop_event  # type: ignore[attr-defined]
    t.start()
    return t
