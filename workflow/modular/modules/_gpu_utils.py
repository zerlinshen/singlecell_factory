"""Shared GPU availability detection for all modules.

Result is cached after first probe so import/device checks run only once.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)
_gpu_ok: bool | None = None


def _probe_gpu() -> bool:
    global _gpu_ok
    if _gpu_ok is not None:
        return _gpu_ok
    try:
        import cupy
        cupy.cuda.runtime.getDeviceCount()
        import rapids_singlecell  # noqa: F401
        _gpu_ok = True
        logger.info("GPU backend available (rapids-singlecell + cupy)")
    except Exception as exc:
        logger.warning("GPU backend probe failed; CPU fallback will be used (%s)", exc)
        _gpu_ok = False
    return _gpu_ok


def gpu_available(mode: str = "auto") -> bool:
    """Return whether GPU should be used under the requested mode."""
    mode = (mode or "auto").lower()
    if mode == "off":
        return False
    ok = _probe_gpu()
    if mode == "force" and not ok:
        raise RuntimeError(
            "GPU mode 'force' requested, but the GPU backend is unavailable or incompatible."
        )
    return ok
