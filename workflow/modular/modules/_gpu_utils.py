"""Shared GPU availability detection for all modules.

Result is cached after first check so the import/device probe runs only once.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)
_gpu_ok: bool | None = None


def gpu_available() -> bool:
    """Return True if rapids-singlecell + cupy + a CUDA device are all usable."""
    global _gpu_ok
    if _gpu_ok is not None:
        return _gpu_ok
    try:
        import cupy  # noqa: F401
        cupy.cuda.runtime.getDeviceCount()
        import rapids_singlecell  # noqa: F401
        _gpu_ok = True
        logger.info("GPU backend available (rapids-singlecell + cupy)")
    except Exception:
        _gpu_ok = False
    return _gpu_ok
