"""Shared GPU availability detection for all modules.

Result is cached after first probe so import/device checks run only once.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)
_gpu_ok: bool | None = None
_cuda_warmed: bool = False


def bind_cuda_context(device: int = 0) -> bool:
    """Pre-warm cupy's cuBLAS/cuSOLVER handles so a later rapids import can't poison them.

    Root cause (probe at projects/gpu-probe-20260527/probe_results.md): importing
    ``rapids_singlecell`` leaves the CUDA device context in a state where cuBLAS/
    cuSOLVER handles *created afterwards* fail on use — ``CUBLAS_STATUS_NOT_INITIALIZED``
    (gemm) / ``CUSOLVER_STATUS_INTERNAL_ERROR`` (eigh/PCA) — even though the handle
    pointer is non-NULL. It is an init-ordering bug, not a Blackwell/CUDA-13 ABI gap.

    EMPIRICALLY-PROVEN FIX (p1_smoke recovery_strategies{2,3}.py + p1_smoke_3lane.py,
    2026-05-27): a post-import ``cupy.cuda.Device(0).use()`` does NOT recover a
    poisoned context in a clean process (the probe's "recovery" was an interactive-
    session artifact where a handle had already been warmed). What DOES work is the
    reverse: create+bind the cuBLAS/cuSOLVER handles BEFORE rapids_singlecell is
    imported. Once warmed, the handles survive the rapids import and all subsequent
    GPU ops pass.

    BOTH backends must be warmed: the cupy handle (rsc PCA/neighbors/scrublet use
    it via cuBLAS/cuSOLVER) AND torch's cuBLAS/cublasLt handle (scVI's encoder
    matmuls use it) — the 3-lane smoke showed scVI failing with
    CUBLAS_STATUS_NOT_INITIALIZED from torch's cublasLt when only cupy was warmed.
    torch warming is guarded so CPU-only / torch-absent installs are unaffected.

    Called inside ``_probe_gpu`` between the cupy import and the rapids import, and
    ``gpu_available`` runs at the top of every GPU module's lane before its lazy
    ``import rapids_singlecell``. Idempotent (warms once, cached); safe to call from
    GPU lanes after a rapids import too (no-op once warmed). Silently no-ops when
    cupy is unavailable or no CUDA device is present, so CPU-only paths are
    unaffected. Returns True if the context is warmed.
    """
    global _cuda_warmed
    if _cuda_warmed:
        return True
    try:
        import cupy
        import numpy as _np
        if cupy.cuda.runtime.getDeviceCount() <= device:
            return False
        cupy.cuda.Device(device).use()
        _w = cupy.asarray(_np.empty((8, 8), dtype=_np.float32))
        _w[...] = 1.0
        _ = _w @ _w                          # creates+binds the cuBLAS handle
        _ = cupy.linalg.eigh(_w @ _w.T)[0]   # creates+binds the cuSOLVER handle
        cupy.cuda.Stream.null.synchronize()
    except Exception as exc:  # pragma: no cover - depends on local GPU/runtime
        logger.debug("bind_cuda_context (cupy) skipped (%s)", exc)
        return False

    # Warm torch's cuBLAS/cublasLt too, for the scVI (torch) lane. Guarded: torch
    # may be absent and that must not break the cupy-only GPU lanes.
    try:
        import torch
        if torch.cuda.is_available():
            _t = torch.ones((32, 32), device=f"cuda:{device}")
            _ = _t @ _t                       # torch cuBLAS
            _ = torch.linalg.matmul(_t, _t)   # torch cublasLt heuristic path
            torch.cuda.synchronize()
    except Exception as exc:  # pragma: no cover - torch optional / GPU-dependent
        logger.debug("bind_cuda_context (torch) skipped (%s)", exc)

    _cuda_warmed = True
    logger.info("CUDA cuBLAS/cuSOLVER context pre-warmed (cupy+torch, pre-rapids, P1 fix)")
    return True


def _probe_gpu() -> bool:
    global _gpu_ok
    if _gpu_ok is not None:
        return _gpu_ok
    try:
        import cupy
        cupy.cuda.runtime.getDeviceCount()
        # P1 fix: warm cuBLAS/cuSOLVER on a clean context BEFORE the rapids import
        # below, otherwise the rapids import poisons handles created afterwards.
        bind_cuda_context()
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
