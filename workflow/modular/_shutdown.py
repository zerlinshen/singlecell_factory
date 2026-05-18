"""Explicit interpreter-shutdown cleanup for C extensions.

Called at the end of cli.main() to drain C-extension finalizers before the
interpreter begins its own GC teardown.  Segfaults from cupy / torch / zarr
occur because their C destructors run in undefined order during interpreter
shutdown when objects are still reachable from sys.modules.  Explicit cleanup
here lets them teardown cleanly while the runtime is still functional.
"""
from __future__ import annotations

import gc
import logging
import sys

logger = logging.getLogger(__name__)


def run_shutdown_cleanup(adata=None) -> None:
    """Release heavy C-extension objects before interpreter shutdown."""

    # 1. Drop the caller's adata reference so the next gc.collect() can free it.
    if adata is not None:
        try:
            del adata
        except Exception:
            pass

    # 2. First GC pass — collect cycles referencing numpy/sparse arrays.
    gc.collect()

    # 3. Torch: free CPU tensors and any cached CUDA state.
    #    Some optional backends (scVI, etc.) import torch.
    if "torch" in sys.modules:
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
                torch.cuda.synchronize()
        except Exception as exc:
            logger.debug("torch shutdown cleanup skipped: %s", exc)

    # 4. CuPy: synchronize device and release memory pool before CUDA context
    #    is torn down.  Without this the CUDA driver destructor races with the
    #    Python GC and segfaults.
    if "cupy" in sys.modules:
        try:
            import cupy
            cupy.get_default_memory_pool().free_all_blocks()
            cupy.get_default_pinned_memory_pool().free_all_blocks()
        except Exception as exc:
            logger.debug("cupy shutdown cleanup skipped: %s", exc)

    # 5. Zarr: consolidated_metadata and remote stores keep file handles open.
    #    Close any zarr.Group / zarr.Array objects still alive.
    if "zarr" in sys.modules:
        try:
            import zarr
            for obj in gc.get_objects():
                if isinstance(obj, (zarr.Group, zarr.Array)):
                    try:
                        store = getattr(obj, "store", None)
                        if store is not None and hasattr(store, "close"):
                            store.close()
                    except Exception:
                        pass
        except Exception as exc:
            logger.debug("zarr shutdown cleanup skipped: %s", exc)

    # 6. Neighbors cache: clear in-process LRU cache before GC teardown.
    try:
        from . import _neighbors_cache
        _neighbors_cache.clear()
    except Exception as exc:
        logger.debug("neighbors_cache clear skipped: %s", exc)

    # 7. Second GC pass — collect anything freed by steps 3-6.
    gc.collect()
    gc.collect()
