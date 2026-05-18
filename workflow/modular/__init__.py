"""Modular single-cell workflow package."""

from typing import TYPE_CHECKING, Any


def _install_series_nonzero_shim() -> None:
    """Shim pd.Series.nonzero for pandas 2.0+.

    pandas removed Series.nonzero in 2.0; upstream scanpy / rapids-singlecell
    still call it on large-cohort DE paths (observed on 800k cells × 33
    clusters during the Wave 3 NC2024 run). Restoring the method as a thin
    wrapper around np.flatnonzero preserves the upstream call shape exactly.
    Idempotent — only installs when missing.
    """
    try:
        import pandas as pd
        import numpy as np
    except Exception:
        return
    if not hasattr(pd.Series, "nonzero"):
        pd.Series.nonzero = lambda self: (np.flatnonzero(self.to_numpy()),)


_install_series_nonzero_shim()




if TYPE_CHECKING:
    from .config import PipelineConfig

__all__ = ["PipelineConfig", "run_pipeline", "main"]


def __getattr__(name: str) -> Any:
    if name == "PipelineConfig":
        from .config import PipelineConfig
        return PipelineConfig
    if name == "run_pipeline":
        from .pipeline import run_pipeline
        return run_pipeline
    if name == "main":
        from .cli import main
        return main
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
