"""Modular single-cell workflow package."""

from typing import TYPE_CHECKING, Any

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
