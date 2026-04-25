"""Workflow package for standard and velocity single-cell pipelines."""

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .modular import PipelineConfig
    from .standard import StandardWorkflowConfig
    from .velocity import VelocityWorkflowConfig

__all__ = [
    "PipelineConfig",
    "StandardWorkflowConfig",
    "VelocityWorkflowConfig",
    "run_pipeline",
    "run_standard_workflow",
    "run_velocity_workflow",
]


def __getattr__(name: str) -> Any:
    if name in {"PipelineConfig", "run_pipeline"}:
        from .modular import PipelineConfig, run_pipeline
        return {"PipelineConfig": PipelineConfig, "run_pipeline": run_pipeline}[name]
    if name in {"StandardWorkflowConfig", "run_standard_workflow"}:
        from .standard import StandardWorkflowConfig, run_standard_workflow
        return {
            "StandardWorkflowConfig": StandardWorkflowConfig,
            "run_standard_workflow": run_standard_workflow,
        }[name]
    if name in {"VelocityWorkflowConfig", "run_velocity_workflow"}:
        from .velocity import VelocityWorkflowConfig, run_velocity_workflow
        return {
            "VelocityWorkflowConfig": VelocityWorkflowConfig,
            "run_velocity_workflow": run_velocity_workflow,
        }[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
