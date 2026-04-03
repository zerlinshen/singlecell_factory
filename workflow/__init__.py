"""Workflow package for standard and velocity single-cell pipelines."""

from .modular import PipelineConfig, run_pipeline
from .standard import StandardWorkflowConfig, run_standard_workflow
from .velocity import VelocityWorkflowConfig, run_velocity_workflow

__all__ = [
    "PipelineConfig",
    "StandardWorkflowConfig",
    "VelocityWorkflowConfig",
    "run_pipeline",
    "run_standard_workflow",
    "run_velocity_workflow",
]
