"""Workflow package for standard and velocity single-cell pipelines."""

from .standard import StandardWorkflowConfig, run_standard_workflow
from .velocity import VelocityWorkflowConfig, run_velocity_workflow

__all__ = [
    "StandardWorkflowConfig",
    "VelocityWorkflowConfig",
    "run_standard_workflow",
    "run_velocity_workflow",
]
