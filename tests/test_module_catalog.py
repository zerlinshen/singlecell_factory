from __future__ import annotations

from workflow.modular.module_catalog import (
    DEFAULT_OPTIONAL_MODULES,
    MANDATORY_MODULES,
    MODULE_SPECS,
    bridge_ready_modules,
    module_dependencies,
    module_help_list,
    modules_by_layer,
    optional_module_names,
)
from workflow.modular.pipeline import MODULE_DEPENDENCIES


def test_module_catalog_is_pipeline_dependency_source() -> None:
    assert module_dependencies() == MODULE_DEPENDENCIES
    assert set(MODULE_SPECS) == set(MODULE_DEPENDENCIES)


def test_mandatory_chain_and_default_optional_modules_are_stable() -> None:
    assert MANDATORY_MODULES == ("cellranger", "qc", "doublet_detection")
    assert DEFAULT_OPTIONAL_MODULES == (
        "clustering",
        "differential_expression",
        "annotation",
        "trajectory",
        "pseudo_velocity",
    )
    assert set(MANDATORY_MODULES).isdisjoint(optional_module_names())


def test_cli_help_uses_optional_modules_only() -> None:
    help_text = module_help_list()
    assert "cellranger" not in help_text
    assert "qc" not in help_text
    assert "paper_repro" in help_text
    assert "pseudobulk_de" in help_text


def test_layers_cover_bridge_ready_outputs() -> None:
    layers = modules_by_layer()
    assert "ingest" in layers
    assert "quality_control" in layers
    assert "annotation" in layers
    assert {"annotation", "composition", "paper_repro"}.issubset(bridge_ready_modules())
