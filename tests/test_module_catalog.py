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
    # 2026-05-20: ambient_correction added between qc and doublet_detection
    # (Phase 2 of ambient_correction_policy.md). It is MANDATORY because the
    # conditional skip is data-driven inside the module, not a DAG-level
    # opt-out -- doublet_detection therefore depends on ambient_correction
    # being present in the DAG.
    assert MANDATORY_MODULES == (
        "cellranger",
        "qc",
        "ambient_correction",
        "doublet_detection",
    )
    assert DEFAULT_OPTIONAL_MODULES == (
        "clustering",
        "differential_expression",
        "annotation",
    )
    assert "trajectory" in MODULE_SPECS
    assert "pseudo_velocity" in MODULE_SPECS
    assert set(MANDATORY_MODULES).isdisjoint(optional_module_names())


def test_cli_help_uses_optional_modules_only() -> None:
    # Check at module-name token level, not substring, because the optional
    # module set legitimately includes names that *contain* mandatory tokens
    # (e.g. ``atac_qc``, ``cross_modality_qc`` contain ``qc``).
    help_tokens = [name.strip() for name in module_help_list().split(",")]
    assert "cellranger" not in help_tokens
    assert "qc" not in help_tokens
    assert "ambient_correction" not in help_tokens
    assert "doublet_detection" not in help_tokens
    assert "paper_repro" in help_tokens
    assert "pseudobulk_de" in help_tokens


def test_layers_cover_bridge_ready_outputs() -> None:
    layers = modules_by_layer()
    assert "ingest" in layers
    assert "quality_control" in layers
    assert "annotation" in layers
    assert {"annotation", "composition", "paper_repro"}.issubset(bridge_ready_modules())
