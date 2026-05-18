"""Extended module_catalog tests — edge cases not in test_module_catalog.py."""
from __future__ import annotations

import pytest

from workflow.modular.module_catalog import (
    MODULE_SPECS,
    ModuleSpec,
    bridge_ready_modules,
    module_dependencies,
    modules_by_layer,
    optional_module_names,
)


def test_module_spec_is_frozen():
    spec = MODULE_SPECS["qc"]
    with pytest.raises((AttributeError, TypeError)):
        spec.name = "changed"  # type: ignore[misc]


def test_all_depends_on_refer_to_known_modules():
    known = set(MODULE_SPECS)
    for name, spec in MODULE_SPECS.items():
        for dep in spec.depends_on:
            assert dep in known, f"{name} depends on unknown module: {dep}"


def test_module_dependencies_returns_sets():
    deps = module_dependencies()
    for name, dep_set in deps.items():
        assert isinstance(dep_set, set), f"{name}: expected set, got {type(dep_set)}"


def test_optional_module_names_excludes_mandatory():
    from workflow.modular.module_catalog import MANDATORY_MODULES
    optional = set(optional_module_names())
    mandatory = set(MANDATORY_MODULES)
    assert optional.isdisjoint(mandatory)


def test_optional_module_names_covers_all_non_mandatory():
    from workflow.modular.module_catalog import MANDATORY_MODULES
    optional = set(optional_module_names())
    expected = set(MODULE_SPECS) - set(MANDATORY_MODULES)
    assert optional == expected


def test_modules_by_layer_contains_all_modules():
    layers = modules_by_layer()
    all_from_layers = {name for names in layers.values() for name in names}
    assert all_from_layers == set(MODULE_SPECS)


def test_bridge_ready_modules_is_subset_of_all():
    bridge = set(bridge_ready_modules())
    assert bridge.issubset(set(MODULE_SPECS))


def test_bridge_ready_modules_matches_spec_flag():
    bridge = set(bridge_ready_modules())
    for name, spec in MODULE_SPECS.items():
        if spec.bridge_ready:
            assert name in bridge, f"{name} has bridge_ready=True but not in bridge_ready_modules()"
        else:
            assert name not in bridge, f"{name} has bridge_ready=False but appears in bridge_ready_modules()"


def test_evolution_depends_on_cnv_and_trajectory():
    spec = MODULE_SPECS["evolution"]
    assert "cnv_inference" in spec.depends_on
    assert "trajectory" in spec.depends_on


def test_spatial_neighborhoods_depends_on_spatial_ingest():
    spec = MODULE_SPECS["spatial_neighborhoods"]
    assert "spatial_ingest" in spec.depends_on
