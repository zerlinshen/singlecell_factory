"""US-016 / D.1 / P-CC.S27 — full registry topo-sort guard.

Asserts that the FULL module registry (mandatories + every optional module
declared in module_catalog.MODULE_SPECS) topo-sorts without cycles. This is
the safety net against accidental cyclic depends_on additions when Wave 2B
vertical slices land (vdj/hic/ribo modules).
"""
from __future__ import annotations

import pytest


def test_full_registry_resolves_no_cycle():
    """Kahn's topo (workflow/modular/pipeline.py:130-163) succeeds on the full DAG."""
    from workflow.modular.pipeline import _build_registry, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES

    registry = _build_registry()
    optional = [name for name in registry if name not in MANDATORY_MODULES]

    # Should NOT raise ValueError("Cyclic dependency detected: ...")
    order = _resolve_execution_order(list(MANDATORY_MODULES), optional)

    # All requested modules end up in order
    requested = set(MANDATORY_MODULES) | set(optional)
    assert set(order) >= requested, f"missing from topo order: {requested - set(order)}"


def test_mandatory_modules_respect_their_own_deps():
    """Each mandatory module precedes every module that lists it as a dependency.

    Note: this is the real topological invariant. Mandatory modules do NOT have to
    precede every optional module — only those that depend on them. Several
    optionals (marker_db_loader, modality_registry) have empty
    depends_on and may topologically appear at the same level as cellranger.
    """
    from workflow.modular.pipeline import _build_registry, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES, MODULE_SPECS

    registry = _build_registry()
    optional = [n for n in registry if n not in MANDATORY_MODULES]
    order = _resolve_execution_order(list(MANDATORY_MODULES), optional)

    for mandatory in MANDATORY_MODULES:
        if mandatory not in order:
            continue
        m_idx = order.index(mandatory)
        # Every module that depends on `mandatory` must appear AFTER it.
        for name, spec in MODULE_SPECS.items():
            if mandatory in spec.depends_on and name in order:
                dependent_idx = order.index(name)
                assert dependent_idx > m_idx, (
                    f"{name} (idx={dependent_idx}) must come after its dep "
                    f"{mandatory} (idx={m_idx})"
                )


def test_registry_count_meets_baseline():
    """Registry should contain at minimum the Wave 1 + Wave 2A modules."""
    from workflow.modular.pipeline import _build_registry

    registry = _build_registry()
    # Wave 1 baseline: 31 modules; Wave 2A added 3 (modality_registry, cross_modality_qc, atac_ingest)
    # Wave 2B target (US-015): ≥58 after adding vdj/hic/ribo + atac_qc + peak_to_gene.
    # For now assert the Wave 2A floor.
    assert len(registry) >= 34, f"registry size shrank below Wave 2A floor of 34 (got {len(registry)})"


def test_no_self_loops_in_module_dependencies():
    """No module declares itself as a dependency (would create a trivial cycle)."""
    from workflow.modular.module_catalog import MODULE_SPECS

    for name, spec in MODULE_SPECS.items():
        assert name not in spec.depends_on, f"{name} has self-loop in depends_on"


def test_dependency_targets_exist_in_registry():
    """Every declared depends_on target must actually exist as a module."""
    from workflow.modular.pipeline import _build_registry
    from workflow.modular.module_catalog import MODULE_SPECS

    registry = _build_registry()
    for name, spec in MODULE_SPECS.items():
        for dep in spec.depends_on:
            assert dep in registry, f"{name}.depends_on references unknown module: {dep}"
