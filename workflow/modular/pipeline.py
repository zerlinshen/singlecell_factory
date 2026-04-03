from __future__ import annotations

import copy
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import json

import pandas as pd

from .config import PipelineConfig
from .context import PipelineContext

logger = logging.getLogger(__name__)

# Module dependency DAG: module_name -> set of modules it requires
MODULE_DEPENDENCIES: dict[str, set[str]] = {
    "cellranger": set(),
    "qc": {"cellranger"},
    "doublet_detection": {"qc"},
    "clustering": {"doublet_detection"},
    "cell_cycle": {"clustering"},
    "batch_correction": {"clustering"},
    "differential_expression": {"clustering"},
    "annotation": {"clustering"},
    "trajectory": {"clustering"},
    "pseudo_velocity": {"trajectory"},
    "rna_velocity": {"clustering"},
    "cnv_inference": {"clustering"},
    "pathway_analysis": {"differential_expression"},
    "cell_communication": {"annotation"},
    "gene_regulatory_network": {"clustering"},
    "compare_10x": {"clustering"},
    "validate_cbioportal": {"differential_expression"},
    "immune_phenotyping": {"annotation"},
    "tumor_microenvironment": {"annotation"},
    "gene_signature_scoring": {"clustering"},
    "evolution": {"cnv_inference", "trajectory"},
}

# Modules that mutate adata structurally (embeddings, X, layers).
# These MUST run sequentially, not in parallel branches.
MUTATING_MODULES = {"batch_correction", "rna_velocity"}


def _resolve_execution_order(mandatory: list[str], optional: list[str]) -> list[str]:
    """Topologically sort modules respecting dependencies.

    Mandatory modules always run first. For optional modules, any missing
    dependencies that are themselves optional are auto-included.
    """
    all_requested = set(mandatory) | set(optional)

    # Auto-include transitive dependencies
    to_process = list(all_requested)
    while to_process:
        mod = to_process.pop()
        for dep in MODULE_DEPENDENCIES.get(mod, set()):
            if dep not in all_requested:
                all_requested.add(dep)
                to_process.append(dep)

    # Topological sort (Kahn's algorithm)
    in_degree: dict[str, int] = {m: 0 for m in all_requested}
    for mod in all_requested:
        for dep in MODULE_DEPENDENCIES.get(mod, set()):
            if dep in all_requested:
                in_degree[mod] = in_degree.get(mod, 0) + 1

    queue = sorted([m for m, d in in_degree.items() if d == 0])
    order = []
    while queue:
        node = queue.pop(0)
        order.append(node)
        for mod in sorted(all_requested):
            if node in MODULE_DEPENDENCIES.get(mod, set()):
                in_degree[mod] -= 1
                if in_degree[mod] == 0:
                    queue.append(mod)

    # Detect cycles: if not all requested modules were emitted, a dependency
    # cycle exists in the induced subgraph.
    remaining = sorted(m for m in all_requested if m not in order)
    if remaining:
        raise ValueError(f"循环依赖检测: {remaining}")

    return order


def _build_registry() -> dict[str, object]:
    """Lazy-import all modules and build the registry.

    Importing inside this function avoids pulling heavy dependencies
    (scrublet, scvelo, liana, etc.) at package-level import time.
    """
    from .modules.annotation import AnnotationModule
    from .modules.batch_correction import BatchCorrectionModule
    from .modules.evolution import EvolutionModule
    from .modules.cell_communication import CellCommunicationModule
    from .modules.cell_cycle import CellCycleModule
    from .modules.cellranger import CellRangerModule
    from .modules.clustering import ClusteringModule
    from .modules.cnv_inference import CNVInferenceModule
    from .modules.compare_10x import Compare10XModule
    from .modules.differential_expression import DifferentialExpressionModule
    from .modules.doublet_detection import DoubletDetectionModule
    from .modules.gene_regulatory_network import GeneRegulatoryNetworkModule
    from .modules.gene_signature_scoring import GeneSignatureScoringModule
    from .modules.immune_phenotyping import ImmunePhenotypingModule
    from .modules.pathway_analysis import PathwayAnalysisModule
    from .modules.pseudo_velocity import PseudoVelocityModule
    from .modules.qc import QCModule
    from .modules.rna_velocity import RNAVelocityModule
    from .modules.trajectory import TrajectoryModule
    from .modules.tumor_microenvironment import TumorMicroenvironmentModule
    from .modules.validate_cbioportal import ValidateCbioPortalModule

    return {
        "cellranger": CellRangerModule(),
        "qc": QCModule(),
        "doublet_detection": DoubletDetectionModule(),
        "clustering": ClusteringModule(),
        "cell_cycle": CellCycleModule(),
        "batch_correction": BatchCorrectionModule(),
        "differential_expression": DifferentialExpressionModule(),
        "annotation": AnnotationModule(),
        "trajectory": TrajectoryModule(),
        "pseudo_velocity": PseudoVelocityModule(),
        "rna_velocity": RNAVelocityModule(),
        "cnv_inference": CNVInferenceModule(),
        "pathway_analysis": PathwayAnalysisModule(),
        "cell_communication": CellCommunicationModule(),
        "gene_regulatory_network": GeneRegulatoryNetworkModule(),
        "compare_10x": Compare10XModule(),
        "validate_cbioportal": ValidateCbioPortalModule(),
        "immune_phenotyping": ImmunePhenotypingModule(),
        "tumor_microenvironment": TumorMicroenvironmentModule(),
        "gene_signature_scoring": GeneSignatureScoringModule(),
        "evolution": EvolutionModule(),
    }


def _prepare_output(cfg: PipelineConfig) -> PipelineContext:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = cfg.output_dir / f"{cfg.project}_{ts}"
    run_dir.mkdir(parents=True, exist_ok=True)
    # figure_dir and table_dir will be set per-module via ctx.set_module_dir()
    return PipelineContext(cfg=cfg, run_dir=run_dir, figure_dir=run_dir, table_dir=run_dir)


def _save_manifest(ctx: PipelineContext) -> Path:
    ctx.flush_figures()
    if ctx.adata is not None:
        ctx.adata.write(ctx.run_dir / "final_adata.h5ad")
    manifest = {
        "project": ctx.cfg.project,
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "run_dir": str(ctx.run_dir),
        "optional_modules": ctx.cfg.optional_modules,
        "module_status": ctx.module_status,
        "metadata": ctx.metadata,
    }
    manifest_path = ctx.run_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    pd.DataFrame(ctx.module_status).to_csv(ctx.run_dir / "module_status.csv", index=False)
    return manifest_path


# ---------------------------------------------------------------------------
# Parallel execution helpers
# ---------------------------------------------------------------------------


def _compute_tiers(execution_order: list[str], completed: set[str]) -> list[list[str]]:
    """Group remaining modules into tiers for parallel execution.

    A tier is a set of modules whose dependencies are all satisfied by
    previously completed tiers plus the *completed* set.
    """
    remaining = [m for m in execution_order if m not in completed]
    done = set(completed)
    tiers: list[list[str]] = []
    while remaining:
        tier = [
            m for m in remaining
            if MODULE_DEPENDENCIES.get(m, set()).issubset(done)
        ]
        if not tier:
            # Shouldn't happen with a valid topo sort, but be safe.
            tier = [remaining[0]]
        tiers.append(tier)
        done.update(tier)
        remaining = [m for m in remaining if m not in done]
    return tiers


def _run_module(mod, ctx: PipelineContext) -> None:
    """Run a single module (used as ThreadPoolExecutor target)."""
    mod.run(ctx)


def _execute_tier(
    tier: list[str],
    registry: dict,
    ctx: PipelineContext,
    mandatory: set[str],
    max_workers: int,
) -> None:
    """Execute a tier of modules, potentially in parallel."""
    # Split into mutating (must run sequentially) and appending (can parallelize).
    mutating = [m for m in tier if m in MUTATING_MODULES]
    appending = [m for m in tier if m not in MUTATING_MODULES]

    # Run mutating modules sequentially first.
    for stage in mutating:
        mod = registry.get(stage)
        if mod is None:
            ctx.status(stage, False, "unknown module")
            continue
        ctx.set_module_dir(stage)
        try:
            mod.run(ctx)
            ctx.status(stage, True, "completed")
            ctx.save_checkpoint(stage)
        except Exception as exc:
            ctx.status(stage, False, str(exc))
            if stage in mandatory:
                raise

    # Run appending modules.
    if max_workers <= 1 or len(appending) <= 1:
        # Sequential execution.
        for stage in appending:
            mod = registry.get(stage)
            if mod is None:
                ctx.status(stage, False, "unknown module")
                continue
            ctx.set_module_dir(stage)
            try:
                mod.run(ctx)
                ctx.status(stage, True, "completed")
                ctx.save_checkpoint(stage)
            except Exception as exc:
                ctx.status(stage, False, str(exc))
                if stage in mandatory:
                    raise
    else:
        _run_parallel_appending(appending, registry, ctx, mandatory, max_workers)


def _run_parallel_appending(
    modules: list[str],
    registry: dict,
    ctx: PipelineContext,
    mandatory: set[str],
    max_workers: int,
) -> None:
    """Run appending-only modules in parallel with copy-on-branch, merge-back."""
    # Create per-module branch contexts with independent adata copies.
    branch_contexts: dict[str, PipelineContext] = {}
    for mod_name in modules:
        mod = registry.get(mod_name)
        if mod is None:
            ctx.status(mod_name, False, "unknown module")
            continue
        branch_ctx = copy.copy(ctx)
        branch_ctx.adata = ctx.adata.copy()
        branch_ctx.module_status = []
        branch_ctx.metadata = dict(ctx.metadata)
        branch_ctx._figure_futures = []
        branch_ctx.set_module_dir(mod_name)
        branch_contexts[mod_name] = branch_ctx

    results: dict[str, tuple[str, str]] = {}
    with ThreadPoolExecutor(max_workers=min(max_workers, len(branch_contexts))) as pool:
        futures = {
            pool.submit(_run_module, registry[name], branch_ctx): name
            for name, branch_ctx in branch_contexts.items()
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
                results[name] = ("ok", "completed")
            except Exception as exc:
                results[name] = ("failed", str(exc))

    # Merge results back into the main context.
    main_obs_cols = set(ctx.adata.obs.columns)
    main_obsm_keys = set(ctx.adata.obsm.keys())

    for mod_name in modules:
        if mod_name not in branch_contexts:
            continue
        branch_ctx = branch_contexts[mod_name]
        ok, msg = results.get(mod_name, ("failed", "not run"))

        if ok == "ok":
            # Merge new obs columns.
            new_cols = set(branch_ctx.adata.obs.columns) - main_obs_cols
            for col in new_cols:
                ctx.adata.obs[col] = branch_ctx.adata.obs[col].values
            # Merge new obsm entries.
            new_obsm = set(branch_ctx.adata.obsm.keys()) - main_obsm_keys
            for key in new_obsm:
                ctx.adata.obsm[key] = branch_ctx.adata.obsm[key]
            # Merge new uns entries.
            for key in branch_ctx.adata.uns:
                if key not in ctx.adata.uns:
                    ctx.adata.uns[key] = branch_ctx.adata.uns[key]
            # Merge metadata.
            ctx.metadata.update(branch_ctx.metadata)
            # Flush any figures the branch produced.
            branch_ctx.flush_figures()

        ctx.status(mod_name, ok == "ok", msg)
        ctx.save_checkpoint(mod_name)

        if ok != "ok" and mod_name in mandatory:
            raise RuntimeError(f"Mandatory module {mod_name} failed: {msg}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run_pipeline(cfg: PipelineConfig) -> Path:
    """Run the modular workflow with mandatory and optional stages."""

    ctx = _prepare_output(cfg)
    registry = _build_registry()
    try:
        # Enable async figure pool if parallel workers > 1.
        if cfg.parallel_workers > 1:
            ctx._figure_pool = ThreadPoolExecutor(max_workers=2)

        mandatory = ["cellranger", "qc", "doublet_detection"]
        mandatory_set = set(mandatory)
        execution_order = _resolve_execution_order(mandatory, cfg.optional_modules)

        # --- Resume from checkpoint ---
        if cfg.resume_from:
            # Find the module just before resume_from in execution order.
            try:
                resume_idx = execution_order.index(cfg.resume_from)
            except ValueError:
                raise ValueError(
                    f"Cannot resume from '{cfg.resume_from}': not in execution order."
                )
            if resume_idx > 0:
                prev_module = execution_order[resume_idx - 1]
                if not ctx.load_checkpoint(prev_module):
                    raise FileNotFoundError(
                        f"No checkpoint found after '{prev_module}'. "
                        f"Run the pipeline with --checkpoint first."
                    )
                logger.info("Resumed from checkpoint after '%s'", prev_module)
            execution_order = execution_order[resume_idx:]

        # --- Execute ---
        if cfg.parallel_workers > 1:
            all_modules = _resolve_execution_order(mandatory, cfg.optional_modules)
            completed = {m for m in all_modules if m not in execution_order}
            tiers = _compute_tiers(execution_order, completed)
            for tier in tiers:
                _execute_tier(tier, registry, ctx, mandatory_set, cfg.parallel_workers)
        else:
            # Original sequential execution.
            for stage in execution_order:
                mod = registry.get(stage)
                if mod is None:
                    ctx.status(stage, False, "unknown module")
                    continue
                ctx.set_module_dir(stage)
                is_mandatory = stage in mandatory_set
                try:
                    mod.run(ctx)
                    ctx.status(stage, True, "completed")
                    ctx.save_checkpoint(stage)
                except Exception as exc:
                    ctx.status(stage, False, str(exc))
                    if is_mandatory:
                        raise

        return _save_manifest(ctx)
    finally:
        if ctx._figure_pool is not None:
            try:
                ctx.flush_figures()
            finally:
                ctx._figure_pool.shutdown(wait=True)
                ctx._figure_pool = None
