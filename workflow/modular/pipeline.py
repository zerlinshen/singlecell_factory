from __future__ import annotations

import copy
import logging
import os
from collections import deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import json
from time import perf_counter

import pandas as pd
from scipy import sparse

from .config import PipelineConfig
from .context import PipelineContext
from .module_catalog import MANDATORY_MODULES, module_dependencies, module_runs_after

logger = logging.getLogger(__name__)

# Defaults tuned for high-memory workstations (~96 GB RAM):
# - allow several large AnnData branch copies when parallel append-only modules run
# - still keep a safety reserve to avoid swapping/OOM on larger cohorts
PARALLEL_COPY_BUDGET_BYTES = 24 * 1024 * 1024 * 1024
MEMORY_RESERVE_BYTES = 8 * 1024 * 1024 * 1024

# Module dependency DAG in the legacy shape expected by older tests/importers.
# The canonical module hierarchy lives in module_catalog.py so CLI help, docs,
# validation scripts, and orchestration all share the same layer contract.
MODULE_DEPENDENCIES: dict[str, set[str]] = module_dependencies()

# Ordering-only hints (e.g. annotation runs_after batch_correction). These are
# deliberately SEPARATE from MODULE_DEPENDENCIES: they never auto-include a
# module (no inclusion auto-pull) and never seed _compute_tiers in-degree.
# They only refine sequencing among modules already in the requested set.
MODULE_RUNS_AFTER: dict[str, set[str]] = module_runs_after()

# Static fallback for modules that mutate adata structurally.
# Prefer the class-level `mutates_structure = True` attribute on modules;
# this set is a safety net for modules that forget to declare it.
_MUTATING_MODULES_FALLBACK = {"batch_correction"}
# Backward-compatibility alias used by older tests/importers.
MUTATING_MODULES = set(_MUTATING_MODULES_FALLBACK)


def _discover_mutating(registry: dict) -> set[str]:
    """Build the set of mutating modules from class attributes + static fallback."""
    discovered = {name for name, mod in registry.items()
                  if getattr(mod, "mutates_structure", False)}
    return discovered | MUTATING_MODULES


def _normalize_status(value: str) -> str:
    token = (value or "").strip().lower()
    if token in {"ok", "completed", "success"}:
        return "ok"
    if token in {"skipped", "skip"}:
        return "skipped"
    return "failed"


def _consume_module_report(
    ctx: PipelineContext,
    module_name: str,
    start_idx: int,
) -> tuple[str, str] | None:
    """Return module-reported status emitted during this stage, if any."""
    for entry in reversed(ctx.module_status[start_idx:]):
        if entry.get("module") != module_name:
            continue
        status = _normalize_status(str(entry.get("status", "failed")))
        message = str(entry.get("message", ""))
        entry["status"] = status
        return status, message
    return None


def _infer_status_from_metadata(ctx: PipelineContext, module_name: str) -> tuple[str, str] | None:
    """Infer skipped status from `<module>_status` metadata conventions."""
    key = f"{module_name}_status"
    raw = ctx.metadata.get(key)
    if not isinstance(raw, str):
        return None
    token = raw.strip().lower()
    if token.startswith("skip") or token.startswith("no_"):
        return "skipped", raw
    return None


def _finalize_stage_success(
    ctx: PipelineContext,
    module_name: str,
    start_idx: int,
) -> tuple[str, str]:
    """Finalize stage status after a successful module return."""
    reported = _consume_module_report(ctx, module_name, start_idx)
    if reported is not None:
        return reported
    inferred = _infer_status_from_metadata(ctx, module_name)
    if inferred is not None:
        status, message = inferred
        ctx.status(module_name, status, message)
        return status, message
    ctx.status(module_name, True, "completed")
    return "ok", "completed"


def _check_requires(mod, ctx: PipelineContext) -> list[str]:
    """Return list of missing keys declared in mod.requires_keys.

    Supports two obsm primitives:
    - ``obsm``: all listed keys must be present.
    - ``obsm_any_of``: at least one listed key must be present.
    """
    missing = []
    adata = ctx.adata
    if adata is None:
        return missing
    reqs = getattr(mod, "requires_keys", {})
    for key in reqs.get("obs", []):
        if key not in adata.obs.columns:
            missing.append(f"obs.{key}")
    for key in reqs.get("obsm", []):
        if key not in adata.obsm:
            missing.append(f"obsm.{key}")
    any_of_keys = reqs.get("obsm_any_of", [])
    if any_of_keys and not any(key in adata.obsm for key in any_of_keys):
        missing.append(f"obsm_any_of({','.join(any_of_keys)})")
    for key in reqs.get("uns", []):
        if key not in adata.uns:
            missing.append(f"uns.{key}")
    return missing


def _kahn_order(all_requested: set[str], adjacency: dict[str, set[str]]) -> list[str] | None:
    """Deterministic Kahn topo sort over `adjacency` restricted to `all_requested`.

    `adjacency[mod]` is the set of nodes `mod` must run after. Returns the
    ordered list, or ``None`` if a cycle/stall prevents emitting every node.
    """
    in_degree: dict[str, int] = {m: 0 for m in all_requested}
    for mod in all_requested:
        for dep in adjacency.get(mod, set()):
            if dep in all_requested:
                in_degree[mod] += 1

    queue = deque(sorted(m for m, d in in_degree.items() if d == 0))
    order: list[str] = []
    while queue:
        node = queue.popleft()
        order.append(node)
        for mod in sorted(all_requested):
            if node in adjacency.get(mod, set()):
                in_degree[mod] -= 1
                if in_degree[mod] == 0:
                    queue.append(mod)

    if len(order) != len(all_requested):
        return None
    return order


def _resolve_execution_order(
    mandatory: list[str],
    optional: list[str],
    dropped_hints_sink: list[dict] | None = None,
) -> list[str]:
    """Topologically sort modules respecting dependencies.

    Mandatory modules always run first. For optional modules, any missing
    dependencies that are themselves optional are auto-included.

    A SECOND, ordering-only pass folds in ``MODULE_RUNS_AFTER`` hints (e.g.
    ``annotation`` runs_after ``batch_correction``). These hints:
      - are restricted to pairs where BOTH endpoints are already in
        ``all_requested`` (they never auto-include — that is the exclusive
        job of ``MODULE_DEPENDENCIES`` above), and
      - never feed ``_compute_tiers`` (which keys off MODULE_DEPENDENCIES).

    The hard cycle ``raise`` stays bound to the ``depends_on``-only graph. If
    the COMBINED (depends_on + runs_after) graph stalls, the runs_after edges
    are DROPPED (never raise) and a loud breadcrumb is recorded so the run
    manifest + log surface that the ordering hint could not be honored.
    """
    all_requested = set(mandatory) | set(optional)

    # Auto-include transitive dependencies (depends_on ONLY — runs_after must
    # never force-pull a module into the run).
    to_process = list(all_requested)
    while to_process:
        mod = to_process.pop()
        for dep in MODULE_DEPENDENCIES.get(mod, set()):
            if dep not in all_requested:
                all_requested.add(dep)
                to_process.append(dep)

    # First pass: depends_on-only graph. A failure here is a real dependency
    # cycle and must raise (preserves the historical contract + tests).
    base_order = _kahn_order(all_requested, MODULE_DEPENDENCIES)
    if base_order is None:
        # Recompute remaining the same way the legacy code did for the message.
        in_degree: dict[str, int] = {m: 0 for m in all_requested}
        for mod in all_requested:
            for dep in MODULE_DEPENDENCIES.get(mod, set()):
                if dep in all_requested:
                    in_degree[mod] += 1
        queue = deque(sorted(m for m, d in in_degree.items() if d == 0))
        emitted: set[str] = set()
        while queue:
            node = queue.popleft()
            emitted.add(node)
            for mod in sorted(all_requested):
                if node in MODULE_DEPENDENCIES.get(mod, set()):
                    in_degree[mod] -= 1
                    if in_degree[mod] == 0:
                        queue.append(mod)
        remaining = sorted(m for m in all_requested if m not in emitted)
        raise ValueError(f"Cyclic dependency detected: {remaining}")

    # Second pass (ordering-only): combine depends_on with runs_after hints,
    # restricted to pairs already present in all_requested. Never auto-include.
    active_hints: dict[str, set[str]] = {}
    for mod in all_requested:
        hints = {
            after for after in MODULE_RUNS_AFTER.get(mod, set())
            if after in all_requested
        }
        if hints:
            active_hints[mod] = hints

    if not active_hints:
        return base_order

    combined: dict[str, set[str]] = {
        m: set(MODULE_DEPENDENCIES.get(m, set())) for m in all_requested
    }
    for mod, hints in active_hints.items():
        combined[mod] = combined.get(mod, set()) | hints

    combined_order = _kahn_order(all_requested, combined)
    if combined_order is not None:
        return combined_order

    # Combined-graph STALL: drop runs_after edges (NEVER raise) and emit a
    # loud breadcrumb. The base depends_on order is still valid and used.
    dropped = {mod: sorted(hints) for mod, hints in sorted(active_hints.items())}
    logger.warning(
        "runs_after ordering hints could not be honored without stalling the "
        "combined graph; DROPPING ordering hints and proceeding on the "
        "depends_on order. Dropped hints: %s",
        dropped,
    )
    if dropped_hints_sink is not None:
        dropped_hints_sink.append({
            "reason": "combined_graph_stall",
            "dropped_runs_after": dropped,
        })
    return base_order


def _build_registry() -> dict[str, object]:
    """Lazy-import all modules and build the registry.

    Importing inside this function avoids pulling heavy dependencies
    (scrublet, scvelo, liana, etc.) at package-level import time.
    """
    from .modules.ambient_correction import AmbientCorrectionModule
    from .modules.annotation import AnnotationModule
    from .modules.batch_correction import BatchCorrectionModule
    from .modules.evolution import EvolutionModule
    from .modules.cell_communication import CellCommunicationModule
    from .modules.cell_cycle import CellCycleModule
    from .modules.cellranger import CellRangerModule
    from .modules.clustering import ClusteringModule
    from .modules.cnv_inference import CNVInferenceModule
    from .modules.differential_expression import DifferentialExpressionModule
    from .modules.doublet_detection import DoubletDetectionModule
    from .modules.gene_regulatory_network import GeneRegulatoryNetworkModule
    from .modules.gene_signature_scoring import GeneSignatureScoringModule
    from .modules.immune_phenotyping import ImmunePhenotypingModule
    from .modules.integration_select import IntegrationSelectModule
    from .modules.pathway_analysis import PathwayAnalysisModule
    from .modules.pseudo_velocity import PseudoVelocityModule
    from .modules.qc import QCModule
    from .modules.rna_velocity import RNAVelocityModule
    from .modules.trajectory import TrajectoryModule
    from .modules.tumor_microenvironment import TumorMicroenvironmentModule
    from .modules.validate_cbioportal import ValidateCbioPortalModule
    from .modules.pseudobulk_de import PseudobulkDEModule
    from .modules.cell_fate import CellFateModule
    from .modules.composition import CompositionModule
    from .modules.metacell import MetacellModule
    from .modules.paper_repro import PaperReproModule
    from .modules.protein_adt import ProteinADTModule
    from .modules.spatial_ingest import SpatialIngestModule
    from .modules.spatial_neighborhoods import SpatialNeighborhoodsModule
    from .modules.multimodal_integration import MultimodalIntegrationModule
    from .modules.marker_db_loader import MarkerDbLoaderModule
    from .modules.context_aware_annotation import ContextAwareAnnotationModule
    from .modules.modality_registry import ModalityRegistryModule
    from .modules.cross_modality_qc import CrossModalityQCModule
    from .modules.atac_ingest import ATACIngestModule
    from .modules.vdj_ingest import VDJIngestModule
    from .modules.vdj_metrics import VDJMetricsModule
    from .modules.atac_qc import ATACQCModule
    from .modules.atac_lsi import ATACLSIModule
    from .modules.peak_to_gene import PeakToGeneModule
    from .modules.hic_ingest import HiCIngestModule
    from .modules.hic_tad import HiCTADModule
    from .modules.ribo_ingest import RiboIngestModule

    return {
        "cellranger": CellRangerModule(),
        "qc": QCModule(),
        "ambient_correction": AmbientCorrectionModule(),
        "doublet_detection": DoubletDetectionModule(),
        "clustering": ClusteringModule(),
        "cell_cycle": CellCycleModule(),
        "integration_select": IntegrationSelectModule(),
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
        "validate_cbioportal": ValidateCbioPortalModule(),
        "immune_phenotyping": ImmunePhenotypingModule(),
        "tumor_microenvironment": TumorMicroenvironmentModule(),
        "gene_signature_scoring": GeneSignatureScoringModule(),
        "evolution": EvolutionModule(),
        "pseudobulk_de": PseudobulkDEModule(),
        "cell_fate": CellFateModule(),
        "composition": CompositionModule(),
        "metacell": MetacellModule(),
        "paper_repro": PaperReproModule(),
        "protein_adt": ProteinADTModule(),
        "spatial_ingest": SpatialIngestModule(),
        "spatial_neighborhoods": SpatialNeighborhoodsModule(),
        "atac_ingest": ATACIngestModule(),
        "atac_qc": ATACQCModule(),
        "atac_lsi": ATACLSIModule(),
        "multimodal_integration": MultimodalIntegrationModule(),
        "marker_db_loader": MarkerDbLoaderModule(),
        "context_aware_annotation": ContextAwareAnnotationModule(),
        "modality_registry": ModalityRegistryModule(),
        "cross_modality_qc": CrossModalityQCModule(),
        "vdj_ingest": VDJIngestModule(),
        "vdj_metrics": VDJMetricsModule(),
        "peak_to_gene": PeakToGeneModule(),
        "hic_ingest": HiCIngestModule(),
        "hic_tad": HiCTADModule(),
        "ribo_ingest": RiboIngestModule(),
    }


def _find_latest_resume_run_dir(cfg: PipelineConfig) -> Path | None:
    """Return the latest run dir for `project` that contains checkpoints."""
    if not cfg.output_dir.exists():
        return None
    candidates = [
        p for p in cfg.output_dir.glob(f"{cfg.project}_*")
        if p.is_dir() and (p / ".checkpoints").exists()
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.stat().st_mtime)
    return candidates[-1]


def _prepare_output(cfg: PipelineConfig) -> PipelineContext:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    run_dir: Path
    if cfg.resume_from:
        run_dir = _find_latest_resume_run_dir(cfg) or (  # fallback keeps old behavior/error path
            cfg.output_dir / f"{cfg.project}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
    else:
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
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
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


# Estimated relative cost per module (higher = slower). Used for longest-job-first
# scheduling within parallel tiers so heavy modules start first.
_MODULE_COST: dict[str, int] = {
    "rna_velocity": 10,
    "cnv_inference": 5,
    "evolution": 5,
    "trajectory": 4,
    "differential_expression": 3,
    "clustering": 3,
    "immune_phenotyping": 2,
    "tumor_microenvironment": 2,
    "pathway_analysis": 2,
    "cell_communication": 2,
    "pseudobulk_de": 2,
    "metacell": 2,
    "cell_fate": 2,
}


def _compute_tiers(execution_order: list[str], completed: set[str]) -> list[list[str]]:
    """Group remaining modules into tiers for parallel execution.

    A tier is a set of modules whose dependencies are all satisfied by
    previously completed tiers plus the *completed* set.  Within each tier,
    modules are sorted by estimated cost (heaviest first) so that the
    longest-running job starts first in the thread pool.
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
        # Sort by cost descending so heavy modules start first in parallel
        tier.sort(key=lambda m: _MODULE_COST.get(m, 1), reverse=True)
        tiers.append(tier)
        done.update(tier)
        remaining = [m for m in remaining if m not in done]
    return tiers


class _SkipModule(Exception):
    """Raised when a non-mandatory module should be skipped due to missing keys or memory."""


def _run_module(mod, ctx: PipelineContext, *, mandatory: bool = False) -> None:
    """Validate requires_keys then run a module. Raises _SkipModule for optional modules."""
    # Wave 3 US-W3-2 — orchestrator-level poison guard. If a prior module
    # raised ClusteringContractViolation and a caller accidentally swallowed
    # it, the next module entry MUST abort here. See _contract_violation.py.
    from ._contract_violation import assert_not_corrupted
    name = getattr(mod, "name", type(mod).__name__)
    assert_not_corrupted(ctx.adata, name)
    missing = _check_requires(mod, ctx)
    if missing:
        name = getattr(mod, "name", type(mod).__name__)
        msg = f"Module '{name}' missing required keys: {', '.join(missing)}"
        if mandatory:
            raise ValueError(msg)
        logger.warning("%s — skipping.", msg)
        raise _SkipModule(msg)
    try:
        mod.run(ctx)
    except Exception as exc:
        # Convert cooperative memory-abort signals to _SkipModule so the
        # pipeline records status="skipped_memory" and continues.
        try:
            from ._mem_guard import MemoryAbortError
        except ImportError:
            raise
        if isinstance(exc, MemoryAbortError):
            name = getattr(mod, "name", type(mod).__name__)
            msg = f"memory abort: {exc}"
            logger.warning("Module '%s' cooperatively aborted: %s", name, exc)
            raise _SkipModule(msg) from exc
        raise


def _record_module_runtime(ctx: PipelineContext, module_name: str, elapsed_seconds: float) -> None:
    """Record per-module wall-time (seconds) in run metadata."""
    runtimes = ctx.metadata.setdefault("module_runtime_sec", {})
    runtimes[module_name] = round(float(elapsed_seconds), 3)


def _estimate_adata_copy_bytes(adata) -> int:
    """Estimate the memory footprint of one AnnData copy."""
    if adata is None:
        return 0

    def _array_bytes(value) -> int:
        if sparse.issparse(value):
            return value.data.nbytes + value.indices.nbytes + value.indptr.nbytes
        return int(getattr(value, "nbytes", 0))

    total = _array_bytes(adata.X)
    total += sum(_array_bytes(layer) for layer in adata.layers.values())
    total += sum(_array_bytes(obsm) for obsm in adata.obsm.values())
    total += sum(_array_bytes(varm) for varm in adata.varm.values())
    # Include obs/var DataFrame memory for more accurate estimates
    if hasattr(adata.obs, "memory_usage"):
        total += int(adata.obs.memory_usage(deep=True).sum())
    if hasattr(adata.var, "memory_usage"):
        total += int(adata.var.memory_usage(deep=True).sum())
    return int(total)


def _get_available_memory_bytes() -> int | None:
    """Return currently available system memory in bytes, if detectable."""
    try:
        import psutil  # type: ignore

        return int(psutil.virtual_memory().available)
    except Exception:
        pass

    try:
        pages = os.sysconf("SC_AVPHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return int(pages * page_size)
    except (AttributeError, ValueError, OSError):
        return None


def _safe_parallel_worker_count(
    ctx: PipelineContext,
    modules: list[str],
    requested_workers: int,
) -> int:
    """Choose a safe parallel worker count based on copy size and free memory."""
    if ctx.adata is None or not modules:
        return 0

    est_copy_bytes = max(_estimate_adata_copy_bytes(ctx.adata), 1)
    hard_budget_workers = max(1, PARALLEL_COPY_BUDGET_BYTES // est_copy_bytes)

    available = _get_available_memory_bytes()
    if available is None:
        memory_workers = hard_budget_workers
    else:
        reserve = min(MEMORY_RESERVE_BYTES, max(available // 4, 0))
        usable = max(0, available - reserve)
        memory_workers = max(1, usable // est_copy_bytes)

    return max(
        0,
        min(requested_workers, len(modules), int(hard_budget_workers), int(memory_workers)),
    )


def _ledger_record_module(ctx: PipelineContext, name: str, status: str, message: str, elapsed: float) -> None:
    """Call ledger.record_module if a ledger is attached to ctx. Never raises."""
    ledger = getattr(ctx, "_ledger", None)
    if ledger is None:
        return
    try:
        import psutil  # type: ignore
        rss = psutil.Process().memory_info().rss
    except Exception:
        rss = 0
    try:
        ledger.record_module(name, status, message, elapsed, rss)
    except Exception as exc:
        logger.warning("ledger.record_module failed: %s", exc)


def _run_sequential(
    modules: list[str],
    registry: dict,
    ctx: PipelineContext,
    mandatory: set[str],
) -> None:
    """Run a list of modules sequentially with checkpoint, status, and timing."""
    for stage in modules:
        mod = registry.get(stage)
        if mod is None:
            ctx.status(stage, False, "unknown module")
            _ledger_record_module(ctx, stage, "failed", "unknown module", 0.0)
            continue
        ctx.set_module_dir(stage)
        pre_status_len = len(ctx.module_status)
        t0 = perf_counter()
        status, message = "failed", ""
        try:
            _run_module(mod, ctx, mandatory=stage in mandatory)
            status, message = _finalize_stage_success(ctx, stage, pre_status_len)
            if status != "ok" and stage in mandatory:
                raise RuntimeError(f"Mandatory module {stage} reported {status}: {message}")
            ctx.save_checkpoint(stage)
        except _SkipModule as exc:
            msg_str = str(exc)
            if "memory" in msg_str.lower() or "abort" in msg_str.lower():
                status, message = "skipped_memory", msg_str
            else:
                status, message = "skipped", msg_str
            ctx.status(stage, status, msg_str)
            ctx.save_checkpoint(stage)
        except Exception as exc:
            status, message = "failed", str(exc)
            ctx.status(stage, False, str(exc))
            if stage in mandatory:
                raise
        finally:
            elapsed = perf_counter() - t0
            _record_module_runtime(ctx, stage, elapsed)
            _ledger_record_module(ctx, stage, status, message, elapsed)


def _execute_tier(
    tier: list[str],
    registry: dict,
    ctx: PipelineContext,
    mandatory: set[str],
    max_workers: int,
    mutating_modules: set[str] | None = None,
) -> None:
    """Execute a tier of modules, potentially in parallel."""
    _mutating_set = mutating_modules if mutating_modules is not None else _MUTATING_MODULES_FALLBACK
    # Split into mutating (must run sequentially) and appending (can parallelize).
    mutating = [m for m in tier if m in _mutating_set]
    appending = [m for m in tier if m not in _mutating_set]

    # Run mutating modules sequentially first.
    _run_sequential(mutating, registry, ctx, mandatory)

    # Run appending modules.
    if max_workers <= 1 or len(appending) <= 1:
        _run_sequential(appending, registry, ctx, mandatory)
    else:
        safe_workers = _safe_parallel_worker_count(ctx, appending, max_workers)
        if safe_workers <= 1:
            logger.warning(
                "Falling back to sequential appending execution due to memory safety guard."
            )
            _run_sequential(appending, registry, ctx, mandatory)
            return
        _run_parallel_appending(appending, registry, ctx, mandatory, safe_workers)


def _warn_dropped_changes(mod_name: str, branch_ad, main_ad) -> None:
    """Log warnings when a parallel branch modified data that merge-back cannot capture."""
    import numpy as np
    from scipy.sparse import issparse

    # --- X mutation check: shape, dtype, and sampled rows ---
    if branch_ad.X is not None and main_ad.X is not None and branch_ad.n_obs > 0:
        x_mutated = False
        if branch_ad.X.shape != main_ad.X.shape:
            x_mutated = True
        elif hasattr(branch_ad.X, "dtype") and hasattr(main_ad.X, "dtype"):
            if branch_ad.X.dtype != main_ad.X.dtype:
                x_mutated = True
        if not x_mutated:
            try:
                # Sample first, middle, and last rows for coverage without full comparison
                n = branch_ad.n_obs
                indices = sorted({0, n // 2, n - 1})
                for idx in indices:
                    b_row = branch_ad.X[idx]
                    m_row = main_ad.X[idx]
                    b_arr = b_row.toarray().ravel() if issparse(b_row) else np.asarray(b_row).ravel()
                    m_arr = m_row.toarray().ravel() if issparse(m_row) else np.asarray(m_row).ravel()
                    if not np.array_equal(b_arr, m_arr):
                        x_mutated = True
                        break
            except Exception:
                pass  # shape mismatch or other edge case — skip check
        if x_mutated:
            logger.warning(
                "Module '%s' modified adata.X but runs as appending "
                "— X changes dropped. Consider adding mutates_structure = True.",
                mod_name,
            )

    # --- Layer checks: new layers AND modifications to existing layers ---
    for layer in set(branch_ad.layers) - set(main_ad.layers):
        logger.warning(
            "Module '%s' added layer '%s' — dropped in parallel merge-back.", mod_name, layer,
        )
    for layer in set(branch_ad.layers) & set(main_ad.layers):
        try:
            b_val = branch_ad.layers[layer]
            m_val = main_ad.layers[layer]
            b_arr = b_val.toarray().ravel()[:100] if issparse(b_val) else np.asarray(b_val).ravel()[:100]
            m_arr = m_val.toarray().ravel()[:100] if issparse(m_val) else np.asarray(m_val).ravel()[:100]
            if not np.array_equal(b_arr, m_arr):
                logger.warning(
                    "Module '%s' modified existing layer '%s' — changes dropped in parallel merge-back.",
                    mod_name, layer,
                )
        except Exception:
            pass

    # --- varm / obsp: new keys ---
    if hasattr(branch_ad, "varm"):
        for key in set(branch_ad.varm) - set(main_ad.varm):
            logger.warning(
                "Module '%s' added varm['%s'] — dropped in parallel merge-back.", mod_name, key,
            )
    if hasattr(branch_ad, "obsp"):
        for key in set(branch_ad.obsp) - set(main_ad.obsp):
            logger.warning(
                "Module '%s' added obsp['%s'] — dropped in parallel merge-back.", mod_name, key,
            )


def _run_parallel_appending(
    modules: list[str],
    registry: dict,
    ctx: PipelineContext,
    mandatory: set[str],
    max_workers: int,
) -> None:
    """Run appending-only modules in parallel with copy-on-branch, merge-back."""
    def _run_module_in_branch(mod_name: str) -> tuple[str, str, str, PipelineContext]:
        mod = registry.get(mod_name)
        if mod is None:
            raise ValueError(f"unknown module: {mod_name}")
        branch_ctx = copy.copy(ctx)
        branch_ctx.adata = ctx.adata.copy()
        branch_ctx.module_status = []
        branch_ctx.metadata = dict(ctx.metadata)
        branch_ctx._figure_futures = []
        branch_ctx._module_dirs = dict(ctx._module_dirs)
        branch_ctx._figure_pool = None
        branch_ctx.set_module_dir(mod_name)
        pre_status_len = len(branch_ctx.module_status)
        try:
            _run_module(mod, branch_ctx, mandatory=mod_name in mandatory)
        except _SkipModule as exc:
            branch_ctx.status(mod_name, "skipped", str(exc))
            return mod_name, "skipped", str(exc), branch_ctx
        status, message = _finalize_stage_success(branch_ctx, mod_name, pre_status_len)
        return mod_name, status, message, branch_ctx

    results: dict[str, tuple[str, str, float]] = {}
    branch_contexts: dict[str, PipelineContext] = {}
    with ThreadPoolExecutor(max_workers=min(max_workers, len(modules))) as pool:
        futures = {pool.submit(_run_module_in_branch, name): name for name in modules}
        starts = {future: perf_counter() for future in futures}
        for future in as_completed(futures):
            name = futures[future]
            elapsed = perf_counter() - starts[future]
            try:
                returned_name, status, message, branch_ctx = future.result()
                branch_contexts[returned_name] = branch_ctx
                results[name] = (status, message, elapsed)
            except Exception as exc:
                results[name] = ("failed", str(exc), elapsed)

    # Merge results back into the main context.
    main_obs_cols = set(ctx.adata.obs.columns)
    main_obsm_keys = set(ctx.adata.obsm.keys())

    for mod_name in modules:
        ok, msg, elapsed = results.get(mod_name, ("failed", "not run", 0.0))

        if ok == "ok" and mod_name in branch_contexts:
            branch_ctx = branch_contexts[mod_name]
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
            # Warn about structural changes that cannot be merged back.
            _warn_dropped_changes(mod_name, branch_ctx.adata, ctx.adata)
            # Merge metadata and module directory registrations.
            # MEDIUM fix (2026-05-19): merge deterministically and warn on
            # cross-module metadata key conflicts so a silent overwrite in
            # parallel-appending tiers cannot produce non-reproducible
            # run_manifest.json values. The outer iteration already runs
            # in fixed `modules` order, so the *last writer* is
            # deterministic, but we still want triage visibility.
            conflicts = sorted(
                k for k in branch_ctx.metadata
                if k in ctx.metadata and ctx.metadata[k] != branch_ctx.metadata[k]
            )
            if conflicts:
                logger.warning(
                    "parallel-appending metadata key conflict in module %s: %s "
                    "(last writer in modules-order wins)",
                    mod_name,
                    conflicts[:5],
                )
                ctx.metadata.setdefault("_parallel_metadata_conflicts", {})[mod_name] = conflicts
            ctx.metadata.update(branch_ctx.metadata)
            ctx._module_dirs.update(branch_ctx._module_dirs)
            # Flush any figures the branch produced.
            branch_ctx.flush_figures()

        ctx.status(mod_name, ok, msg)
        _record_module_runtime(ctx, mod_name, elapsed)
        _ledger_record_module(ctx, mod_name, ok, msg, elapsed)
        ctx.save_checkpoint(mod_name)

        if ok != "ok" and mod_name in mandatory:
            raise RuntimeError(f"Mandatory module {mod_name} failed: {msg}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run_pipeline(cfg: PipelineConfig, ledger=None) -> Path:
    """Run the modular workflow with mandatory and optional stages."""

    pipeline_t0 = perf_counter()
    ctx = _prepare_output(cfg)
    if ledger is not None:
        ctx._ledger = ledger
    registry = _build_registry()
    mutating_set = _discover_mutating(registry)
    _watchdog_thread = None
    try:
        # Enable async figure pool if parallel workers > 1.
        if cfg.parallel_workers > 1:
            ctx._figure_pool = ThreadPoolExecutor(max_workers=1)

        # Start memory watchdog when SC_MEM_GUARD is active.
        if os.environ.get("SC_MEM_GUARD", "").lower() == "on":
            from . import _mem_watchdog
            from ._mem_guard import MemoryGuard
            MemoryGuard.clear_abort()
            _watchdog_thread = _mem_watchdog.start(ctx)

        mandatory = list(MANDATORY_MODULES)
        mandatory_set = set(mandatory)
        dropped_hints_sink: list[dict] = []
        execution_order = _resolve_execution_order(
            mandatory, cfg.optional_modules, dropped_hints_sink
        )
        if dropped_hints_sink:
            # Loud breadcrumb in run_manifest.json: ordering hints were dropped
            # to avoid a combined-graph stall (annotation/batch_correction
            # leiden race fix). The depends_on order was used instead.
            ctx.metadata["runs_after_hints_dropped"] = dropped_hints_sink

        # --- Resume from checkpoint ---
        if cfg.resume_from:
            if not ctx._checkpoint_dir.exists():
                raise FileNotFoundError(
                    f"No checkpoint directory found for project '{cfg.project}'. "
                    "Run once with --checkpoint before using --resume-from."
                )
            # Find the module just before resume_from in execution order.
            try:
                resume_idx = execution_order.index(cfg.resume_from)
            except ValueError:
                raise ValueError(
                    f"Cannot resume from '{cfg.resume_from}': not in execution order."
                )
            if resume_idx > 0:
                # Search backwards for the nearest available checkpoint.
                loaded = False
                for search_idx in range(resume_idx - 1, -1, -1):
                    candidate = execution_order[search_idx]
                    if ctx.load_checkpoint(candidate):
                        logger.info("Resumed from checkpoint after '%s'", candidate)
                        loaded = True
                        break
                if not loaded:
                    raise FileNotFoundError(
                        f"No checkpoint found before '{cfg.resume_from}'. "
                        f"Run the pipeline with --checkpoint first."
                    )
            execution_order = execution_order[resume_idx:]

        # --- Execute ---
        if cfg.parallel_workers > 1:
            all_modules = _resolve_execution_order(mandatory, cfg.optional_modules)
            completed = {m for m in all_modules if m not in execution_order}
            tiers = _compute_tiers(execution_order, completed)
            for tier in tiers:
                _execute_tier(tier, registry, ctx, mandatory_set, cfg.parallel_workers, mutating_set)
        else:
            # Original sequential execution.
            _run_sequential(execution_order, registry, ctx, mandatory_set)

        ctx.metadata["pipeline_wall_seconds"] = round(perf_counter() - pipeline_t0, 3)
        manifest_path = _save_manifest(ctx)
        try:
            ledger = getattr(ctx, "_ledger", None)
            if ledger is not None:
                final_adata_path = ctx.run_dir / "final_adata.h5ad"
                ledger.record_end(final_adata_path if final_adata_path.exists() else None)
                ledger.write()
        except Exception as _ledger_exc:
            logger.warning("RunLedger finalization failed: %s", _ledger_exc)
        return manifest_path
    finally:
        if _watchdog_thread is not None:
            try:
                from . import _mem_watchdog
                _mem_watchdog.stop(_watchdog_thread)
            except Exception as _wd_exc:
                logger.warning("Watchdog stop failed: %s", _wd_exc)
        if ctx._figure_pool is not None:
            try:
                ctx.flush_figures()
            finally:
                ctx._figure_pool.shutdown(wait=True)
                ctx._figure_pool = None
