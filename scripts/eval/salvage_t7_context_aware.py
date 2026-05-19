"""Salvage T7 candidate annotation from existing after_annotation.h5ad checkpoint.

The full pipeline OOM'd at clustering's GPU post-processing step (post-Harmony
neighbors/UMAP/Leiden recomputation) because the legacy context_aware_annotation
module called X.todense() on the full 795k x 30k sparse matrix.

This script:
  1. Loads the after_annotation.h5ad checkpoint (21GB on disk, sparse in RAM).
  2. Runs the PATCHED ContextAwareAnnotationModule (now sparse-aware, per-cluster
     slicing instead of full densification).
  3. Writes marker_intelligence_candidate.json to a metrics/ directory next to
     the source checkpoint.

Output: ready for T8 (compare_marker_intelligence.py) to compute AC-1/AC-2.
"""
from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace

import anndata as ad

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("salvage_t7")


@dataclass
class StubContext:
    adata: object
    random_state: int
    run_dir: Path
    cfg: object
    metadata: dict = field(default_factory=dict)
    _status: list = field(default_factory=list)

    def status(self, name: str, status: str, msg: str) -> None:
        self._status.append((name, status, msg))
        logger.info("status[%s] = %s (%s)", name, status, msg)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--checkpoint",
        default="/home/zerlinshen/projects/nc-reproduction/runs/2026-05-14T1738Z-d192836/python/NC2024_NSCLC_CANDIDATE_20260515_013812/.checkpoints/after_annotation.h5ad",
        help="Path to after_annotation.h5ad checkpoint",
    )
    parser.add_argument(
        "--run-dir",
        default="/home/zerlinshen/projects/nc-reproduction/runs/2026-05-14T1738Z-d192836/python/NC2024_NSCLC_CANDIDATE_20260515_013812",
        help="Run dir where metrics/ should be written",
    )
    parser.add_argument("--tissue", default="lung")
    parser.add_argument("--condition", default="NSCLC")
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--validate-context", action="store_true", default=True)
    args = parser.parse_args()

    cp = Path(args.checkpoint)
    if not cp.exists():
        logger.error("checkpoint not found: %s", cp)
        return 2

    logger.info("loading %s (this may take ~1 min for 21GB)", cp)
    adata = ad.read_h5ad(cp)
    logger.info("loaded: shape=%s, leiden in obs: %s, marker_db_index in uns: %s",
                adata.shape, "leiden" in adata.obs.columns, "marker_db_index" in adata.uns)

    if "leiden" not in adata.obs.columns:
        logger.error("missing required key: obs[leiden]")
        return 2
    sys.path.insert(0, "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory")
    from workflow.modular.modules.marker_db_loader import MarkerDbLoaderModule  # noqa: E402
    from workflow.modular.modules.context_aware_annotation import ContextAwareAnnotationModule  # noqa: E402

    cfg = SimpleNamespace(
        tissue=args.tissue,
        condition=args.condition,
        validate_context=args.validate_context,
        context_mismatch_threshold=0.3,
        context_min_cells=20,
    )
    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    ctx = StubContext(
        adata=adata,
        random_state=args.random_state,
        run_dir=run_dir,
        cfg=cfg,
    )

    # Re-run marker_db_loader to pick up the curated overlay (the checkpoint
    # h5ad has a stale marker_db_index from before the curated YAML existed).
    logger.info("re-running MarkerDbLoader to refresh marker_db_index (incl. curated overlay)...")
    MarkerDbLoaderModule().run(ctx)
    refreshed = adata.uns.get("marker_db_index", {})
    logger.info("refreshed marker_db_index: dbs=%s, total_entries=%d",
                list(refreshed.keys()),
                sum(len(df) for df in refreshed.values()))

    if not refreshed:
        logger.error("marker_db_index still empty after refresh — check ref/markers/")
        return 2

    module = ContextAwareAnnotationModule()
    logger.info("running ContextAwareAnnotationModule with refreshed index...")
    module.run(ctx)
    logger.info("done. metadata: %s", ctx.metadata)
    logger.info("statuses: %s", ctx._status)

    candidate_json = run_dir / "metrics" / "marker_intelligence_candidate.json"
    if candidate_json.exists():
        logger.info("candidate JSON written: %s (%d bytes)", candidate_json, candidate_json.stat().st_size)
    else:
        logger.error("candidate JSON NOT written — module may have early-exited")
        return 1

    # Save updated adata to a smaller sidecar (just obs columns we added)
    obs_out = run_dir / "context_aware_annotation"
    obs_out.mkdir(parents=True, exist_ok=True)
    cols_to_save = ["leiden", "context_aware_celltype"]
    if "context_aware_substate" in adata.obs.columns:
        cols_to_save.append("context_aware_substate")
    adata.obs[cols_to_save].to_csv(obs_out / "cell_type_annotation.csv", index=True)
    logger.info("obs columns CSV: %s", obs_out / "cell_type_annotation.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main())
