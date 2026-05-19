"""Context-aware annotation module: (tissue, condition)-driven cell type annotation.

Consumes adata.uns["marker_db_index"] from marker_db_loader and cluster labels
from clustering. Produces adata.obs["context_aware_celltype"] and optionally
adata.obs["context_aware_substate"]. Leaves annotation.py DEFAULT_MARKERS
untouched (AC-10 zero-regression).
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "CellMarker2": {
        "title": "CellMarker 2.0: an updated database of manually curated cell markers in human/mouse",
        "authors": "Hu et al.",
        "journal": "Nucleic Acids Research",
        "year": "2023",
        "doi": "10.1093/nar/gkac947",
        "description": "Primary marker DB consumed by marker_db_loader; used here only as a DB source, not as a scoring method.",
    },
}

# Implementation note (2026-05-19 scientific code review):
# - Cluster -> cell-type scoring is `mean(cluster_X[:, marker_idx])` followed
#   by argmax with a tolerance tie-break. This is an in-house cluster-voting
#   heuristic in the spirit of Tirosh et al. 2016, NOT the scType
#   specificity-weighted scoring (Ianevski 2022). The previous reference to
#   scTypeDB was removed because the implementation does not match that
#   method.
# - Substate detection uses `startswith()` string matching against cell-type
#   names — this is a pragmatic in-house convention with no published
#   precedent. It is intentionally simple; treat results as exploratory.

# Default threshold for context validation warnings (AC-4)
_DEFAULT_MISMATCH_THRESHOLD = 0.3
# Minimum cells per cluster for validation (AC-4, R5)
_DEFAULT_MIN_CELLS = 20


class ContextAwareAnnotationModule:
    """Opt-in context-aware cell type annotation driven by (tissue, condition) marker DBs.

    Requires:
      - adata.obs["leiden"] from clustering
      - adata.uns["marker_db_index"] from marker_db_loader

    Produces:
      - adata.obs["context_aware_celltype"]
      - adata.obs["context_aware_substate"] (when condition-specific sub-states exist)
      - runs/<run-id>/metrics/marker_intelligence_candidate.json
      - runs/<run-id>/metrics/context_validation.json (when --validate-context)
    """

    name = "context_aware_annotation"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "obs": ["leiden"],
        "uns": ["marker_db_index"],
    }
    provides_keys: dict[str, list[str]] = {
        "obs": ["context_aware_celltype"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata

        if "leiden" not in adata.obs:
            raise ValueError(f"{self.name} requires clustering (leiden column absent).")

        marker_db_index = adata.uns.get("marker_db_index", {})
        if not marker_db_index:
            logger.warning(
                "%s: marker_db_index is empty — run marker_db_loader first. Skipping.", self.name
            )
            ctx.status(self.name, "skipped", "marker_db_index empty")
            return

        # Seed for any tie-breaking (AC-10 reproducibility)
        rng = np.random.default_rng(ctx.random_state)

        # Build unified marker table from all loaded DBs
        db_frames = list(marker_db_index.values())
        all_markers = pd.concat(db_frames, ignore_index=True)

        # Score each cluster against each cell type
        cluster_labels = adata.obs["leiden"].astype(str)
        clusters = sorted(cluster_labels.unique())
        cluster_celltype: dict[str, str] = {}
        cluster_substate: dict[str, str | None] = {}
        cluster_scores: dict[str, float] = {}

        # Build per-cell-type gene sets (union across DBs, positive markers only)
        cell_types = all_markers["cell_type"].unique()
        ct_genes: dict[str, list[str]] = {}
        for ct in cell_types:
            ct_df = all_markers[all_markers["cell_type"] == ct]
            pos = ct_df[ct_df["marker_type"] == "positive"]["marker_gene"].tolist()
            ct_genes[ct] = [g for g in pos if g in adata.var_names]

        ct_genes = {ct: genes for ct, genes in ct_genes.items() if genes}
        if not ct_genes:
            logger.warning("%s: no valid marker genes found in dataset var_names. Skipping.", self.name)
            ctx.status(self.name, "skipped", "no marker genes in dataset")
            return

        # Memory-bounded per-cluster mean expression:
        # avoid densifying the full X matrix (would be ~93GB at 795k cells × 30k genes).
        # Sparse-aware loop: slice per cluster only, mean across cells natively.
        import scipy.sparse as sp
        X = adata.X
        is_sparse = sp.issparse(X)

        gene_idx = {g: i for i, g in enumerate(adata.var_names)}

        for cluster in clusters:
            mask = cluster_labels.values == cluster
            cell_count = mask.sum()
            if cell_count == 0:
                cluster_celltype[cluster] = "Unknown"
                cluster_substate[cluster] = None
                cluster_scores[cluster] = 0.0
                continue

            if is_sparse:
                cluster_X = X[mask]
                # sparse .mean returns np.matrix shape (1, n_genes); flatten to 1-D
                cluster_mean = np.asarray(cluster_X.mean(axis=0)).flatten().astype(np.float32)
            else:
                cluster_X = np.asarray(X[mask], dtype=np.float32)
                cluster_mean = cluster_X.mean(axis=0)

            # Score each cell type: mean expression of positive markers
            scores: dict[str, float] = {}
            for ct, genes in ct_genes.items():
                idxs = [gene_idx[g] for g in genes if g in gene_idx]
                if not idxs:
                    scores[ct] = 0.0
                    continue
                scores[ct] = float(cluster_mean[idxs].mean())

            if not scores or max(scores.values()) == 0.0:
                cluster_celltype[cluster] = "Unknown"
                cluster_substate[cluster] = None
                cluster_scores[cluster] = 0.0
                continue

            # Tie-break with rng for determinism (AC-10)
            best_score = max(scores.values())
            top_cts = [ct for ct, s in scores.items() if s >= best_score * 0.999]
            if len(top_cts) > 1:
                best_ct = top_cts[int(rng.integers(0, len(top_cts)))]
            else:
                best_ct = top_cts[0]

            cluster_celltype[cluster] = best_ct
            cluster_scores[cluster] = best_score

            # Sub-state detection: look for condition-specific sub-types within best_ct
            subtype_rows = all_markers[
                (all_markers["cell_type"].str.startswith(best_ct)) &
                (all_markers["cell_type"] != best_ct)
            ]
            substate_scores: dict[str, float] = {}
            for sub_ct in subtype_rows["cell_type"].unique():
                sub_genes = subtype_rows[
                    (subtype_rows["cell_type"] == sub_ct) &
                    (subtype_rows["marker_type"] == "positive")
                ]["marker_gene"].tolist()
                sub_genes = [g for g in sub_genes if g in gene_idx]
                if not sub_genes:
                    continue
                idxs = [gene_idx[g] for g in sub_genes]
                substate_scores[sub_ct] = float(cluster_mean[idxs].mean())

            if substate_scores and max(substate_scores.values()) > 0.0:
                cluster_substate[cluster] = max(substate_scores, key=lambda k: substate_scores[k])
            else:
                cluster_substate[cluster] = None

        # Map cluster labels to cells
        adata.obs["context_aware_celltype"] = cluster_labels.map(cluster_celltype).astype(str)

        has_substate = any(v is not None for v in cluster_substate.values())
        if has_substate:
            adata.obs["context_aware_substate"] = cluster_labels.map(
                lambda c: cluster_substate.get(c)
            )
            logger.info("%s: condition-specific sub-states detected for some clusters.", self.name)

        # Write candidate metrics JSON (AC-1/AC-2 eval input)
        metrics_dir = ctx.run_dir / "metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)
        candidate_json = {
            "run_id": ctx.run_dir.name,
            "tissue": getattr(ctx.cfg, "tissue", "lung"),
            "condition": getattr(ctx.cfg, "condition", "NSCLC"),
            "n_clusters": len(clusters),
            "per_cluster": [
                {
                    "cluster_id": c,
                    "cell_type": cluster_celltype[c],
                    "substate": cluster_substate.get(c),
                    "score": cluster_scores[c],
                }
                for c in clusters
            ],
            "dbs_used": list(marker_db_index.keys()),
        }
        (metrics_dir / "marker_intelligence_candidate.json").write_text(
            json.dumps(candidate_json, indent=2, ensure_ascii=False), encoding="utf-8"
        )

        # Optional --validate-context validation (AC-4)
        validate = getattr(ctx.cfg, "validate_context", False)
        if validate:
            self._run_context_validation(
                ctx, adata, clusters, cluster_celltype, cluster_scores,
                X, is_sparse, gene_idx, ct_genes, metrics_dir,
            )

        ctx.metadata["context_aware_n_cell_types"] = len(set(cluster_celltype.values()))
        ctx.metadata["context_aware_has_substates"] = has_substate
        ctx.metadata["context_aware_random_state"] = ctx.random_state
        logger.info(
            "%s: annotated %d clusters → %d cell types (substates: %s)",
            self.name, len(clusters), ctx.metadata["context_aware_n_cell_types"], has_substate,
        )

    def _run_context_validation(
        self,
        ctx: PipelineContext,
        adata,
        clusters: list[str],
        cluster_celltype: dict[str, str],
        cluster_scores: dict[str, float],
        X,
        is_sparse: bool,
        gene_idx: dict[str, int],
        ct_genes: dict[str, list[str]],
        metrics_dir: Path,
    ) -> None:
        """Compute per-cluster context validation scores and emit context_validation.json.

        Warn-only on mismatch (AC-4): never raises, never hard-errors.
        """
        import numpy as np
        declared_ct = getattr(ctx.cfg, "condition", "NSCLC")
        threshold = getattr(ctx.cfg, "context_mismatch_threshold", _DEFAULT_MISMATCH_THRESHOLD)
        min_cells = getattr(ctx.cfg, "context_min_cells", _DEFAULT_MIN_CELLS)

        cluster_labels = adata.obs["leiden"].astype(str)
        validation: dict[str, dict] = {}
        mismatches: list[str] = []

        for cluster in clusters:
            mask = cluster_labels.values == cluster
            cell_count = int(mask.sum())
            if cell_count < min_cells:
                continue

            assigned_ct = cluster_celltype[cluster]
            score = cluster_scores[cluster]

            # Compute top expressed cell type via raw score (sparse-aware per-cluster mean)
            if is_sparse:
                cluster_X = X[mask]
                cluster_mean = np.asarray(cluster_X.mean(axis=0)).flatten().astype(np.float32)
            else:
                cluster_X = np.asarray(X[mask], dtype=np.float32)
                cluster_mean = cluster_X.mean(axis=0)
            raw_scores = {}
            for ct, genes in ct_genes.items():
                idxs = [gene_idx[g] for g in genes if g in gene_idx]
                if not idxs:
                    continue
                raw_scores[ct] = float(cluster_mean[idxs].mean())

            observed_top = max(raw_scores, key=lambda k: raw_scores[k]) if raw_scores else "Unknown"
            normalized_score = score / (max(raw_scores.values()) + 1e-9) if raw_scores else 0.0

            validation[cluster] = {
                "declared_condition": declared_ct,
                "assigned_cell_type": assigned_ct,
                "observed_top": observed_top,
                "score": round(float(normalized_score), 4),
                "n_cells": cell_count,
            }

            if normalized_score < threshold:
                mismatches.append(cluster)
                import sys
                logger.warning(
                    "context mismatch: cluster %s assigned as '%s' but score %.3f < threshold %.3f",
                    cluster, assigned_ct, normalized_score, threshold,
                )

        (metrics_dir / "context_validation.json").write_text(
            json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        ctx.metadata["context_validation_mismatches"] = len(mismatches)
        ctx.metadata["context_validation_threshold"] = threshold
