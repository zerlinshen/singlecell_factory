from __future__ import annotations

import json
import logging
import matplotlib
from pathlib import Path

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext


__references__ = {
    "Tirosh_marker_scoring_2016": {
        "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
        "authors": "Tirosh et al.",
        "journal": "Science",
        "year": "2016",
        "doi": "10.1126/science.aad0501",
        "description": "Cluster-vote marker-mean scoring approach for cell type assignment.",
    },
    "Stuart_label_transfer_2019": {
        "title": "Comprehensive Integration of Single-Cell Data",
        "authors": "Stuart et al.",
        "journal": "Cell",
        "year": "2019",
        "doi": "10.1016/j.cell.2019.05.031",
        "description": "Reference-based label transfer principles; this module uses a simpler sklearn.NearestNeighbors KNN majority vote on a labeled reference.",
    },
}



# Curated marker catalogue used ONLY when neither ctx.cfg.markers nor the
# marker_db_loader provided a marker map. The authoritative path for any
# scientific run is to populate `adata.uns["marker_db_index"]` via
# `marker_db_loader` (which pins CellMarker 2.0 / PanglaoDB / CellTypist /
# scTypeDB version strings into the run manifest). This dict is a tumor /
# immune / stromal starter set assembled by hand from common single-cell
# atlas conventions; it is NOT a versioned DB.
#
# When tuning markers for a new tissue or condition, prefer (a) configuring
# `marker_db_loader` for that (tissue, condition) pair or (b) passing
# `--markers-json` to the CLI rather than editing this dict.
DEFAULT_MARKERS = {
    "Tumor epithelial": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1"],
    "T cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "IL7R"],
    "NK cell": ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "B cell": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Myeloid/Macro": ["LYZ", "CD68", "CST3", "FCGR3A", "AIF1"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL3A1", "FAP"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "Plasma cell": ["IGHG1", "IGKC", "MZB1", "XBP1"],
    "Mast cell": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "Dendritic cell": ["CD1C", "CLEC9A", "FCER1A", "IRF8"],
}

EPITHELIAL_QC_MARKERS = ("EPCAM", "KRT8", "KRT18")

logger = logging.getLogger(__name__)


class AnnotationModule:
    """Optional module: marker-scoring based cell-type annotation with confidence scores."""

    name = "annotation"
    requires_keys = {"obs": ["leiden"]}
    provides_keys = {"obs": ["cell_type"]}

    @staticmethod
    def _resolve_leiden_key(adata) -> str:
        """Pick the cluster column to bind cell_type against.

        Prefer ``leiden_corrected`` (the post-batch-correction snapshot written
        by batch_correction) so cell_type is a deterministic function of the
        FINAL post-correction clustering regardless of execution mode. Fall
        back to ``leiden`` when batch_correction did not run.
        """
        if "leiden_corrected" in adata.obs:
            return "leiden_corrected"
        return "leiden"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "leiden" not in adata.obs:
            raise ValueError("Annotation requires clustered AnnData.")
        leiden_key = self._resolve_leiden_key(adata)
        ctx.metadata["annotation_leiden_source"] = leiden_key

        marker_map = ctx.cfg.markers or DEFAULT_MARKERS
        available = {
            cell_type: [gene for gene in genes if gene in adata.var_names]
            for cell_type, genes in marker_map.items()
        }
        available = {cell_type: genes for cell_type, genes in available.items() if genes}
        if not available:
            raise ValueError("No valid marker genes found in dataset.")

        strategy = getattr(ctx.cfg, "annotation_strategy", "cluster_voting")
        ctx.metadata["annotation_strategy"] = strategy

        # Score genes on full adata (per cell) — populates score_* obs columns and
        # provides per-cell confidence used by reference mapping conservative gate.
        for cell_type, genes in available.items():
            sc.tl.score_genes(adata, genes, score_name=f"score_{cell_type}", use_raw=False)

        score_cols = [f"score_{cell_type}" for cell_type in available]
        score_mat = adata.obs[score_cols].copy()
        score_mat.columns = list(available.keys())

        if strategy == "cluster_voting":
            # Cluster-level assignment uses raw mean expression per marker gene set,
            # not aggregated score_genes values. score_genes background subtraction is
            # unreliable when absolute marker expression is low (background genes in the
            # same expression bin dominate the signal). Raw mean expression correctly
            # captures relative marker enrichment per cluster.
            cluster_score_matrix = self._cluster_mean_expression(adata, available, leiden_key)
            cluster_label_map = cluster_score_matrix.idxmax(axis=1).to_dict()
            adata.obs["cell_type"] = adata.obs[leiden_key].astype(str).map(cluster_label_map).astype(str)
            # Confidence: per-cell score_genes margin (max - second), for reference mapping gate
            max_scores = score_mat.max(axis=1)
            vals = score_mat.values
            if vals.shape[1] >= 2:
                second_scores = pd.Series(
                    np.partition(vals, -2, axis=1)[:, -2], index=score_mat.index,
                )
            else:
                second_scores = pd.Series(np.zeros(vals.shape[0], dtype=float), index=score_mat.index)
            adata.obs["annotation_confidence"] = (max_scores - second_scores).values
            # Unknown gate: cluster-level raw mean max score
            cluster_max_score = cluster_score_matrix.max(axis=1)
            low_conf_clusters = set(
                cluster_max_score[cluster_max_score < ctx.cfg.annotation_confidence_threshold].index.astype(str)
            )
            low_conf = adata.obs[leiden_key].astype(str).isin(low_conf_clusters)
            adata.obs.loc[low_conf, "cell_type"] = "Unknown"
            # Write audit trail
            cluster_score_matrix.to_csv(ctx.table_dir / "cluster_score_matrix.csv")
        else:
            # cell_argmax: original per-cell path
            max_scores = score_mat.max(axis=1)
            vals = score_mat.values
            if vals.shape[1] >= 2:
                second_scores = pd.Series(
                    np.partition(vals, -2, axis=1)[:, -2], index=score_mat.index,
                )
            else:
                second_scores = pd.Series(np.zeros(vals.shape[0], dtype=float), index=score_mat.index)
            adata.obs["cell_type"] = score_mat.idxmax(axis=1).values
            adata.obs["annotation_confidence"] = (max_scores - second_scores).values
            low_conf = max_scores < ctx.cfg.annotation_confidence_threshold
            adata.obs.loc[low_conf, "cell_type"] = "Unknown"
            # Write cluster_score_matrix.csv for consistency
            cluster_score_matrix = score_mat.groupby(
                adata.obs[leiden_key].astype(str), observed=True
            ).mean()
            cluster_score_matrix.index.name = leiden_key
            cluster_score_matrix.to_csv(ctx.table_dir / "cluster_score_matrix.csv")

        adata.obs["cell_type_marker"] = adata.obs["cell_type"].astype(str)

        self._try_reference_mapping(adata, ctx)

        low_conf_mask = adata.obs["cell_type"] == "Unknown"
        ctx.metadata["annotation_unknown_pct"] = round(
            float(low_conf_mask.sum()) / max(len(low_conf_mask), 1) * 100, 2
        )

        score_cols = [f"score_{ct}" for ct in available]
        base_cols = ["leiden", "cell_type", "annotation_confidence", "cell_type_marker"] + score_cols
        extra_cols = [c for c in ["reference_cell_type", "reference_confidence"] if c in adata.obs.columns]
        adata.obs[base_cols + extra_cols].to_csv(ctx.table_dir / "cell_type_annotation.csv")

        cluster_summary = (
            adata.obs.groupby("leiden", observed=True)["cell_type"]
            .agg(lambda x: x.value_counts().idxmax())
        ).rename("majority_cell_type")
        cluster_summary.to_csv(ctx.table_dir / "cluster_majority_cell_type.csv")
        self._write_epithelial_marker_qc(adata, ctx)

        # --- Visualizations ---
        sc.pl.umap(adata, color=["cell_type"], show=False, legend_loc="on data")
        plt.savefig(ctx.figure_dir / "umap_cell_type.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Confidence UMAP
        sc.pl.umap(adata, color=["annotation_confidence"], show=False, cmap="viridis")
        plt.savefig(ctx.figure_dir / "umap_annotation_confidence.png", dpi=160, bbox_inches="tight")
        plt.close()

        if "reference_cell_type" in adata.obs:
            sc.pl.umap(adata, color=["reference_cell_type"], show=False, legend_loc="on data")
            plt.savefig(ctx.figure_dir / "umap_reference_cell_type.png", dpi=160, bbox_inches="tight")
            plt.close()
        if "reference_confidence" in adata.obs:
            sc.pl.umap(adata, color=["reference_confidence"], show=False, cmap="magma")
            plt.savefig(ctx.figure_dir / "umap_reference_confidence.png", dpi=160, bbox_inches="tight")
            plt.close()

        # Stacked bar: cell type composition per cluster
        self._plot_composition(adata, ctx)

    @staticmethod
    def _cluster_mean_expression(
        adata, available: dict[str, list[str]], leiden_key: str = "leiden"
    ) -> pd.DataFrame:
        """Compute mean expression of each marker gene set per leiden cluster.

        Returns DataFrame with rows=cluster, cols=cell_type, values=mean log-expr of gene set.
        Uses raw matrix slicing to avoid score_genes background subtraction artefacts.
        """
        X = adata.X
        leiden_vals = adata.obs[leiden_key].astype(str).values
        clusters = sorted(set(leiden_vals), key=lambda c: (not c.isdigit(), int(c) if c.isdigit() else c))
        gene_index = {g: i for i, g in enumerate(adata.var_names)}

        rows: dict[str, dict[str, float]] = {}
        for cluster in clusters:
            mask = leiden_vals == cluster
            x_clust = X[mask]
            if sparse.issparse(x_clust):
                mean_expr = np.asarray(x_clust.mean(axis=0)).flatten()
            else:
                mean_expr = np.asarray(x_clust, dtype=np.float64).mean(axis=0)
            row: dict[str, float] = {}
            for cell_type, genes in available.items():
                idxs = [gene_index[g] for g in genes if g in gene_index]
                row[cell_type] = float(mean_expr[idxs].mean()) if idxs else 0.0
            rows[cluster] = row

        df = pd.DataFrame(rows).T
        df.index.name = leiden_key
        return df

    @staticmethod
    def _try_reference_mapping(adata, ctx: PipelineContext) -> None:
        """Optional KNN label transfer from a reference h5ad."""
        ref_path = ctx.cfg.reference_adata
        if ref_path is None:
            ctx.metadata["reference_mapping_status"] = "skipped_no_reference"
            return

        ref_path = Path(ref_path)
        if not ref_path.exists():
            ctx.metadata["reference_mapping_status"] = "skipped_missing_reference_file"
            ctx.metadata["reference_mapping_skip_reason"] = str(ref_path)
            return

        try:
            ref = ad.read_h5ad(ref_path)
        except Exception as exc:
            ctx.metadata["reference_mapping_status"] = "skipped_reference_read_error"
            ctx.metadata["reference_mapping_skip_reason"] = str(exc)
            return

        label_key = ctx.cfg.reference_label_key
        if label_key not in ref.obs.columns:
            ctx.metadata["reference_mapping_status"] = "skipped_missing_reference_label_key"
            ctx.metadata["reference_mapping_skip_reason"] = label_key
            return

        shared = adata.var_names.intersection(ref.var_names)
        if len(shared) < 50:
            ctx.metadata["reference_mapping_status"] = "skipped_too_few_shared_genes"
            ctx.metadata["reference_mapping_shared_genes"] = int(len(shared))
            return

        query_x = adata[:, shared].X
        ref_x = ref[:, shared].X
        query_x = query_x if sparse.issparse(query_x) else np.asarray(query_x, dtype=np.float32)
        ref_x = ref_x if sparse.issparse(ref_x) else np.asarray(ref_x, dtype=np.float32)
        query_x = normalize(query_x, norm="l2", axis=1, copy=False)
        ref_x = normalize(ref_x, norm="l2", axis=1, copy=False)

        k = min(int(ctx.cfg.reference_k), int(ref.n_obs))
        if k < 1:
            ctx.metadata["reference_mapping_status"] = "skipped_invalid_k"
            return

        nn = NearestNeighbors(n_neighbors=k, metric="cosine", algorithm="brute")
        nn.fit(ref_x)
        _, idx = nn.kneighbors(query_x)
        ref_labels = ref.obs[label_key].astype(str).to_numpy()

        pred_labels: list[str] = []
        pred_conf: list[float] = []
        for neighbors in idx:
            votes = ref_labels[neighbors]
            values, counts = np.unique(votes, return_counts=True)
            top = int(np.argmax(counts))
            pred_labels.append(str(values[top]))
            pred_conf.append(float(counts[top]) / float(k))

        adata.obs["reference_cell_type"] = pd.Series(pred_labels, index=adata.obs_names, dtype="object")
        adata.obs["reference_confidence"] = pd.Series(pred_conf, index=adata.obs_names, dtype=float)

        mode = str(ctx.cfg.reference_override_mode).strip().lower()
        base_mask = adata.obs["reference_confidence"] >= float(ctx.cfg.reference_min_confidence)
        if mode == "all":
            override_mask = base_mask
        else:
            marker_unknown = adata.obs["cell_type_marker"].astype(str) == "Unknown"
            marker_low = adata.obs["annotation_confidence"] < float(ctx.cfg.annotation_confidence_threshold)
            override_mask = base_mask & (marker_unknown | marker_low)
        adata.obs.loc[override_mask, "cell_type"] = adata.obs.loc[override_mask, "reference_cell_type"].values

        applied = int(override_mask.sum())
        total = int(adata.n_obs)
        ctx.metadata["reference_mapping_status"] = "completed"
        ctx.metadata["reference_mapping_source"] = str(ref_path)
        ctx.metadata["reference_mapping_label_key"] = label_key
        ctx.metadata["reference_mapping_shared_genes"] = int(len(shared))
        ctx.metadata["reference_mapping_k"] = int(k)
        ctx.metadata["reference_mapping_override_mode"] = mode if mode in {"all", "conservative"} else "conservative"
        ctx.metadata["reference_mapping_overridden_cells"] = applied
        ctx.metadata["reference_mapping_overridden_pct"] = round(100.0 * applied / max(total, 1), 2)
        logger.info(
            "Reference mapping applied: %d/%d cells overridden (mode=%s, k=%d, min_conf=%.3f).",
            applied,
            total,
            ctx.metadata["reference_mapping_override_mode"],
            k,
            float(ctx.cfg.reference_min_confidence),
        )

    @staticmethod
    def _plot_composition(adata, ctx: PipelineContext) -> None:
        ct_counts = adata.obs.groupby(["leiden", "cell_type"], observed=True).size().unstack(fill_value=0)
        ct_frac = ct_counts.div(ct_counts.sum(axis=1), axis=0)
        ct_frac.plot(kind="bar", stacked=True, figsize=(12, 5), colormap="tab20")
        plt.ylabel("Fraction")
        plt.xlabel("Leiden cluster")
        plt.title("Cell type composition per cluster")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "cell_type_composition.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _marker_expression_frame(adata, markers: list[str]) -> pd.DataFrame:
        expr = adata[:, markers].X
        if sparse.issparse(expr):
            # densify-allowed: marker panel is O(n_cells x <=3 genes)
            expr = expr.toarray()
        else:
            expr = np.asarray(expr)
        return pd.DataFrame(expr, index=adata.obs_names, columns=markers)

    @staticmethod
    def _summarize_marker_distribution(expr: pd.DataFrame, labels, group_name: str) -> pd.DataFrame:
        work = expr.copy()
        work[group_name] = pd.Series(labels, index=expr.index).astype(str)
        rows = []
        for group, sub in work.groupby(group_name, observed=True):
            marker_values = sub.drop(columns=[group_name])
            for marker in marker_values.columns:
                values = marker_values[marker].to_numpy(dtype=float)
                rows.append(
                    {
                        group_name: str(group),
                        "marker": marker,
                        "n_cells": int(values.shape[0]),
                        "mean": float(np.mean(values)),
                        "median": float(np.median(values)),
                        "pct_positive": round(float(np.mean(values > 0.0)) * 100.0, 3),
                    }
                )
        return pd.DataFrame(rows)

    @classmethod
    def _write_epithelial_marker_qc(cls, adata, ctx: PipelineContext) -> None:
        present = [marker for marker in EPITHELIAL_QC_MARKERS if marker in adata.var_names]
        ctx.metadata["epithelial_marker_qc_present_markers"] = present
        if not present:
            ctx.metadata["epithelial_marker_qc_status"] = "skipped_no_epithelial_markers"
            return

        expr = cls._marker_expression_frame(adata, present)
        if "leiden" in adata.obs:
            cls._summarize_marker_distribution(
                expr,
                adata.obs["leiden"],
                "leiden",
            ).to_csv(ctx.table_dir / "epithelial_marker_summary_by_leiden.csv", index=False)
        if "cell_type" in adata.obs:
            cls._summarize_marker_distribution(
                expr,
                adata.obs["cell_type"],
                "cell_type",
            ).to_csv(ctx.table_dir / "epithelial_marker_summary_by_cell_type.csv", index=False)

        epithelial_mask = adata.obs["cell_type"].astype(str).str.contains(
            "epithelial",
            case=False,
            regex=False,
        ) if "cell_type" in adata.obs else pd.Series(False, index=adata.obs_names)
        marker_score = expr.mean(axis=1)
        qc_payload = {
            "present_markers": present,
            "epithelial_cells": int(epithelial_mask.sum()),
            "non_epithelial_cells": int((~epithelial_mask).sum()),
            "epithelial_marker_score_median": (
                float(marker_score[epithelial_mask].median()) if epithelial_mask.any() else None
            ),
            "non_epithelial_marker_score_median": (
                float(marker_score[~epithelial_mask].median()) if (~epithelial_mask).any() else None
            ),
            "status": "completed",
        }
        ctx.metadata["epithelial_marker_qc_status"] = "completed"
        ctx.metadata["epithelial_marker_score_median"] = qc_payload[
            "epithelial_marker_score_median"
        ]
        (ctx.table_dir / "epithelial_marker_qc.json").write_text(
            json.dumps(qc_payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        if "X_umap" in adata.obsm:
            sc.pl.umap(adata, color=present, show=False, cmap="viridis")
            plt.savefig(
                ctx.figure_dir / "umap_epithelial_markers.png",
                dpi=160,
                bbox_inches="tight",
            )
            plt.close()
