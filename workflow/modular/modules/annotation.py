from __future__ import annotations

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

logger = logging.getLogger(__name__)


class AnnotationModule:
    """Optional module: marker-scoring based cell-type annotation with confidence scores."""

    name = "annotation"
    requires_keys = {"obs": ["leiden"]}
    provides_keys = {"obs": ["cell_type"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "leiden" not in adata.obs:
            raise ValueError("Annotation requires clustered AnnData.")

        marker_map = ctx.cfg.markers or DEFAULT_MARKERS
        available = {
            cell_type: [gene for gene in genes if gene in adata.var_names]
            for cell_type, genes in marker_map.items()
        }
        available = {cell_type: genes for cell_type, genes in available.items() if genes}
        if not available:
            raise ValueError("No valid marker genes found in dataset.")

        for cell_type, genes in available.items():
            sc.tl.score_genes(adata, genes, score_name=f"score_{cell_type}", use_raw=False)

        score_cols = [f"score_{cell_type}" for cell_type in available]
        score_mat = adata.obs[score_cols].copy()
        score_mat.columns = list(available.keys())

        # Assign cell type with confidence scoring
        max_scores = score_mat.max(axis=1)
        vals = score_mat.values
        if vals.shape[1] >= 2:
            second_scores = pd.Series(
                np.partition(vals, -2, axis=1)[:, -2], index=score_mat.index,
            )
        else:
            second_scores = pd.Series(np.zeros(vals.shape[0], dtype=float), index=score_mat.index)
        confidence = max_scores - second_scores

        adata.obs["cell_type"] = score_mat.idxmax(axis=1).values
        adata.obs["annotation_confidence"] = confidence.values
        # Mark low-confidence assignments as Unknown
        low_conf = max_scores < ctx.cfg.annotation_confidence_threshold
        adata.obs.loc[low_conf, "cell_type"] = "Unknown"
        adata.obs["cell_type_marker"] = adata.obs["cell_type"].astype(str)

        self._try_reference_mapping(adata, ctx)

        ctx.metadata["annotation_unknown_pct"] = round(
            float(low_conf.sum()) / len(low_conf) * 100, 2
        )

        base_cols = ["leiden", "cell_type", "annotation_confidence", "cell_type_marker"] + score_cols
        extra_cols = [c for c in ["reference_cell_type", "reference_confidence"] if c in adata.obs.columns]
        adata.obs[base_cols + extra_cols].to_csv(ctx.table_dir / "cell_type_annotation.csv")

        cluster_summary = (
            adata.obs.groupby("leiden", observed=True)["cell_type"]
            .agg(lambda x: x.value_counts().idxmax())
        ).rename("majority_cell_type")
        cluster_summary.to_csv(ctx.table_dir / "cluster_majority_cell_type.csv")

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
