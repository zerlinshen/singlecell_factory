"""TF → gene regulatory network from chromVAR + peak-to-gene linkage (Wave 5.1 / MV3).

# STATUS: promote-on-reuse

Builds an edge list of TF → target_gene relationships by combining:
- per-cell TF motif accessibility scores (adata.obsm['chromvar_scores'] from
  the chromvar module)
- peak-to-gene linkages (adata.uns['peak_to_gene_linkages'] from the
  peak_to_gene module)
- per-cell gene expression (adata.X)

Edge weight = Pearson correlation between TF activity per cell and target
gene expression per cell, restricted to (TF, gene) pairs where at least one
peak in the gene's linkage set carries the TF motif annotation. p-value via
two-sided Fisher transform; BH-FDR per TF.

Inputs (via cfg):
  --tf-network-top-edges INT      Max edges to retain by |weight| (default 5000)
  --tf-network-fdr-alpha FLOAT    BH-FDR threshold per TF (default 0.05)
  --tf-network-min-cells INT      Skip TF-gene pairs with fewer expressing cells (default 20)

Outputs:
  <run-dir>/tf_network/edges.csv      tf, target_gene, weight, p_value, fdr, n_supporting_peaks
  <run-dir>/tf_network/summary.json   n_tfs, n_genes, n_edges
  adata.uns['tf_network_edges']       full edge DataFrame
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse
from scipy import stats

from ..context import PipelineContext
from .._contract_violation import ModuleContractError

logger = logging.getLogger(__name__)


__references__ = {
    "Aibar_SCENIC_2017": {
        "title": "SCENIC: single-cell regulatory network inference and clustering",
        "authors": "Aibar S, González-Blas CB, Moerman T, et al.",
        "journal": "Nature Methods",
        "year": "2017",
        "doi": "10.1038/nmeth.4463",
        "description": "TF → gene network inference framework; this module implements a chromVAR-anchored variant that requires both motif accessibility and peak-to-gene support.",
    },
    "Trevino_2021": {
        "title": "Chromatin and gene-regulatory dynamics of the developing human neocortex at single-cell resolution",
        "authors": "Trevino AE, Müller F, Andersen J, et al.",
        "journal": "Cell",
        "year": "2021",
        "doi": "10.1016/j.cell.2021.07.039",
        "description": "Source paper for F-6 TF-gene network reproduction (Wave-5.6 MV3).",
    },
}


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    n = len(pvalues)
    if n == 0:
        return np.array([])
    order = np.argsort(pvalues)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    raw_q = pvalues * n / ranks
    q = np.minimum.accumulate(raw_q[order][::-1])[::-1]
    result = np.empty(n)
    result[order] = q
    return np.minimum(result, 1.0)


def _pearson_dense_vs_sparse_block(
    tf_scores: np.ndarray,
    X_block: scipy.sparse.spmatrix,
) -> np.ndarray:
    """Pearson r between a 1D dense vector (n_cells,) and each column of
    a sparse matrix (n_cells, k). Cell axis stays sparse on the gene side
    (the gene block is column-subset of adata.X, not a full densification).
    Returns (k,) array of r values.
    """
    n = tf_scores.shape[0]
    if X_block.shape[0] != n:
        raise ValueError("row count mismatch tf_scores vs X_block")
    tf_mean = float(tf_scores.mean())
    tf_centered = tf_scores - tf_mean
    tf_std = float(tf_scores.std())
    if tf_std <= 0:
        return np.zeros(X_block.shape[1])
    sum_x = np.asarray(X_block.sum(axis=0)).ravel()
    mean_x = sum_x / n
    sum_x2 = np.asarray(X_block.power(2).sum(axis=0)).ravel()
    var_x = sum_x2 / n - mean_x ** 2
    std_x = np.sqrt(np.maximum(var_x, 0.0))
    # cross = tf_centered^T . X_block  → (1, k) dense
    cross_dense = tf_centered.reshape(1, n) @ X_block
    cross = np.asarray(cross_dense).ravel() / n
    denom = tf_std * std_x
    r = np.zeros_like(mean_x)
    valid = denom > 0
    r[valid] = cross[valid] / denom[valid]
    return np.clip(r, -1.0, 1.0)


def _fisher_pvalue(r: np.ndarray, n: int) -> np.ndarray:
    """Two-sided p-value from Pearson r via Fisher z-transform."""
    # Avoid r == ±1 exact
    r_safe = np.clip(r, -0.999999, 0.999999)
    z = 0.5 * np.log((1 + r_safe) / (1 - r_safe))
    se = 1.0 / np.sqrt(max(n - 3, 1))
    p = 2.0 * (1.0 - stats.norm.cdf(np.abs(z) / se))
    return p


class TFNetworkModule:
    """Build TF → target gene network from chromVAR + peak-to-gene linkages."""

    name = "tf_network"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {
        "obsm": ["chromvar_scores"],
        "uns": ["chromvar_tfs", "peak_to_gene_linkages"],
    }
    provides_keys: dict[str, list[str]] = {"uns": ["tf_network_edges"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")
        adata = ctx.adata
        cfg = ctx.cfg

        if "chromvar_scores" not in adata.obsm:
            ctx.status(self.name, "skipped", "chromvar_scores absent — run chromvar first")
            ctx.metadata["tf_network_status"] = "skipped_no_chromvar"
            return
        tfs = adata.uns.get("chromvar_tfs")
        if tfs is None or len(tfs) == 0:
            ctx.status(self.name, "skipped", "chromvar_tfs empty")
            ctx.metadata["tf_network_status"] = "skipped_no_tfs"
            return

        linkages = adata.uns.get("peak_to_gene_linkages")
        if linkages is None or (isinstance(linkages, pd.DataFrame) and linkages.empty):
            ctx.status(self.name, "skipped",
                       "peak_to_gene_linkages empty — run peak_to_gene module first")
            ctx.metadata["tf_network_status"] = "skipped_no_linkages"
            return
        if not isinstance(linkages, pd.DataFrame):
            linkages = pd.DataFrame(linkages)

        tf_peak_ann_path = getattr(cfg, "tf_peak_annotation_path", None)
        if not tf_peak_ann_path or not Path(tf_peak_ann_path).exists():
            ctx.status(self.name, "skipped",
                       "tf_peak_annotation_path required for TF→peak filter")
            ctx.metadata["tf_network_status"] = "skipped_no_annotation"
            return

        ann = pd.read_csv(tf_peak_ann_path, sep="\t", dtype={"tf_name": str, "peak_id": str})

        # Build peak_id → peak_idx
        atac_var = adata.uns.get("atac_var")
        if atac_var is None:
            raise ModuleContractError("tf_network: adata.uns['atac_var'] missing.")
        if "peak_id" in atac_var.columns:
            atac_var = atac_var.reset_index(drop=True)
            peak_id_to_idx = {pid: i for i, pid in enumerate(atac_var["peak_id"].astype(str))}
        else:
            peak_id_to_idx = {
                f"{c}:{s}-{e}": i for i, (c, s, e) in enumerate(zip(
                    atac_var["seqnames"].astype(str),
                    atac_var["start"].astype(int),
                    atac_var["end"].astype(int),
                ))
            }
        ann["peak_idx"] = ann["peak_id"].map(peak_id_to_idx)
        ann = ann.dropna(subset=["peak_idx"])
        ann["peak_idx"] = ann["peak_idx"].astype(int)
        tf_to_peaks: dict[str, set[int]] = {
            str(tf): set(grp["peak_idx"].tolist())
            for tf, grp in ann.groupby("tf_name")
        }

        # linkages: peak_idx + gene; build gene → set(peak_idx)
        if "peak_idx" not in linkages.columns or "gene" not in linkages.columns:
            raise ModuleContractError(
                "tf_network: peak_to_gene_linkages missing 'peak_idx' or 'gene' column"
            )
        gene_to_peaks: dict[str, set[int]] = {}
        for _, row in linkages.iterrows():
            gene_to_peaks.setdefault(str(row["gene"]), set()).add(int(row["peak_idx"]))

        # Build candidate (tf, gene) pairs with supporting peak count
        candidate: list[tuple[str, str, int]] = []
        for tf, tf_peaks in tf_to_peaks.items():
            if tf not in tfs:
                continue
            for gene, gene_peaks in gene_to_peaks.items():
                overlap = tf_peaks & gene_peaks
                if overlap:
                    candidate.append((tf, gene, len(overlap)))
        if not candidate:
            ctx.status(self.name, "skipped",
                       "no TF-gene pairs with shared supporting peaks")
            ctx.metadata["tf_network_status"] = "skipped_no_candidate"
            adata.uns["tf_network_edges"] = pd.DataFrame(
                columns=["tf", "target_gene", "weight", "p_value", "fdr", "n_supporting_peaks"]
            )
            return

        # Pre-index TF scores by name
        tf_to_col = {tf: i for i, tf in enumerate(tfs)}
        Z = adata.obsm["chromvar_scores"]
        var_index = list(adata.var.index)
        gene_to_var = {g: i for i, g in enumerate(var_index)}

        # Group candidates by TF to amortize per-TF score extraction + per-gene block
        records: list[dict] = []
        X = adata.X
        if not scipy.sparse.issparse(X):
            raise ModuleContractError("tf_network: adata.X must be sparse.")
        X_csc = X.tocsc()
        n_cells = X.shape[0]
        fdr_alpha = float(getattr(cfg, "tf_network_fdr_alpha", 0.05))
        min_cells = int(getattr(cfg, "tf_network_min_cells", 20))

        by_tf: dict[str, list[tuple[str, int]]] = {}
        for tf, gene, n_sup in candidate:
            by_tf.setdefault(tf, []).append((gene, n_sup))

        for tf, gene_list in by_tf.items():
            tf_scores = np.asarray(Z[:, tf_to_col[tf]], dtype=np.float64)
            gene_idx_list = []
            gene_name_list = []
            sup_list = []
            for gene, n_sup in gene_list:
                gi = gene_to_var.get(gene)
                if gi is None:
                    continue
                gene_idx_list.append(gi)
                gene_name_list.append(gene)
                sup_list.append(n_sup)
            if not gene_idx_list:
                continue
            X_block = X_csc[:, gene_idx_list]
            # min-cells filter
            nonzero_per_gene = np.asarray((X_block != 0).sum(axis=0)).ravel()
            keep = nonzero_per_gene >= min_cells
            if not keep.any():
                continue
            X_block_keep = X_block[:, keep]
            kept_names = [gene_name_list[i] for i, k in enumerate(keep) if k]
            kept_sup = [sup_list[i] for i, k in enumerate(keep) if k]
            r = _pearson_dense_vs_sparse_block(tf_scores, X_block_keep)
            p = _fisher_pvalue(r, n_cells)
            q = _bh_fdr(p)
            for gene, weight, pval, qval, n_sup in zip(kept_names, r, p, q, kept_sup):
                if qval < fdr_alpha:
                    records.append({
                        "tf": tf,
                        "target_gene": gene,
                        "weight": float(weight),
                        "p_value": float(pval),
                        "fdr": float(qval),
                        "n_supporting_peaks": int(n_sup),
                    })

        edges_df = pd.DataFrame(records) if records else pd.DataFrame(
            columns=["tf", "target_gene", "weight", "p_value", "fdr", "n_supporting_peaks"]
        )
        if not edges_df.empty:
            edges_df = edges_df.reindex(
                edges_df["weight"].abs().sort_values(ascending=False).index
            )
            top = int(getattr(cfg, "tf_network_top_edges", 5000))
            edges_df = edges_df.head(top).reset_index(drop=True)

        adata.uns["tf_network_edges"] = edges_df

        out_dir = ctx.run_dir / "tf_network"
        out_dir.mkdir(parents=True, exist_ok=True)
        edges_df.to_csv(out_dir / "edges.csv", index=False)
        (out_dir / "summary.json").write_text(
            json.dumps({
                "n_edges": int(len(edges_df)),
                "n_unique_tfs": int(edges_df["tf"].nunique()) if not edges_df.empty else 0,
                "n_unique_genes": int(edges_df["target_gene"].nunique()) if not edges_df.empty else 0,
                "fdr_alpha": fdr_alpha,
                "min_cells_per_gene": min_cells,
            }, indent=2),
            encoding="utf-8",
        )

        ctx.metadata["tf_network_status"] = "ok"
        ctx.metadata["tf_network_n_edges"] = int(len(edges_df))
        logger.info("tf_network: %d edges at FDR<%.2f", len(edges_df), fdr_alpha)
