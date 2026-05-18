"""Ribo-seq (ribosome profiling) ingest (Wave 2B / P2.S21).

Memory-safe pandas-based ingest. Reads a per-gene ribosome footprint table
and pairs it with RNA expression to compute translation efficiency (TE) per
gene per sample.

Inputs (via cfg):
  --ribo-footprints-path PATH       TSV with columns: gene, sample, footprint_count
  --ribo-rna-counts-path PATH       (optional) TSV with columns: gene, sample, rna_count
                                     when omitted, TE is computed from the RNA assay in adata
  --ribo-min-footprint-count INT    minimum footprint count to keep a gene (default 10)

Outputs:
  adata.uns["ribo_translation_efficiency"]   DataFrame: gene, sample, footprint_count,
                                              rna_count, te (footprint / rna with pseudocount)
  adata.uns["ribo_ingest_metadata"]          { n_genes, n_samples, te_mean, te_median }
  runs/<run-id>/ribo_ingest/ribo_summary.json
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "Ingolia_RiboSeq_2009": {
        "title": "Genome-Wide Analysis in Vivo of Translation with Nucleotide Resolution Using Ribosome Profiling",
        "authors": "Ingolia et al.",
        "journal": "Science",
        "year": "2009",
        "doi": "10.1126/science.1168978",
        "description": "Original ribosome profiling (Ribo-seq) methodology defining the footprint count semantics consumed here.",
    },
    "Brar_Weissman_RiboSeq_review_2015": {
        "title": "Ribosome profiling reveals the what, when, where and how of protein synthesis",
        "authors": "Brar, Weissman",
        "journal": "Nature Reviews Molecular Cell Biology",
        "year": "2015",
        "doi": "10.1038/nrm3950",
        "description": "Translation efficiency (footprint / RNA) interpretation framework adopted here.",
    },
}


def _aggregate_rna_per_gene_sample(adata) -> pd.DataFrame:
    """Sum RNA counts per (gene, sample) from adata.X. Sparse-safe; no densification."""
    if "sample" not in adata.obs.columns:
        return pd.DataFrame(columns=["gene", "sample", "rna_count"])
    sample_col = adata.obs["sample"].astype(str)
    gene_names = adata.var_names.to_numpy()
    out_rows: list[dict] = []
    for sample_id, mask in sample_col.groupby(sample_col).groups.items():
        # mask is a list of obs index labels; convert to row indices
        row_idx = adata.obs_names.get_indexer(mask)
        X_sub = adata.X[row_idx]
        if sp.issparse(X_sub):
            gene_sums = np.asarray(X_sub.sum(axis=0)).flatten()
        else:
            gene_sums = np.asarray(X_sub).sum(axis=0)
        for g, c in zip(gene_names, gene_sums):
            if c > 0:
                out_rows.append({"gene": str(g), "sample": str(sample_id), "rna_count": float(c)})
    return pd.DataFrame(out_rows)


class RiboIngestModule:
    """Ribosome profiling ingest: footprint counts + translation efficiency."""

    name = "ribo_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "uns": ["ribo_translation_efficiency", "ribo_ingest_metadata"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        footprints_path = getattr(ctx.cfg, "ribo_footprints_path", None)
        if not footprints_path:
            ctx.status(self.name, "skipped", "ribo_footprints_path not set")
            ctx.metadata["ribo_ingest_status"] = "skipped_no_input"
            return

        footprints_path = Path(footprints_path)
        if not footprints_path.exists():
            ctx.status(self.name, "skipped", f"footprints missing: {footprints_path}")
            ctx.metadata["ribo_ingest_status"] = "skipped_path_missing"
            return

        min_fp = int(getattr(ctx.cfg, "ribo_min_footprint_count", 10))
        logger.info("%s: loading ribo footprints from %s", self.name, footprints_path)
        footprints = pd.read_csv(footprints_path, sep="\t")
        required = {"gene", "sample", "footprint_count"}
        if not required.issubset(footprints.columns):
            ctx.status(self.name, "skipped",
                       f"footprints schema missing required cols: {required - set(footprints.columns)}")
            ctx.metadata["ribo_ingest_status"] = "skipped_bad_schema"
            return

        footprints["gene"] = footprints["gene"].astype(str)
        footprints["sample"] = footprints["sample"].astype(str)
        footprints["footprint_count"] = pd.to_numeric(footprints["footprint_count"], errors="coerce").fillna(0)
        footprints = footprints[footprints["footprint_count"] >= min_fp]

        # RNA counts source
        rna_path = getattr(ctx.cfg, "ribo_rna_counts_path", None)
        if rna_path and Path(rna_path).exists():
            rna_df = pd.read_csv(rna_path, sep="\t")
            rna_df["gene"] = rna_df["gene"].astype(str)
            rna_df["sample"] = rna_df["sample"].astype(str)
            rna_df["rna_count"] = pd.to_numeric(rna_df["rna_count"], errors="coerce").fillna(0)
        else:
            logger.info("%s: deriving RNA counts from adata.X (sparse-safe per-sample sum)", self.name)
            rna_df = _aggregate_rna_per_gene_sample(adata)

        # Outer-join footprints + rna on (gene, sample); fill missing with 0
        te_df = footprints.merge(rna_df, on=["gene", "sample"], how="outer").fillna({"footprint_count": 0.0, "rna_count": 0.0})
        # Translation efficiency = footprint / (rna + 1) — pseudocount for stability
        te_df["te"] = te_df["footprint_count"] / (te_df["rna_count"] + 1.0)

        adata.uns["ribo_translation_efficiency"] = te_df
        meta = {
            "n_genes_with_footprints": int(footprints["gene"].nunique()),
            "n_samples": int(te_df["sample"].nunique()),
            "te_mean": float(te_df["te"].mean()) if not te_df.empty else 0.0,
            "te_median": float(te_df["te"].median()) if not te_df.empty else 0.0,
            "min_footprint_threshold": min_fp,
        }
        adata.uns["ribo_ingest_metadata"] = meta

        out_dir = ctx.run_dir / "ribo_ingest"
        out_dir.mkdir(parents=True, exist_ok=True)
        te_df.to_csv(out_dir / "translation_efficiency.csv", index=False)
        (out_dir / "ribo_summary.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

        ctx.metadata["ribo_ingest_status"] = "ok"
        ctx.metadata["ribo_n_genes_with_footprints"] = meta["n_genes_with_footprints"]
        logger.info("%s: %d genes x %d samples; TE mean=%.3f median=%.3f",
                    self.name, meta["n_genes_with_footprints"], meta["n_samples"],
                    meta["te_mean"], meta["te_median"])
