"""Spatial neighborhoods analysis module.

Optional analysis layer that runs on top of ``SpatialIngestModule``. When
``squidpy`` is available we compute:

- Moran's I (spatial autocorrelation) for the top-N highly variable genes.
- Neighborhood enrichment by ``cfg.group_by`` (typically a cluster or
  cell-type annotation).
- Co-occurrence probability across the same grouping.

When ``squidpy`` is NOT installed this module logs a clean message and exits
without raising. That is the documented contract for this lane: spatial
neighborhood analysis is an optional analytical add-on.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "palla_squidpy_2022": {
        "title": "Squidpy: a scalable framework for spatial omics analysis",
        "authors": "Palla et al.",
        "journal": "Nature Methods",
        "year": "2022",
        "doi": "10.1038/s41592-021-01358-2",
    },
}


@dataclass(frozen=True)
class SpatialNeighborhoodsConfig:
    """Local config for the spatial neighborhoods module."""

    group_by: str = "leiden"
    n_top_genes: int = 50
    n_neighs: int = 6
    coord_type: str = "generic"  # squidpy spatial_neighbors coord_type
    obsm_key: str = "spatial"


class SpatialNeighborhoodsModule:
    """Optional module: Moran's I + neighborhood enrichment + co-occurrence.

    Depends on ``SpatialIngestModule`` having populated
    ``adata.obsm["spatial"]``. Cleanly skips when squidpy is missing.
    """

    name = "spatial_neighborhoods"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {"obsm": ["spatial"]}
    provides_keys: dict[str, list[str]] = {
        "uns": ["spatial_morans_i", "spatial_nhood_enrichment", "spatial_co_occurrence"],
    }

    def __init__(self, config: SpatialNeighborhoodsConfig | None = None) -> None:
        self.config = config or SpatialNeighborhoodsConfig()

    # ------------------------------------------------------------------
    # Public entrypoint
    # ------------------------------------------------------------------

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        try:
            import squidpy as sq  # type: ignore
        except ImportError:
            logger.info("squidpy not available; spatial_neighborhoods skipped")
            ctx.metadata["spatial_neighborhoods_status"] = "skipped_no_squidpy"
            ctx.status(self.name, "skipped", "squidpy not installed")
            return

        cfg = self._merge_runtime_config(ctx)
        if cfg.obsm_key not in adata.obsm:
            logger.info(
                "%s: adata.obsm['%s'] missing; skipping spatial neighborhood analysis.",
                self.name, cfg.obsm_key,
            )
            ctx.metadata["spatial_neighborhoods_status"] = "skipped_no_coords"
            ctx.status(self.name, "skipped", f"obsm['{cfg.obsm_key}'] missing")
            return

        # Build the spatial neighbor graph. squidpy mutates adata.obsp in place.
        sq.gr.spatial_neighbors(
            adata, coord_type=cfg.coord_type, n_neighs=cfg.n_neighs,
            spatial_key=cfg.obsm_key,
        )

        # ---- Moran's I on top-N variable genes -------------------------------
        morans_path = None
        try:
            genes = self._select_morans_genes(adata, cfg.n_top_genes)
            if genes:
                sq.gr.spatial_autocorr(adata, mode="moran", genes=genes)
                morans_df = adata.uns.get("moranI")
                if isinstance(morans_df, pd.DataFrame):
                    morans_path = ctx.table_dir / "spatial_morans_i.parquet"
                    morans_df.to_parquet(morans_path)
                    adata.uns["spatial_morans_i"] = {
                        "n_genes": int(len(morans_df)),
                        "table_path": str(morans_path),
                    }
        except Exception as exc:
            logger.warning("%s: Moran's I failed: %s", self.name, exc)

        # ---- Neighborhood enrichment ----------------------------------------
        enrich_path = None
        if cfg.group_by in adata.obs.columns:
            try:
                sq.gr.nhood_enrichment(adata, cluster_key=cfg.group_by)
                enrich = adata.uns.get(f"{cfg.group_by}_nhood_enrichment")
                if isinstance(enrich, dict) and "zscore" in enrich:
                    z = np.asarray(enrich["zscore"])
                    enrich_path = ctx.table_dir / "spatial_nhood_enrichment_zscore.parquet"
                    pd.DataFrame(z).to_parquet(enrich_path)
                    adata.uns["spatial_nhood_enrichment"] = {
                        "group_by": cfg.group_by,
                        "shape": list(z.shape),
                        "table_path": str(enrich_path),
                    }
                    self._plot_heatmap(z, ctx, "spatial_nhood_enrichment.png",
                                       title=f"Neighborhood enrichment ({cfg.group_by})")
            except Exception as exc:
                logger.warning("%s: nhood_enrichment failed: %s", self.name, exc)
        else:
            logger.info(
                "%s: group_by='%s' not in adata.obs; skipping nhood_enrichment.",
                self.name, cfg.group_by,
            )

        # ---- Co-occurrence ---------------------------------------------------
        cooc_path = None
        if cfg.group_by in adata.obs.columns:
            try:
                sq.gr.co_occurrence(adata, cluster_key=cfg.group_by)
                cooc = adata.uns.get(f"{cfg.group_by}_co_occurrence")
                if isinstance(cooc, dict) and "occ" in cooc:
                    occ = np.asarray(cooc["occ"])
                    cooc_path = ctx.table_dir / "spatial_co_occurrence.parquet"
                    # Flatten to a long-form table for parquet portability.
                    df = pd.DataFrame(occ.reshape(occ.shape[0], -1))
                    df.to_parquet(cooc_path)
                    adata.uns["spatial_co_occurrence"] = {
                        "group_by": cfg.group_by,
                        "shape": list(occ.shape),
                        "table_path": str(cooc_path),
                    }
            except Exception as exc:
                logger.warning("%s: co_occurrence failed: %s", self.name, exc)

        ctx.metadata["spatial_neighborhoods_status"] = "ok"
        ctx.metadata["spatial_neighborhoods_outputs"] = {
            "morans_i": str(morans_path) if morans_path else None,
            "nhood_enrichment": str(enrich_path) if enrich_path else None,
            "co_occurrence": str(cooc_path) if cooc_path else None,
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _merge_runtime_config(self, ctx: PipelineContext) -> SpatialNeighborhoodsConfig:
        runtime = getattr(ctx.cfg, "spatial_neighborhoods", None)
        if runtime is None:
            return self.config
        if isinstance(runtime, SpatialNeighborhoodsConfig):
            return runtime
        if isinstance(runtime, dict):
            base = self.config
            return SpatialNeighborhoodsConfig(
                group_by=str(runtime.get("group_by", base.group_by)),
                n_top_genes=int(runtime.get("n_top_genes", base.n_top_genes)),
                n_neighs=int(runtime.get("n_neighs", base.n_neighs)),
                coord_type=str(runtime.get("coord_type", base.coord_type)),
                obsm_key=str(runtime.get("obsm_key", base.obsm_key)),
            )
        return self.config

    @staticmethod
    def _select_morans_genes(adata, n_top_genes: int) -> list[str]:
        """Return up to n_top_genes gene names for Moran's I.

        Prefer ``var['highly_variable']`` when available, fall back to genes
        with the highest mean expression. Returns at most ``n_top_genes``
        names; never raises.
        """
        if "highly_variable" in adata.var.columns:
            hv = adata.var.index[adata.var["highly_variable"].astype(bool)].tolist()
            if hv:
                return hv[:n_top_genes]
        # Fallback: top genes by mean expression. Avoid densifying X by working
        # via a per-column sum where possible.
        try:
            from scipy import sparse

            if sparse.issparse(adata.X):
                # densify-allowed: row mean across all cells produces a 1D vector of size n_genes,
                # which is small (10s of thousands) regardless of cell count.
                means = np.asarray(adata.X.mean(axis=0)).ravel()
            else:
                means = np.asarray(adata.X).mean(axis=0)
        except Exception:
            return list(adata.var_names[:n_top_genes])
        order = np.argsort(-means)[:n_top_genes]
        return [str(adata.var_names[i]) for i in order]

    @staticmethod
    def _plot_heatmap(matrix: np.ndarray, ctx: PipelineContext, filename: str, *, title: str) -> None:
        if matrix.size == 0:
            return
        fig, ax = plt.subplots(figsize=(5, 4))
        im = ax.imshow(matrix, cmap="viridis", aspect="auto")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(ctx.figure_dir / filename, dpi=160, bbox_inches="tight")
        plt.close(fig)
