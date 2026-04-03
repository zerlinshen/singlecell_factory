from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.ndimage import uniform_filter1d

from ..context import PipelineContext


class CNVInferenceModule:
    """Optional module: infer copy number variations from scRNA-seq expression data.

    Uses a sliding-window approach over chromosomal gene positions to detect
    large-scale gains/losses — critical for distinguishing malignant vs normal cells
    in tumor samples. Inspired by InferCNV/CopyKAT methodology.
    """

    name = "cnv_inference"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("CNV inference requires AnnData.")

        cfg = ctx.cfg.cnv
        reference_group = cfg.reference_group  # e.g., "Fibroblast" or None

        # Attempt to load chromosome position annotations
        gene_pos = self._get_gene_positions(adata)
        if gene_pos is None or gene_pos.empty:
            # Auto-annotate from Ensembl/BioMart
            self._annotate_gene_positions(adata)
            gene_pos = self._get_gene_positions(adata)

        if gene_pos is not None and not gene_pos.empty:
            self._run_window_cnv(adata, ctx, gene_pos, reference_group)
        else:
            # Fallback: use infercnvpy if available
            self._run_infercnvpy(adata, ctx, reference_group)

    def _run_infercnvpy(self, adata, ctx: PipelineContext, reference_group: str | None) -> None:
        """Use infercnvpy package if available."""
        try:
            import infercnvpy as cnv
        except ImportError:
            raise ValueError("infercnvpy not installed and no gene position data for window-based CNV")

        if reference_group and "cell_type" in adata.obs:
            ref_key = "cell_type"
            ref_cat = [reference_group]
        else:
            ref_key = None
            ref_cat = None

        cnv.tl.infercnv(
            adata,
            reference_key=ref_key,
            reference_cat=ref_cat,
            window_size=ctx.cfg.cnv.window_size,
        )
        cnv.tl.cnv_score(adata)

        adata.obs["cnv_score"] = adata.obs.get("cnv_score", 0.0)
        adata.obs[["cnv_score"]].to_csv(ctx.table_dir / "cnv_scores.csv")

        # CNV heatmap
        try:
            cnv.pl.chromosome_heatmap(adata, groupby="leiden" if "leiden" in adata.obs else None)
            plt.savefig(ctx.figure_dir / "cnv_heatmap.png", dpi=160, bbox_inches="tight")
            plt.close()
        except Exception:
            plt.close("all")

        ctx.metadata["cnv_mean_score"] = float(adata.obs["cnv_score"].mean())

    def _run_window_cnv(
        self, adata, ctx: PipelineContext, gene_pos: pd.DataFrame,
        reference_group: str | None,
    ) -> None:
        """Expression-based sliding window CNV inference."""
        cfg = ctx.cfg.cnv

        # Align genes with positional data
        common = adata.var_names.intersection(gene_pos.index)
        if len(common) < 100:
            raise ValueError(f"Only {len(common)} genes with position data — too few for CNV.")

        gene_pos = gene_pos.loc[common].sort_values(["chromosome", "start"])
        expr = adata[:, gene_pos.index].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        expr = np.array(expr, dtype=np.float32)

        # Center expression per gene
        if reference_group and "cell_type" in adata.obs:
            ref_mask = adata.obs["cell_type"] == reference_group
            if ref_mask.sum() > 10:
                ref_mean = expr[ref_mask].mean(axis=0)
            else:
                ref_mean = expr.mean(axis=0)
        else:
            ref_mean = expr.mean(axis=0)

        centered = expr - ref_mean

        # Sliding window smoothing per chromosome (vectorized)
        chromosomes = gene_pos["chromosome"].values
        smoothed = np.zeros_like(centered)
        window = cfg.window_size

        for chrom in np.unique(chromosomes):
            chrom_idx = np.where(chromosomes == chrom)[0]
            if len(chrom_idx) < 3:
                smoothed[:, chrom_idx] = centered[:, chrom_idx]
                continue
            smoothed[:, chrom_idx] = uniform_filter1d(
                centered[:, chrom_idx],
                size=min(window, len(chrom_idx)),
                axis=1,
                mode="nearest",
            )

        # CNV score per cell: variance of smoothed signal
        cnv_score = np.var(smoothed, axis=1)
        adata.obs["cnv_score"] = cnv_score
        adata.obsm["X_cnv"] = smoothed

        adata.obs[["cnv_score"]].to_csv(ctx.table_dir / "cnv_scores.csv")
        ctx.metadata["cnv_mean_score"] = float(cnv_score.mean())

        # Visualizations
        if "X_umap" in adata.obsm:
            sc.pl.umap(adata, color="cnv_score", show=False, cmap="Reds")
            plt.savefig(ctx.figure_dir / "cnv_score_umap.png", dpi=160, bbox_inches="tight")
            plt.close()

        # Heatmap: top cells by CNV score
        self._plot_cnv_heatmap(adata, smoothed, gene_pos, ctx)

        # Malignant vs normal classification
        threshold = np.percentile(cnv_score, cfg.malignant_percentile)
        adata.obs["cnv_class"] = np.where(cnv_score > threshold, "Malignant", "Normal")
        ctx.metadata["cnv_malignant_cells"] = int((adata.obs["cnv_class"] == "Malignant").sum())

        (ctx.table_dir / "cnv_classification.json").write_text(
            json.dumps({
                "threshold": float(threshold),
                "malignant_count": ctx.metadata["cnv_malignant_cells"],
                "normal_count": int((adata.obs["cnv_class"] == "Normal").sum()),
            }, indent=2),
            encoding="utf-8",
        )

    @staticmethod
    def _plot_cnv_heatmap(adata, smoothed, gene_pos, ctx: PipelineContext) -> None:
        """Plot CNV heatmap sorted by chromosome position."""
        fig, ax = plt.subplots(figsize=(16, 8))

        # Sort cells by cluster if available
        if "leiden" in adata.obs:
            cell_order = adata.obs.sort_values("leiden").index
            cell_idx = [list(adata.obs_names).index(c) for c in cell_order]
        else:
            cell_idx = list(range(adata.n_obs))

        im = ax.imshow(
            smoothed[cell_idx],
            aspect="auto", cmap="RdBu_r", vmin=-0.5, vmax=0.5,
            interpolation="none",
        )

        # Chromosome boundaries
        chroms = gene_pos["chromosome"].values
        boundaries = []
        for i in range(1, len(chroms)):
            if chroms[i] != chroms[i - 1]:
                boundaries.append(i)
                ax.axvline(i, color="black", linewidth=0.3, alpha=0.5)

        plt.colorbar(im, ax=ax, label="Relative expression")
        ax.set_xlabel("Genes (ordered by chromosome)")
        ax.set_ylabel("Cells")
        ax.set_title("Inferred CNV heatmap")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "cnv_heatmap.png", dpi=120, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _annotate_gene_positions(adata) -> None:
        """Auto-annotate adata.var with chromosome and start position from Ensembl."""
        try:
            from pybiomart import Server
        except ImportError:
            return

        try:
            server = Server(host="http://www.ensembl.org")
            dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[
                "hsapiens_gene_ensembl"
            ]
            result = dataset.query(
                attributes=[
                    "hgnc_symbol",
                    "chromosome_name",
                    "start_position",
                    "end_position",
                ],
            )
        except Exception:
            return

        # Identify the gene name column (varies by Ensembl version)
        gene_col = None
        for col in ["HGNC symbol", "Gene name", "hgnc_symbol"]:
            if col in result.columns:
                gene_col = col
                break
        if gene_col is None:
            return

        # Clean: keep only standard chromosomes
        valid_chroms = {str(c) for c in range(1, 23)} | {"X", "Y"}
        result = result[result["Chromosome/scaffold name"].isin(valid_chroms)]
        result = result.dropna(subset=[gene_col])
        result = result[result[gene_col] != ""]
        result = result.drop_duplicates(subset=[gene_col], keep="first")
        result = result.set_index(gene_col)

        # Map to adata.var
        matched = adata.var_names.intersection(result.index)
        if len(matched) < 100:
            return

        adata.var["chromosome"] = result.loc[matched, "Chromosome/scaffold name"].reindex(
            adata.var_names
        )
        adata.var["start"] = result.loc[matched, "Gene start (bp)"].reindex(
            adata.var_names
        )
        adata.var["end"] = result.loc[matched, "Gene end (bp)"].reindex(
            adata.var_names
        )

    @staticmethod
    def _get_gene_positions(adata) -> pd.DataFrame | None:
        """Try to extract chromosome/position from adata.var or return None."""
        if "chromosome" in adata.var.columns and "start" in adata.var.columns:
            return adata.var[["chromosome", "start"]].dropna()

        # Try common alternative column names
        for chrom_col in ["chr", "chrom", "Chromosome"]:
            for pos_col in ["start", "Start", "gene_start", "start_position"]:
                if chrom_col in adata.var.columns and pos_col in adata.var.columns:
                    df = adata.var[[chrom_col, pos_col]].rename(
                        columns={chrom_col: "chromosome", pos_col: "start"}
                    )
                    return df.dropna()
        return None
