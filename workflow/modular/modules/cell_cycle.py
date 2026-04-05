from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ..context import PipelineContext

# Tirosh et al. 2016 cell cycle gene sets (commonly used in Scanpy/Seurat)
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


class CellCycleModule:
    """Optional module: score cell cycle phases (S/G2M) and optionally regress out."""

    name = "cell_cycle"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Cell cycle scoring requires AnnData.")

        # Filter gene lists to those present in the dataset (case-insensitive)
        var_upper_to_actual = {v.upper(): v for v in adata.var_names}
        s_genes = [var_upper_to_actual[g.upper()] for g in S_GENES if g.upper() in var_upper_to_actual]
        g2m_genes = [var_upper_to_actual[g.upper()] for g in G2M_GENES if g.upper() in var_upper_to_actual]

        if not s_genes or not g2m_genes:
            ctx.status("cell_cycle", False, "Insufficient cell cycle genes found")
            return

        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

        # Record phase distribution
        phase_counts = adata.obs["phase"].value_counts().to_dict()
        ctx.metadata["cell_cycle_phases"] = phase_counts

        # Optionally regress out cell cycle effects
        if ctx.cfg.regress_cell_cycle:
            sc.pp.regress_out(adata, ["S_score", "G2M_score"])

        # Save per-cell scores
        adata.obs[["S_score", "G2M_score", "phase"]].to_csv(
            ctx.table_dir / "cell_cycle_scores.csv"
        )

        # Visualization: UMAP colored by phase (if UMAP exists)
        if "X_umap" in adata.obsm:
            fig, axes = plt.subplots(1, 3, figsize=(18, 5))
            sc.pl.umap(adata, color="phase", ax=axes[0], show=False)
            axes[0].set_title("Cell cycle phase")
            sc.pl.umap(adata, color="S_score", ax=axes[1], show=False, cmap="YlOrRd")
            axes[1].set_title("S phase score")
            sc.pl.umap(adata, color="G2M_score", ax=axes[2], show=False, cmap="YlOrRd")
            axes[2].set_title("G2M phase score")
            plt.tight_layout()
            plt.savefig(ctx.figure_dir / "cell_cycle_umap.png", dpi=160, bbox_inches="tight")
            plt.close()
