from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext


__references__ = {
    "Tirosh_cellcycle_2016": {
        "title": "Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma",
        "authors": "Tirosh et al.",
        "journal": "Nature",
        "year": "2016",
        "doi": "10.1038/nature20123",
        "description": "Source of canonical S and G2/M cell-cycle gene sets; sc.tl.score_genes_cell_cycle implements this.",
    },
}


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
            ctx.status("cell_cycle", "skipped", "Insufficient cell cycle genes found")
            return

        # random_state: score_genes_cell_cycle's control-bin sampling is
        # RNG-dependent (Tirosh 2016); propagate the canonical seed.
        sc.tl.score_genes_cell_cycle(
            adata, s_genes=s_genes, g2m_genes=g2m_genes,
            random_state=getattr(ctx, "random_state", 0),
        )

        # Record phase distribution
        phase_counts = adata.obs["phase"].value_counts().to_dict()
        ctx.metadata["cell_cycle_phases"] = phase_counts

        # Optionally regress out cell cycle effects.
        if ctx.cfg.regress_cell_cycle:
            # Guard (Luecken & Theis 2019) + Wave-2 design C: this module runs
            # AFTER clustering (module_catalog depends_on=("clustering",)), which
            # already computed PCA/neighbors/UMAP/Leiden. regress_out only
            # rewrites adata.X, so the existing embedding/clustering will NOT
            # change. Stamp that layout is exploratory for cell-cycle-corrected
            # biology rather than reordering the DAG (Wave-3 product option A).
            msg = (
                "regress_cell_cycle=True but cell_cycle runs after clustering; "
                "regress_out rewrites adata.X only and CANNOT affect the already-"
                "computed PCA/UMAP/Leiden embedding. Layout remains pre-regression; "
                "cell_cycle_regressed_after_layout=True and layout claim is "
                "exploratory. To regress cell cycle into the embedding, run "
                "regression before clustering (Wave-3 option A)."
            )
            ctx.status("cell_cycle", "warning", msg)
            ctx.metadata["cell_cycle_regress_out_embedding_effect"] = (
                "none_runs_after_clustering"
            )
            ctx.metadata["cell_cycle_regress_out_warning"] = msg
            ctx.metadata["cell_cycle_regressed_after_layout"] = True
            ctx.metadata["cell_cycle_layout_claimable"] = False
            ctx.metadata["cell_cycle_layout_claim_reason"] = (
                "regress_after_layout_exploratory"
            )
            # Downgrade clustering claim on the authoritative batch_risk surface
            # (manifest consumers read batch_risk.clustering_claim_status).
            risk = ctx.metadata.get("batch_risk")
            if not isinstance(risk, dict):
                risk = {}
                ctx.metadata["batch_risk"] = risk
            prev = risk.get("clustering_claim_status")
            if prev not in (None, "exploratory"):
                risk["clustering_claim_status_before_cell_cycle"] = prev
            risk["clustering_claim_status"] = "exploratory"
            risk["clustering_claim_reason"] = (
                "exploratory_cell_cycle_regressed_after_layout"
            )
            # Mirror top-level keys for modules that only read flat metadata.
            ctx.metadata["clustering_claim_status"] = "exploratory"
            ctx.metadata["clustering_claim_reason"] = risk["clustering_claim_reason"]
            if isinstance(getattr(adata, "uns", None), dict):
                cc = dict(adata.uns.get("cell_cycle") or {})
                cc.update(
                    {
                        "regressed_after_layout": True,
                        "layout_claimable": False,
                        "layout_claim_reason": "regress_after_layout_exploratory",
                        "embedding_effect": "none_runs_after_clustering",
                    }
                )
                adata.uns["cell_cycle"] = cc
            sc.pp.regress_out(adata, ["S_score", "G2M_score"])
        else:
            ctx.metadata["cell_cycle_regressed_after_layout"] = False

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
            plt.savefig(ctx.figure_dir / "cell_cycle_umap.png", bbox_inches="tight")
            plt.close()
