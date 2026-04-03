from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from ..context import PipelineContext
from . import score_gene_sets

# Immune subtype markers based on Zheng et al. (Science 2017),
# Zhang et al. (Nature 2018), and Tirosh et al. (Science 2016).
IMMUNE_SUBTYPES = {
    "CD4_naive": ["CCR7", "LEF1", "SELL", "TCF7", "IL7R"],
    "CD4_memory": ["IL7R", "S100A4", "ANXA1", "CCR6", "AQP3"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
    "Th1": ["TBX21", "IFNG", "CXCR3", "IL12RB2", "STAT4"],
    "Th2": ["GATA3", "IL4", "IL5", "IL13", "CCR4"],
    "Th17": ["RORC", "IL17A", "CCR6", "IL23R", "IL22"],
    "CD8_effector": ["CD8A", "GZMB", "PRF1", "NKG7", "GNLY", "IFNG"],
    "CD8_memory": ["CD8A", "IL7R", "GPR183", "SELL", "TCF7"],
    "CD8_exhausted": ["CD8A", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"],
    "NK_cytotoxic": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "FCGR3A"],
    "Macro_M1": ["CD68", "NOS2", "IL1B", "TNF", "CXCL10", "CD80"],
    "Macro_M2": ["CD68", "CD163", "MRC1", "MSR1", "TGFB1", "IL10"],
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "IRF8"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "CD1E"],
    "pDC": ["LILRA4", "IRF7", "TCF4", "CLEC4C"],
}

EXHAUSTION_GENES = ["LAG3", "HAVCR2", "TIGIT", "PDCD1", "CTLA4", "TOX", "TOX2", "ENTPD1"]
CYTOTOXICITY_GENES = ["GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "GNLY", "FASLG"]
ACTIVATION_GENES = ["CD69", "CD38", "HLA-DRA", "ICOS", "TNFRSF9", "TNFRSF4"]

# Parent cell types from annotation module that are considered immune
IMMUNE_PARENT_TYPES = {
    "T cell", "NK cell", "Myeloid/Macro", "Dendritic cell",
    "Mast cell", "B cell", "Plasma cell",
}


class ImmunePhenotypingModule:
    """Optional module: deep immune cell subtyping with functional scoring.

    Assigns fine-grained immune subtypes (CD4 naive/memory/Treg, CD8 effector/
    exhausted, M1/M2 macrophages, DC subsets) and computes continuous exhaustion,
    cytotoxicity, and activation scores per cell.

    Based on marker panels from:
    - Zheng et al. (Science 2017) — pan-cancer T cell landscape
    - Zhang et al. (Nature 2018) — single-cell landscape of immune cells
    - Tirosh et al. (Science 2016) — melanoma TME dissection

    Requires cell type annotations from the annotation module.
    """

    name = "immune_phenotyping"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "cell_type" not in adata.obs.columns:
            raise ValueError("Immune phenotyping requires cell type annotations.")

        var_names = set(adata.var_names if adata.raw is None else adata.raw.var_names)

        # --- Score immune subtypes ---
        scored_subtypes = score_gene_sets(adata, IMMUNE_SUBTYPES, "immune")
        available_subtypes = {s: IMMUNE_SUBTYPES[s] for s in scored_subtypes}

        if not available_subtypes:
            raise ValueError("Insufficient immune subtype marker genes found in dataset.")

        # Assign fine-grained immune subtype (only to immune cells)
        score_cols = [f"immune_{s}" for s in available_subtypes]
        score_mat = adata.obs[score_cols].values
        subtype_names = list(available_subtypes.keys())

        best_idx = np.argmax(score_mat, axis=1)
        best_scores = np.max(score_mat, axis=1)
        assigned = np.array([subtype_names[i] for i in best_idx])

        # Only assign to cells that are immune in the parent annotation
        is_immune = adata.obs["cell_type"].isin(IMMUNE_PARENT_TYPES).values
        assigned[~is_immune] = "Non-immune"
        assigned[is_immune & (best_scores < 0.05)] = "Immune_unassigned"

        adata.obs["immune_subtype"] = assigned

        # --- Functional signature scores ---
        functional_sets = {
            "exhaustion": EXHAUSTION_GENES,
            "cytotoxicity": CYTOTOXICITY_GENES,
            "activation": ACTIVATION_GENES,
        }
        score_gene_sets(adata, functional_sets, "immune")

        # --- Save tables ---
        out_cols = ["immune_subtype"] + [
            c for c in adata.obs.columns
            if c.startswith("immune_") and c != "immune_subtype"
        ]
        adata.obs[out_cols].to_csv(ctx.table_dir / "immune_phenotyping.csv")

        # Per-cluster summary
        if "leiden" in adata.obs:
            summary = (
                adata.obs[adata.obs["immune_subtype"] != "Non-immune"]
                .groupby(["leiden", "immune_subtype"], observed=True)
                .size()
                .reset_index(name="count")
            )
            summary.to_csv(ctx.table_dir / "immune_subtype_summary.csv", index=False)

        # Metadata
        subtype_counts = adata.obs["immune_subtype"].value_counts().to_dict()
        ctx.metadata["immune_subtypes_detected"] = {
            k: v for k, v in subtype_counts.items()
            if k not in ("Non-immune", "Immune_unassigned")
        }

        # --- Visualizations ---
        self._plot_visualizations(adata, ctx)

    def _plot_visualizations(self, adata, ctx: PipelineContext) -> None:
        if "X_umap" not in adata.obsm:
            return

        # UMAP by immune subtype
        sc.pl.umap(adata, color="immune_subtype", show=False)
        plt.savefig(ctx.figure_dir / "umap_immune_subtype.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Exhaustion and cytotoxicity UMAPs
        for score_name, fname in [
            ("immune_exhaustion", "immune_exhaustion_umap.png"),
            ("immune_cytotoxicity", "immune_cytotoxicity_umap.png"),
        ]:
            if score_name in adata.obs:
                sc.pl.umap(adata, color=score_name, show=False, cmap="YlOrRd")
                plt.savefig(ctx.figure_dir / fname, dpi=160, bbox_inches="tight")
                plt.close()

        # Immune subtype composition per cluster
        immune_cells = adata.obs[
            ~adata.obs["immune_subtype"].isin(["Non-immune", "Immune_unassigned"])
        ]
        if len(immune_cells) > 0 and "leiden" in immune_cells.columns:
            ct_counts = immune_cells.groupby(["leiden", "immune_subtype"], observed=True).size().unstack(fill_value=0)
            ct_frac = ct_counts.div(ct_counts.sum(axis=1), axis=0)
            ct_frac.plot(kind="bar", stacked=True, figsize=(14, 5), colormap="tab20")
            plt.ylabel("Fraction")
            plt.xlabel("Leiden cluster")
            plt.title("Immune subtype composition per cluster")
            plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
            plt.tight_layout()
            plt.savefig(
                ctx.figure_dir / "immune_subtype_composition.png", dpi=160, bbox_inches="tight",
            )
            plt.close()

        # Signature heatmap: exhaustion/cytotoxicity/activation per immune subtype
        sig_cols = [c for c in ["immune_exhaustion", "immune_cytotoxicity", "immune_activation"]
                    if c in adata.obs.columns]
        if sig_cols and len(immune_cells) > 0:
            sig_data = adata.obs.loc[immune_cells.index, sig_cols + ["immune_subtype"]]
            mean_sigs = sig_data.groupby("immune_subtype", observed=True)[sig_cols].mean()
            if len(mean_sigs) >= 2:
                fig, ax = plt.subplots(figsize=(8, max(4, len(mean_sigs) * 0.4)))
                im = ax.imshow(mean_sigs.values, aspect="auto", cmap="YlOrRd")
                ax.set_xticks(range(len(sig_cols)))
                ax.set_xticklabels([c.replace("immune_", "") for c in sig_cols], rotation=45)
                ax.set_yticks(range(len(mean_sigs)))
                ax.set_yticklabels(mean_sigs.index, fontsize=8)
                plt.colorbar(im, ax=ax, label="Mean score")
                ax.set_title("Immune functional signatures per subtype")
                plt.tight_layout()
                plt.savefig(
                    ctx.figure_dir / "immune_signature_heatmap.png", dpi=160, bbox_inches="tight",
                )
                plt.close()
