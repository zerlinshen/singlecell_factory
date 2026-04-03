from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..context import PipelineContext


class CellCommunicationModule:
    """Optional module: cell-cell communication analysis via ligand-receptor interactions.

    Supports multiple backends:
    1. LIANA (preferred) — multi-method consensus scoring
    2. Fallback — manual ligand-receptor database scoring

    Requires cell type annotations in adata.obs['cell_type'].
    """

    name = "cell_communication"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Cell communication requires AnnData.")
        if "cell_type" not in adata.obs.columns:
            raise ValueError("Cell communication requires cell type annotations.")

        # Try LIANA first
        try:
            self._run_liana(adata, ctx)
            return
        except ImportError:
            pass

        # Fallback: manual L-R scoring
        self._run_manual_lr(adata, ctx)

    def _run_liana(self, adata, ctx: PipelineContext) -> None:
        """Run LIANA multi-method consensus."""
        import liana as li

        li.mt.rank_aggregate(
            adata,
            groupby="cell_type",
            resource_name="consensus",
            use_raw=True,
            verbose=False,
        )

        if hasattr(adata, "uns") and "liana_res" in adata.uns:
            results = adata.uns["liana_res"]
            if isinstance(results, pd.DataFrame) and not results.empty:
                results.to_csv(ctx.table_dir / "cell_communication_liana.csv", index=False)

                top_interactions = results.head(50)
                ctx.metadata["cc_top_interactions"] = len(top_interactions)

                # Dot plot
                try:
                    li.pl.dotplot(
                        adata=adata,
                        colour="magnitude_rank",
                        size="specificity_rank",
                        inverse_size=True,
                        inverse_colour=True,
                        source_labels=adata.obs["cell_type"].unique()[:5],
                        target_labels=adata.obs["cell_type"].unique()[:5],
                        top_n=20,
                        figure_size=(12, 8),
                    )
                    plt.savefig(
                        ctx.figure_dir / "cell_communication_dotplot.png",
                        dpi=160, bbox_inches="tight",
                    )
                    plt.close()
                except Exception:
                    plt.close("all")

    def _run_manual_lr(self, adata, ctx: PipelineContext) -> None:
        """Fallback: score known ligand-receptor pairs manually."""
        import numpy as np

        # Curated L-R pairs relevant to tumor microenvironment
        lr_pairs = [
            ("CD274", "PDCD1"),      # PD-L1 / PD-1
            ("CD80", "CTLA4"),       # CD80 / CTLA-4
            ("CXCL12", "CXCR4"),    # SDF-1 / CXCR4
            ("CCL2", "CCR2"),        # CCL2 / CCR2
            ("VEGFA", "FLT1"),       # VEGF-A / VEGFR1
            ("VEGFA", "KDR"),        # VEGF-A / VEGFR2
            ("TGFB1", "TGFBR2"),    # TGF-b1 / TGF-bR2
            ("TNF", "TNFRSF1A"),    # TNF-a / TNFR1
            ("FAS", "FASLG"),        # Fas / FasL
            ("HGF", "MET"),          # HGF / c-Met
            ("PDGFA", "PDGFRA"),     # PDGF-A / PDGFR-a
            ("DLL1", "NOTCH1"),      # Delta / Notch
            ("WNT5A", "FZD5"),       # Wnt5a / Frizzled-5
            ("SPP1", "CD44"),        # Osteopontin / CD44
            ("CDH1", "CDH1"),        # E-cadherin homophilic
            ("CXCL10", "CXCR3"),    # CXCL10 / CXCR3
            ("IL6", "IL6R"),         # IL-6 / IL-6R
            ("CCL5", "CCR5"),        # RANTES / CCR5
        ]

        var_names = set(adata.var_names)
        valid_pairs = [(l, r) for l, r in lr_pairs if l in var_names and r in var_names]

        if not valid_pairs:
            raise ValueError("No valid ligand-receptor pairs found in dataset.")

        # Collect all genes needed for valid L-R pairs
        needed_genes = sorted({g for pair in valid_pairs for g in pair})
        cell_types = adata.obs["cell_type"].unique()
        expr = adata.raw.to_adata() if adata.raw else adata

        gene_to_idx = {g: i for i, g in enumerate(expr.var_names)}
        gene_indices = [gene_to_idx[g] for g in needed_genes]
        gene_col_map = {g: i for i, g in enumerate(needed_genes)}

        # Pre-compute mean expression matrix: (n_cell_types, n_needed_genes)
        ct_list = []
        mean_mat_rows = []
        for ct in cell_types:
            mask = (adata.obs["cell_type"] == ct).values
            if mask.sum() < 5:
                continue
            ct_list.append(ct)
            chunk = expr.X[mask][:, gene_indices]
            if hasattr(chunk, "toarray"):
                chunk = chunk.toarray()
            mean_mat_rows.append(np.asarray(chunk, dtype=np.float32).mean(axis=0))

        if not ct_list:
            raise ValueError("No cell types with >= 5 cells.")

        mean_mat = np.stack(mean_mat_rows)  # (n_cell_types, n_needed_genes)

        # Score all L-R pairs via pre-computed matrix lookups
        results = []
        for ligand, receptor in valid_pairs:
            lig_col = gene_col_map[ligand]
            rec_col = gene_col_map[receptor]
            lig_vec = mean_mat[:, lig_col]   # (n_cell_types,)
            rec_vec = mean_mat[:, rec_col]   # (n_cell_types,)
            score_mat = lig_vec[:, None] * rec_vec[None, :]  # outer product
            for i, source in enumerate(ct_list):
                for j, target in enumerate(ct_list):
                    score = score_mat[i, j]
                    if score > 0:
                        results.append({
                            "ligand": ligand,
                            "receptor": receptor,
                            "source": source,
                            "target": target,
                            "ligand_expr": round(float(lig_vec[i]), 4),
                            "receptor_expr": round(float(rec_vec[j]), 4),
                            "lr_score": round(float(score), 4),
                        })

        if not results:
            raise ValueError("No significant L-R interactions detected.")

        df = pd.DataFrame(results).sort_values("lr_score", ascending=False)
        df.to_csv(ctx.table_dir / "cell_communication_lr.csv", index=False)
        ctx.metadata["cc_top_interactions"] = min(50, len(df))

        # Heatmap of top interactions
        self._plot_lr_heatmap(df, ctx)

    @staticmethod
    def _plot_lr_heatmap(df: pd.DataFrame, ctx: PipelineContext) -> None:
        """Plot heatmap of top L-R interaction scores."""
        import numpy as np

        top = df.head(30).copy()
        top["pair"] = top["ligand"] + " → " + top["receptor"]
        top["interaction"] = top["source"] + " → " + top["target"]

        pivot = top.pivot_table(
            index="pair", columns="interaction", values="lr_score", aggfunc="max",
        ).fillna(0)

        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 0.8), max(6, len(pivot) * 0.4)))
        im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=7)
        plt.colorbar(im, ax=ax, label="L-R score")
        ax.set_title("Top ligand-receptor interactions")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "cell_communication_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()


