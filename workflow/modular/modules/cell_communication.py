from __future__ import annotations

import logging

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..context import PipelineContext


logger = logging.getLogger(__name__)


__references__ = {
    "Dimitrov_LIANA_2022": {
        "title": "Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data",
        "authors": "Dimitrov et al.",
        "journal": "Nature Communications",
        "year": "2022",
        "doi": "10.1038/s41467-022-30755-0",
        "description": "LIANA benchmark + framework. Used here via liana-py.",
    },
    "Efremova_CellPhoneDB_2020": {
        "title": "CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes",
        "authors": "Efremova et al.",
        "journal": "Nature Protocols",
        "year": "2020",
        "doi": "10.1038/s41596-020-0292-x",
        "description": "Source ligand-receptor resource consumed by LIANA.",
    },
}



class CellCommunicationModule:
    """Optional module: cell-cell communication analysis via ligand-receptor interactions.

    Supports multiple backends:
    1. LIANA (preferred) — multi-method consensus scoring
    2. Fallback — manual ligand-receptor database scoring

    Requires cell type annotations in adata.obs['cell_type'].
    """

    name = "cell_communication"
    requires_keys = {"obs": ["cell_type"]}

    def run(self, ctx: PipelineContext) -> None:
        import os
        if os.environ.get("SC_MEM_GUARD", "").lower() == "on":
            import logging as _logging
            from .._mem_guard import MemoryGuard, MemoryGuardError
            _logger = _logging.getLogger(__name__)
            try:
                with MemoryGuard(ctx, self.name) as mg:
                    mg.check("entry")
                    return self._run_impl(ctx)
            except MemoryGuardError as exc:
                _logger.warning("MemoryGuard: %s", exc)
                ctx.metadata.setdefault("mem_warnings", []).append(
                    {"module": self.name, "error": str(exc)}
                )
        return self._run_impl(ctx)

    def _run_impl(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Cell communication requires AnnData.")
        if "cell_type" not in adata.obs.columns:
            raise ValueError("Cell communication requires cell type annotations.")

        # LIANA (multi-method consensus over several L-R scoring functions) and the
        # built-in manual L-R scorer are NOT equivalent: the manual path scores a single
        # heuristic against a bundled ligand-receptor list, with no consensus across
        # methods and no permutation-based specificity. Substituting one for the other
        # silently makes the two indistinguishable downstream.
        #
        # This previously read `except ImportError: pass` -- no log line at all, and the
        # module recorded no engine or claim status anywhere. liana is declared in
        # environment.yml but installed only in sc_gpu, so every default CPU-env run took
        # the manual path invisibly. Found 2026-08-02 by the engine-fallback marking gate,
        # which caught what an earlier hand-written sweep missed.
        try:
            self._run_liana(adata, ctx)
            engine = "liana_consensus"
            inference_class = "multi_method_consensus"
            inference_status = "supported_confirmatory"
            claimable = True
        except ImportError:
            self._run_manual_lr(adata, ctx)
            engine = "manual_lr_fallback"
            inference_class = "single_heuristic_lr_fallback"
            inference_status = "exploratory_nonclaimable_lr_fallback"
            claimable = False
            logger.warning(
                "LIANA is NOT INSTALLED -- fell back to manual ligand-receptor scoring: a "
                "single heuristic over a bundled L-R list, with no multi-method consensus "
                "and no permutation-based specificity. Result marked non-claimable "
                "(inference_status=%s). Install liana (declared in environment.yml) to "
                "obtain a claimable result.",
                inference_status,
            )

        ctx.metadata["cell_communication_engine"] = engine
        ctx.metadata["cell_communication_inference_class"] = inference_class
        ctx.metadata["cell_communication_inference_status"] = inference_status
        ctx.metadata["cell_communication_claimable"] = claimable

    def _run_liana(self, adata, ctx: PipelineContext) -> None:
        """Run LIANA multi-method consensus."""
        import liana as li

        li.mt.rank_aggregate(
            adata,
            groupby="cell_type",
            resource_name="consensus",
            use_raw=False,
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
                        ctx.figure_dir / "cell_communication_dotplot.png", bbox_inches="tight",
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

        # Pre-compute mean expression matrix: (n_cell_types, n_needed_genes).
        # The sparse path (single-pass indicator-matrix groupby) and the dense
        # path produce numerically identical means; pick sparse whenever X is
        # sparse (the post-QC default) to avoid per-cell-type densification.
        # SC_CELLCOMM_ENGINE overrides: "sparse" forces sparse, "dense" forces dense.
        import os
        import scipy.sparse as _sp
        engine = os.environ.get("SC_CELLCOMM_ENGINE", "auto").lower()
        x_is_sparse = _sp.issparse(expr.X) or hasattr(expr.X, "tocsr")
        if engine == "sparse":
            use_sparse = True
        elif engine == "dense":
            use_sparse = False
        else:  # "auto" (default): sparse when X is sparse
            use_sparse = x_is_sparse
        ct_list = []
        if use_sparse and hasattr(expr.X, "tocsr"):
            from .._sparse_utils import sparse_groupby_mean
            X_sub = expr.X[:, gene_indices]
            valid_cts = [ct for ct in cell_types
                         if (adata.obs["cell_type"] == ct).values.sum() >= 5]
            if not valid_cts:
                raise ValueError("No cell types with >= 5 cells.")
            ct_array = adata.obs["cell_type"].values
            keep_mask = np.isin(ct_array, valid_cts)
            X_kept = X_sub[keep_mask]
            labels_kept = ct_array[keep_mask]
            mean_dict = sparse_groupby_mean(X_kept, labels_kept)
            mean_mat_rows = []
            for ct in valid_cts:
                ct_list.append(ct)
                mean_mat_rows.append(np.asarray(mean_dict[ct], dtype=np.float32))
        else:
            from .._mem_guard import MemoryGuard, MemoryAbortError
            mean_mat_rows = []
            for ct in cell_types:
                if MemoryGuard.abort_requested():
                    raise MemoryAbortError("watchdog abort during cell communication cell-type loop")
                mask = (adata.obs["cell_type"] == ct).values
                if mask.sum() < 5:
                    continue
                ct_list.append(ct)
                chunk = expr.X[mask][:, gene_indices]
                if hasattr(chunk, "toarray"):
                    # densify-allowed: per-cell-type row-slice (≥5 cells) × n_needed_genes; LIANA fallback path, size bounded by gene list
                    chunk = chunk.toarray()
                mean_mat_rows.append(np.asarray(chunk, dtype=np.float32).mean(axis=0))

        if not ct_list:
            raise ValueError("No cell types with >= 5 cells.")

        mean_mat = np.stack(mean_mat_rows)  # (n_cell_types, n_needed_genes)

        # Score all L-R pairs via vectorized matrix lookups
        result_frames = []
        for ligand, receptor in valid_pairs:
            lig_col = gene_col_map[ligand]
            rec_col = gene_col_map[receptor]
            lig_vec = mean_mat[:, lig_col]   # (n_cell_types,)
            rec_vec = mean_mat[:, rec_col]   # (n_cell_types,)
            score_mat = lig_vec[:, None] * rec_vec[None, :]  # outer product
            ii, jj = np.where(score_mat > 0)
            if len(ii) == 0:
                continue
            result_frames.append(pd.DataFrame({
                "ligand": ligand,
                "receptor": receptor,
                "source": [ct_list[i] for i in ii],
                "target": [ct_list[j] for j in jj],
                "ligand_expr": np.round(lig_vec[ii], 4),
                "receptor_expr": np.round(rec_vec[jj], 4),
                "lr_score": np.round(score_mat[ii, jj], 4),
            }))
        results = result_frames

        if not results:
            raise ValueError("No significant L-R interactions detected.")

        df = pd.concat(results, ignore_index=True).sort_values("lr_score", ascending=False)
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
        plt.savefig(ctx.figure_dir / "cell_communication_heatmap.png", bbox_inches="tight")
        plt.close()


