from __future__ import annotations

import json
import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


class PathwayAnalysisModule:
    """Optional module: gene set enrichment analysis on DE results.

    Supports multiple backends:
    1. gseapy (preferred) — uses MSigDB gene sets
    2. decoupler — uses PROGENy pathway activity scoring
    3. Fallback — simple overlap-based enrichment with bundled gene sets
    """

    name = "pathway_analysis"

    # Hallmark gene sets subset for fallback mode
    HALLMARK_SETS = {
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": [
            "JUNB", "CXCL2", "ATF3", "NFKBIA", "TNFAIP3", "PTGS2", "CXCL1",
            "IER3", "CD83", "CCL20", "CXCL3", "MAFF", "NFKB2", "TNFAIP2",
        ],
        "HALLMARK_HYPOXIA": [
            "PGK1", "LDHA", "ALDOA", "ENO1", "GAPDH", "TPI1", "PGAM1",
            "PKM", "GPI", "VEGFA", "SLC2A1", "HK2", "BNIP3", "CA9",
        ],
        "HALLMARK_P53_PATHWAY": [
            "CDKN1A", "MDM2", "BAX", "FAS", "GADD45A", "SESN1", "TP53I3",
            "RRM2B", "DDB2", "FDXR", "PMAIP1", "BBC3", "ZMAT3", "BTG2",
        ],
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": [
            "VIM", "CDH2", "FN1", "SNAI2", "MMP2", "MMP3", "SPARC",
            "COL1A1", "COL3A1", "ACTA2", "TAGLN", "THY1", "LOX", "TWIST1",
        ],
        "HALLMARK_INFLAMMATORY_RESPONSE": [
            "IL1B", "CXCL8", "CCL2", "CXCL10", "ICAM1", "IRF1", "SOCS3",
            "STAT1", "IL6", "CCL5", "CXCL11", "ISG15", "GBP1", "TAP1",
        ],
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": [
            "IRF1", "STAT1", "GBP1", "ISG15", "CXCL10", "CXCL11", "TAP1",
            "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "CIITA", "IDO1",
        ],
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION": [
            "ATP5F1A", "ATP5F1B", "COX5A", "COX7C", "NDUFA1", "NDUFB8",
            "UQCRC1", "SDHA", "SDHB", "CYC1", "ATP5MC1", "NDUFS3",
        ],
        "HALLMARK_MYC_TARGETS_V1": [
            "NCL", "NPM1", "LDHA", "ENO1", "CCT5", "RAN", "SRM", "NOP56",
            "HSPD1", "PHB", "PA2G4", "CDK4", "PCNA", "RRP1",
        ],
        "HALLMARK_APOPTOSIS": [
            "CASP3", "CASP8", "BAX", "BCL2L11", "BID", "DIABLO", "CYCS",
            "TNFRSF10B", "BIRC3", "CFLAR", "FAS", "LMNA", "PMAIP1",
        ],
        "HALLMARK_ANGIOGENESIS": [
            "VEGFA", "FLT1", "KDR", "PECAM1", "CDH5", "ANGPT2", "NRP1",
            "FGF2", "PDGFRB", "TEK", "THBS1", "COL4A1",
        ],
    }

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Pathway analysis requires AnnData.")

        # Load DE results from the differential_expression module's output dir
        de_dir = ctx.module_output_dir("differential_expression")
        de_csv = de_dir / "marker_genes.csv" if de_dir else None
        if de_csv is None or not de_csv.exists():
            raise ValueError("Pathway analysis requires DE results (run differential_expression first).")

        de_df = pd.read_csv(de_csv)

        # Try gseapy first, then decoupler, then fallback
        results = None
        try:
            results = self._run_gseapy(de_df, adata, ctx)
        except ImportError:
            logger.info("gseapy not available, trying decoupler backend")
        except Exception as exc:
            logger.warning("gseapy enrichment failed: %s", exc)

        if results is None:
            try:
                results = self._run_decoupler(adata, ctx)
            except ImportError:
                logger.info("decoupler not available, using built-in fallback")
            except Exception as exc:
                logger.warning("decoupler pathway analysis failed: %s", exc)

        if results is None:
            results = self._run_fallback(de_df, adata, ctx)

        if results is not None and not results.empty:
            results.to_csv(ctx.table_dir / "pathway_enrichment.csv", index=False)
            ctx.metadata["pathway_top_terms"] = (
                results.head(10)["term"].tolist() if "term" in results.columns else []
            )
            self._plot_enrichment(results, ctx)

    def _run_gseapy(self, de_df: pd.DataFrame, adata, ctx) -> pd.DataFrame | None:
        """Run enrichment using gseapy with MSigDB Hallmark gene sets."""
        import gseapy as gp

        all_results = []
        for cluster in de_df["group"].unique():
            cluster_genes = (
                de_df[de_df["group"] == cluster]
                .sort_values("scores", ascending=False)["names"]
                .tolist()
            )
            if len(cluster_genes) < 5:
                continue
            try:
                enr = gp.enrich(
                    gene_list=cluster_genes[:200],
                    gene_sets="MSigDB_Hallmark_2020",
                    organism="human",
                    outdir=None,
                    no_plot=True,
                )
                if enr.results is not None and not enr.results.empty:
                    res = enr.results.copy()
                    res["cluster"] = cluster
                    all_results.append(res)
            except Exception:
                continue

        if not all_results:
            return None

        combined = pd.concat(all_results, ignore_index=True)
        combined = combined.rename(columns={"Term": "term", "Adjusted P-value": "padj"})
        return combined

    def _run_decoupler(self, adata, ctx) -> pd.DataFrame | None:
        """Run pathway activity scoring using decoupler + PROGENy."""
        import decoupler as dc

        progeny = dc.get_progeny(organism="human", top=300)
        dc.run_mlm(
            mat=adata,
            net=progeny,
            source="source",
            target="target",
            weight="weight",
            use_raw=False,
        )
        if "mlm_estimate" not in adata.obsm:
            return None

        # Per-cluster mean pathway activity
        acts = pd.DataFrame(
            adata.obsm["mlm_estimate"],
            index=adata.obs_names,
            columns=progeny["source"].unique(),
        )
        acts["leiden"] = adata.obs["leiden"].values

        cluster_acts = acts.groupby("leiden", observed=True).mean()
        cluster_acts.to_csv(ctx.table_dir / "pathway_activity_per_cluster.csv")

        # Heatmap
        fig, ax = plt.subplots(figsize=(10, max(4, len(cluster_acts) * 0.4)))
        im = ax.imshow(cluster_acts.values, aspect="auto", cmap="RdBu_r")
        ax.set_xticks(range(cluster_acts.shape[1]))
        ax.set_xticklabels(cluster_acts.columns, rotation=90, fontsize=8)
        ax.set_yticks(range(cluster_acts.shape[0]))
        ax.set_yticklabels(cluster_acts.index)
        ax.set_xlabel("Pathway")
        ax.set_ylabel("Cluster")
        plt.colorbar(im, ax=ax, label="Activity score")
        plt.title("PROGENy pathway activity per cluster")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pathway_activity_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()

        # Flatten to result format
        rows = []
        for cluster in cluster_acts.index:
            for pathway in cluster_acts.columns:
                rows.append({
                    "cluster": cluster,
                    "term": pathway,
                    "score": cluster_acts.loc[cluster, pathway],
                })
        return pd.DataFrame(rows)

    def _run_fallback(self, de_df: pd.DataFrame, adata, ctx) -> pd.DataFrame:
        """Simple overlap-based enrichment with built-in hallmark gene sets."""
        from scipy.stats import hypergeom

        background = set(adata.var_names)
        all_results = []

        for cluster in de_df["group"].unique():
            gene_list = set(
                de_df[de_df["group"] == cluster]
                .sort_values("scores", ascending=False)["names"]
                .head(100)
                .tolist()
            )
            for term, gs_genes in self.HALLMARK_SETS.items():
                gs_in_bg = set(gs_genes) & background
                overlap = gene_list & gs_in_bg
                if not gs_in_bg or not overlap:
                    continue
                M = len(background)
                n = len(gs_in_bg)
                N = len(gene_list)
                k = len(overlap)
                pval = hypergeom.sf(k - 1, M, n, N)
                all_results.append({
                    "cluster": cluster,
                    "term": term,
                    "overlap": k,
                    "gene_set_size": n,
                    "pval": pval,
                    "genes": ",".join(sorted(overlap)),
                })

        if not all_results:
            return pd.DataFrame()

        df = pd.DataFrame(all_results)
        # Multiple testing correction (Benjamini-Hochberg FDR)
        n_tests = len(df)
        df = df.sort_values("pval")
        ranks = np.arange(1, n_tests + 1)
        raw_adj = df["pval"].values * n_tests / ranks
        df["padj"] = np.minimum.accumulate(raw_adj[::-1])[::-1]
        df["padj"] = df["padj"].clip(upper=1.0)
        return df.sort_values("padj")

    @staticmethod
    def _plot_enrichment(results: pd.DataFrame, ctx) -> None:
        """Bar plot of top enriched terms."""
        if "padj" in results.columns:
            top = results.sort_values("padj").drop_duplicates("term").head(15)
            y_val = -np.log10(top["padj"].clip(lower=1e-50))
            xlabel = "-log10(adj p-value)"
        elif "score" in results.columns:
            top = results.groupby("term", observed=True)["score"].mean().abs().nlargest(15).reset_index()
            y_val = top["score"]
            xlabel = "Mean activity score"
        else:
            return

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(range(len(top)), y_val, color="steelblue", edgecolor="black", linewidth=0.5)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["term"].values, fontsize=8)
        ax.set_xlabel(xlabel)
        ax.set_title("Top enriched pathways")
        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pathway_enrichment_bar.png", dpi=160, bbox_inches="tight")
        plt.close()
