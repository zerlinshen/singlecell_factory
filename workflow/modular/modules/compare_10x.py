from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from ..context import PipelineContext


class Compare10XModule:
    """Optional module: compare our outputs against 10X Cell Ranger outputs."""

    name = "compare_10x"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("10X comparison requires AnnData.")

        outs = Path(ctx.cfg.cellranger.outs_dir).parent
        metrics_csv = outs / "metrics_summary.csv"
        cluster_csv = (
            outs / "analysis" / "clustering" / "gene_expression_graphclust" / "clusters.csv"
        )

        rows = []
        if metrics_csv.exists():
            metrics = pd.read_csv(metrics_csv).iloc[0]
            rows.extend(
                [
                    (
                        "Estimated Number of Cells",
                        _pnum(metrics.get("Estimated Number of Cells")),
                        float(adata.n_obs),
                    ),
                    (
                        "Median Genes per Cell",
                        _pnum(metrics.get("Median Genes per Cell")),
                        float(adata.obs["n_genes_by_counts"].median()),
                    ),
                    (
                        "Median UMI Counts per Cell",
                        _pnum(metrics.get("Median UMI Counts per Cell")),
                        float(adata.obs["total_counts"].median()),
                    ),
                    (
                        "Total Genes Detected",
                        _pnum(metrics.get("Total Genes Detected")),
                        float(ctx.metadata.get("genes_after_qc", adata.n_vars)),
                    ),
                ]
            )
        df = pd.DataFrame(rows, columns=["metric", "tenx", "ours"])
        if not df.empty:
            df["diff_pct"] = (df["ours"] - df["tenx"]) / df["tenx"] * 100
            df.to_csv(ctx.table_dir / "metrics_vs_10x.csv", index=False)

        summary: dict[str, float] = {}
        if cluster_csv.exists() and "leiden" in adata.obs:
            cl = pd.read_csv(cluster_csv).rename(columns={"Barcode": "barcode", "Cluster": "cluster_10x"})
            cl = cl.set_index("barcode")
            obs = adata.obs[["leiden"]].join(cl, how="inner").dropna()
            if len(obs) > 0:
                y_true = obs["cluster_10x"].astype(str).values
                y_pred = obs["leiden"].astype(str).values
                summary["cluster_ari"] = float(adjusted_rand_score(y_true, y_pred))
                summary["cluster_nmi"] = float(normalized_mutual_info_score(y_true, y_pred))
                summary["cluster_overlap_cells"] = int(len(obs))

        (ctx.table_dir / "compare_10x_summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        ctx.metadata["compare_10x"] = summary


def _pnum(value: object) -> float:
    if value is None:
        return float("nan")
    raw = str(value).replace(",", "").replace("%", "").strip()
    try:
        return float(raw)
    except ValueError:
        return float("nan")
