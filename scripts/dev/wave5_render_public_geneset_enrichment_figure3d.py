#!/usr/bin/env python3
"""Render Trevino Figure 3D public-matrix gene-set enrichment evidence.

Figure 3D in the paper reports GO / curated gene-set enrichments for the five
interaction clusters from Figure 3C. This script uses the public-matrix Figure
3C derived-k5 link summary and local author-distributed gene-set files. It does
not claim exact topGO/GO database-version parity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    paths = {
        "png": out_dir / f"{stem}.png",
        "svg": out_dir / f"{stem}.svg",
        "pdf": out_dir / f"{stem}.pdf",
    }
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    }


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    order = np.argsort(pvals)
    ranked = pvals[order]
    n = len(pvals)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out


def clean_gene(x: object) -> str:
    return str(x).strip().upper()


def parse_simple_geneset(path: Path) -> tuple[str, set[str], str]:
    lines = [line.strip() for line in path.read_text(errors="replace").splitlines() if line.strip()]
    if not lines:
        return path.stem, set(), "empty"
    # Standard GMT-like author text: first line set name, second line description.
    if len(lines) >= 2 and lines[1].startswith(">"):
        name = lines[0].strip()
        genes = {clean_gene(x) for x in lines[2:] if not x.startswith(">")}
        return name, genes, "named_text"
    # Some disease-prioritization files are plain one-gene-per-line lists.
    return path.stem, {clean_gene(x) for x in lines if not x.startswith(">")}, "plain_list"


def parse_go_summary(path: Path, min_genes: int) -> list[tuple[str, set[str], str]]:
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return []
    if "Symbol" not in df.columns or "Annotated Term" not in df.columns:
        return []
    out: list[tuple[str, set[str], str]] = []
    for term, sub in df.groupby("Annotated Term", observed=True):
        genes = {clean_gene(x) for x in sub["Symbol"].dropna()}
        if len(genes) >= min_genes:
            out.append((f"{path.stem}:{term}", genes, "go_summary_tsv"))
    return out


def load_gene_sets(gene_set_dir: Path, min_genes: int) -> pd.DataFrame:
    rows = []
    seen: set[tuple[str, tuple[str, ...]]] = set()
    for path in sorted(gene_set_dir.glob("*.txt")):
        if path.name.startswith("._"):
            continue
        parsed: list[tuple[str, set[str], str]]
        if path.name.startswith("GO_term_summary"):
            parsed = parse_go_summary(path, min_genes)
        else:
            parsed = [parse_simple_geneset(path)]
        for name, genes, source_type in parsed:
            genes = {g for g in genes if g and g != "NA"}
            if len(genes) < min_genes:
                continue
            key = (name, tuple(sorted(genes)))
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "gene_set": name,
                    "source_file": str(path),
                    "source_type": source_type,
                    "n_genes_in_set": len(genes),
                    "genes": genes,
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--link-summary", required=True, type=Path)
    p.add_argument("--gene-set-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--cluster-col", default="derived_k5_cluster")
    p.add_argument("--gene-col", default="Gene symbol")
    p.add_argument("--min-gene-set-size", type=int, default=3)
    p.add_argument("--top-n-per-cluster", type=int, default=8)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    links = pd.read_csv(args.link_summary, sep="\t")
    if args.cluster_col not in links.columns:
        raise KeyError(f"Missing cluster column {args.cluster_col}")
    if args.gene_col not in links.columns:
        raise KeyError(f"Missing gene column {args.gene_col}")
    links["_gene_clean"] = links[args.gene_col].map(clean_gene)
    cluster_genes = {
        str(cluster): set(sub["_gene_clean"].dropna())
        for cluster, sub in links.groupby(args.cluster_col, observed=True)
    }
    universe = set(links["_gene_clean"].dropna())
    gene_sets = load_gene_sets(args.gene_set_dir, args.min_gene_set_size)
    if gene_sets.empty:
        raise ValueError(f"No gene sets loaded from {args.gene_set_dir}")

    rows = []
    for cluster, genes in cluster_genes.items():
        not_cluster = universe.difference(genes)
        for _, gs in gene_sets.iterrows():
            gs_genes = set(gs["genes"]).intersection(universe)
            if len(gs_genes) < args.min_gene_set_size:
                continue
            a = len(genes.intersection(gs_genes))
            b = len(genes.difference(gs_genes))
            c = len(not_cluster.intersection(gs_genes))
            d = len(not_cluster.difference(gs_genes))
            if a == 0:
                continue
            odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append(
                {
                    "cluster": cluster,
                    "gene_set": gs["gene_set"],
                    "source_file": gs["source_file"],
                    "source_type": gs["source_type"],
                    "overlap_genes": ",".join(sorted(genes.intersection(gs_genes))),
                    "overlap_n": a,
                    "cluster_gene_n": len(genes),
                    "gene_set_n_in_universe": len(gs_genes),
                    "universe_gene_n": len(universe),
                    "odds_ratio": float(odds) if np.isfinite(odds) else None,
                    "p_value": float(pval),
                }
            )
    if not rows:
        raise ValueError("No nonzero enrichment overlaps found")
    enrich = pd.DataFrame(rows)
    enrich["fdr_bh"] = bh_adjust(enrich["p_value"].to_numpy())
    enrich["neg_log10_fdr"] = -np.log10(np.maximum(enrich["fdr_bh"], np.finfo(float).tiny))
    enrich = enrich.sort_values(["cluster", "fdr_bh", "p_value", "gene_set"]).reset_index(drop=True)
    enrichment_tsv = args.out_dir / "figure3d_geneset_enrichment.tsv"
    enrich.to_csv(enrichment_tsv, sep="\t", index=False)

    top = (
        enrich.sort_values(["cluster", "fdr_bh", "p_value", "gene_set"])
        .groupby("cluster", observed=True)
        .head(args.top_n_per_cluster)
        .copy()
    )
    top_tsv = args.out_dir / "figure3d_top_genesets_by_cluster.tsv"
    top.to_csv(top_tsv, sep="\t", index=False)

    cluster_order = sorted(cluster_genes)
    # Keep the plot readable by taking union of top terms; truncate labels only in the plot.
    plot_terms = list(dict.fromkeys(top["gene_set"].tolist()))
    plot_df = top[top["gene_set"].isin(plot_terms)].copy()
    x_map = {c: i for i, c in enumerate(cluster_order)}
    y_map = {term: i for i, term in enumerate(reversed(plot_terms))}
    plot_df["x"] = plot_df["cluster"].map(x_map)
    plot_df["y"] = plot_df["gene_set"].map(y_map)
    plot_df["size"] = 30 + 25 * np.minimum(plot_df["neg_log10_fdr"], 8)
    plot_df["color"] = np.log2(np.array([x if x and np.isfinite(x) else 1.0 for x in plot_df["odds_ratio"]], dtype=float))

    fig_h = max(4.5, 0.24 * len(plot_terms))
    fig, ax = plt.subplots(figsize=(8, fig_h))
    sc = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["size"],
        c=plot_df["color"],
        cmap="viridis",
        edgecolor="black",
        linewidth=0.25,
    )
    ax.set_xticks(range(len(cluster_order)))
    ax.set_xticklabels(cluster_order)
    ax.set_yticks(range(len(plot_terms)))
    labels = [term if len(term) <= 80 else term[:77] + "..." for term in reversed(plot_terms)]
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Figure 3C derived k=5 interaction cluster")
    ax.set_ylabel("local gene set")
    ax.set_title("Trevino Fig. 3D public-matrix local gene-set enrichment")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.75)
    cbar.set_label("log2 odds ratio")
    ax.grid(axis="x", alpha=0.2)
    figure_outputs = {"geneset_dotplot": save_triple(fig, args.out_dir, "figure3d_geneset_enrichment_dotplot")}

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "GENERATED_PUBLIC_MATRIX_LOCAL_GENESET_ENRICHMENT_NOT_TOPGO_PARITY",
        "truth_boundary": (
            "Public processed-matrix Figure 3C derived-k5 clusters plus local author gene_sets. "
            "This is not exact topGO v2.36.0 / GO database-version parity."
        ),
        "link_summary": str(args.link_summary),
        "gene_set_dir": str(args.gene_set_dir),
        "out_dir": str(args.out_dir),
        "parameters": {
            "cluster_col": args.cluster_col,
            "gene_col": args.gene_col,
            "min_gene_set_size": args.min_gene_set_size,
            "top_n_per_cluster": args.top_n_per_cluster,
        },
        "cluster_gene_counts": {k: int(len(v)) for k, v in cluster_genes.items()},
        "universe_gene_n": int(len(universe)),
        "gene_sets_loaded": int(len(gene_sets)),
        "tested_enrichment_rows": int(len(enrich)),
        "significant_rows_fdr_0_05": int((enrich["fdr_bh"] <= 0.05).sum()),
        "top_terms_by_cluster": {
            cluster: sub[["gene_set", "overlap_n", "odds_ratio", "p_value", "fdr_bh"]].head(args.top_n_per_cluster).to_dict("records")
            for cluster, sub in enrich.groupby("cluster", observed=True)
        },
        "outputs": {
            "enrichment_tsv": str(enrichment_tsv),
            "top_tsv": str(top_tsv),
            "figure_outputs": figure_outputs,
        },
        "checks": {
            "all_clusters_have_genes": all(len(v) > 0 for v in cluster_genes.values()),
            "gene_sets_loaded_gt_0": len(gene_sets) > 0,
            "enrichment_rows_gt_0": len(enrich) > 0,
            "significant_rows_present": bool((enrich["fdr_bh"] <= 0.05).any()),
        },
    }
    summary_path = args.out_dir / "figure3d_geneset_enrichment_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
