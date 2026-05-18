#!/usr/bin/env python3
"""Render Trevino Figure 2E public/local GPC gene-set enrichment.

The paper defines GPCs as top-decile gene activity-expression correlations
linked to >10 CREs, yielding 185 genes after restricting to the paper's 1,999
variable dorsal-forebrain genes.  The exact Table S2 / variable-gene universe is
not present in the current public resource set, so this script applies the
published rule to the author GA/RNA correlation table and then uses a
count-aligned top-185 subset for Figure 2E-style enrichment.

Enrichment is local Fisher exact testing against gene sets shipped in the public
author code checkout.  This is not exact topGO v2.36.0 / GO database-version
parity.
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


PAPER_GPC_COUNT = 185


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    paths = {ext: out_dir / f"{stem}.{ext}" for ext in ("png", "svg", "pdf")}
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    }


def bh_adjust(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    out = np.empty(n, dtype=float)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = n - rank + 1
        val = min(prev, p[idx] * n / i)
        out[idx] = val
        prev = val
    return out.tolist()


def parse_gene_set(path: Path) -> tuple[str, str, set[str]]:
    text = path.read_text(errors="replace").splitlines()
    name = path.stem
    description = ""
    genes: set[str] = set()
    if not text:
        return name, description, genes
    if text[0].startswith("MGI Gene/Marker ID") and "\tSymbol\t" in text[0]:
        df = pd.read_csv(path, sep="\t", dtype=str)
        genes = {str(x).upper() for x in df["Symbol"].dropna() if str(x).strip()}
        terms = df["Annotated Term"].dropna().astype(str).value_counts()
        description = terms.index[0] if not terms.empty else ""
        name = path.stem
        return name, description, genes
    name = text[0].strip() or path.stem
    start = 1
    if len(text) > 1 and text[1].startswith(">"):
        description = text[1].lstrip("> ").strip()
        start = 2
    for line in text[start:]:
        line = line.strip()
        if not line or line.startswith(">") or "\t" in line:
            continue
        genes.add(line.upper())
    return name, description, genes


def load_gene_sets(gene_set_dir: Path, universe: set[str]) -> list[dict]:
    out: list[dict] = []
    for path in sorted(gene_set_dir.glob("*.txt")):
        if path.name.startswith("._"):
            continue
        name, desc, genes = parse_gene_set(path)
        genes = {g for g in genes if g in universe}
        if len(genes) < 3:
            continue
        out.append({"name": name, "description": desc, "path": str(path), "genes": genes})
    return out


def fisher_enrichment(gpc: set[str], universe: set[str], gene_sets: list[dict]) -> pd.DataFrame:
    rows: list[dict] = []
    bg_non_gpc = universe - gpc
    for gs in gene_sets:
        genes = set(gs["genes"])
        a = len(gpc & genes)
        b = len(gpc - genes)
        c = len(bg_non_gpc & genes)
        d = len(bg_non_gpc - genes)
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "gene_set": gs["name"],
                "description": gs["description"],
                "source_path": gs["path"],
                "overlap": a,
                "gpc_size": len(gpc),
                "gene_set_size_in_universe": len(genes),
                "universe_size": len(universe),
                "odds_ratio": float(odds) if np.isfinite(odds) else np.inf,
                "p_value": float(pval),
                "overlap_genes": ",".join(sorted(gpc & genes)),
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_value_bh"] = bh_adjust(df["p_value"].tolist())
        df = df.sort_values(["q_value_bh", "p_value", "odds_ratio"], ascending=[True, True, False])
    return df


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ga-rna-link-table", required=True, type=Path)
    p.add_argument("--gene-set-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.ga_rna_link_table, sep="\t")
    df = df.dropna(subset=["gene.symbol", "spearman.cor"]).copy()
    df["gene.symbol"] = df["gene.symbol"].astype(str)
    df["gene_upper"] = df["gene.symbol"].str.upper()
    df["spearman.cor"] = pd.to_numeric(df["spearman.cor"], errors="coerce")
    df["linked_cre_count"] = pd.to_numeric(df["linked_cre_count"], errors="coerce").fillna(0).astype(int)
    df = df[np.isfinite(df["spearman.cor"])].copy()
    df = (
        df.sort_values(["linked_cre_count", "spearman.cor", "gene.symbol"], ascending=[False, False, True])
        .drop_duplicates("gene_upper", keep="first")
        .copy()
    )
    universe = set(df["gene_upper"])
    corr_cutoff = float(df["spearman.cor"].quantile(0.9))
    rule = df[(df["spearman.cor"] >= corr_cutoff) & (df["linked_cre_count"] > 10)].copy()
    rule = rule.sort_values(["spearman.cor", "linked_cre_count", "gene.symbol"], ascending=[False, False, True])
    count_aligned = rule.head(PAPER_GPC_COUNT).copy()
    gpc = set(count_aligned["gene_upper"])
    gene_sets = load_gene_sets(args.gene_set_dir, universe)
    enrich = fisher_enrichment(gpc, universe, gene_sets)

    gpc_path = args.out_dir / "figure2e_gpc_count_aligned_top185.tsv"
    rule_path = args.out_dir / "figure2e_gpc_rule_candidates.tsv"
    enrich_path = args.out_dir / "figure2e_local_geneset_enrichment.tsv"
    count_aligned.drop(columns=["gene_upper"]).to_csv(gpc_path, sep="\t", index=False)
    rule.drop(columns=["gene_upper"]).to_csv(rule_path, sep="\t", index=False)
    enrich.to_csv(enrich_path, sep="\t", index=False)

    figure_outputs: dict[str, dict[str, str]] = {}
    top = enrich.head(12).copy()
    if not top.empty:
        top["minus_log10_q"] = -np.log10(np.maximum(top["q_value_bh"].astype(float), 1e-300))
        fig_h = max(4.0, 0.42 * len(top) + 1.4)
        fig, ax = plt.subplots(figsize=(8.5, fig_h), constrained_layout=True)
        labels = top["gene_set"].astype(str).str.replace("_", " ", regex=False)
        ax.barh(labels[::-1], top["minus_log10_q"].to_numpy()[::-1], color="#4C78A8")
        ax.set_xlabel("-log10(BH q-value)")
        ax.set_title("Figure 2E proxy — local gene-set enrichment of derived GPCs")
        for i, (_, row) in enumerate(top.iloc[::-1].iterrows()):
            ax.text(
                row["minus_log10_q"] + 0.02,
                i,
                f"{int(row['overlap'])}/{int(row['gene_set_size_in_universe'])}",
                va="center",
                fontsize=7,
            )
        figure_outputs["figure2e_local_geneset_enrichment_barplot"] = save_triple(
            fig, args.out_dir, "figure2e_local_geneset_enrichment_barplot"
        )

    tf_rows = enrich[enrich["gene_set"].astype(str).str.contains("TRANSCRIPTION_FACTOR|TF", case=False, regex=True)]
    tf_best = tf_rows.iloc[0].to_dict() if not tf_rows.empty else None
    checks = {
        "count_aligned_gpc_size_185": int(count_aligned.shape[0]) == PAPER_GPC_COUNT,
        "published_rule_candidates_ge_185": int(rule.shape[0]) >= PAPER_GPC_COUNT,
        "gene_sets_tested_ge_5": len(gene_sets) >= 5,
        "tf_gene_set_tested": tf_best is not None,
        "tf_enrichment_nominal_p_lt_0_05": bool(tf_best is not None and float(tf_best["p_value"]) < 0.05),
        "all_output_hashes_recorded": all(
            "png_sha256" in entry and "svg_sha256" in entry and "pdf_sha256" in entry
            for entry in figure_outputs.values()
        ),
    }
    status = (
        "GENERATED_LOCAL_GPC_GENESET_ENRICHMENT_NOT_TOPGO_PARITY"
        if checks["count_aligned_gpc_size_185"] and checks["tf_gene_set_tested"]
        else "GENERATED_LOCAL_GPC_GENESET_ENRICHMENT_WITH_VALIDATION_WARNINGS"
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "panel": "Figure 2E",
        "status": status,
        "truth_boundary": (
            "GPCs were derived from the public author GA/RNA correlation table using the paper rule "
            "(top decile correlation and >10 linked CREs), then count-aligned to 185 genes because the exact "
            "paper Table S2 / 1,999 variable-gene universe is not present in the current resource set. "
            "Enrichment uses local author-code gene sets with Fisher exact tests, not topGO v2.36.0 parity."
        ),
        "inputs": {
            "ga_rna_link_table": str(args.ga_rna_link_table),
            "gene_set_dir": str(args.gene_set_dir),
        },
        "outputs": {
            "out_dir": str(args.out_dir),
            "gpc_count_aligned_top185": str(gpc_path),
            "gpc_rule_candidates": str(rule_path),
            "local_geneset_enrichment": str(enrich_path),
            "figure_outputs": figure_outputs,
        },
        "metrics": {
            "universe_genes": len(universe),
            "correlation_top_decile_cutoff": corr_cutoff,
            "published_rule_candidate_count": int(rule.shape[0]),
            "count_aligned_gpc_count": int(count_aligned.shape[0]),
            "paper_gpc_count": PAPER_GPC_COUNT,
            "gene_sets_tested": len(gene_sets),
            "top_enrichment_terms": enrich.head(10).drop(columns=["overlap_genes"]).to_dict(orient="records")
            if not enrich.empty
            else [],
            "tf_best": tf_best,
        },
        "checks": checks,
    }
    summary_path = args.out_dir / "figure2e_public_local_gpc_enrichment_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
