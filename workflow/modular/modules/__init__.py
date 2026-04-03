"""Optional and mandatory modules for modular workflow."""
from __future__ import annotations

import scanpy as sc


def score_gene_sets(
    adata,
    gene_sets: dict[str, list[str]],
    prefix: str,
    *,
    use_raw: bool = True,
    min_genes: int = 2,
) -> list[str]:
    """Score multiple gene sets against an AnnData object.

    Returns the list of gene set names that were successfully scored
    (had at least *min_genes* present in the dataset).
    """
    var_names = set(adata.var_names if adata.raw is None else adata.raw.var_names)
    scored: list[str] = []
    for name, genes in gene_sets.items():
        valid = [g for g in genes if g in var_names]
        if len(valid) >= min_genes:
            sc.tl.score_genes(adata, valid, score_name=f"{prefix}_{name}", use_raw=use_raw)
            scored.append(name)
    return scored
