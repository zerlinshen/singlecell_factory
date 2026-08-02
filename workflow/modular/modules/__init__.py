"""Optional and mandatory modules for modular workflow."""
from __future__ import annotations

from .._gene_symbols import resolve_expression_axis


def score_gene_sets(
    adata,
    gene_sets: dict[str, list[str]],
    prefix: str,
    *,
    use_raw: bool = False,
    min_genes: int = 2,
) -> list[str]:
    """Score multiple gene sets against an AnnData object.

    Returns the list of gene set names that were successfully scored
    (had at least *min_genes* present in the dataset).
    """
    import scanpy as sc

    # Membership and Scanpy must address the same matrix. ``use_raw=False`` is
    # the default because factory gene sets are symbol-keyed and ingest
    # normalises the live axis while preserving an Ensembl-indexed raw snapshot.
    axis = resolve_expression_axis(adata, use_raw=use_raw)
    var_names = axis.gene_set
    scored: list[str] = []
    for name, genes in gene_sets.items():
        valid = [g for g in genes if g in var_names]
        if len(valid) >= min_genes:
            sc.tl.score_genes(adata, valid, score_name=f"{prefix}_{name}", use_raw=use_raw)
            scored.append(name)
    return scored
