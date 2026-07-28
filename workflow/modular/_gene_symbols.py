"""Gene-identifier namespace resolution for the ingest boundary.

Why this exists
---------------
The entire downstream factory addresses genes by HGNC-style **symbol**: QC flags
mitochondrial/ribosomal/haemoglobin genes with ``startswith("MT-"/"RPS"/"RPL"/
"HBA"/"HBB")``, the annotation marker database is symbol-keyed, and cell-cycle,
signature-scoring and immune-phenotyping gene sets are all symbol-based. The
Cell Ranger ingest path establishes that contract explicitly via
``sc.read_10x_mtx(..., var_names="gene_symbols")``.

The ``prepared_input.h5ad`` / ``prepared_input.zarr`` ingest path did **not**
enforce it. That matters because the CZ CELLxGENE Discover schema — the most
common public distribution format for published atlases — *mandates* Ensembl
gene IDs as ``var`` index and carries the symbol in ``var["feature_name"]``.
Feeding such a file in produced a silent scientific failure rather than an
error: every ``startswith("MT-")`` test returned False, so ``pct_counts_mt``,
``pct_counts_ribo`` and ``pct_counts_hb`` were 0.0 for every cell and the
corresponding QC thresholds filtered nothing while still being reported as
applied. Mitochondrial fraction is a primary QC covariate for detecting stressed
and dying cells (Luecken & Theis 2019, *Mol Syst Biol* 15:e8746; Heumos et al.
2023, *Nat Rev Genet* 24:550-572), so a silently inert mitochondrial filter
admits low-quality cells into every downstream result.

This module restores namespace parity between the two ingest paths: when the
``var`` index is Ensembl-like and a symbol column is present, ``var_names``
become symbols (uniquified, exactly as ``read_10x_mtx`` does) and the original
Ensembl IDs are preserved in ``var["ensembl_id"]``. Nothing is guessed: if no
usable symbol column exists the index is left untouched and the caller is told,
so QC can fail loudly instead of silently measuring nothing.
"""

from __future__ import annotations

import logging
import re
from typing import Any

logger = logging.getLogger(__name__)

# ENSG00000121410 / ENSMUSG00000051951 / ENSG00000121410.5 (versioned)
_ENSEMBL_RE = re.compile(r"^ENS[A-Z]{0,6}[GTP]\d{6,}(\.\d+)?$", re.IGNORECASE)

# Ordered by how authoritative the column is. "feature_name" is the CELLxGENE
# schema field; the others cover Cell Ranger h5 and common community exports.
_SYMBOL_COLUMN_CANDIDATES: tuple[str, ...] = (
    "feature_name",
    "gene_symbol",
    "gene_symbols",
    "gene_name",
    "gene_names",
    "symbol",
    "SYMBOL",
    "hgnc_symbol",
)

# Below this fraction of Ensembl-shaped entries we treat the axis as symbols.
_ENSEMBL_FRACTION_THRESHOLD = 0.5


def ensembl_fraction(values) -> float:
    """Fraction of ``values`` that look like Ensembl gene/transcript IDs."""
    import pandas as pd

    s = pd.Index(values).astype(str)
    if len(s) == 0:
        return 0.0
    return float(s.to_series().str.match(_ENSEMBL_RE).mean())


def looks_like_ensembl(values) -> bool:
    """True when the majority of ``values`` are Ensembl-shaped identifiers."""
    return ensembl_fraction(values) >= _ENSEMBL_FRACTION_THRESHOLD


def find_symbol_column(var) -> str | None:
    """Return the first ``var`` column that plausibly holds gene symbols."""
    import pandas as pd

    for col in _SYMBOL_COLUMN_CANDIDATES:
        if col not in var.columns:
            continue
        values = pd.Series(var[col]).astype(str)
        if values.empty:
            continue
        # Reject placeholder columns and columns that are themselves Ensembl IDs.
        missing = values.str.lower().isin({"nan", "none", "", "na"}).mean()
        if missing > 0.5:
            continue
        if ensembl_fraction(values) >= _ENSEMBL_FRACTION_THRESHOLD:
            continue
        return col
    return None


def normalize_var_to_symbols(adata) -> dict[str, Any]:
    """Make ``adata.var_names`` gene symbols when the input is Ensembl-indexed.

    Mutates ``adata`` in place. Returns a provenance dict suitable for
    ``ctx.metadata`` so the manifest records exactly what happened:

    ``status``
        ``symbols_already`` — index was already symbol-like, nothing changed.
        ``converted``       — index was Ensembl, symbols applied from ``source_column``.
        ``ensembl_unresolved`` — index is Ensembl but no symbol column exists.
    """
    import pandas as pd

    provenance: dict[str, Any] = {
        "n_genes": int(adata.n_vars),
        "ensembl_fraction": round(ensembl_fraction(adata.var_names), 4),
    }

    if not looks_like_ensembl(adata.var_names):
        provenance["status"] = "symbols_already"
        provenance["source_column"] = None
        return provenance

    col = find_symbol_column(adata.var)
    if col is None:
        provenance["status"] = "ensembl_unresolved"
        provenance["source_column"] = None
        provenance["candidate_columns_checked"] = list(_SYMBOL_COLUMN_CANDIDATES)
        provenance["available_var_columns"] = [str(c) for c in adata.var.columns]
        logger.warning(
            "GENE_NAMESPACE: var_names look like Ensembl IDs (%.1f%%) but no symbol "
            "column was found among %s. Symbol-keyed logic (mitochondrial/ribosomal "
            "QC flags, marker annotation, signature scoring) CANNOT match and will "
            "be reported as inactive. Available var columns: %s",
            100 * provenance["ensembl_fraction"],
            list(_SYMBOL_COLUMN_CANDIDATES),
            [str(c) for c in adata.var.columns],
        )
        return provenance

    # Preserve the original identifiers before overwriting the index.
    if "ensembl_id" not in adata.var.columns:
        adata.var["ensembl_id"] = pd.Index(adata.var_names).astype(str)

    symbols = pd.Series(adata.var[col]).astype(str)
    # Fall back to the Ensembl ID for genes with no symbol, so the axis stays
    # complete and addressable rather than collapsing to "nan".
    blank = symbols.str.lower().isin({"nan", "none", "", "na"}).to_numpy()
    symbols = symbols.mask(blank, pd.Series(adata.var["ensembl_id"]).astype(str).to_numpy())

    adata.var_names = pd.Index(symbols.to_numpy(), dtype=object)
    adata.var_names_make_unique()

    provenance["status"] = "converted"
    provenance["source_column"] = col
    provenance["n_without_symbol"] = int(blank.sum())
    logger.info(
        "GENE_NAMESPACE: converted %d Ensembl-indexed genes to symbols from "
        "var[%r] (%d had no symbol and kept their Ensembl ID); original IDs "
        "preserved in var['ensembl_id'].",
        adata.n_vars, col, int(blank.sum()),
    )
    return provenance
