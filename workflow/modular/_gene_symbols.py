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


def describe_raw_axis(adata) -> dict[str, Any]:
    """Describe ``adata.raw``'s gene axis relative to ``adata.var_names``.

    ``normalize_var_to_symbols`` rewrites ``adata.var_names`` but deliberately
    leaves ``adata.raw`` alone: ``.raw`` is an immutable snapshot of the
    pre-processing state, and rewriting it would change ingest semantics for
    every module and every already-prepared artifact. The consequence is that
    on canonically-ingested CELLxGENE data the two axes address genes in
    *different namespaces*, so a module that tests membership against ``adata``
    and then reads values out of ``adata.raw`` passes its own guard and raises
    on the lookup.

    That divergence is legitimate but must not be invisible. This records it so
    the manifest, the modules, and the reviewer can all see it.

    ``raw_axis_status``
        ``absent``        — no ``.raw``; ``expr`` falls through to ``adata``.
        ``matches_var``   — same namespace as ``var_names``; safe to cross-index.
        ``ensembl_while_var_symbols`` — the divergence described above.
        ``diverged_other`` — axes differ for some other reason (e.g. ``.raw``
        retains genes dropped by HVG selection).
    """
    provenance: dict[str, Any] = {}
    raw = getattr(adata, "raw", None)
    if raw is None:
        provenance["raw_axis_status"] = "absent"
        provenance["raw_axis_diverged"] = False
        return provenance

    raw_names = list(raw.var_names)
    provenance["raw_n_genes"] = len(raw_names)

    if list(adata.var_names) == raw_names:
        provenance["raw_axis_status"] = "matches_var"
        provenance["raw_axis_diverged"] = False
        return provenance

    raw_ensembl = round(ensembl_fraction(raw_names), 4)
    provenance["raw_ensembl_fraction"] = raw_ensembl
    provenance["raw_axis_diverged"] = True
    if raw_ensembl >= _ENSEMBL_FRACTION_THRESHOLD and not looks_like_ensembl(adata.var_names):
        provenance["raw_axis_status"] = "ensembl_while_var_symbols"
        logger.warning(
            "GENE_NAMESPACE: adata.raw is Ensembl-indexed (%.1f%%) while adata.var_names "
            "are symbols. Symbol-keyed lookups against .raw WILL fail. Modules must test "
            "membership against the axis they index (see governance/"
            "raw_axis_namespace_divergence_2026-08-02.md).",
            100 * raw_ensembl,
        )
    else:
        provenance["raw_axis_status"] = "diverged_other"
    return provenance


def _backfill_ensembl_id_from_raw(adata, provenance: dict[str, Any]) -> None:
    """Recover ``var["ensembl_id"]`` when ``var`` was symbol-ified outside this factory.

    ``normalize_var_to_symbols`` writes ``ensembl_id`` only when IT does the conversion.
    But the commonest way an h5ad arrives already symbol-indexed is the one-liner
    ``adata.var_names = adata.var["feature_name"]``, which touches neither ``.raw`` nor
    ``ensembl_id``. That file then takes the ``symbols_already`` branch, and
    ``describe_raw_axis`` correctly reports ``ensembl_while_var_symbols`` — so the
    pipeline KNOWS the axes are split, is holding the Ensembl IDs in
    ``adata.raw.var_names``, and yet every downstream re-keying attempt fails for want of
    the column. Measured on the real LUSC file prepared that way: the trajectory aligner
    refused and the HVG restriction was lost permanently, with the fix in place.

    The mapping is POSITIONAL, so it is only taken when ``.raw`` has exactly as many
    genes as ``var`` — the shape produced by the standard ``adata.raw = adata.copy()`` at
    ingest, before any subsetting. That assumption is recorded in the provenance rather
    than hidden, because a reordered ``.raw`` would make it wrong and nothing here can
    detect that.
    """
    import pandas as pd

    if provenance.get("raw_axis_status") != "ensembl_while_var_symbols":
        return
    if "ensembl_id" in adata.var.columns:
        return
    if provenance.get("raw_n_genes") != int(adata.n_vars):
        provenance["ensembl_id_backfill"] = "skipped_raw_length_mismatch"
        return

    adata.var["ensembl_id"] = pd.Index(adata.raw.var_names).astype(str)
    provenance["ensembl_id_backfill"] = "from_raw_positional"
    logger.warning(
        "GENE_NAMESPACE: var_names were already symbols but .raw is Ensembl-indexed and "
        "var['ensembl_id'] was absent, so downstream re-keying had no join key. Backfilled "
        "ensembl_id from adata.raw.var_names POSITIONALLY (%d genes, lengths match). This "
        "assumes .raw preserves var's gene order, which holds for the standard "
        "`adata.raw = adata.copy()` ingest; recorded as ensembl_id_backfill in provenance.",
        int(adata.n_vars),
    )


def normalize_var_to_symbols(adata) -> dict[str, Any]:
    """Make ``adata.var_names`` gene symbols when the input is Ensembl-indexed.

    Mutates ``adata`` in place. Returns a provenance dict suitable for
    ``ctx.metadata`` so the manifest records exactly what happened:

    ``status``
        ``symbols_already`` — index was already symbol-like, nothing changed.
        ``converted``       — index was Ensembl, symbols applied from ``source_column``.
        ``ensembl_unresolved`` — index is Ensembl but no symbol column exists.

    The returned dict always also carries the ``raw_axis_*`` keys from
    ``describe_raw_axis`` — see that function for why ``.raw`` is not rewritten.
    """
    import pandas as pd

    provenance: dict[str, Any] = {
        "n_genes": int(adata.n_vars),
        "ensembl_fraction": round(ensembl_fraction(adata.var_names), 4),
    }

    if not looks_like_ensembl(adata.var_names):
        provenance["status"] = "symbols_already"
        provenance["source_column"] = None
        provenance.update(describe_raw_axis(adata))
        _backfill_ensembl_id_from_raw(adata, provenance)
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
        provenance.update(describe_raw_axis(adata))
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
    provenance.update(describe_raw_axis(adata))
    logger.info(
        "GENE_NAMESPACE: converted %d Ensembl-indexed genes to symbols from "
        "var[%r] (%d had no symbol and kept their Ensembl ID); original IDs "
        "preserved in var['ensembl_id'].",
        adata.n_vars, col, int(blank.sum()),
    )
    return provenance
