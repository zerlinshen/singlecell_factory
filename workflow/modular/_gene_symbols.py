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
from dataclasses import dataclass
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


@dataclass(frozen=True)
class ExpressionAxis:
    """An expression object coupled to the immutable gene keys used to index it."""

    expression: Any
    gene_names: tuple[str, ...]
    source: str

    @property
    def gene_set(self) -> frozenset[str]:
        return frozenset(self.gene_names)


def resolve_expression_axis(adata, *, use_raw: bool = True) -> ExpressionAxis:
    """Return the selected expression matrix and the exact gene index it owns.

    Downstream modules historically selected ``adata.raw`` and then independently read
    ``adata.var_names``. Those axes legitimately diverge after ingest normalises the live
    axis to symbols. This is the single runtime boundary for that choice: callers receive
    the matrix and an immutable tuple captured from *that matrix*, never a free-floating
    name set that could have come from another AnnData object.
    """
    raw = getattr(adata, "raw", None)
    expression = raw.to_adata() if use_raw and raw is not None else adata
    source = "adata.raw" if use_raw and raw is not None else "adata"

    gene_names = tuple(str(name) for name in expression.var_names)
    n_vars = int(expression.n_vars)
    matrix_shape = getattr(expression.X, "shape", None)
    if len(gene_names) != n_vars or matrix_shape is None or len(matrix_shape) != 2 \
            or int(matrix_shape[1]) != n_vars:
        raise RuntimeError(
            f"Expression-axis contract violated for {source}: "
            f"n_vars={n_vars}, names={len(gene_names)}, X.shape={matrix_shape}."
        )

    seen: set[str] = set()
    duplicates: set[str] = set()
    for name in gene_names:
        if name in seen:
            duplicates.add(name)
        else:
            seen.add(name)
    if duplicates:
        raise ValueError(
            f"Expression-axis contract violated for {source}: duplicate gene names "
            f"cannot form an unambiguous lookup index (examples: {sorted(duplicates)[:5]})."
        )

    if int(expression.n_obs) != int(adata.n_obs) or not expression.obs_names.equals(adata.obs_names):
        raise RuntimeError(
            f"Expression-axis contract violated for {source}: observation rows are not "
            "aligned to adata, so cell masks cannot be applied safely."
        )

    return ExpressionAxis(expression=expression, gene_names=gene_names, source=source)


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


# `var["ensembl_id"]` is the join key every cross-namespace repair depends on, and it is
# PERSISTED into final_adata.h5ad. Only these origins are load-bearing evidence that it is
# correct; anything else is a string column of unknown provenance that happens to have the
# right name.
TRUSTED_ENSEMBL_ID_SOURCES = frozenset({
    "converted_in_factory",            # this module wrote it during Ensembl->symbol conversion
    "backfilled_from_raw_positional",  # this module wrote it from an order-verified `.raw`
    "verified_against_raw",            # shipped by the preparer, and it matches an order-verified `.raw`
})


def classify_ensembl_id(adata, provenance: dict[str, Any]) -> None:
    """Record WHERE ``var["ensembl_id"]`` came from, and whether it can be trusted.

    A review refuted the assumption that the backfill is the only origin of a wrong join
    key. A pre-existing column was trusted with no verification at all, by both the
    conversion path and the backfill's early return — and an upstream preparer doing the
    same unsafe thing one step earlier (symbol-ify, reorder, then write
    ``var["ensembl_id"] = raw.var_names``) reproduced the corruption EXACTLY: 489 of
    17,764 correct on the real file, 290 of 2,000 HVGs correct, aligner reporting
    ``recovered_via_ensembl_id``. Worse than the bug it mirrors, because that one at
    least left ``ensembl_id_backfill`` behind; this path left no trace whatsoever.

    Seurat->h5ad conversions and hand-prepared CELLxGENE files routinely ship an
    ``ensembl_id`` / ``gene_ids`` column, so the unverified path is the LIKELIER one.
    """
    import pandas as pd

    if "ensembl_id" not in adata.var.columns:
        provenance["ensembl_id_source"] = "absent"
        return
    if provenance.get("ensembl_id_source"):
        return  # this module wrote it on this call and already said so

    if (provenance.get("raw_axis_status") != "ensembl_while_var_symbols"
            or provenance.get("raw_n_genes") != int(adata.n_vars)):
        # Nothing to check it against.
        provenance["ensembl_id_source"] = "preexisting_unverified"
        return

    witness = _positional_order_witness(adata)
    if witness is None:
        provenance["ensembl_id_source"] = "preexisting_unverified"
        return

    shipped = pd.Series(adata.var["ensembl_id"]).astype(str).to_numpy()
    from_raw = pd.Index(adata.raw.var_names).astype(str).to_numpy()
    if shipped.shape == from_raw.shape and (shipped == from_raw).all():
        provenance["ensembl_id_source"] = "verified_against_raw"
        provenance["ensembl_id_verified_via"] = witness
        return

    provenance["ensembl_id_source"] = "preexisting_contradicts_raw"
    logger.warning(
        "GENE_NAMESPACE: var['ensembl_id'] was shipped with this file, but the gene ORDER "
        "of var and raw.var is confirmed identical (via %r) and the column does NOT match "
        "adata.raw.var_names. The column is therefore wrong. It is left in place but "
        "marked untrusted, so cross-namespace repairs will refuse it rather than join "
        "through it.", witness,
    )


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

    The mapping is POSITIONAL, and equal lengths are NOT sufficient to justify it.
    **AnnData does not reorder ``.raw`` when the gene axis is sliced or reordered**, so
    ``adata[:, sorted_order]`` leaves ``var`` permuted and ``.raw`` untouched, with the
    lengths still equal. An earlier version of this function trusted the length alone.
    Measured on the real LUSC file prepared with two ubiquitous one-liners
    (``var_names = var["feature_name"]``, then sort the gene axis by symbol): the
    backfilled column was correct for **489 of 17,764 genes — 2.75%** — and, because
    ``ensembl_id`` is PERSISTED into ``final_adata.h5ad``, that 97%-wrong mapping would
    have reached every downstream consumer with nothing marking it.

    So the order is verified, not assumed: a column carried by BOTH ``var`` and
    ``raw.var`` must agree positionally. Without such a column, or when it disagrees,
    the backfill is refused and the reason recorded. Refusing costs the HVG restriction;
    guessing corrupts the gene identity of the whole run.
    """
    import pandas as pd

    if provenance.get("raw_axis_status") != "ensembl_while_var_symbols":
        return
    if "ensembl_id" in adata.var.columns:
        return
    if provenance.get("raw_n_genes") != int(adata.n_vars):
        provenance["ensembl_id_backfill"] = "skipped_raw_length_mismatch"
        return

    witness = _positional_order_witness(adata)
    if witness is None:
        provenance["ensembl_id_backfill"] = "skipped_raw_order_unverifiable"
        logger.warning(
            "GENE_NAMESPACE: .raw is Ensembl-indexed and var['ensembl_id'] is absent, but "
            "no column shared by var and raw.var confirms they are in the same gene "
            "ORDER. Equal lengths do not establish that -- AnnData leaves .raw untouched "
            "when the gene axis is reordered. Refusing to backfill rather than persist a "
            "possibly-wrong ensembl_id into the run artifact."
        )
        return

    adata.var["ensembl_id"] = pd.Index(adata.raw.var_names).astype(str)
    provenance["ensembl_id_backfill"] = "from_raw_positional"
    provenance["ensembl_id_source"] = "backfilled_from_raw_positional"
    provenance["ensembl_id_backfill_witness"] = witness
    logger.info(
        "GENE_NAMESPACE: backfilled var['ensembl_id'] from adata.raw.var_names "
        "positionally (%d genes); gene ORDER verified via the shared column %r.",
        int(adata.n_vars), witness,
    )


# A witness must distinguish EVERY ordering. Any tied composite key leaves a permutation
# within that tie group invisible, even if 99.99% of the remaining genes are unique.
# Positional stable-ID assignment is all-or-nothing, so majority distinctness is not a
# proof threshold.


def _positional_order_witness(adata) -> str | None:
    """Proof that ``var`` and ``raw.var`` are in the same gene order, or ``None``.

    Compares every column carried by BOTH frames, as one composite key, positionally.
    On the real file ``raw.var`` carries ``feature_name`` — the detector that catches the
    permuted-order corruption instantly (positional agreement 489/17,764 when permuted,
    17,764/17,764 when intact).

    Returns the composite's description so the provenance records WHAT was checked, not
    merely that something was.
    """
    import pandas as pd

    raw = getattr(adata, "raw", None)
    if raw is None:
        return None
    try:
        raw_var = raw.var
    except (AttributeError, ValueError):  # pragma: no cover - defensive
        return None

    shared = [c for c in raw_var.columns
              if c in adata.var.columns and len(adata.var) == len(raw_var)]
    if not shared:
        return None

    agreeing = []
    for col in shared:
        left = pd.Series(adata.var[col]).astype(str).to_numpy()
        right = pd.Series(raw_var[col]).astype(str).to_numpy()
        if left.shape == right.shape and (left == right).all():
            agreeing.append(str(col))
    if not agreeing:
        return None

    composite = pd.Series(
        ["\x1f".join(vals) for vals in
         zip(*(pd.Series(adata.var[c]).astype(str).tolist() for c in agreeing))]
    )
    if not len(composite):
        return None
    tied = composite.duplicated(keep=False)
    if bool(tied.any()):
        n_tied_positions = int(tied.sum())
        n_tied_keys = int(composite[tied].nunique())
        logger.warning(
            "GENE_NAMESPACE: columns %s agree positionally between var and raw.var, but "
            "their composite key has %d tied values covering %d of %d genes. Any "
            "permutation within a tied group preserves the witness, so it cannot prove "
            "the two frames are in the same gene ORDER. Treating the order as "
            "unverified.",
            agreeing, n_tied_keys, n_tied_positions, len(composite),
        )
        return None
    return "+".join(agreeing)


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
        classify_ensembl_id(adata, provenance)
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
        classify_ensembl_id(adata, provenance)
        return provenance

    # Preserve the original identifiers before overwriting the index.
    if "ensembl_id" not in adata.var.columns:
        adata.var["ensembl_id"] = pd.Index(adata.var_names).astype(str)
        provenance["ensembl_id_source"] = "converted_in_factory"

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
    classify_ensembl_id(adata, provenance)
    logger.info(
        "GENE_NAMESPACE: converted %d Ensembl-indexed genes to symbols from "
        "var[%r] (%d had no symbol and kept their Ensembl ID); original IDs "
        "preserved in var['ensembl_id'].",
        adata.n_vars, col, int(blank.sum()),
    )
    return provenance
