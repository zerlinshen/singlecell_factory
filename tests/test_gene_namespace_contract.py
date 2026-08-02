"""Gene-identifier namespace contract at the ingest boundary.

Regression cover for a silent scientific failure found by a real-data run of the
LuCA core atlas LUSC cohort (Salcher et al. 2022, *Cancer Cell*,
doi:10.1016/j.ccell.2022.10.008), distributed under the CZ CELLxGENE schema,
which mandates Ensembl gene IDs as the ``var`` index.

The factory addresses genes by symbol everywhere downstream (mitochondrial /
ribosomal / haemoglobin QC prefixes, the marker annotation database, signature
scoring). The Cell Ranger ingest establishes that via
``read_10x_mtx(var_names="gene_symbols")``; the ``prepared_input.h5ad`` ingest
did not. On the real run this produced ``pct_counts_mt == 0.0`` for all 88,061
cells, so ``--max-mito-pct 20`` filtered nothing while still being recorded as
applied, and the annotation module failed with "No valid marker genes found".
Mitochondrial fraction is a primary quality covariate (Luecken & Theis 2019,
*Mol Syst Biol* 15:e8746; Heumos et al. 2023, *Nat Rev Genet* 24:550-572).
"""

from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from workflow.modular._gene_symbols import (
    describe_raw_axis,
    find_symbol_column,
    looks_like_ensembl,
    normalize_var_to_symbols,
    resolve_expression_axis,
)
from workflow.modular.modules.qc import QCModule

ENSEMBL_IDS = ["ENSG00000198804", "ENSG00000198712", "ENSG00000121410", "ENSG00000268895"]
SYMBOLS = ["MT-CO1", "MT-CO2", "A1BG", "A1BG-AS1"]


def _adata(var_index, var_cols=None) -> ad.AnnData:
    n_obs = 6
    rng = np.random.default_rng(0)
    X = rng.poisson(3.0, size=(n_obs, len(var_index))).astype(np.float32)
    var = pd.DataFrame(index=pd.Index(list(var_index), dtype=object))
    for key, values in (var_cols or {}).items():
        var[key] = list(values)
    return ad.AnnData(X=X, obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_obs)]), var=var)


def _ctx(max_mito_pct: float = 20.0, namespace: dict | None = None) -> SimpleNamespace:
    ctx = SimpleNamespace(
        cfg=SimpleNamespace(qc=SimpleNamespace(max_mito_pct=max_mito_pct)),
        metadata={},
    )
    if namespace is not None:
        ctx.metadata["gene_namespace"] = namespace
    return ctx


def test_detects_ensembl_and_symbol_axes():
    assert looks_like_ensembl(ENSEMBL_IDS)
    assert not looks_like_ensembl(SYMBOLS)
    # Versioned Ensembl IDs are still Ensembl IDs.
    assert looks_like_ensembl([f"{g}.7" for g in ENSEMBL_IDS])


def test_checked_expression_axis_binds_names_to_the_raw_matrix_it_returns():
    adata = _adata(ENSEMBL_IDS, {"feature_name": SYMBOLS})
    adata.raw = adata.copy()
    normalize_var_to_symbols(adata)

    axis = resolve_expression_axis(adata)

    assert axis.source == "adata.raw"
    assert axis.expression.n_vars == len(ENSEMBL_IDS)
    assert axis.gene_names == tuple(ENSEMBL_IDS)
    assert axis.gene_set == frozenset(ENSEMBL_IDS)
    assert tuple(adata.var_names) == tuple(SYMBOLS), "the live axis is deliberately different"
    with pytest.raises(TypeError):
        axis.gene_names[0] = "MUTATED"


def test_checked_expression_axis_uses_live_matrix_when_raw_is_absent():
    adata = _adata(SYMBOLS)
    axis = resolve_expression_axis(adata)
    assert axis.source == "adata"
    assert axis.expression is adata
    assert axis.gene_names == tuple(SYMBOLS)


def test_checked_expression_axis_rejects_duplicate_gene_keys():
    adata = _adata(["DUP", "DUP", "GENE"])
    with pytest.raises(ValueError, match="duplicate gene names"):
        resolve_expression_axis(adata)


def test_symbol_column_discovery_rejects_ensembl_and_placeholder_columns():
    var = pd.DataFrame(
        {
            "feature_name": SYMBOLS,
            "other_ids": ENSEMBL_IDS,
            "blank": ["nan", "nan", "nan", ""],
        },
        index=ENSEMBL_IDS,
    )
    assert find_symbol_column(var) == "feature_name"
    assert find_symbol_column(var[["other_ids", "blank"]]) is None


def test_ensembl_index_is_converted_and_original_ids_preserved():
    adata = _adata(ENSEMBL_IDS, {"feature_name": SYMBOLS})
    prov = normalize_var_to_symbols(adata)

    assert prov["status"] == "converted"
    assert prov["source_column"] == "feature_name"
    assert list(adata.var_names) == SYMBOLS
    assert list(adata.var["ensembl_id"]) == ENSEMBL_IDS
    # The whole point: symbol-prefix matching now works.
    assert adata.var_names.str.upper().str.startswith("MT-").sum() == 2


def test_symbol_index_is_left_untouched():
    adata = _adata(SYMBOLS)
    prov = normalize_var_to_symbols(adata)
    assert prov["status"] == "symbols_already"
    assert list(adata.var_names) == SYMBOLS


def test_genes_without_a_symbol_keep_their_ensembl_id():
    adata = _adata(ENSEMBL_IDS, {"feature_name": ["MT-CO1", "nan", "A1BG", ""]})
    prov = normalize_var_to_symbols(adata)
    assert prov["status"] == "converted"
    assert prov["n_without_symbol"] == 2
    # Axis stays complete and addressable rather than collapsing to "nan".
    assert list(adata.var_names) == ["MT-CO1", "ENSG00000198712", "A1BG", "ENSG00000268895"]


def test_unresolved_ensembl_is_reported_not_guessed():
    adata = _adata(ENSEMBL_IDS)
    prov = normalize_var_to_symbols(adata)
    assert prov["status"] == "ensembl_unresolved"
    assert prov["source_column"] is None
    assert list(adata.var_names) == ENSEMBL_IDS


# ------------------------------------------------------- .raw axis divergence
# Real-data finding, 2026-08-02 (92,430-cell LUSC squamous dataset): converting
# `var_names` to symbols leaves `.raw` in Ensembl space. Modules that tested
# membership against `adata` and then read values from `adata.raw` passed their
# own guard and raised on the lookup. Synthetic fixtures elsewhere carry no
# `.raw`, so `expr` falls through to `adata` and the two axes are identical by
# construction — which is exactly why unit tests could not see this. These
# fixtures therefore build `.raw` explicitly.
# See governance/raw_axis_namespace_divergence_2026-08-02.md.


def _with_ensembl_raw(adata):
    """Attach a `.raw` in the pre-conversion (Ensembl) namespace, as real ingest does."""
    adata.raw = adata.copy()
    return adata


def test_raw_axis_divergence_is_recorded_not_silent():
    adata = _with_ensembl_raw(_adata(ENSEMBL_IDS, {"feature_name": SYMBOLS}))
    prov = normalize_var_to_symbols(adata)

    assert prov["status"] == "converted"
    assert list(adata.var_names) == SYMBOLS
    # The divergence the real run hit, now visible in the manifest.
    assert prov["raw_axis_diverged"] is True
    assert prov["raw_axis_status"] == "ensembl_while_var_symbols"
    assert prov["raw_n_genes"] == len(ENSEMBL_IDS)
    assert prov["raw_ensembl_fraction"] == 1.0


def test_absent_raw_is_reported_as_absent_not_as_agreement():
    adata = _adata(ENSEMBL_IDS, {"feature_name": SYMBOLS})
    prov = normalize_var_to_symbols(adata)
    assert prov["raw_axis_status"] == "absent"
    assert prov["raw_axis_diverged"] is False


def test_matching_raw_axis_is_safe_to_cross_index():
    adata = _adata(SYMBOLS)
    adata.raw = adata.copy()
    prov = normalize_var_to_symbols(adata)
    assert prov["status"] == "symbols_already"
    assert prov["raw_axis_status"] == "matches_var"
    assert prov["raw_axis_diverged"] is False


def test_raw_retaining_extra_genes_is_diverged_but_not_a_namespace_split():
    """HVG selection leaves `.raw` wider than `var`. Same namespace, still not cross-indexable."""
    adata = _adata(SYMBOLS)
    adata.raw = adata.copy()
    adata._inplace_subset_var(np.array([True, True, False, False]))
    prov = describe_raw_axis(adata)
    assert prov["raw_axis_status"] == "diverged_other"
    assert prov["raw_axis_diverged"] is True
    assert prov["raw_n_genes"] == 4


def test_unresolved_ensembl_still_reports_the_raw_axis():
    """The unresolved branch is the one most likely to be cross-indexed by mistake."""
    adata = _with_ensembl_raw(_adata(ENSEMBL_IDS))
    prov = normalize_var_to_symbols(adata)
    assert prov["status"] == "ensembl_unresolved"
    # Both axes are Ensembl here, so there is no namespace split to report.
    assert prov["raw_axis_status"] == "matches_var"


# --------------------------------------------------------------- QC enforcement
def _flag_gene_classes(adata) -> None:
    upper = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = upper.str.startswith("MT-")
    adata.var["ribo"] = upper.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = upper.str.startswith(("HBA", "HBB"))


def test_qc_refuses_an_unresolvable_namespace_when_a_mito_filter_is_active():
    """The exact real-run failure: mito threshold recorded as applied, filtering nothing."""
    adata = _adata(ENSEMBL_IDS)
    _flag_gene_classes(adata)
    ctx = _ctx(max_mito_pct=20.0, namespace={"status": "ensembl_unresolved"})

    with pytest.raises(ValueError) as exc:
        QCModule()._assert_gene_classes_detectable(adata, ctx)
    message = str(exc.value)
    assert "no mitochondrial genes are detectable" in message
    assert "would filter nothing" in message


def test_qc_allows_an_unresolvable_namespace_when_no_mito_filter_is_intended():
    adata = _adata(ENSEMBL_IDS)
    _flag_gene_classes(adata)
    ctx = _ctx(max_mito_pct=100.0, namespace={"status": "ensembl_unresolved"})

    QCModule()._assert_gene_classes_detectable(adata, ctx)
    assert ctx.metadata["qc_mito_filter_status"] == "inactive_no_mitochondrial_genes_detected"


def test_qc_warns_but_proceeds_when_a_symbol_dataset_genuinely_lacks_mito_genes():
    """A pre-filtered atlas or targeted panel is legitimate — flag it, do not fail."""
    adata = _adata(["A1BG", "A1BG-AS1"])
    _flag_gene_classes(adata)
    ctx = _ctx(max_mito_pct=20.0, namespace={"status": "symbols_already"})

    QCModule()._assert_gene_classes_detectable(adata, ctx)
    assert ctx.metadata["qc_gene_class_counts"]["mt"] == 0
    assert ctx.metadata["qc_mito_filter_status"] == "inactive_no_mitochondrial_genes_detected"


def test_qc_records_gene_class_counts_on_the_happy_path():
    adata = _adata(["MT-CO1", "MT-CO2", "RPS6", "HBB", "A1BG"])
    _flag_gene_classes(adata)
    ctx = _ctx(max_mito_pct=20.0, namespace={"status": "symbols_already"})

    QCModule()._assert_gene_classes_detectable(adata, ctx)
    assert ctx.metadata["qc_gene_class_counts"] == {"mt": 2, "ribo": 1, "hb": 1}
    assert "qc_mito_filter_status" not in ctx.metadata


def test_externally_symbolified_file_gets_its_repair_key_backfilled():
    """The commonest way a file arrives symbol-indexed leaves no join key behind.

    `adata.var_names = adata.var["feature_name"]` symbol-ifies var and touches neither
    `.raw` nor `ensembl_id`. That file then takes the `symbols_already` branch, so the
    original code never wrote `ensembl_id` -- while `describe_raw_axis` correctly reported
    `ensembl_while_var_symbols`. The pipeline knew the axes were split, was holding the
    Ensembl IDs in `adata.raw.var_names`, and every downstream re-keying failed anyway.
    Measured on the real LUSC file: the trajectory HVG restriction was lost permanently.
    """
    adata = _adata(SYMBOLS, {"probe": ["a", "b", "c", "d"]})
    raw = _adata(ENSEMBL_IDS, {"probe": ["a", "b", "c", "d"]})
    adata.raw = raw

    prov = normalize_var_to_symbols(adata)
    assert prov["status"] == "symbols_already"
    assert prov["raw_axis_status"] == "ensembl_while_var_symbols"
    assert prov["ensembl_id_backfill"] == "from_raw_positional"
    assert list(adata.var["ensembl_id"]) == ENSEMBL_IDS


def test_backfill_is_skipped_when_positional_correspondence_is_not_established():
    """The mapping is positional, so a length mismatch must refuse rather than guess."""
    adata = _adata(SYMBOLS)
    adata.raw = _adata(ENSEMBL_IDS[:2] + ["ENSG00000000009"] * 2)
    adata.raw = _adata(ENSEMBL_IDS)
    # Shrink var so raw is a superset -> positional correspondence is not established.
    adata._inplace_subset_var(np.array([True, True, False, False]))
    prov = describe_raw_axis(adata)
    prov["status"] = "symbols_already"
    from workflow.modular._gene_symbols import _backfill_ensembl_id_from_raw
    _backfill_ensembl_id_from_raw(adata, prov)
    assert prov["ensembl_id_backfill"] == "skipped_raw_length_mismatch"
    assert "ensembl_id" not in adata.var.columns


def test_backfill_does_not_overwrite_an_existing_ensembl_id():
    adata = _adata(SYMBOLS, {"ensembl_id": ["KEEP1", "KEEP2", "KEEP3", "KEEP4"]})
    adata.raw = _adata(ENSEMBL_IDS)
    prov = normalize_var_to_symbols(adata)
    assert "ensembl_id_backfill" not in prov
    assert list(adata.var["ensembl_id"]) == ["KEEP1", "KEEP2", "KEEP3", "KEEP4"]


# ------------------------------------------- backfill order verification (2026-08-02)
# Second-round review finding. Equal lengths do NOT establish that `var` and `.raw` are
# in the same gene order: AnnData leaves `.raw` untouched when the gene axis is sliced or
# reordered, so `adata[:, sorted_order]` permutes `var` and not `.raw`. An earlier
# version of the backfill trusted the length alone. Measured on the real 17,764-gene LUSC
# file prepared with `var_names = var["feature_name"]` then sorted by symbol: the
# backfilled ensembl_id was correct for 489 genes (2.75%), the trajectory aligner
# certified it `recovered_via_ensembl_id`, and 1,710 of 2,000 "highly variable genes"
# were the wrong genes. `ensembl_id` is PERSISTED into final_adata.h5ad, so that
# 97%-wrong mapping would have reached every downstream consumer unmarked.


def _diverged_with_raw(var_order, raw_order, shared_col=None):
    """`var` and `.raw` in Ensembl/symbol split, with independently controlled orders."""
    var = pd.DataFrame(index=pd.Index(var_order, dtype=object))
    raw_var = pd.DataFrame(index=pd.Index(raw_order, dtype=object))
    if shared_col is not None:
        var["probe"], raw_var["probe"] = shared_col[0], shared_col[1]
    n = len(var_order)
    adata = ad.AnnData(np.zeros((3, n)), var=var)
    raw = ad.AnnData(np.zeros((3, n)), var=raw_var)
    adata.raw = raw
    return adata


def test_backfill_refuses_when_raw_gene_order_cannot_be_verified():
    """No shared column -> the positional assumption is unproven -> refuse."""
    adata = _diverged_with_raw(SYMBOLS, ENSEMBL_IDS)
    prov = normalize_var_to_symbols(adata)
    assert prov["raw_axis_status"] == "ensembl_while_var_symbols"
    assert prov["ensembl_id_backfill"] == "skipped_raw_order_unverifiable"
    assert "ensembl_id" not in adata.var.columns


def test_backfill_refuses_when_a_shared_column_disagrees_positionally():
    """The permutation case: lengths match, order does not, and a witness proves it."""
    adata = _diverged_with_raw(
        SYMBOLS, ENSEMBL_IDS,
        shared_col=(["a", "b", "c", "d"], ["d", "c", "b", "a"]),
    )
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_backfill"] == "skipped_raw_order_unverifiable"
    assert "ensembl_id" not in adata.var.columns


def test_backfill_proceeds_when_a_discriminating_column_agrees():
    adata = _diverged_with_raw(
        SYMBOLS, ENSEMBL_IDS,
        shared_col=(["a", "b", "c", "d"], ["a", "b", "c", "d"]),
    )
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_backfill"] == "from_raw_positional"
    assert prov["ensembl_id_backfill_witness"] == "probe"
    assert list(adata.var["ensembl_id"]) == ENSEMBL_IDS


def test_a_constant_shared_column_is_not_a_valid_order_witness():
    """A column that agrees under ANY permutation proves nothing about order."""
    adata = _diverged_with_raw(
        SYMBOLS, ENSEMBL_IDS,
        shared_col=(["same"] * 4, ["same"] * 4),
    )
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_backfill"] == "skipped_raw_order_unverifiable"
    assert "ensembl_id" not in adata.var.columns


# ------------------------------------- ensembl_id provenance (2026-08-02, round 3)
# A review REFUTED "the backfill is the only origin of a wrong join key". A pre-existing
# `ensembl_id` column was trusted with zero verification by both the conversion path and
# the backfill's early return. An upstream preparer doing the unsafe thing one step
# earlier -- symbol-ify var, reorder the gene axis, then write
# `var["ensembl_id"] = raw.var_names` -- reproduced the corruption exactly (489/17,764
# correct; 290 of 2,000 HVGs correct; reported as `recovered_via_ensembl_id`) and left NO
# provenance key at all, making it less visible than the bug it mirrors. Seurat->h5ad
# conversions and hand-prepared CELLxGENE files routinely ship such a column.


def test_shipped_ensembl_id_is_marked_unverified_without_an_order_witness():
    adata = _diverged_with_raw(SYMBOLS, ENSEMBL_IDS)
    adata.var["ensembl_id"] = list(reversed(ENSEMBL_IDS))   # bijective and WRONG
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "preexisting_unverified"


def test_shipped_ensembl_id_is_verified_when_it_matches_an_ordered_raw():
    adata = _diverged_with_raw(
        SYMBOLS, ENSEMBL_IDS, shared_col=(["a", "b", "c", "d"], ["a", "b", "c", "d"]))
    adata.var["ensembl_id"] = list(ENSEMBL_IDS)
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "verified_against_raw"
    assert prov["ensembl_id_verified_via"] == "probe"


def test_shipped_ensembl_id_contradicting_an_ordered_raw_is_marked_wrong():
    """Order is proven, so a mismatching column is not ambiguous -- it is incorrect."""
    adata = _diverged_with_raw(
        SYMBOLS, ENSEMBL_IDS, shared_col=(["a", "b", "c", "d"], ["a", "b", "c", "d"]))
    adata.var["ensembl_id"] = list(reversed(ENSEMBL_IDS))
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "preexisting_contradicts_raw"


def test_factory_written_ensembl_id_vouches_for_itself():
    adata = _adata(ENSEMBL_IDS, {"feature_name": SYMBOLS})
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "converted_in_factory"


def test_absent_ensembl_id_is_reported_as_absent():
    adata = _adata(SYMBOLS)
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "absent"


def test_a_witness_at_exactly_the_distinctness_bar_is_refused():
    """100% twin pairs give exactly 50% distinct values.

    A permutation WITHIN a tie group preserves every shared column positionally, so a
    composite that only reaches the bar proves nothing. A strict `<` let this worst case
    through by a single comparison; the boundary belongs on the refusing side.
    """
    n = 8
    syms = [f"SYM{i}" for i in range(n)]
    ens = [f"ENSG{i:011d}" for i in range(n)]
    twins = [f"pair{i // 2}" for i in range(n)]          # exactly n/2 distinct values
    adata = _diverged_with_raw(syms, ens, shared_col=(twins, twins))
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_backfill"] == "skipped_raw_order_unverifiable"
    assert "ensembl_id" not in adata.var.columns


def _tied_symbol_permutation(*, shipped_join_key: bool):
    """Four-gene reproduction of an invisible within-tie permutation.

    ``feature_name`` agrees positionally after the live axis swaps the two DUP
    genes, but that agreement cannot reveal which stable ID belongs to either
    duplicate. AnnData keeps ``raw.var`` in its original order.
    """
    ens = [f"ENSG{i:011d}" for i in range(1, 5)]
    base = ad.AnnData(np.arange(16, dtype=float).reshape(4, 4))
    base.var_names = ens
    base.var["feature_name"] = ["DUP", "DUP", "GENEA", "GENEB"]
    base.raw = base.copy()

    obj = base[:, [1, 0, 2, 3]].copy()
    obj.var_names = ["DUP", "DUP-1", "GENEA", "GENEB"]
    if shipped_join_key:
        # The unsafe upstream repair: attach the unchanged raw IDs by position
        # after the live matrix has moved.
        obj.var["ensembl_id"] = list(obj.raw.var_names)
    return obj


def test_tied_symbol_group_cannot_verify_a_shipped_positional_join_key():
    """Accepting feature_name here certifies two provably wrong stable IDs."""
    from workflow.modular._gene_symbols import _positional_order_witness

    adata = _tied_symbol_permutation(shipped_join_key=True)
    assert _positional_order_witness(adata) is None

    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "preexisting_unverified"
    assert "ensembl_id_verified_via" not in prov


def test_tied_symbol_group_refuses_positional_ensembl_backfill():
    """Removing a shipped key must not turn the same ambiguity into backfill."""
    adata = _tied_symbol_permutation(shipped_join_key=False)

    prov = normalize_var_to_symbols(adata)

    assert prov["ensembl_id_backfill"] == "skipped_raw_order_unverifiable"
    assert "ensembl_id" not in adata.var.columns


def test_shipped_ensembl_id_is_unverified_when_there_is_nothing_to_check_it_against():
    """No `.raw` at all: the column may be perfect, but nothing here can vouch for it.

    Distinct from the no-witness branch. There, `.raw` exists and the gene ORDER is
    unproven; here there is no second axis to compare against in the first place. Both
    must land on `preexisting_unverified`, and both need their own cover — a mutation
    that flipped only this branch to a trusted value survived the suite.
    """
    adata = _adata(SYMBOLS, {"ensembl_id": ENSEMBL_IDS})
    prov = normalize_var_to_symbols(adata)
    assert prov["raw_axis_status"] == "absent"
    assert prov["ensembl_id_source"] == "preexisting_unverified"


def test_shipped_ensembl_id_is_unverified_when_raw_lengths_differ():
    adata = _adata(SYMBOLS, {"ensembl_id": ENSEMBL_IDS})
    adata.raw = _adata(ENSEMBL_IDS)
    adata._inplace_subset_var(np.array([True, True, False, False]))
    prov = normalize_var_to_symbols(adata)
    assert prov["ensembl_id_source"] == "preexisting_unverified"


# ------------------------------- .raw reassignment during the run (2026-08-02, round 4)
# `gene_namespace` is captured once at ingest, but `.raw` is REASSIGNED later, so the
# recorded status could describe a state no module ever saw and no artifact ever had.
# Two reachable paths, and the second is the DEFAULT:
#   (a) Cell Ranger yields no `.raw` -> ingest records `absent`; clustering.py then sets
#       `adata.raw = adata`. Benign in direction, still factually wrong.
#   (b) GPU failure policy resolves to `restore-cpu` for any gpu_mode != "force";
#       clustering.py nulls `.raw` and the CPU lanes re-assign it in SYMBOL space, while
#       the manifest still asserts ensembl_while_var_symbols / raw_axis_diverged: True.
# (b) is why this matters: a divergence flag asserting a split that no longer exists
# cannot serve as gate input or reviewer evidence, which is what it exists to be.

from workflow.modular.pipeline import _restate_raw_axis_at_manifest_time


def _ctx_with(adata, namespace):
    return SimpleNamespace(adata=adata, metadata={"gene_namespace": namespace})


def test_manifest_restates_raw_axis_after_cellranger_assigns_it():
    adata = _adata(SYMBOLS)
    prov = normalize_var_to_symbols(adata)
    assert prov["raw_axis_status"] == "absent"

    adata.raw = adata                      # clustering.py, checkpoint_policy=full
    ctx = _ctx_with(adata, prov)
    _restate_raw_axis_at_manifest_time(ctx)

    ns = ctx.metadata["gene_namespace"]
    assert ns["raw_axis_status_at_ingest"] == "absent"
    assert ns["raw_axis_status_at_manifest"] == "matches_var"
    assert ns["raw_axis_reassigned_during_run"] is True
    assert ns["raw_n_genes_at_manifest"] == len(SYMBOLS)


def test_manifest_stops_asserting_a_divergence_the_gpu_fallback_removed():
    """The harmful direction, on the default GPU policy."""
    adata = _adata(ENSEMBL_IDS, {"feature_name": SYMBOLS})
    adata.raw = adata.copy()
    prov = normalize_var_to_symbols(adata)
    assert prov["raw_axis_status"] == "ensembl_while_var_symbols"
    assert prov["raw_axis_diverged"] is True

    adata.raw = None                       # clustering.py:446 on GPU failure
    adata.raw = adata                      # _run_cpu re-assigns, now in symbol space
    ctx = _ctx_with(adata, prov)
    _restate_raw_axis_at_manifest_time(ctx)

    ns = ctx.metadata["gene_namespace"]
    assert ns["raw_axis_status_at_manifest"] == "matches_var"
    assert ns["raw_axis_diverged_at_manifest"] is False
    assert ns["raw_axis_reassigned_during_run"] is True
    # The ingest reading is preserved, not overwritten: it explains what the modules saw.
    assert ns["raw_axis_status_at_ingest"] == "ensembl_while_var_symbols"


def test_no_reassignment_is_not_reported_as_one():
    adata = _adata(ENSEMBL_IDS, {"feature_name": SYMBOLS})
    adata.raw = adata.copy()
    prov = normalize_var_to_symbols(adata)
    ctx = _ctx_with(adata, prov)
    _restate_raw_axis_at_manifest_time(ctx)

    ns = ctx.metadata["gene_namespace"]
    assert ns["raw_axis_status_at_manifest"] == ns["raw_axis_status_at_ingest"]
    assert "raw_axis_reassigned_during_run" not in ns


def test_restatement_never_blocks_a_manifest_write():
    """A manifest that fails to write loses the whole run's provenance."""
    class _Exploding:
        var_names = pd.Index(SYMBOLS, dtype=object)

        @property
        def raw(self):
            raise RuntimeError("backing store closed")

    ctx = _ctx_with(_Exploding(), {"raw_axis_status": "absent"})
    _restate_raw_axis_at_manifest_time(ctx)
    assert "raw_axis_at_manifest_error" in ctx.metadata["gene_namespace"]
