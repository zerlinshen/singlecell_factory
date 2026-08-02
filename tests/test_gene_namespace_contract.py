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
