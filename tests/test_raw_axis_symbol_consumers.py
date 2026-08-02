"""Execution regressions for symbol-keyed consumers on split live/raw axes."""

from __future__ import annotations

import json
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from workflow.modular._gene_symbols import resolve_expression_axis
from workflow.modular.modules import score_gene_sets
from workflow.modular.modules.cell_communication import CellCommunicationModule
from workflow.modular.modules.gene_signature_scoring import GeneSignatureScoringModule
from workflow.modular.modules.immune_phenotyping import ImmunePhenotypingModule
from workflow.modular.modules.tumor_microenvironment import TumorMicroenvironmentModule


SYMBOLS = [
    "CD274",
    "PDCD1",
    "GZMA",
    "PRF1",
    "CD8A",
    "GZMB",
    "NKG7",
    "GNLY",
    "IFNG",
    "LAG3",
    "HAVCR2",
    "TIGIT",
    "CTLA4",
    "TOX",
]
ENSEMBL = [f"ENSG{i:011d}" for i in range(1, len(SYMBOLS) + 1)]


def _split_axis_adata() -> ad.AnnData:
    n_cells = 12
    obs = pd.DataFrame(
        {"cell_type": ["T cell"] * 6 + ["Tumor epithelial"] * 6},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    raw = ad.AnnData(
        X=np.full((n_cells, len(SYMBOLS)), 2.0, dtype=np.float32),
        obs=obs.copy(),
        var=pd.DataFrame(index=ENSEMBL),
    )
    adata = raw.copy()
    adata.raw = raw
    adata.var_names = SYMBOLS
    adata.var["ensembl_id"] = ENSEMBL
    return adata


def _ctx(tmp_path, adata: ad.AnnData, **cfg) -> SimpleNamespace:
    return SimpleNamespace(
        adata=adata,
        table_dir=tmp_path,
        figure_dir=tmp_path,
        metadata={},
        cfg=SimpleNamespace(**cfg),
    )


def _patch_score_genes(monkeypatch) -> None:
    import scanpy as sc

    def fake_score_genes(adata, genes, *, score_name, use_raw):
        assert use_raw is False
        assert set(genes) <= set(adata.var_names)
        adata.obs[score_name] = float(len(genes))

    monkeypatch.setattr(sc.tl, "score_genes", fake_score_genes)


def test_expression_axis_can_select_the_live_symbol_matrix() -> None:
    adata = _split_axis_adata()

    axis = resolve_expression_axis(adata, use_raw=False)

    assert axis.expression is adata
    assert axis.source == "adata"
    assert axis.gene_names == tuple(SYMBOLS)


def test_score_gene_sets_filters_the_same_live_axis_it_scores(monkeypatch) -> None:
    adata = _split_axis_adata()
    _patch_score_genes(monkeypatch)

    scored = score_gene_sets(
        adata,
        {"cyt": ["GZMA", "PRF1"]},
        "sig",
        use_raw=False,
    )

    assert scored == ["cyt"]
    assert "sig_cyt" in adata.obs


def test_manual_cell_communication_executes_on_live_symbols(tmp_path, monkeypatch) -> None:
    adata = _split_axis_adata()
    ctx = _ctx(tmp_path, adata)
    monkeypatch.setattr(CellCommunicationModule, "_plot_lr_heatmap", lambda *a: None)

    CellCommunicationModule()._run_manual_lr(adata, ctx)

    result = pd.read_csv(tmp_path / "cell_communication_lr.csv")
    assert {"CD274", "PDCD1"} <= set(result[["ligand", "receptor"]].stack())


def test_tme_and_immune_consumers_execute_on_live_symbols(tmp_path, monkeypatch) -> None:
    adata = _split_axis_adata()
    _patch_score_genes(monkeypatch)

    tme = TumorMicroenvironmentModule()
    monkeypatch.setattr(tme, "_plot_visualizations", lambda *a: None)
    tme.run(_ctx(tmp_path, adata))
    assert "cyt_score" in adata.obs
    assert (tmp_path / "checkpoint_expression.csv").is_file()

    immune = ImmunePhenotypingModule()
    monkeypatch.setattr(immune, "_plot_visualizations", lambda *a: None)
    immune.run(_ctx(tmp_path, adata))
    assert "immune_subtype" in adata.obs
    assert (tmp_path / "immune_phenotyping.csv").is_file()


def test_user_gene_signature_executes_on_live_symbols(tmp_path, monkeypatch) -> None:
    adata = _split_axis_adata()
    _patch_score_genes(monkeypatch)
    contract = tmp_path / "signatures.json"
    contract.write_text(json.dumps({"cyt": ["GZMA", "PRF1"]}), encoding="utf-8")
    cfg = SimpleNamespace(
        use_builtin=False,
        signature_json=contract,
    )
    module = GeneSignatureScoringModule()
    monkeypatch.setattr(module, "_plot_visualizations", lambda *a: None)

    module.run(_ctx(tmp_path, adata, gene_signature=cfg))

    assert "sig_cyt" in adata.obs
    assert (tmp_path / "gene_signature_scores.csv").is_file()
