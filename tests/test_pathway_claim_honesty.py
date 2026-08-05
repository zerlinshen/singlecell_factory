"""Claim honesty: ORA must not be stamped as rank-based GSEA."""
from __future__ import annotations

import pandas as pd
import pytest

from workflow.modular.modules.pathway_analysis import PathwayAnalysisModule


class _Ctx:
    random_state = 0


def test_run_gseapy_ora_engine_tag_when_enrich_only(monkeypatch):
    """If prerank is unavailable/fails and enrich succeeds, engine is gseapy_ora."""
    mod = PathwayAnalysisModule()

    class _FakeEnr:
        def __init__(self):
            self.results = pd.DataFrame(
                {"Term": ["HALLMARK_HYPOXIA"], "Adjusted P-value": [0.01]}
            )

    class _FakeGP:
        def prerank(self, *a, **k):
            raise RuntimeError("force ORA path")

        def enrich(self, *a, **k):
            return _FakeEnr()

    import sys
    import types

    fake = types.ModuleType("gseapy")
    fake.prerank = _FakeGP().prerank
    fake.enrich = _FakeGP().enrich
    monkeypatch.setitem(sys.modules, "gseapy", fake)

    de = pd.DataFrame(
        {
            "group": ["0"] * 20,
            "names": [f"G{i}" for i in range(20)],
            "scores": list(range(20, 0, -1)),
        }
    )
    results, engine = mod._run_gseapy(de, adata=None, ctx=_Ctx())
    assert engine == "gseapy_ora"
    assert results is not None and not results.empty
    assert "term" in results.columns or "Term" in results.columns


def test_inference_map_has_ora_distinct_from_gsea():
    """Source contract: ORA and GSEA must carry distinct inference classes."""
    from pathlib import Path

    mod_path = Path(
        "/home/zerlinshen/Bioinformatics Research Pipeline/"
        "singlecell_factory/workflow/modular/modules/pathway_analysis.py"
    )
    text = mod_path.read_text()
    assert '"gseapy_ora"' in text or "'gseapy_ora'" in text
    assert "over_representation_analysis" in text
    assert "supported_ora_not_gsea" in text
    assert "gseapy_mixed" in text
    # Must not map gseapy_ora to rank_based_gsea
    assert 'gseapy_ora": ("rank_based_gsea"' not in text


def test_run_gseapy_prerank_engine_when_scores_ok(monkeypatch):
    """Successful prerank path stamps gseapy_gsea (not ORA)."""
    mod = PathwayAnalysisModule()

    class _FakePrerank:
        def __init__(self):
            self.res2d = pd.DataFrame(
                {"Term": ["HALLMARK_HYPOXIA"], "FDR q-val": [0.01]}
            )

    class _FakeGP:
        def prerank(self, *a, **k):
            return _FakePrerank()

        def enrich(self, *a, **k):
            raise AssertionError("enrich must not run when prerank succeeds")

    import sys
    import types

    fake = types.ModuleType("gseapy")
    fake.prerank = _FakeGP().prerank
    fake.enrich = _FakeGP().enrich
    monkeypatch.setitem(sys.modules, "gseapy", fake)

    de = pd.DataFrame(
        {
            "group": ["0"] * 20,
            "names": [f"G{i}" for i in range(20)],
            "scores": list(range(20, 0, -1)),
        }
    )
    results, engine = mod._run_gseapy(de, adata=None, ctx=_Ctx())
    assert engine == "gseapy_gsea"
    assert results is not None and not results.empty


def test_run_gseapy_mixed_when_one_cluster_prerank_one_ora(monkeypatch):
    """Mixed cluster backends must not advertise pure GSEA."""
    mod = PathwayAnalysisModule()
    calls = {"prerank": 0, "enrich": 0}

    class _FakePrerank:
        res2d = pd.DataFrame({"Term": ["HALLMARK_MYC"], "FDR q-val": [0.02]})

    class _FakeEnr:
        results = pd.DataFrame(
            {"Term": ["HALLMARK_HYPOXIA"], "Adjusted P-value": [0.03]}
        )

    class _FakeGP:
        def prerank(self, *a, **k):
            calls["prerank"] += 1
            # Fail on second cluster so ORA is used for group "1"
            if calls["prerank"] > 1:
                raise RuntimeError("prerank fail cluster1")
            return _FakePrerank()

        def enrich(self, *a, **k):
            calls["enrich"] += 1
            return _FakeEnr()

    import sys
    import types

    fake = types.ModuleType("gseapy")
    fake.prerank = _FakeGP().prerank
    fake.enrich = _FakeGP().enrich
    monkeypatch.setitem(sys.modules, "gseapy", fake)

    de = pd.DataFrame(
        {
            "group": ["0"] * 20 + ["1"] * 20,
            "names": [f"G{i}" for i in range(20)] * 2,
            "scores": list(range(20, 0, -1)) * 2,
        }
    )
    results, engine = mod._run_gseapy(de, adata=None, ctx=_Ctx())
    assert engine == "gseapy_mixed"
    assert results is not None
    assert set(results["pathway_method"]) == {"prerank_gsea", "ora_enrich"}


def test_run_metadata_stamps_ora_not_confirmatory_gsea(monkeypatch, tmp_path):
    """Full run() path: gseapy_ora must not set rank_based_gsea claim."""
    from pathlib import Path
    import types

    mod = PathwayAnalysisModule()

    class _FakeEnr:
        results = pd.DataFrame(
            {"Term": ["HALLMARK_HYPOXIA"], "Adjusted P-value": [0.01]}
        )

    class _FakeGP:
        def prerank(self, *a, **k):
            raise RuntimeError("force ora")

        def enrich(self, *a, **k):
            return _FakeEnr()

    import sys

    fake = types.ModuleType("gseapy")
    fake.prerank = _FakeGP().prerank
    fake.enrich = _FakeGP().enrich
    monkeypatch.setitem(sys.modules, "gseapy", fake)

    de_dir = tmp_path / "de"
    de_dir.mkdir()
    pd.DataFrame(
        {
            "group": ["0"] * 20,
            "names": [f"G{i}" for i in range(20)],
            "scores": list(range(20, 0, -1)),
        }
    ).to_csv(de_dir / "marker_genes.csv", index=False)

    class _FullCtx:
        def __init__(self):
            self.adata = types.SimpleNamespace()  # truthy, unused by gseapy path
            self.metadata = {}
            self.random_state = 0
            self.table_dir = tmp_path / "tables"
            self.table_dir.mkdir()
            self.figure_dir = tmp_path / "figures"
            self.figure_dir.mkdir()

        def module_output_dir(self, name):
            return de_dir if name == "differential_expression" else None

    monkeypatch.setattr(
        PathwayAnalysisModule, "_plot_enrichment", lambda self, *a, **k: None
    )
    ctx = _FullCtx()
    mod.run(ctx)
    assert ctx.metadata["pathway_engine"] == "gseapy_ora"
    assert ctx.metadata["pathway_inference_class"] == "over_representation_analysis"
    assert ctx.metadata["pathway_inference_status"] == "supported_ora_not_gsea"
    assert ctx.metadata["pathway_claimable"] is True
