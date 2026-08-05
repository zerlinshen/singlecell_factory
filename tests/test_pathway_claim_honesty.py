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
