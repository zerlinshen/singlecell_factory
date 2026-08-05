"""Tests for stage LIANA contrast summarizer (Wave-5)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.liana_stage_contrast_summary import jaccard, main, pair_key


def test_jaccard_and_pair_key():
    df = pd.DataFrame(
        {
            "ligand_complex": ["A", "B"],
            "receptor_complex": ["X", "Y"],
            "source": ["T cell", "NK cell"],
            "target": ["Tumor epithelial", "T cell"],
        }
    )
    s = pair_key(df)
    assert len(s) == 2
    assert jaccard(s, s) == 1.0
    assert jaccard(s, set()) == 0.0


def test_main_writes_summary(tmp_path: Path):
    rows = pd.DataFrame(
        {
            "source": ["T cell", "Myeloid/Macro"],
            "target": ["NK cell", "T cell"],
            "ligand_complex": ["HLA-DRA", "CD80"],
            "receptor_complex": ["LAG3", "CTLA4"],
            "magnitude_rank": [0.1, 0.2],
        }
    )
    p1 = tmp_path / "s1.csv"
    p2 = tmp_path / "s2.csv"
    rows.to_csv(p1, index=False)
    rows.iloc[[0]].to_csv(p2, index=False)
    out = tmp_path / "out.json"
    rc = main(["--stage-dir", f"I={p1}", "--stage-dir", f"II={p2}", "--output", str(out)])
    assert rc == 0
    assert out.is_file()
