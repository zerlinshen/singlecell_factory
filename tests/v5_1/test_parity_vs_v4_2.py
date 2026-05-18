from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
PARITY = REPO_ROOT / "scripts" / "v5_1" / "parity_vs_v4_2.py"


def _load_parity():
    spec = importlib.util.spec_from_file_location("wave5_parity_vs_v4_2", PARITY)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _adata_with_pairs(pairs: list[tuple[int, str]]):
    return SimpleNamespace(
        uns={
            "peak_to_gene_top1000": pd.DataFrame(
                {"peak_idx": [p for p, _ in pairs], "gene": [g for _, g in pairs]}
            )
        }
    )


def test_peak_gene_replay_baseline_is_jaccard_one_not_trevino_overlap():
    parity = _load_parity()

    assert parity.V42_BASELINES["peak_gene_overlap_top1000_strict"] == 1.0


def test_compute_peak_gene_overlap_returns_one_for_identical_rankings():
    parity = _load_parity()
    pairs = [(i, f"gene{i}") for i in range(1000)]

    overlap = parity._compute_peak_gene_overlap(_adata_with_pairs(pairs), _adata_with_pairs(pairs))

    assert overlap == 1.0


def test_compute_peak_gene_overlap_uses_true_jaccard_for_partial_overlap():
    parity = _load_parity()
    left = [(i, f"gene{i}") for i in range(1000)]
    right = [(i, f"gene{i}") for i in range(500)] + [
        (10_000 + i, f"other{i}") for i in range(500)
    ]

    overlap = parity._compute_peak_gene_overlap(_adata_with_pairs(left), _adata_with_pairs(right))

    assert overlap == round(500 / 1500, 4)
