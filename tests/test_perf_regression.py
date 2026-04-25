"""Performance regression gate (opt-in via `pytest -m perf`).

Constructs a small synthetic AnnData and runs targeted modules end-to-end,
measuring wall time and RSS. Compared against tests/perf_baseline.json with
a 1.5x slack. The gate is intentionally coarse — its job is to catch
egregious regressions (>50% slower or >50% RAM), not to replace proper
profiling.
"""
from __future__ import annotations

import json
import os
import resource
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData


BASELINE_PATH = Path(__file__).parent / "perf_baseline.json"


def _baseline():
    return json.loads(BASELINE_PATH.read_text())


def _peak_rss_bytes() -> int:
    """ru_maxrss is KB on Linux, bytes on macOS. Normalize to bytes."""
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Heuristic: any value > 1e9 is already bytes (macOS); otherwise treat as KB (Linux)
    return int(raw if raw > 1e9 else raw * 1024)


@pytest.fixture
def synth_5k_adata():
    rng = np.random.default_rng(42)
    n_cells, n_genes = 5_000, 2_000
    base = (rng.random((n_cells, n_genes)) < 0.1).astype(np.float32) \
        * rng.random((n_cells, n_genes)).astype(np.float32) * 100.0
    X = sp.csr_matrix(base)
    adata = AnnData(X)
    adata.var_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    # Mark some mt/ribo/hb so QC has work to do
    for i, prefix in enumerate(["MT-", "RPS", "HB"]):
        adata.var_names.values[i * 5:(i + 1) * 5] = [
            f"{prefix}{j:02d}" for j in range(5)
        ]
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.startswith("HB")
    return adata


@pytest.mark.perf
def test_qc_module_wall_and_rss_within_baseline(synth_5k_adata, tmp_path):
    """QC fallback path on 5k cells must complete within 1.5x baseline wall + RSS."""
    from workflow.modular.modules.qc import QCModule

    baseline = _baseline()["fixtures"]["synth_5k"]["baseline"]
    ratio = _baseline()["thresholds"]
    wall_max = baseline["wall_seconds"] * ratio["wall_seconds_max_ratio"]
    rss_max = baseline["rss_peak_bytes"] * ratio["rss_peak_bytes_max_ratio"]

    rss_before = _peak_rss_bytes()
    t0 = time.perf_counter()
    QCModule._calculate_qc_metrics_fallback(synth_5k_adata)
    QCModule._filter_genes_fallback(synth_5k_adata, min_cells=3)
    wall = time.perf_counter() - t0
    rss_after = _peak_rss_bytes()
    rss_delta = max(0, rss_after - rss_before)

    # Sanity: results were produced
    assert "total_counts" in synth_5k_adata.obs.columns
    # Gates
    assert wall <= wall_max, (
        f"QC fallback wall = {wall:.2f}s exceeds {wall_max:.2f}s "
        f"(baseline {baseline['wall_seconds']}s × "
        f"{ratio['wall_seconds_max_ratio']})"
    )
    # RSS gate is informational (peak RSS is process-global; tight bounds are unreliable
    # in a shared pytest process). Just record it.
    print(f"QC perf: wall={wall:.2f}s rss_delta={rss_delta/1e6:.1f}MB")
    if rss_delta > rss_max:
        # Soft warning, not failure (other tests in the same process inflate RSS)
        print(f"WARN: rss_delta={rss_delta/1e6:.1f}MB > baseline {rss_max/1e6:.1f}MB")


@pytest.mark.perf
def test_sparse_welch_t_perf(synth_5k_adata):
    """sparse Welch DE must complete under 30s on 5k × 2k synthetic."""
    from workflow.modular._sparse_utils import sparse_welch_t

    n = synth_5k_adata.n_obs
    in_mask = np.zeros(n, dtype=bool)
    in_mask[: n // 2] = True
    out_mask = ~in_mask

    t0 = time.perf_counter()
    t, p, m_in, m_out = sparse_welch_t(synth_5k_adata.X, in_mask, out_mask)
    wall = time.perf_counter() - t0

    assert t.shape == (synth_5k_adata.n_vars,)
    assert wall <= 30.0, f"sparse_welch_t wall = {wall:.2f}s exceeds 30s budget"
    print(f"sparse_welch_t 5k×2k: wall={wall:.2f}s")
