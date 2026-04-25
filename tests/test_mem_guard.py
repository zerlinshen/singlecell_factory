import numpy as np
import pytest
from unittest.mock import patch

from workflow.modular._mem_guard import MemoryGuard


def _make_ctx(n_cells=100, n_genes=50):
    from anndata import AnnData
    import scipy.sparse as sp

    class _FakeCtx:
        pass

    ctx = _FakeCtx()
    ctx.adata = AnnData(sp.random(n_cells, n_genes, density=0.3, format="csr"))
    ctx.metadata = {}
    return ctx


# --- Decision table: 4 outcomes ---

def test_decision_go():
    mg = MemoryGuard()
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=10 * 1024**3):
        decision = mg.estimate_and_decide(100 * 1024**2)  # 100 MB vs 10 GB
    assert decision == "go"


def test_decision_chunk():
    mg = MemoryGuard()
    # planned = 60% of available → chunk (>= 50% threshold)
    avail = 1 * 1024**3
    planned = int(avail * 0.6)
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=avail):
        decision = mg.estimate_and_decide(planned)
    assert decision == "chunk"


def test_decision_abort():
    mg = MemoryGuard()
    avail = 1 * 1024**3
    planned = int(avail * 0.95)  # > 90% threshold
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=avail):
        decision = mg.estimate_and_decide(planned)
    assert decision == "abort"


def test_decision_backed_when_no_memory_info():
    mg = MemoryGuard()
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=None):
        decision = mg.estimate_and_decide(1024**3)
    assert decision == "backed"


# --- check() uses adata copy size ---

def test_check_go_with_small_adata():
    ctx = _make_ctx(10, 10)
    mg = MemoryGuard(ctx=ctx, label="test")
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=10 * 1024**3):
        result = mg.check()
    assert result == "go"


def test_check_no_adata():
    class _FakeCtx:
        adata = None
        metadata = {}

    mg = MemoryGuard(ctx=_FakeCtx(), label="nodata")
    assert mg.check() == "go"


def test_check_no_ctx():
    mg = MemoryGuard(ctx=None)
    assert mg.check() == "go"


# --- _quick_check static method ---

def test_quick_check_go():
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=4 * 1024**3):
        result = MemoryGuard._quick_check(10 * 1024**2)
    assert result == "go"


def test_quick_check_chunk():
    avail = 2 * 1024**3
    planned = int(avail * 0.7)
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=avail):
        result = MemoryGuard._quick_check(planned)
    assert result == "chunk"


def test_quick_check_abort():
    avail = 2 * 1024**3
    planned = int(avail * 0.95)
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=avail):
        result = MemoryGuard._quick_check(planned)
    assert result == "abort"


def test_quick_check_no_info():
    with patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=None):
        result = MemoryGuard._quick_check(1024**3)
    assert result == "go"
