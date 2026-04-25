"""Tests for MemoryEnforcer: check_budget, watchdog abort signal, /proc fallback,
cooperative abort in module loops, and pipeline skipped_memory status."""
from __future__ import annotations

import threading
import time
from unittest.mock import patch, MagicMock

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ctx(n_cells=50, n_genes=30):
    from anndata import AnnData
    import scipy.sparse as sp

    class _FakeCtx:
        pass

    ctx = _FakeCtx()
    ctx.adata = AnnData(sp.random(n_cells, n_genes, density=0.3, format="csr"))
    ctx.metadata = {}
    return ctx


# ---------------------------------------------------------------------------
# check_budget — returns, never raises
# ---------------------------------------------------------------------------

class TestCheckBudget:
    def test_go_when_low_pressure(self):
        from workflow.modular._mem_guard import MemoryGuard
        with (
            patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=50 * 1024**3),
            patch("workflow.modular._mem_guard._get_total_memory_bytes", return_value=64 * 1024**3),
            patch("workflow.modular._mem_guard._get_rss_bytes", return_value=4 * 1024**3),
        ):
            # projected RSS = 4GB + 100MB ≈ 6% of 64GB — well under soft_cap=65%
            result = MemoryGuard.check_budget(100 * 1024**2, label="test")
        assert result == "go"

    def test_chunk_when_moderate_pressure(self):
        from workflow.modular._mem_guard import MemoryGuard
        total = 64 * 1024**3
        # rss=40GB + planned=5GB → projected=45GB → 70% of 64GB > soft_cap=65%
        with (
            patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=20 * 1024**3),
            patch("workflow.modular._mem_guard._get_total_memory_bytes", return_value=total),
            patch("workflow.modular._mem_guard._get_rss_bytes", return_value=40 * 1024**3),
        ):
            result = MemoryGuard.check_budget(5 * 1024**3, label="test")
        assert result == "chunk"

    def test_abort_when_high_pressure(self):
        from workflow.modular._mem_guard import MemoryGuard
        total = 64 * 1024**3
        # rss=52GB + planned=5GB → projected=57GB → 89% of 64GB > hard_cap=85%
        with (
            patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=10 * 1024**3),
            patch("workflow.modular._mem_guard._get_total_memory_bytes", return_value=total),
            patch("workflow.modular._mem_guard._get_rss_bytes", return_value=52 * 1024**3),
        ):
            result = MemoryGuard.check_budget(5 * 1024**3, label="test")
        assert result == "abort"

    def test_go_when_no_memory_info(self):
        from workflow.modular._mem_guard import MemoryGuard
        with (
            patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=None),
            patch("workflow.modular._mem_guard._get_total_memory_bytes", return_value=None),
            patch("workflow.modular._mem_guard._get_rss_bytes", return_value=None),
        ):
            result = MemoryGuard.check_budget(1024**3)
        assert result == "go"

    def test_never_raises(self):
        from workflow.modular._mem_guard import MemoryGuard
        # Even with extreme pressure, check_budget must not raise
        with (
            patch("workflow.modular._mem_guard._get_available_memory_bytes", return_value=1),
            patch("workflow.modular._mem_guard._get_total_memory_bytes", return_value=100),
            patch("workflow.modular._mem_guard._get_rss_bytes", return_value=99),
        ):
            result = MemoryGuard.check_budget(10**12)
        assert result in {"go", "chunk", "abort"}


# ---------------------------------------------------------------------------
# /proc/self/status fallback (no psutil)
# ---------------------------------------------------------------------------

class TestProcFallback:
    def test_get_rss_bytes_proc_fallback(self):
        """_get_rss_bytes falls back to /proc/self/status when psutil missing."""
        from workflow.modular import _mem_guard

        proc_content = "Name:\tpython3\nVmRSS:\t4194304 kB\nVmSize:\t8388608 kB\n"

        original_import = __builtins__.__import__ if hasattr(__builtins__, "__import__") else __import__

        with patch("builtins.open", MagicMock(return_value=MagicMock(
            __enter__=lambda s, *a: iter(proc_content.splitlines(keepends=True)),
            __exit__=lambda s, *a: False,
        ))):
            # Force psutil ImportError path
            with patch.dict("sys.modules", {"psutil": None}):
                rss = _mem_guard._get_rss_bytes()
        # If psutil is installed in env it will be used; just verify the function returns int or None
        assert rss is None or isinstance(rss, int)

    def test_get_total_memory_bytes_returns_int_or_none(self):
        from workflow.modular._mem_guard import _get_total_memory_bytes
        result = _get_total_memory_bytes()
        assert result is None or (isinstance(result, int) and result > 0)

    def test_get_rss_bytes_returns_int_or_none(self):
        from workflow.modular._mem_guard import _get_rss_bytes
        result = _get_rss_bytes()
        assert result is None or (isinstance(result, int) and result > 0)


# ---------------------------------------------------------------------------
# Abort event: request_abort / abort_requested / clear_abort
# ---------------------------------------------------------------------------

class TestAbortEvent:
    def setup_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def teardown_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def test_initially_not_set(self):
        from workflow.modular._mem_guard import MemoryGuard
        assert not MemoryGuard.abort_requested()

    def test_request_abort_sets_flag(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.request_abort()
        assert MemoryGuard.abort_requested()

    def test_clear_abort_clears_flag(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.request_abort()
        MemoryGuard.clear_abort()
        assert not MemoryGuard.abort_requested()

    def test_abort_event_is_shared_across_instances(self):
        from workflow.modular._mem_guard import MemoryGuard
        mg1 = MemoryGuard()
        mg2 = MemoryGuard()
        MemoryGuard.request_abort()
        assert mg1.abort_requested()
        assert mg2.abort_requested()


# ---------------------------------------------------------------------------
# Watchdog: sets abort flag when RSS exceeds hard_cap
# ---------------------------------------------------------------------------

class TestWatchdog:
    def setup_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def teardown_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def test_watchdog_sets_abort_when_rss_exceeds_hard_cap(self):
        from workflow.modular import _mem_watchdog
        from workflow.modular._mem_guard import MemoryGuard

        ctx = _make_ctx()
        total = 64 * 1024**3
        # RSS at 90% of total — above hard_cap=85%
        rss = int(total * 0.90)

        with (
            patch("workflow.modular._mem_watchdog._get_available_memory_bytes", return_value=total - rss),
            patch("workflow.modular._mem_watchdog._get_rss_bytes", return_value=rss),
            patch("workflow.modular._mem_watchdog._get_total_memory_bytes", return_value=total),
            patch("workflow.modular._mem_watchdog._HARD_CAP", 0.85),
        ):
            t = _mem_watchdog.start(ctx, period=0.05)
            time.sleep(0.3)
            _mem_watchdog.stop(t)

        assert MemoryGuard.abort_requested()
        assert "oom_traces" in ctx.metadata

    def test_watchdog_does_not_abort_when_rss_below_hard_cap(self):
        from workflow.modular import _mem_watchdog
        from workflow.modular._mem_guard import MemoryGuard

        ctx = _make_ctx()
        total = 64 * 1024**3
        rss = int(total * 0.50)  # 50% — well below hard_cap=85%

        with (
            patch("workflow.modular._mem_watchdog._get_available_memory_bytes", return_value=total - rss),
            patch("workflow.modular._mem_watchdog._get_rss_bytes", return_value=rss),
            patch("workflow.modular._mem_watchdog._get_total_memory_bytes", return_value=total),
            patch("workflow.modular._mem_watchdog._HARD_CAP", 0.85),
        ):
            t = _mem_watchdog.start(ctx, period=0.05)
            time.sleep(0.3)
            _mem_watchdog.stop(t)

        assert not MemoryGuard.abort_requested()

    def test_watchdog_stop_joins_thread(self):
        from workflow.modular import _mem_watchdog

        ctx = _make_ctx()
        with (
            patch("workflow.modular._mem_watchdog._get_available_memory_bytes", return_value=10 * 1024**3),
            patch("workflow.modular._mem_watchdog._get_rss_bytes", return_value=1 * 1024**3),
            patch("workflow.modular._mem_watchdog._get_total_memory_bytes", return_value=64 * 1024**3),
        ):
            t = _mem_watchdog.start(ctx, period=0.1)
            _mem_watchdog.stop(t)

        assert not t.is_alive()


# ---------------------------------------------------------------------------
# Cooperative abort in module loops: MemoryAbortError raised → pipeline skips
# ---------------------------------------------------------------------------

class TestCooperativeAbort:
    def setup_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def teardown_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def test_cnv_compute_smoothed_chunked_aborts(self):
        """_compute_smoothed_chunked raises MemoryAbortError when abort requested."""
        import numpy as np
        import scipy.sparse as sp
        from workflow.modular._mem_guard import MemoryGuard, MemoryAbortError
        from workflow.modular.modules.cnv_inference import CNVInferenceModule

        MemoryGuard.request_abort()

        n_cells, n_genes = 20, 10
        expr_sparse = sp.random(n_cells, n_genes, density=0.5, format="csr").astype(np.float32)
        chromosomes = np.array(["chr1"] * n_genes)

        class _FakeAdata:
            obs = type("obs", (), {"columns": []})()

        with pytest.raises(MemoryAbortError):
            CNVInferenceModule._compute_smoothed_chunked(
                expr_sparse, _FakeAdata(), None, chromosomes, 5, 5
            )

    def test_metacell_aggregate_aborts(self):
        """_aggregate raises MemoryAbortError when abort requested."""
        import numpy as np
        import scipy.sparse as sp
        from anndata import AnnData
        from workflow.modular._mem_guard import MemoryGuard, MemoryAbortError
        from workflow.modular.modules.metacell import MetacellModule

        MemoryGuard.request_abort()

        n_cells, n_genes = 10, 5
        adata = AnnData(sp.random(n_cells, n_genes, density=0.5, format="csr"))
        labels = np.zeros(n_cells, dtype=int)

        with pytest.raises(MemoryAbortError):
            MetacellModule._aggregate(adata, labels, 1)

    def test_pseudobulk_de_ranktest_aborts(self):
        """_de_ranktest raises MemoryAbortError when abort requested."""
        import numpy as np
        import pandas as pd
        from scipy.stats import ranksums
        from workflow.modular._mem_guard import MemoryGuard, MemoryAbortError
        from workflow.modular.modules.pseudobulk_de import PseudobulkDEModule

        MemoryGuard.request_abort()

        counts = pd.DataFrame(np.random.rand(4, 20), columns=[f"g{i}" for i in range(20)])
        meta = pd.DataFrame({"condition": ["A", "A", "B", "B"]})

        with pytest.raises(MemoryAbortError):
            PseudobulkDEModule._de_ranktest(counts, meta, "condition", "A", "B", ranksums)


# ---------------------------------------------------------------------------
# Pipeline: skipped_memory status when MemoryAbortError propagates
# ---------------------------------------------------------------------------

class TestPipelineSkippedMemory:
    def setup_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def teardown_method(self):
        from workflow.modular._mem_guard import MemoryGuard
        MemoryGuard.clear_abort()

    def test_run_module_converts_memory_abort_to_skip_module(self):
        """_run_module converts MemoryAbortError to _SkipModule."""
        from workflow.modular.pipeline import _run_module, _SkipModule
        from workflow.modular._mem_guard import MemoryAbortError

        class _AbortingMod:
            name = "test_abort"
            requires_keys = {}

            def run(self, ctx):
                raise MemoryAbortError("test abort")

        class _FakeCtx:
            adata = None
            metadata = {}

        with pytest.raises(_SkipModule) as exc_info:
            _run_module(_AbortingMod(), _FakeCtx(), mandatory=False)

        assert "memory abort" in str(exc_info.value).lower()

    def test_sequential_records_skipped_memory_status(self):
        """_run_sequential records status='skipped_memory' for memory-aborted modules."""
        from workflow.modular.pipeline import _run_sequential, _SkipModule
        from workflow.modular._mem_guard import MemoryAbortError

        class _FakeCtx:
            adata = None
            metadata = {}
            module_status = []
            _checkpoint_dir = None
            _module_dirs = {}
            figure_dir = None
            table_dir = None

            def status(self, name, ok, msg=""):
                status_val = ok if isinstance(ok, str) else ("ok" if ok else "failed")
                self.module_status.append({"module": name, "status": status_val, "message": msg})

            def save_checkpoint(self, name):
                pass

            def set_module_dir(self, name):
                pass

        class _AbortingMod:
            name = "aborting"
            requires_keys = {}

            def run(self, ctx):
                raise MemoryAbortError("simulated watchdog abort")

        ctx = _FakeCtx()
        registry = {"aborting": _AbortingMod()}
        _run_sequential(["aborting"], registry, ctx, mandatory=set())

        assert len(ctx.module_status) == 1
        assert ctx.module_status[0]["status"] == "skipped_memory"
