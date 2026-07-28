"""Clustering thread count is a recorded, explicit part of the run record.

Why this file exists: ``sc.settings.n_jobs`` is not a speed knob for clustering,
it is an INPUT TO THE NEIGHBOUR GRAPH. scanpy hands it to pynndescent's
transformer, which wraps its NN-descent fit in ``numba.set_num_threads()``, and
NN-descent seeds one RNG state per thread. Measured on this suite's own LUSC PCA
(20,000 cells x 40 PCs, n_neighbors=15, random_state=42, one host, one input):

    n_jobs= 1 -> 29 clusters
    n_jobs= 8 -> 31 clusters
    n_jobs=32 -> 29 clusters, ARI 0.931 against n_jobs=1

Clustering used to set ``sc.settings.n_jobs = os.cpu_count()``, so the cluster
count was a property of the machine while ``run_manifest.json`` recorded only
``clustering_random_state`` — the run was not reproducible from its own record.

These tests pin the two properties that fix it: the resolved count is recorded,
and the default does not track the host's core count.
"""
from __future__ import annotations

import os

import numpy as np
import pytest
from anndata import AnnData

from workflow.modular.config import CellRangerConfig, PipelineConfig
from workflow.modular.context import PipelineContext
import workflow.modular.modules.clustering as clustering_mod
from workflow.modular.modules.clustering import (
    _DEFAULT_CLUSTERING_N_JOBS,
    _N_JOBS_ENV,
    ClusteringModule,
)


def _fake_host(monkeypatch, cores: int) -> None:
    """Pretend the run is on a `cores`-core machine.

    Both knobs matter: the old code read ``os.cpu_count()``, and the clamp reads
    numba's launch-time cap, so a test host is only convincing if it moves both.
    """
    monkeypatch.setattr(os, "cpu_count", lambda: cores)
    try:
        import numba

        monkeypatch.setattr(numba.config, "NUMBA_NUM_THREADS", cores, raising=False)
    except ImportError:  # pragma: no cover — numba ships with scanpy in this env
        pass


def _run_clustering(monkeypatch, tmp_path, cores: int) -> PipelineContext:
    """Run ClusteringModule on a fake host with the compute lanes stubbed out.

    Only the thread-count resolution is under test, so _run_cpu just labels the
    cells; the point is what lands in ctx.metadata and in sc.settings.n_jobs.
    """
    _fake_host(monkeypatch, cores)

    def fake_run_cpu(self, adata, cfg, ctx):
        adata.obs["leiden"] = "0"

    monkeypatch.setattr(ClusteringModule, "_run_cpu", fake_run_cpu)
    monkeypatch.setattr(ClusteringModule, "_plot_umap_clusters", lambda self, a, c: None)

    cfg = PipelineConfig(
        project="threads",
        output_dir=tmp_path / "out",
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        gpu_mode="off",  # keep the assertion on the CPU lane the measurement used
    )
    ctx = PipelineContext(
        cfg=cfg,
        run_dir=tmp_path / "run",
        figure_dir=tmp_path / "run",
        table_dir=tmp_path / "run",
        adata=AnnData(np.ones((4, 4), dtype=float)),
    )
    ClusteringModule().run(ctx)
    return ctx


def test_resolved_thread_count_is_recorded(monkeypatch, tmp_path):
    """The count the graph actually ran at must be in the manifest metadata."""
    ctx = _run_clustering(monkeypatch, tmp_path, cores=16)

    assert ctx.metadata["clustering_n_jobs"] == _DEFAULT_CLUSTERING_N_JOBS
    assert ctx.metadata["clustering_n_jobs_source"] == "default"
    assert ctx.metadata["clustering_n_jobs_requested"] == _DEFAULT_CLUSTERING_N_JOBS
    assert ctx.metadata["clustering_n_jobs_available"] == 16
    # Recorded value must be the value scanpy was actually driven with, not a
    # parallel bookkeeping number that could drift from it.
    assert clustering_mod.sc.settings.n_jobs == ctx.metadata["clustering_n_jobs"]


def test_default_thread_count_does_not_track_host_core_count(monkeypatch, tmp_path):
    """Same seed + same input on two hosts must resolve the same thread count.

    This is the reproducibility property. Under the old
    ``sc.settings.n_jobs = os.cpu_count()`` these two arms resolve 7 and 31,
    which is a different kNN graph and a different cluster count.
    """
    small = _run_clustering(monkeypatch, tmp_path / "a", cores=7)
    small_setting = clustering_mod.sc.settings.n_jobs
    large = _run_clustering(monkeypatch, tmp_path / "b", cores=31)
    large_setting = clustering_mod.sc.settings.n_jobs

    # Assert on what scanpy was driven with, not only on the record, so this
    # cannot pass by recording a number the graph never used.
    assert small_setting == large_setting  # pre-fix: 7 != 31
    assert small.metadata["clustering_n_jobs"] == large.metadata["clustering_n_jobs"]
    assert small.metadata["clustering_n_jobs"] == _DEFAULT_CLUSTERING_N_JOBS


def test_env_var_sets_thread_count_explicitly(monkeypatch, tmp_path):
    """A run can pin the graph's thread count, which is how an archive is replayed."""
    monkeypatch.setenv(_N_JOBS_ENV, "3")
    ctx = _run_clustering(monkeypatch, tmp_path, cores=16)

    assert ctx.metadata["clustering_n_jobs"] == 3
    assert ctx.metadata["clustering_n_jobs_source"] == f"env:{_N_JOBS_ENV}"
    assert clustering_mod.sc.settings.n_jobs == 3
    assert "clustering_n_jobs_clamped" not in ctx.metadata


def test_thread_count_clamped_to_host_capacity_is_recorded(monkeypatch, tmp_path):
    """Replaying a 64-thread record on an 8-thread host must say so, not crash.

    ``numba.set_num_threads`` raises above the launch-time cap, so the request
    has to be clamped; a clamp means the archived graph was NOT reproduced, so
    it belongs in the record.
    """
    monkeypatch.setenv(_N_JOBS_ENV, "64")
    ctx = _run_clustering(monkeypatch, tmp_path, cores=8)

    assert ctx.metadata["clustering_n_jobs"] == 8
    assert ctx.metadata["clustering_n_jobs_requested"] == 64
    assert "clustering_n_jobs_clamped" in ctx.metadata
    assert "64" in ctx.metadata["clustering_n_jobs_clamped"]


@pytest.mark.parametrize("value", ["-1", "0", "all", ""])
def test_non_reproducible_thread_request_is_rejected(monkeypatch, tmp_path, value):
    """"-1 = every core" is the host-dependence being removed; refuse it loudly.

    Empty string falls through to the default (that is how the env var is
    switched off), so it is the one value here that must NOT raise.
    """
    monkeypatch.setenv(_N_JOBS_ENV, value)
    if value == "":
        ctx = _run_clustering(monkeypatch, tmp_path, cores=16)
        assert ctx.metadata["clustering_n_jobs_source"] == "default"
        return
    with pytest.raises(ValueError, match="thread count"):
        _run_clustering(monkeypatch, tmp_path, cores=16)


# Child script for the empirical canary below. It has to run out-of-process
# because tests/conftest.py sets NUMBA_DISABLE_JIT=1, which routes pynndescent
# around the parallel NN-descent kernel that the thread dependence lives in.
_CANARY_SOURCE = '''
import numpy as np, scanpy as sc, numba
from anndata import AnnData

rng = np.random.default_rng(0)
n_obs = 12_000                      # >= 8192: see the calling test's docstring
pcs = rng.normal(size=(n_obs, 40)).astype(np.float32)
pcs[: n_obs // 2] += 3.0            # two blobs, so the graph is not pure noise

graphs = {}
for n_jobs in (1, int(numba.config.NUMBA_NUM_THREADS)):
    adata = AnnData(np.zeros((n_obs, 2), dtype=np.float32))
    adata.obsm["X_pca"] = pcs.copy()
    sc.settings.n_jobs = n_jobs
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep="X_pca",
                    method="umap", random_state=42)
    graphs[n_jobs] = adata.obsp["distances"].tocsr().indices.copy()

single, multi = graphs[1], graphs[max(graphs)]
differing = int((single != multi).sum())
print(f"CANARY threads={sorted(graphs)} differing_edges={differing} "
      f"identical={np.array_equal(single, multi)}")
'''


@pytest.mark.skipif(
    os.environ.get("SC_THREAD_DETERMINISM_EMPIRICAL") != "1",
    reason="Empirical canary: ~30 s and needs a JIT-enabled subprocess. "
           "Opt in with SC_THREAD_DETERMINISM_EMPIRICAL=1.",
)
def test_neighbour_graph_actually_depends_on_thread_count(tmp_path):
    """Canary for the premise: the kNN graph moves with the thread count.

    This is the measurement the fix rests on, kept executable so the rationale
    can be re-checked rather than believed. If it ever FAILS, pynndescent has
    become thread-count invariant and the fixed default in clustering.py can be
    relaxed — a failure here is a green light, not a defect.

    Needs n_obs >= 8192: below that scanpy takes the exact brute-force shortcut
    for euclidean metrics (scanpy/neighbors/__init__.py:733), which is
    thread-independent by construction.
    """
    import subprocess
    import sys

    script = tmp_path / "canary.py"
    script.write_text(_CANARY_SOURCE)
    env = dict(os.environ, NUMBA_DISABLE_JIT="0")
    proc = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True, timeout=900, env=env,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    verdict = [ln for ln in proc.stdout.splitlines() if ln.startswith("CANARY")]
    assert verdict, proc.stdout[-2000:]
    assert "identical=False" in verdict[0], (
        f"kNN graph no longer varies with thread count ({verdict[0]}) — see docstring."
    )
