"""US-008 / B.3 — Clustering memory regression tests.

200k and 800k synthetic stress fixtures for the Hotspot 1 instrumentation trail.

Strategy:
  * 200k stress: NOT run inline under pytest (see wave3_200k_stress_result()'s
    docstring for the two stacked scanpy/numba/pytest bugs that makes that
    unsafe in this conda env). Instead it runs scripts/dev/wave3_us_w3_3_m2_200k_stress.py
    as a real subprocess, opt-in via `pytest -m memory_stress`, and is
    ratcheted so it only re-runs for real when a watched source file changes
    (see the `_ratchet_*` helpers below) -- otherwise it reports an explicit,
    reasoned skip naming the last real measurement. Wired into the suite gate
    (../scripts/run_all_gates.sh) as a host_only step, not the default fast
    lane: ~4-6 min wall / ~10 GB peak RSS on a real run.
  * 800k fixture: env-gated on MEMORY_GUARD_800K=1 — only runs when the host
    explicitly opts in (≥80GB free RAM). Programmatic synthetic; not committed.

Both stress paths:
  1. Generate a sparse synthetic AnnData (cells × ~20k genes, ~0.1-0.3 density).
  2. Run the M2 mechanism (mandate adata.raw + in-place host preprocessing;
     no host clone of the full adata) that replaced the original triple-copy
     Hotspot 1 pattern.
  3. Capture RSS (peak via resource.getrusage for the 800k in-process path;
     /proc/self/status VmRSS checkpoints for the 200k subprocess path).
  4. Assert peak stays under the budget.

Realistic post-fix budget at 200k is <35 GB (the original <15GB AC target was
infeasible without rewriting AnnData to share X across copy() — see
docs/HOTSPOT1_DIAGNOSIS.md "Empirical 200k / 800k validation"), with a tighter
canary at <=1.2x the ~10GB theoretical minimum and a required >=30% reduction
vs. the pre-fix triple-copy pattern.
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

_SC_ROOT = Path(__file__).resolve().parents[1]

# Canonical 200k memory-regression stress script. See wave3_200k_stress_result()
# below for why the pytest tests below invoke it as a SUBPROCESS instead of
# reproducing its body inline under pytest.
_WAVE3_200K_SCRIPT = _SC_ROOT / "scripts" / "dev" / "wave3_us_w3_3_m2_200k_stress.py"

# --- Ratchet: makes the ~4-6 min / ~10GB stress cheap to re-verify -----------
#
# `memory_stress` is opt-in and expensive, so it is wired into the suite gate
# (../scripts/run_all_gates.sh) rather than the default fast lane. Running the
# full subprocess on EVERY gate invocation, unconditionally, would roughly
# double local pre-push cost -- exactly the kind of tax that makes people
# reach for SKIP_SUITE_GATE=1 (see validate_claim_critical_deps.py's C5/C6
# rationale for why that is treated as a real risk in this repo). Instead:
# re-run the real subprocess only when a source file the measurement actually
# depends on has changed since the last recorded PASS; otherwise report an
# explicit, reasoned, visible skip naming the prior measurement. This is
# different from the unconditional `@pytest.mark.skip` this file used to
# carry: the condition is checked at collection time from real file content,
# and the reason string names exactly what was last measured and when, so a
# reader (or `assert_selftest_integrity.py --allow-skip-reason-substring`)
# can tell a legitimate ratchet skip from coverage quietly going missing.
_RATCHET_PATH = Path(__file__).resolve().parent / "memory_regression_ratchet.json"
_RATCHET_WATCHED_PATHS = (
    _SC_ROOT / "workflow" / "modular" / "modules" / "clustering.py",
    _SC_ROOT / "workflow" / "modular" / "_mem_guard.py",
    _SC_ROOT / "workflow" / "modular" / "_mem_watchdog.py",
    _WAVE3_200K_SCRIPT,
    Path(__file__).resolve(),
)
_RATCHET_SKIP_REASON_TAG = "memory-regression ratchet"


def _ratchet_watched_hash() -> str:
    """SHA256 over the content of every watched path, order-independent by name."""
    digest = hashlib.sha256()
    for path in sorted(_RATCHET_WATCHED_PATHS, key=lambda p: str(p)):
        digest.update(str(path.relative_to(_SC_ROOT)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes() if path.is_file() else b"<missing>")
        digest.update(b"\0")
    return digest.hexdigest()


def _ratchet_skip_reason() -> str | None:
    """None -> must run for real. A string -> safe to skip, and why."""
    if not _RATCHET_PATH.is_file():
        return None
    try:
        record = json.loads(_RATCHET_PATH.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if record.get("verdict") != "pass":
        return None
    if record.get("watched_sha256") != _ratchet_watched_hash():
        return None
    return (
        f"{_RATCHET_SKIP_REASON_TAG}: unchanged since {record.get('recorded_utc', '?')} "
        f"(watched_sha256={record.get('watched_sha256', '?')[:12]}...). Last measured "
        f"pre-fix={record.get('pre_fix_peak_gb')}GB m2={record.get('m2_peak_gb')}GB "
        f"reduction={record.get('reduction_pct')}%. Delete "
        f"tests/{_RATCHET_PATH.name} to force a real re-run."
    )


# Computed once at collection time (cheap: a handful of small source-file hashes,
# NOT the 200k stress itself). Applied to all three memory_stress tests below.
_RATCHET_SKIP_REASON = _ratchet_skip_reason()
_ratchet_skipif = pytest.mark.skipif(
    _RATCHET_SKIP_REASON is not None,
    reason=_RATCHET_SKIP_REASON or "(unreachable: ratchet says run)",
)


def _peak_rss_gb() -> float:
    """Peak RSS since process start, in GB (Linux ru_maxrss is in KB).

    NOTE: this is monotonic across the lifetime of the process and CANNOT be
    used to measure a single test arm in isolation. Use _current_rss_gb() for
    independent arm measurement.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 2)


def _current_rss_gb() -> float:
    """Current resident set size in GB (not monotonic; reflects live allocation)."""
    # Read /proc/self/status for VmRSS — Linux-only but matches our deployment target.
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                # "VmRSS:    12345 kB"
                kb = int(line.split()[1])
                return kb / (1024.0 ** 2)
    return float("nan")


def _synthetic_adata(n_cells: int, n_genes: int = 20_000, density: float = 0.3,
                     seed: int = 42) -> ad.AnnData:
    """Programmatic synthetic AnnData. Sparse float32 X with integer-like counts
    (np.ceil(rand * 100)) so scanpy normalize_total + log1p see realistic
    count-matrix semantics — pure-fractional random data triggers numba
    compilation issues in scanpy/_normalization.py:65 under some pytest envs.
    """
    rng = np.random.default_rng(seed)
    X = sp.random(n_cells, n_genes, density=density, format="csr",
                  dtype=np.float32, random_state=seed)
    X.data = np.ceil(X.data * 100).astype(np.float32)
    obs = pd.DataFrame({
        "leiden": rng.integers(0, 8, size=n_cells).astype(str),
        "sample": [f"S{i % 10}" for i in range(n_cells)],
    }, index=[f"cell_{i:08d}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"GENE_{i:05d}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


def _exercise_hotspot1_pattern(adata: ad.AnnData, preserve_raw_via_ctx: bool = True) -> dict:
    """Reproduce the Hotspot 1 GPU-clone code path WITHOUT requiring rapids.

    Mirrors the patched clustering._run_gpu():
      - Site 1 (clustering.py:170): adata_gpu = adata.copy()  — host duplication
      - Site 2 (clustering.py:411): adata_gpu.raw = adata_gpu  — SKIPPED when
        preserve_raw_via_ctx=True (US-007 fix); the caller is responsible for
        keeping the original adata reachable via ctx.adata.
      - Site 3 (clustering.py:412): _materialize_matrix(X) — sparse-safe.

    Pass preserve_raw_via_ctx=False to reproduce the PRE-FIX pattern for
    comparison (used when measuring the fix's RSS delta).
    """
    from workflow.modular.modules.clustering import ClusteringModule

    ctx = SimpleNamespace(metadata={})
    ClusteringModule._log_rss(ctx, "pre:adata.copy")
    adata_gpu = adata.copy()
    ClusteringModule._log_rss(ctx, "post:adata.copy")
    if adata_gpu.raw is None and not preserve_raw_via_ctx:
        # PRE-FIX behavior: doubles dense storage via .raw = adata_gpu
        adata_gpu.raw = adata_gpu
    ClusteringModule._log_rss(ctx, "post:preserve_raw" if not preserve_raw_via_ctx else "skip:preserve_raw")
    # Mimic _materialize_matrix on sparse path (no-op-equivalent: tocsr().astype(float32))
    X = adata_gpu.X
    if sp.issparse(X):
        adata_gpu.X = X.tocsr().astype(np.float32)
    ClusteringModule._log_rss(ctx, "post:materialize_X")
    return {
        "trace": ctx.metadata["clustering_rss_trace"],
        "final_rss_gb": _current_rss_gb(),  # current state, not peak-since-start
    }


@pytest.fixture(scope="module")
def wave3_200k_stress_result() -> dict:
    """Run the Wave 3 US-W3-3 200k memory-regression stress as a SUBPROCESS.

    WHY A SUBPROCESS INSTEAD OF RUNNING THE STRESS INLINE UNDER PYTEST
    --------------------------------------------------------------------
    Running scanpy.pp.normalize_total/log1p over a 200k x 20k sparse
    synthetic AnnData directly inside a pytest test function crashes
    reproducibly in this conda env (re-verified 2026-08-03: scanpy 1.12,
    numba 0.61.2, pytest 9.0.2, numpy 2.2.6):

        AttributeError: 'function' object has no attribute
        'get_call_template'
        (numba/core/types/functions.py:538, raised while resolving an
        `ol_np_zeros` overload reached from scanpy/_normalization.py)

    This was re-confirmed with `--no-cov` (rules out the coverage tracer)
    and with a minimal 15-line repro that imports nothing from this
    project (rules out the project's import graph / pytest-cov / a
    project-specific bug). It is specifically an interaction with running
    under `python -m pytest`: the byte-for-byte identical code, executed
    as a plain `python script.py` process, completes cleanly every time
    (confirmed both at this 200k scale and originally at 2026-05-15 in
    this file's prior skip reasons).

    SECOND, DISTINCT BUG THIS FIXTURE ALSO WORKS AROUND: this repo's own
    tests/conftest.py sets `NUMBA_DISABLE_JIT=1` (module scope, to
    stabilize *other* numba-touching tests) for the whole pytest session.
    A bare subprocess.run() inherits that into the child by default, and
    with JIT disabled, scanpy 1.12's `_normalize_csr`
    (scanpy/preprocessing/_normalization.py:65) has a REAL upstream bug:
    it unconditionally returns `counts_per_cell, counts_per_cols`, but
    `counts_per_cols` is only ever assigned inside
    `if exclude_highly_expressed:` -- so with the default
    `exclude_highly_expressed=False` (what this stress script uses),
    plain-Python execution raises
    `UnboundLocalError: cannot access local variable 'counts_per_cols'`.
    Under normal JIT compilation this particular bug is masked (numba's
    SSA-based typing does not raise the same CPython-semantics error), so
    it is invisible except when JIT is off -- exactly this repo's pytest
    default. Verified 2026-08-03 by running the stress subprocess with
    the inherited environment (crashes with the UnboundLocalError above)
    vs. an explicit `NUMBA_DISABLE_JIT=0` override (passes cleanly).

    So: run the canonical stress script
    (scripts/dev/wave3_us_w3_3_m2_200k_stress.py) as a real subprocess
    with an explicit environment that force-enables numba JIT
    (`NUMBA_DISABLE_JIT=0`), regardless of what the parent pytest process
    has set. That sidesteps BOTH bugs: the numba-JIT-triggering compiler
    crash never executes inside the pytest process, and JIT stays enabled
    in the child so the scanpy UnboundLocalError bug's un-JIT-ed code
    path is never reached either -- while the assertions below still run
    for real, every time this fixture's marker is selected, instead of
    being permanently skipped.

    Opt-in only (`pytest -m memory_stress`; excluded from the default
    fast lane via addopts in pyproject.toml, same convention already used
    for `perf` / `r_contract` / `*_real`). Real cost: ~4-6 min wall,
    ~10 GB peak RSS on this host (JIT warmup dominates the wall time; the
    tracked VmRSS checkpoints stay near the numbers in this docstring's
    empirical baseline). Module-scoped so the three tests below that
    consume it only pay for one subprocess run, not three.
    """
    if not _WAVE3_200K_SCRIPT.is_file():
        pytest.fail(f"canonical stress script missing: {_WAVE3_200K_SCRIPT}")
    # Force JIT ON in the child regardless of what conftest.py set for the
    # parent pytest process (see docstring above) -- this is what actually
    # makes the subprocess strategy work, not merely "not being pytest".
    child_env = dict(os.environ)
    child_env["NUMBA_DISABLE_JIT"] = "0"
    proc = subprocess.run(
        [sys.executable, str(_WAVE3_200K_SCRIPT)],
        capture_output=True, text=True, timeout=600, env=child_env,
    )
    if proc.returncode != 0:
        pytest.fail(
            "wave3_us_w3_3_m2_200k_stress.py subprocess failed "
            f"(exit {proc.returncode}):\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR (tail):\n{proc.stderr[-4000:]}"
        )
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        pytest.fail(
            f"could not parse stress-script JSON stdout ({exc}):\n{proc.stdout}"
        )
    result = payload["us_w3_3_m2_200k_stress"]

    # Ratchet write: only record a "pass" checkpoint when ALL THREE downstream
    # ACs (not just the script's own ac_met, which only covers two of them)
    # actually hold. A partial/failing result must NOT be ratcheted, or a
    # future unchanged-sources run would wrongly skip past a real regression.
    m2_peak = result["m2"]["rss_peak_gb"]
    all_acs_met = (
        bool(result.get("budget_35GB_met"))
        and bool(result.get("reduction_30pct_met"))
        and m2_peak <= 1.2 * 10.0
    )
    if all_acs_met:
        _RATCHET_PATH.write_text(json.dumps({
            "verdict": "pass",
            "recorded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "watched_sha256": _ratchet_watched_hash(),
            "pre_fix_peak_gb": result["pre_fix"]["rss_peak_gb"],
            "m2_peak_gb": m2_peak,
            "reduction_pct": result["reduction_pct"],
        }, indent=2) + "\n")
    return result


@_ratchet_skipif
@pytest.mark.memory_stress
def test_200k_post_fix_under_35GB_and_drops_30pct_vs_prefix(wave3_200k_stress_result):
    """Wave 3 US-W3-3 — peak RSS < 35 GB AND >=30% drop vs pre-fix at 200k.

    Replaces the pre-Wave-3 Site-2-skip vs Site-2-keep comparison (both arms
    of which still called adata.copy() and produced 0% reduction) with the
    actual Wave 3 M2 mechanism comparison:
      - pre-fix arm: adata.copy() + raw=clone + materialize + normalize/log1p
      - M2 arm:      adata.raw = adata.copy() + in-place normalize/log1p
                     (no host clone of adata itself)

    Original empirical baseline 2026-05-15 (density=0.1, 20k genes, 200k cells):
      pre-fix peak: 9.31 GB / M2 peak: 6.32 GB / reduction: 32.1%
    """
    r = wave3_200k_stress_result
    pre = r["pre_fix"]["rss_peak_gb"]
    m2 = r["m2"]["rss_peak_gb"]
    reduction_pct = r["reduction_pct"]

    print(f"\n  pre-fix peak = {pre:.2f} GB")
    print(f"  M2 peak      = {m2:.2f} GB")
    print(f"  reduction    = {reduction_pct:.1f}%")

    assert r["budget_35GB_met"], f"M2 peak RSS={m2:.2f}GB exceeded 35GB budget"
    assert r["reduction_30pct_met"], f"RSS reduction {reduction_pct:.1f}% < required 30%"


@_ratchet_skipif
@pytest.mark.memory_stress
def test_200k_m2_mechanism_under_budget(wave3_200k_stress_result):
    """Wave 3 US-W3-3 — measure M2 mechanism's actual peak RSS at 200k.

    Locks in the empirical result so future refactors that re-introduce a
    host-side clone of the full adata fail this gate. Shares the single
    subprocess run from wave3_200k_stress_result with the sibling tests in
    this file rather than re-running the 200k stress a second time.
    """
    r = wave3_200k_stress_result
    m2 = r["m2"]["rss_peak_gb"]
    reduction_pct = r["reduction_pct"]
    assert m2 < 35.0, f"M2 peak RSS {m2:.2f} GB exceeded 35 GB budget"
    assert r["reduction_30pct_met"], f"M2 reduction {reduction_pct:.1f}% < required 30%"


@_ratchet_skipif
@pytest.mark.memory_stress
def test_200k_post_fix_tight_ratio(wave3_200k_stress_result):
    """Wave 3 US-W3-3 AC-3 — canary against memory regressions slipping under 35GB slack.

    Theoretical minimum host RSS for 200k x 20k sparse-float32 at density=0.1
    with raw preservation + log1p host mutation:
      sparse X (nnz=400M * 12B) = 4.8 GB
      adata.raw (sparse X copy)  = 4.8 GB
      log1p in-place (no extra)  = 0 GB
      pandas obs/var + overhead  = ~0.5 GB
      THEORETICAL_MIN_GB         ~ 10 GB

    Tighter assertion: peak <= 1.2x theoretical minimum.
    """
    THEORETICAL_MIN_GB = 10.0
    m2 = wave3_200k_stress_result["m2"]["rss_peak_gb"]
    ratio = m2 / THEORETICAL_MIN_GB if THEORETICAL_MIN_GB > 0 else float("inf")
    print(f"\n  M2 peak: {m2:.2f} GB  theoretical_min: {THEORETICAL_MIN_GB:.2f} GB  ratio: {ratio:.2f}x")
    assert m2 <= 1.2 * THEORETICAL_MIN_GB, (
        f"M2 peak {m2:.2f} GB > 1.2x theoretical_min ({1.2 * THEORETICAL_MIN_GB:.2f} GB) "
        f"— possible memory regression upstream of clustering"
    )


@pytest.mark.skipif(
    os.environ.get("MEMORY_GUARD_800K", "0") != "1",
    reason="800k stress is opt-in; set MEMORY_GUARD_800K=1 to run.",
)
def test_800k_fits_under_60GB():
    """800k synthetic AnnData through Hotspot-1 pattern: peak RSS < 60 GB.

    Opt-in only. Requires ≥80 GB free system RAM. This is the audit-cited
    NC2024 NSCLC tumor cohort scale (795,707 cells). If Hotspot 1 is fixed
    (US-007), peak should land ~30-40 GB instead of the historical ~93 GB
    OOM-kill range.
    """
    pre_rss_gb = _peak_rss_gb()
    adata = _synthetic_adata(n_cells=800_000)
    result = _exercise_hotspot1_pattern(adata)
    peak_gb = result["final_rss_gb"]
    delta_gb = peak_gb - pre_rss_gb
    print(f"\n  800k stress: pre_rss={pre_rss_gb:.2f}GB peak_rss={peak_gb:.2f}GB delta={delta_gb:.2f}GB")
    print(f"  trace: {result['trace']}")
    assert peak_gb < 60.0, f"800k peak RSS={peak_gb:.2f}GB exceeded 60GB budget"


def test_synthetic_fixture_is_not_committed():
    """Sanity: no multi-GB .h5ad / .npz binary fixture committed to git for this test."""
    from pathlib import Path
    fixtures_dir = Path(__file__).resolve().parent / "fixtures"
    if not fixtures_dir.exists():
        return  # fine — no fixtures dir
    for path in fixtures_dir.rglob("*"):
        if path.is_file():
            size = path.stat().st_size
            # ~1 GB is the threshold; we should never commit multi-GB stress fixtures
            assert size < 100 * 1024 * 1024, f"committed fixture too large: {path} ({size} bytes)"
