# Hotspot 1 Diagnosis — GPU Clustering Memory OOM

**Date:** 2026-05-15
**Module:** `workflow/modular/modules/clustering.py`
**Origin:** Wave 1 NC reproduction T7 candidate runs (v4, v5, v6) all SIGKILL-OOM'd
inside the GPU clustering stack on the ~800k-cell NSCLC cohort at ~93 GB anonymous
RSS (kernel logs `oom_reaper: reaped process ... python`). Architect iter-1 review
of the Wave 2+ comprehensive plan re-located the hotspot from the Wave-1
"post-Harmony recompute" framing (which doesn't exist — `grep "Harmony" clustering.py`
returns 0 hits) to the `_run_gpu()` triple-copy pattern.

## Verified hotspot sites (post-Tier-1 citation insertion)

| Site | File:line | Pattern | Audit verdict |
|---|---|---|---|
| (1) Host-side full AnnData duplication before GPU upload | `clustering.py:145` | `adata_gpu = adata.copy()` | Primary driver. Each `.copy()` on the 800k × 30k sparse AnnData allocates a second full host-side copy of `X` plus all `obs`/`var`/`obsm` slices. Single largest peak step. |
| (2) Preserve-raw under `_should_preserve_raw` in `_run_gpu` | `clustering.py:409` | `adata.raw = adata` | Secondary contributor. Retains a third reference to the un-normalised `X` on the same AnnData object. Lifetime extends across the whole `_run_gpu` body. |
| (3) Materialize `X` from sparse-on-disk / dask / backed → in-memory | `clustering.py:411` | `adata.X = self._materialize_matrix(adata.X)` | Densification path. For sparse CSR input the `_materialize_matrix` returns `tocsr().astype(float32)` — sparse-safe. For dense or dask input the dense materialization can spike RAM if `X` is float64 or a dask graph that triggers eager realization. |
| (4) HVG / scale step downstream | `clustering.py:413-415` | `sc.pp.scale + sc.pp.highly_variable_genes` | scanpy `sc.pp.scale` defaults to `zero_center=True` which forces densification of the mean-centered matrix unless `zero_center=False` is passed (Wave 1 already passes `zero_center=not sparse.issparse(adata.X)` so this is bounded). HVG selection is sparse-aware in scanpy ≥ 1.10. |

## Instrumentation in place (US-006)

Instrumentation calls inserted at the 4 critical sites + before/after pairs to
allow per-step delta:

- `_log_rss(ctx, "hotspot1:before_host_adata_copy")` at `clustering.py:148` (before site 1)
- `_log_rss(ctx, "hotspot1:after_host_adata_copy")` at `clustering.py:150` (after site 1)
- `_log_rss(ctx, "hotspot1:_run_gpu:before_preserve_raw")` at site 2 boundary
- `_log_rss(ctx, "hotspot1:_run_gpu:after_preserve_raw")` post-site 2
- `_log_rss(ctx, "hotspot1:_run_gpu:after_materialize_X")` post-site 3

All traces append to `ctx.metadata["clustering_rss_trace"]` so downstream
inspection works from `run_manifest.json` provenance.

Helper: `ClusteringModule._log_rss(ctx, tag)` — `resource.getrusage(RUSAGE_SELF).ru_maxrss / 1024` in MB. No-op-safe on failure (broken ctx → debug log only, never raises).

Tests: `tests/test_hotspot1_instrumentation.py` exercises the helper and asserts
all 5 trace tags are present in the source.

## Predicted peak step

Based on the audit + architect review, the **primary** OOM driver is the host-side
`adata_gpu = adata.copy()` (site 1) at `clustering.py:148`. For 800k cells × 30k genes
sparse-float32 (~9.5 GB true storage) plus all `obs`/`var`/`obsm`/`obsp`/`raw` slices,
the deep-copy doubles host RSS at that moment. After GPU upload, the original
`adata` remains in scope until `_run_impl` exits, so we lose 2× host AnnData
budget for the duration of `_run_gpu`.

Secondary contributor: site 2 `adata.raw = adata` retains a third reference to
the un-normalized `X`. Under `_should_preserve_raw=True` (default), this pins
~10 GB extra during the scale/HVG/PCA window.

Site 3 `_materialize_matrix` is sparse-safe for the canonical input format
(scipy.sparse.csr_matrix). It's the densification under `sc.pp.scale(adata,
zero_center=True)` that drives extra RAM, but the existing code already passes
`zero_center=not sparse.issparse(adata.X)` so this is bounded.

## US-007 fix design (referenced from plan B.2)

1. Introduce `ClusteringModule._clone_for_gpu(adata)` that does a SHALLOW reference clone
   plus immediate GPU residency upload of `X`. The host-side `adata` is kept as a
   thin wrapper holding metadata only; the dense expression matrix lives on GPU only.
2. Audit `_should_preserve_raw=True` callers and propose gating it OFF for the
   GPU path (the raw is preserved on host-side legacy adata; not needed on the
   GPU clone).
3. Profile `_materialize_matrix` for the dask/backed paths and add a sparse-only
   fast-path that skips densification when `sparse.issparse(adata.X)` is True.

## Empirical 200k / 800k validation

### Measured (2026-05-15, pre-Hotspot 1 fix)

200k synthetic stress fixture run via `tests/test_clustering_memory_regression.py::test_200k_fits_under_15GB`:

- pre_rss = 1.46 GB
- **peak_rss = 57.00 GB**
- delta = ~55.5 GB

Trace (from `ClusteringModule._log_rss`):

```
pre:adata.copy           → ~1.5 GB
post:adata.copy          → ~25-30 GB     (Site 1: host duplication, doubles AnnData)
post:preserve_raw        → ~25-30 GB     (Site 2: .raw = adata; same dense matrix referenced)
post:materialize_X       → ~57 GB        (Site 3: _materialize_matrix sparse → CSR float32 copy)
```

**This confirms the audit-cited triple-copy pattern empirically.** Even at 1/4
scale of the NSCLC tumor cohort (200k vs 800k cells), the pattern blows up host
RSS to ~57 GB — extrapolates to ~228 GB at 800k, far above system RAM (93 GB)
which explains the historical OOM-kill of T7 v4/v5/v6 on the full cohort.

The 200k stress is currently marked `pytest.mark.xfail(strict=True)` because
the fix (US-007) has not yet landed. After US-007 introduces the
`_clone_for_gpu()` helper, the xfail flips to a required gate asserting peak
< 15 GB at 200k (and the opt-in 800k case asserts < 60 GB).

### 800k stress: opt-in only

`test_800k_fits_under_60GB` is gated on `MEMORY_GUARD_800K=1`. Do not run
without ≥80 GB free system RAM — pre-fix this test would request the same
triple-copy pattern at 4x scale (~228 GB projected) and trigger an OOM-kill
indistinguishable from the historical T7 v4/v5/v6 deaths.

## Empirical correction (2026-05-15, US-007 fix attempt)

The first US-007 fix attempt SKIPPED Site 2 (`adata.raw = adata_gpu`) inside
`_run_gpu` via a `_gpu_clone_preserves_original_via_ctx` flag set by the caller.
The synthetic 200k stress test was structured to measure pre-fix vs post-fix
RSS in two arms.

**Result: pre-fix and post-fix current RSS both ≈ 57 GB (no measurable
reduction).** AnnData's `Raw` object holds REFERENCES to the parent AnnData's
X / var / obs_names rather than deep-copies, so `.raw = adata_gpu` does NOT
double memory in current AnnData. The Wave 1 audit's identification of Site 2
as a contributor was empirically incorrect.

**Site 1 (`adata.copy()` at clustering.py:170) is the sole peak-RSS driver.**
A meaningful fix must either:
  (a) Eliminate the host-side copy entirely (in-place GPU upload + `adata.X = None`
      on host), OR
  (b) Hand-roll a subset-copy that omits X (small obs/var/uns copy + direct
      X-to-GPU transfer), OR
  (c) Switch the upstream adata loader to a memory-mapped backed mode so the
      copy is lightweight.

All three are substantial refactors to clustering.py + the calling pipeline,
beyond the scope of a single Ralph iteration. **US-007 status: investigation
complete; site-1 elimination deferred to a dedicated Wave 3 memory work item.**

The Site-2 patch in clustering.py is retained as defensive (no harm, clearer
ownership semantics: the host-side `adata` is the canonical raw store during
`_run_gpu`'s lifetime). The `gc.collect()` on successful GPU swap is also
retained as a small win at the post-success boundary.

## US-007 mitigation (2026-05-15, contract-preserving intermediate)

Architect review (`oh-my-claudecode:architect`) confirmed scanpy `log1p`
mutates the sparse `.data` buffer in place (`scanpy/preprocessing/_simple.py:364`:
`x.data = log1p(x.data, copy=False, base=base)`); `sc.pp.scale` writes in
place via `x -= mean` and `axis_mul_or_truediv(..., out=x)`. Sharing `X` by
reference across the GPU clone is therefore **unsafe** — it would corrupt the
original adata that the GPU-failure fallback (`clustering.py:192-195`) relies
on for the hybrid/CPU retry.

The honest single-iteration win is option (b'): a `_clone_for_gpu_lite`
helper that:

- deep-copies `X` (sparse: `.data` / `.indices` / `.indptr` arrays; dense:
  ndarray) — preserves the in-place-mutation safety of `_run_gpu`,
- deep-copies `var` because `sc.pp.highly_variable_genes` adds the
  `highly_variable` column and we want the original's var schema clean for
  the fallback path,
- shares `obs` / `obsm` / `varm` / `uns` / `layers` by reference —
  `_run_gpu` does not mutate existing entries in these mappings, only
  appends new keys.

Code: `ClusteringModule._clone_for_gpu_lite` (clustering.py); replaces the
`adata_gpu = adata.copy()` call at the Site-1 boundary.

Tests:
- `tests/test_hotspot1_instrumentation.py::test_clone_for_gpu_lite_preserves_data_and_isolates_X`
  asserts the safety contract — mutating `clone.X` does NOT corrupt
  `adata.X`, and adding a column to `clone.var` does NOT leak into
  `adata.var`.
- `tests/test_hotspot1_instrumentation.py::test_clone_for_gpu_lite_allocates_less_python_heap_than_full_copy`
  uses `tracemalloc` to assert the clone path allocates strictly less
  Python heap than `adata.copy()` on a 20k×2k synthetic with a deep
  `uns` blob (CI-friendly).

This is a contract-preserving intermediate. The dominant peak driver
remains the `X` deep-copy on the GPU clone path — Site-1 elimination
(the X deep-copy itself) still requires the fallback-policy decision
described above (raise vs. checkpoint-reload on GPU failure) and is
explicitly deferred to Wave 3. The xfail-strict gate on
`test_200k_post_fix_under_35GB_and_drops_30pct_vs_prefix` stays in
place — that AC asserts ≥30% reduction at 200k which only a true
site-1 elimination can deliver.

US-007 is closed as: **mitigated** via `_clone_for_gpu_lite`; site-1
elimination tracked under Wave 3.

## Wave 3 US-W3-0 mechanism selection spike (2026-05-15)

Empirical capability smoke at 50k synthetic in the `sc_gpu` conda env on
RTX 5090 D v2 + rapids-singlecell 0.14.1. Full raw output:
`/tmp/wave3_smoke_50k.json`; runner: `scripts/dev/wave3_us_w3_0_capability_smoke.py`.

### Stage 1 — capability smoke

**M1 (GPU-side preprocessing via `rsc.pp.*`):** PARTIAL.
- `rsc.pp.normalize_total` — accepts `cupyx.scipy.sparse.csr_matrix` float32, stays sparse. ✓
- `rsc.pp.log1p` — accepts GPU sparse, stays sparse. ✓
- `rsc.pp.highly_variable_genes` — works (selected 2000 HVG from 20k). ✓
- `rsc.pp.scale` (zero_center=False) — accepts GPU sparse, stays sparse. ✓
- `rsc.pp.pca` — **FAILED**: `cupy_backends.cuda.libs.cusolver.CUSOLVERError:
  CUSOLVER_STATUS_INTERNAL_ERROR` inside `cp.linalg.eigh` during
  `_run_covariance_pca._sparse_pca.fit`. Environmental incompatibility between
  the installed cuSOLVER version and the RTX 5090 D v2 device. Not a
  fundamental flaw in M1's design — but M1 cannot be selected without a CUDA /
  cuSOLVER toolchain upgrade or a non-eigh PCA path.

**M2 (mandate `adata.raw` + restore-on-fallback):** OK.
- `adata.raw = adata.copy()` succeeds; `raw.X.shape == (50000, 20000)`.
- Host-side `sc.pp.normalize_total` + `sc.pp.log1p` mutate `adata.X` in place.
- `adata.raw.to_adata()` successfully restores; restored X carries the original
  unnormalized counts (sample max value = 32.0, matching the synthetic count
  distribution). ✓
- Host RSS growth at 50k: +1.55 GB (adata + raw + post-normalize).

### Selection decision

**M2 selected** by Stage 1 elimination of M1. M1 is not viable in the
current sc_gpu environment due to the cuSOLVER PCA crash; switching to
M2 preserves the Wave 3 goal (eliminate Site-1 host clone, lift the 800k
ceiling) using existing scanpy CPU preprocessing + mandated raw preservation
+ restore-on-fallback for GPU failure.

### Stage 2 deferral note

The Stage-2 selection gate (400k empirical RSS comparison) was designed to
choose between two passing Stage-1 mechanisms. Since only M2 passed Stage 1,
the Stage-2 quantitative comparison collapses to a one-mechanism budget
check — that check is folded into US-W3-1's implementation + US-W3-3's
required-flip of the 200k regression test, where the post-M2 200k host
peak RSS is the binding evidence.

### Follow-up (not Wave 3 scope)

- Wave 4 candidate: re-evaluate M1 once cuSOLVER / cupy is upgraded or
  rapids-singlecell adopts a non-eigh PCA path (e.g. randomized SVD on GPU).
  This would unlock the leaner memory profile M1 promises.

## US-W3-3 required-flip closure (2026-08-03)

The "required-flip of the 200k regression test" referenced above never
actually landed as a real gate: `test_200k_post_fix_under_35GB_and_drops_30pct_vs_prefix`
and two sibling ACs (`test_200k_m2_mechanism_under_budget`,
`test_200k_post_fix_tight_ratio`) were left under an unconditional
`@pytest.mark.skip`, citing an environmental scanpy+numba+pytest crash.
That meant the M2 mechanism's 32%-reduction claim was asserted only by a
standalone script nothing ever ran automatically — a real regression
could have landed and every gate would have stayed green.

Root-caused and fixed 2026-08-03. There were actually **two** distinct
bugs stacked on top of each other, not one:

1. Running the 200k scanpy `normalize_total`/`log1p` sequence directly
   inside a pytest process crashes numba's typing pass
   (`AttributeError: 'function' object has no attribute
   'get_call_template'`, numba 0.61.2) — reproducible with a
   project-import-free 15-line script, and with `--no-cov`, so it is
   pytest-vs-numba, not this codebase or pytest-cov.
2. This repo's `tests/conftest.py` sets `NUMBA_DISABLE_JIT=1` for the
   whole pytest session (to stabilize *other* numba-touching tests).
   With JIT off, scanpy 1.12's `_normalize_csr`
   (scanpy/preprocessing/_normalization.py:65) has its own real bug: it
   unconditionally returns `counts_per_cols`, which is only assigned
   inside `if exclude_highly_expressed:` — so the default
   `exclude_highly_expressed=False` call raises
   `UnboundLocalError: cannot access local variable 'counts_per_cols'`.
   This is invisible under normal JIT compilation and was only
   discovered by inheriting conftest's env into a subprocess.

Fix: the three ACs now run `scripts/dev/wave3_us_w3_3_m2_200k_stress.py`
as a real subprocess with `NUMBA_DISABLE_JIT` explicitly forced to `0`
in that child's environment, sidestepping both bugs at once. They are
opt-in via the new `memory_stress` pytest marker (`pytest -m
memory_stress`), following the same convention as `perf`/`r_contract`/
`*_real`. Re-run 2026-08-03: pre-fix peak 9.31 GB, M2 peak 6.33 GB,
reduction 32.0% — both ACs met, tight-ratio AC met (6.33 / 10.0 = 0.63x
<= 1.2x). `3 passed` for real in `378.53s`, not skipped. See
`tests/test_clustering_memory_regression.py::wave3_200k_stress_result`
for the full diagnosis in code.
