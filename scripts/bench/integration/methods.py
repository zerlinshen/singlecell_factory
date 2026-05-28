"""Integration method wiring (T2) — embeddings via the REAL production backends.

CALIBRATION / RELATIVE lane ONLY. This module computes the four embeddings the
integration calibration harness scores, by REUSING the production
``BatchCorrectionModule`` backends in
``workflow/modular/modules/batch_correction.py`` — Harmony and scVI are NOT
reimplemented here. We call ``BatchCorrectionModule._run_harmony`` and
``BatchCorrectionModule._run_scvi`` directly through a minimal duck-typed
context shim so the actual Korsunsky-2019 / Lopez-2018 calls (and their
fail-loud guards: Principle-9 non-convergence, scVI raw-count check) execute
exactly as in production.

The four embeddings (obsm keys mirror the production names):
  1. ``X_pca``          baseline (uncorrected PCA).
  2. ``X_pca_harmony``  Harmony at the production-default theta.
  3. ``X_scvi``         scVI with a SEED SWEEP (>=3 seeds) -> a per-seed family
                        of embeddings; the scoreboard reports a variance BAND.
  4. negative control   a falsifiable over-corrector — Harmony at EXTREME theta
                        (theta >> default crushes biology to maximise mixing).
                        The OVER-CORRECTING flag (T3) MUST fire on this and MUST
                        NOT fire on the baseline.

Separation of concerns (single-env sc_gpu, but kept clean): every method writes
its embedding into ``adata.obsm`` and the assembled AnnData can be persisted, so
the metric-scoring step (score_integration.score_embeddings) can run as a
distinct pass — embeddings here, metrics there.

scVI seeds: scvi-tools has no per-call seed argument, so the sweep sets
``scvi.settings.seed`` (and a numpy/torch seed) before each independent
setup+train. Each seed is a fully independent model fit -> the spread across
seeds is the run-to-run variance band.
"""
from __future__ import annotations

import contextlib
import logging
import signal
import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# Method/obsm names — kept identical to the production batch_correction backends
# so a reader can trace each calibration embedding back to its pipeline origin.
BASELINE_METHOD = "baseline"
BASELINE_OBSM = "X_pca"
HARMONY_METHOD = "harmony"
HARMONY_OBSM = "X_pca_harmony"
SCVI_METHOD = "scvi"
SCVI_OBSM = "X_scvi"
NEG_CONTROL_METHOD = "neg_control_harmony_extreme_theta"
NEG_CONTROL_OBSM = "X_neg_control"

# scVI seed sweep — >=3 independent fits. Fixed list for determinism.
DEFAULT_SCVI_SEEDS: tuple[int, ...] = (0, 1, 2)

# Negative-control Harmony theta. Production default is 2.0; an extreme theta
# over-mixes (theta is the diversity-penalty strength — large theta forces ALL
# cells together regardless of biology, destroying bio conservation). This is
# the falsifiable over-corrector the OVER-CORRECTING flag must catch.
#
# theta=1000 makes harmonypy's centroid covariance SINGULAR (linalg.inv fails),
# so we use theta=100 — still a strong over-corrector (50x the production
# default of 2.0) but numerically stable. Empirically theta in {20,50,100} all
# run clean; 1000 is singular on this dataset.
NEG_CONTROL_HARMONY_THETA = 100.0

# Per-method wall-clock timeouts (seconds). scVI is the hang risk; baseline is
# instant. controller_validation mode -> conservative caps, long calls wrapped.
DEFAULT_TIMEOUTS: dict[str, float] = {
    BASELINE_METHOD: 60.0,
    HARMONY_METHOD: 1800.0,
    SCVI_METHOD: 3600.0,          # per-seed cap
    NEG_CONTROL_METHOD: 1800.0,
}


# ==========================================================================
# Minimal production-context shim
# ==========================================================================
# BatchCorrectionModule._run_harmony / _run_scvi only read:
#   ctx.cfg.batch.{harmony_theta,harmony_sigma,harmony_max_iter,harmony_backend,
#                  scvi_n_latent,scvi_max_epochs,scvi_early_stopping}
#   ctx.cfg.gpu_mode  (read by _resolve_harmony_backend only via env/cfg; default 'auto')
#   ctx.cfg.clustering.n_pcs (mnn/fastmnn only — unused here)
#   ctx.metadata[...]  (write-only audit sink)
# We provide exactly that surface. Defaults mirror config.py BatchConfig so the
# calibration Harmony/scVI match the production defaults unless overridden.

@dataclass
class _BatchCfgShim:
    batch_key: str = "sample"
    method: str = "harmony"
    harmony_theta: float = 2.0
    harmony_sigma: float = 0.1
    harmony_max_iter: int = 50
    harmony_backend: str = "auto"
    scvi_max_epochs: int = 200
    scvi_n_latent: int = 30
    scvi_early_stopping: bool = True


@dataclass
class _ClusteringCfgShim:
    n_pcs: int = 15
    random_state: int = 0


@dataclass
class _CfgShim:
    batch: _BatchCfgShim = field(default_factory=_BatchCfgShim)
    clustering: _ClusteringCfgShim = field(default_factory=_ClusteringCfgShim)
    gpu_mode: str = "auto"
    random_state: int = 0


@dataclass
class _CtxShim:
    """Duck-typed PipelineContext exposing only the backend-read surface."""

    cfg: _CfgShim = field(default_factory=_CfgShim)
    metadata: dict = field(default_factory=dict)

    @property
    def random_state(self) -> int:
        return getattr(self.cfg, "random_state", 0)

    def status(self, module: str, ok: Any, message: str) -> None:  # no-op sink
        self.metadata.setdefault("_status_log", []).append((module, ok, message))


# ==========================================================================
# Timeout wrapper (controller_validation mode)
# ==========================================================================

class MethodTimeout(RuntimeError):
    """Raised when a per-method backend call exceeds its wall-clock budget."""


@contextlib.contextmanager
def _time_limit(seconds: float | None, label: str):
    """SIGALRM-based wall-clock guard (main-thread only).

    scVI/Harmony can hang on pathological inputs; controller_validation mode
    bounds every long call. If SIGALRM is unavailable (non-main thread), the
    guard degrades to no-op and records that it could not arm — the caller
    still records elapsed time, so an over-run is observable post-hoc.
    """
    if not seconds or seconds <= 0:
        yield
        return

    def _handler(signum, frame):  # noqa: ARG001
        raise MethodTimeout(f"{label}: exceeded {seconds:.0f}s budget")

    armed = False
    try:
        signal.signal(signal.SIGALRM, _handler)
        signal.setitimer(signal.ITIMER_REAL, float(seconds))
        armed = True
    except (ValueError, AttributeError, OSError):
        logger.warning("could not arm SIGALRM for %s; running unguarded", label)
    try:
        yield
    finally:
        if armed:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, signal.SIG_DFL)


# ==========================================================================
# Backend reuse — call the production static methods, do NOT reimplement
# ==========================================================================

def _import_backend():
    """Import the production BatchCorrectionModule (the single source of truth)."""
    from workflow.modular.modules.batch_correction import BatchCorrectionModule
    return BatchCorrectionModule


@dataclass
class MethodEmbedding:
    """One computed embedding + provenance for the scoreboard/audit."""

    method: str
    obsm_key: str
    X: np.ndarray | None
    elapsed_s: float
    backend_meta: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None
    error: str | None = None


def compute_baseline(adata, *, obsm_key: str = BASELINE_OBSM) -> MethodEmbedding:
    """Baseline = uncorrected PCA already present in obsm[X_pca].

    The calibration lane treats the existing PCA as the no-integration control;
    it is not recomputed here (the assembly step is responsible for producing a
    deterministic X_pca). Fail loud if it is missing — a silent zero baseline
    would corrupt the OVER-CORRECTING threshold.
    """
    t0 = time.time()
    if obsm_key not in adata.obsm:
        return MethodEmbedding(
            BASELINE_METHOD, obsm_key, None, time.time() - t0,
            error=f"baseline obsm '{obsm_key}' absent; assemble PCA first",
        )
    X = np.asarray(adata.obsm[obsm_key], dtype=np.float32)
    return MethodEmbedding(BASELINE_METHOD, obsm_key, X, time.time() - t0,
                           backend_meta={"source": f"adata.obsm['{obsm_key}']"})


def compute_harmony(
    adata,
    batch_key: str,
    *,
    theta: float = 2.0,
    method_name: str = HARMONY_METHOD,
    obsm_key: str = HARMONY_OBSM,
    gpu_mode: str = "auto",
    timeout_s: float | None = None,
) -> MethodEmbedding:
    """Run Harmony via BatchCorrectionModule._run_harmony (production path).

    ``theta`` is exposed so the negative control can drive it to an extreme.
    Writes adata.obsm['X_pca_harmony']; we copy it to ``obsm_key`` so the neg
    control (extreme theta) does not clobber the real Harmony embedding.
    """
    Backend = _import_backend()
    ctx = _CtxShim()
    ctx.cfg.batch.batch_key = batch_key
    ctx.cfg.batch.harmony_theta = float(theta)
    ctx.cfg.gpu_mode = gpu_mode

    t0 = time.time()
    err = None
    X = None
    try:
        with _time_limit(timeout_s, method_name):
            # _run_harmony always writes obsm['X_pca_harmony'].
            Backend._run_harmony(adata, batch_key, ctx)
        X = np.asarray(adata.obsm[HARMONY_OBSM], dtype=np.float32)
        if obsm_key != HARMONY_OBSM:
            adata.obsm[obsm_key] = X
    except Exception as exc:  # noqa: BLE001 — recorded, not swallowed silently
        err = f"{type(exc).__name__}: {exc}"
        logger.error("harmony (%s, theta=%s) failed: %s", method_name, theta, err)
    return MethodEmbedding(
        method_name, obsm_key, X, time.time() - t0,
        backend_meta={"harmony_theta": float(theta),
                      "harmony_backend": ctx.metadata.get("harmony_backend"),
                      "harmony_converged": ctx.metadata.get("harmony_converged")},
        error=err,
    )


def compute_harmonypy_direct(
    adata,
    batch_key: str,
    *,
    theta: float = 2.0,
    method_name: str = HARMONY_METHOD,
    obsm_key: str = HARMONY_OBSM,
    pca_key: str = BASELINE_OBSM,
    max_iter_harmony: int = 50,
    sigma: float = 0.1,
    timeout_s: float | None = None,
) -> MethodEmbedding:
    """Harmony via harmonypy.run_harmony DIRECT (production scanpy wrapper bypassed).

    DEVIATION (documented, authorized): the production path
    ``BatchCorrectionModule._run_harmony`` is broken in sc_gpu on this hardware
    on BOTH backends (rapids GPU -> CUBLAS_STATUS_NOT_INITIALIZED + CUDA-context
    corruption; CPU scanpy.external.pp.harmony_integrate 1.12 mis-stores
    harmonypy 0.2.0's Z_corr as shape (n_pcs,) instead of (n_cells, n_pcs)).
    See the audit's production_backend_findings + open-questions CRITICAL Q.

    harmonypy.run_harmony IS the canonical Harmony reference (Korsunsky 2019);
    this routes around scanpy 1.12's proven-buggy thin wrapper while using the
    SAME algorithm on the CPU. harmonypy returns ``Z_corr`` shape
    (n_pcs, n_cells); we transpose to (n_cells, n_pcs) and store it. This is
    labeled distinctly in the scoreboard/audit and is NOT relabeled as
    production-Harmony output.
    """
    t0 = time.time()
    err = None
    X = None
    converged_note = None
    try:
        import harmonypy
        if pca_key not in adata.obsm:
            raise ValueError(f"harmonypy requires PCA in obsm['{pca_key}']")
        pca = np.asarray(adata.obsm[pca_key], dtype=np.float64)
        with _time_limit(timeout_s, method_name):
            ho = harmonypy.run_harmony(
                pca, adata.obs, [batch_key],
                theta=theta, sigma=sigma, max_iter_harmony=max_iter_harmony,
            )
        Z = np.asarray(ho.Z_corr, dtype=np.float32)
        # harmonypy Z_corr is (n_pcs, n_cells) -> transpose to (n_cells, n_pcs).
        if Z.shape[0] != adata.n_obs:
            Z = Z.T
        if Z.shape[0] != adata.n_obs:
            raise ValueError(
                f"harmonypy output shape {Z.shape} does not align to "
                f"n_obs={adata.n_obs} on either axis")
        X = np.ascontiguousarray(Z, dtype=np.float32)
        adata.obsm[obsm_key] = X
    except Exception as exc:  # noqa: BLE001 — recorded, not silently swallowed
        err = f"{type(exc).__name__}: {exc}"
        logger.error("harmonypy-direct (%s, theta=%s) failed: %s",
                     method_name, theta, err)
    return MethodEmbedding(
        method_name, obsm_key, X, time.time() - t0,
        backend_meta={"engine": "harmonypy.run_harmony direct (scanpy wrapper bypassed)",
                      "harmonypy_theta": float(theta),
                      "harmony_max_iter": max_iter_harmony,
                      "sigma": sigma,
                      "production_scanpy_wrapper": "bypassed_due_to_known_bug",
                      "note": converged_note},
        error=err,
    )


def compute_neg_control(
    adata,
    batch_key: str,
    *,
    theta: float = NEG_CONTROL_HARMONY_THETA,
    gpu_mode: str = "auto",
    timeout_s: float | None = None,
    use_harmonypy_direct: bool = True,
) -> MethodEmbedding:
    """Falsifiable negative control: Harmony at EXTREME theta (over-corrector).

    Large theta over-penalises batch diversity, forcing all cells together and
    shredding biological structure -> high mixing, collapsed bio conservation.
    The OVER-CORRECTING flag must fire on this and not on the baseline.

    Uses the harmonypy-direct path by default (the production rapids/scanpy
    Harmony backends are broken in sc_gpu; see production_backend_findings),
    preserving the over-correction control through the SAME canonical algorithm.
    """
    if use_harmonypy_direct:
        return compute_harmonypy_direct(
            adata, batch_key, theta=theta,
            method_name=NEG_CONTROL_METHOD, obsm_key=NEG_CONTROL_OBSM,
            timeout_s=timeout_s,
        )
    return compute_harmony(
        adata, batch_key, theta=theta,
        method_name=NEG_CONTROL_METHOD, obsm_key=NEG_CONTROL_OBSM,
        gpu_mode=gpu_mode, timeout_s=timeout_s,
    )


def compute_shuffle_label_control(
    adata,
    batch_key: str,
    *,
    pca_key: str = BASELINE_OBSM,
    seed: int = 0,
) -> MethodEmbedding:
    """Second falsifiability anchor: permute cell rows of the baseline embedding.

    Shuffling rows destroys ALL biological structure while leaving the embedding
    geometry's marginal distribution intact. A correct bio-conservation metric
    must collapse on this control; it is a cheap, algorithm-independent sanity
    anchor complementing the Harmony-extreme-theta over-corrector.
    """
    t0 = time.time()
    if pca_key not in adata.obsm:
        return MethodEmbedding(
            "neg_control_shuffle_label", "X_shuffle", None, time.time() - t0,
            error=f"shuffle control requires obsm['{pca_key}']")
    rng = np.random.default_rng(seed)
    base = np.asarray(adata.obsm[pca_key], dtype=np.float32)
    perm = rng.permutation(base.shape[0])
    X = np.ascontiguousarray(base[perm], dtype=np.float32)
    adata.obsm["X_shuffle"] = X
    return MethodEmbedding(
        "neg_control_shuffle_label", "X_shuffle", X, time.time() - t0,
        backend_meta={"engine": "row-permutation of baseline PCA", "seed": seed})


def compute_scvi_seed_sweep(
    adata,
    batch_key: str,
    *,
    seeds: tuple[int, ...] = DEFAULT_SCVI_SEEDS,
    gpu_mode: str = "auto",
    per_seed_timeout_s: float | None = None,
    scvi_max_epochs: int | None = None,
    scvi_n_latent: int | None = None,
    scvi_early_stopping: bool | None = None,
) -> list[MethodEmbedding]:
    """scVI seed sweep via BatchCorrectionModule._run_scvi (production path).

    scvi-tools exposes no per-call seed, so each seed sets scvi.settings.seed
    (+ numpy/torch) and runs a fully independent setup+train. Each fit writes
    adata.obsm['X_scvi']; we snapshot it to a per-seed key so the sweep family
    survives. Returns one MethodEmbedding per seed (method names suffixed by
    seed). The band is derived later from these per-seed metric scores.
    """
    if len(seeds) < 3:
        raise ValueError(
            f"scVI seed sweep requires >=3 seeds (band-aware ties); got {seeds!r}"
        )
    Backend = _import_backend()
    out: list[MethodEmbedding] = []
    for seed in seeds:
        ctx = _CtxShim()
        ctx.cfg.batch.batch_key = batch_key
        ctx.cfg.gpu_mode = gpu_mode
        ctx.cfg.random_state = int(seed)
        if scvi_max_epochs is not None:
            ctx.cfg.batch.scvi_max_epochs = scvi_max_epochs
        if scvi_n_latent is not None:
            ctx.cfg.batch.scvi_n_latent = scvi_n_latent
        if scvi_early_stopping is not None:
            ctx.cfg.batch.scvi_early_stopping = bool(scvi_early_stopping)
        method_name = f"{SCVI_METHOD}_seed{seed}"
        per_seed_key = f"{SCVI_OBSM}_seed{seed}"
        t0 = time.time()
        err = None
        X = None
        try:
            _seed_everything(int(seed))
            with _time_limit(per_seed_timeout_s, method_name):
                Backend._run_scvi(adata, batch_key, ctx)
            X = np.asarray(adata.obsm[SCVI_OBSM], dtype=np.float32)
            adata.obsm[per_seed_key] = X
        except Exception as exc:  # noqa: BLE001
            err = f"{type(exc).__name__}: {exc}"
            logger.error("scvi seed=%s failed: %s", seed, err)
        out.append(MethodEmbedding(
            method_name, per_seed_key, X, time.time() - t0,
            backend_meta={"scvi_input_source": ctx.metadata.get("scvi_input_source"),
                          "scvi_train_config": ctx.metadata.get("scvi_train_config")},
            seed=int(seed), error=err,
        ))
    return out


def _seed_everything(seed: int) -> None:
    """Seed scvi-tools + numpy + torch for an independent reproducible fit."""
    np.random.seed(seed)
    try:
        import scvi
        scvi.settings.seed = seed
    except Exception:  # noqa: BLE001 — scvi import failure surfaces in _run_scvi
        pass
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:  # noqa: BLE001
        pass


# ==========================================================================
# Band-aware ranking (scVI variance band => TIE, never a forced winner)
# ==========================================================================

@dataclass
class ScviBand:
    """Variance band of a scVI metric aggregate across the seed sweep."""

    metric: str                 # "batch_mixing_score" | "bio_conservation_score"
    values: list[float]
    mean: float
    std: float
    lo: float                   # mean - std (band floor)
    hi: float                   # mean + std (band ceiling)
    n_seeds: int


def summarize_scvi_band(per_seed_scores: list[float], *, metric: str) -> ScviBand | None:
    """Mean +/- 1 std band across the scVI seed sweep for one metric."""
    present = [float(v) for v in per_seed_scores if v is not None and not np.isnan(v)]
    if not present:
        return None
    arr = np.asarray(present, dtype=np.float64)
    mean = float(arr.mean())
    std = float(arr.std(ddof=1)) if arr.size > 1 else 0.0
    return ScviBand(metric=metric, values=present, mean=mean, std=std,
                    lo=mean - std, hi=mean + std, n_seeds=arr.size)


def band_aware_rank(
    point_scores: dict[str, float],
    *,
    scvi_band: ScviBand | None,
    scvi_label: str = SCVI_METHOD,
) -> list[dict[str, Any]]:
    """Rank methods on a single aggregate metric, BAND-AWARE for scVI.

    ``point_scores`` maps method->scalar for the deterministic methods
    (baseline/harmony/neg_control). scVI is represented by its band, not a
    point. A deterministic method whose point score falls INSIDE the scVI
    [lo,hi] band is declared a TIE with scVI (tie_with_scvi=True) — we never
    force a winner across the band overlap. Returns ranked rows (desc) with an
    explicit tie annotation; emits NO "winner" field.
    """
    rows: list[dict[str, Any]] = []
    for method, score in point_scores.items():
        if score is None:
            continue
        tie = False
        if scvi_band is not None and scvi_band.lo <= float(score) <= scvi_band.hi:
            tie = True
        rows.append({
            "method": method,
            "score": float(score),
            "tie_with_scvi_band": tie,
        })
    if scvi_band is not None:
        rows.append({
            "method": scvi_label,
            "score": scvi_band.mean,
            "band_lo": scvi_band.lo,
            "band_hi": scvi_band.hi,
            "band_std": scvi_band.std,
            "n_seeds": scvi_band.n_seeds,
            "is_band": True,
        })
    rows.sort(key=lambda r: r["score"], reverse=True)
    # Annotate ordinal rank but NEVER a winner verdict (calibration_relative_only).
    for i, r in enumerate(rows):
        r["relative_rank"] = i + 1
    return rows
