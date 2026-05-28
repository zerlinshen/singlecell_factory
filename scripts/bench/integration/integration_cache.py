"""Default-ON cache for the discovery integration-selection gate (plan §D7).

"Per-run" means "per dataset VERSION". A literal per-run scVI retrain (3 seeds x
3600s on large LUSC cohorts) is infeasible, so the cache IS the per-run
semantics: a recommendation is keyed on
``(data_hash, batch_key, batch_state_hash, scVI_config, code_version,
gate_params)``. A cache HIT returns the prior recommendation with NO scVI
retrain; a MISS triggers a fresh evaluation.

``gate_params`` (W12) folds in the parameters that change the gate's
RECOMMENDATION — ``margin_mix``, Harmony ``theta``/``sigma``/``max_iter``,
scVI ``max_epochs``/``n_latent``/``early_stopping``, and ``n_neighbors`` — so a
run that differs ONLY in one of those does NOT alias to a stale recommendation.
When none are supplied the digest is identical to the legacy 4-tuple key
(back-compat).

Loud bypasses:
  - ``--force-refresh`` : ignore any cached entry, recompute, overwrite.
  - FAST mode           : single scVI seed, no band -> "NOT FOR FINAL CLAIM".
  - ``--skip-scvi``     : Harmony-vs-baseline only (GPU unavailable) -> loud.

This is a NEW layer; it does not modify the reviewed engine.
"""
from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# Bump when the gate's metric/selection logic changes so stale recommendations
# are invalidated (part of the cache key, plan §D7).
GATE_CODE_VERSION = "discovery-gate-v2"

DEFAULT_SCVI_SEEDS: tuple[int, ...] = (0, 1, 2)


def compute_data_hash(X: np.ndarray) -> str:
    """Stable sha256 (truncated) of the baseline embedding bytes."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(X, dtype=np.float32).tobytes())
    return h.hexdigest()


def compute_batch_state_hash(batch_labels, obs_names=None) -> str:
    """Stable hash of the batch labels and row order used by discovery scoring."""
    labels = [str(x) for x in np.asarray(batch_labels, dtype=object).tolist()]
    names = None if obs_names is None else [str(x) for x in list(obs_names)]
    payload = {
        "batch_labels": labels,
        "obs_names": names,
    }
    raw = json.dumps(payload, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@dataclass
class CacheKey:
    """The discovery-gate cache key (§D7).

    ``gate_params`` (W12 fix) folds the parameters that change the gate's
    RECOMMENDATION into the key: the mixing margin (``margin_mix``), the Harmony
    candidate's ``theta``/``sigma``/``max_iter``, the scVI training params
    (``scvi_max_epochs``, ``scvi_n_latent``, ``scvi_early_stopping``), and the
    per-batch-Leiden / downstream ``n_neighbors``.
    Previously the key hashed only (data_hash, batch_key, scvi_config,
    code_version), so two runs that differed ONLY in (say) ``margin_mix`` or
    ``harmony_theta`` aliased to the SAME cache entry and the second run silently
    replayed a recommendation computed under different gate parameters. Empty
    string preserves the legacy digest for callers that pass none of the new
    params (back-compat for the existing cache-key tests).
    """

    data_hash: str
    batch_key: str
    scvi_config: str
    code_version: str = GATE_CODE_VERSION
    gate_params: str = ""
    batch_state_hash: str = ""

    def digest(self) -> str:
        parts = [self.data_hash, self.batch_key, self.scvi_config, self.code_version]
        # Append optional fields only when present so legacy no-extra calls keep
        # their pre-W12 digest for compatibility tests and old persisted caches.
        if self.batch_state_hash:
            parts.append(f"batch_state={self.batch_state_hash}")
        # Append gate_params only when present so the legacy (no-gate-params)
        # digest is byte-identical to the pre-W12 key (back-compat).
        if self.gate_params:
            parts.append(self.gate_params)
        raw = "|".join(parts)
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]

    def as_dict(self) -> dict[str, str]:
        return {
            "data_hash": self.data_hash,
            "batch_key": self.batch_key,
            "scvi_config": self.scvi_config,
            "code_version": self.code_version,
            "gate_params": self.gate_params,
            "batch_state_hash": self.batch_state_hash,
            "digest": self.digest(),
        }


def make_cache_key(
    baseline_embedding: np.ndarray,
    batch_key: str,
    *,
    scvi_seeds: tuple[int, ...] = DEFAULT_SCVI_SEEDS,
    fast: bool = False,
    skip_scvi: bool = False,
    code_version: str = GATE_CODE_VERSION,
    margin_mix: float | None = None,
    harmony_theta: float | None = None,
    harmony_sigma: float | None = None,
    harmony_max_iter: int | None = None,
    scvi_max_epochs: int | None = None,
    scvi_n_latent: int | None = None,
    scvi_early_stopping: bool | None = None,
    n_neighbors: int | None = None,
    batch_labels=None,
    obs_names=None,
) -> CacheKey:
    """Build a CacheKey. scVI config encodes seeds + FAST/skip flags so a FAST
    or skip-scvi recommendation never aliases a full recommendation.

    W12 fix: the gate-affecting params (``margin_mix``, Harmony
    ``theta``/``sigma``/``max_iter``, scVI training params, ``n_neighbors``) are
    folded into the key via ``CacheKey.gate_params`` so a run that changes ANY
    of them gets a fresh recommendation instead of silently replaying a
    recommendation computed under different parameters. Any param left ``None``
    is omitted from the encoding; when ALL are ``None`` the digest is identical
    to the pre-W12 key (back-compat for callers/tests that pass none)."""
    if skip_scvi:
        scvi_config = "skip_scvi"
    elif fast:
        scvi_config = f"fast_seed{scvi_seeds[0]}"
    else:
        scvi_config = "seeds=" + ",".join(str(s) for s in scvi_seeds)
    gate_params = _encode_gate_params(
        margin_mix=margin_mix,
        harmony_theta=harmony_theta,
        harmony_sigma=harmony_sigma,
        harmony_max_iter=harmony_max_iter,
        scvi_max_epochs=scvi_max_epochs,
        scvi_n_latent=scvi_n_latent,
        scvi_early_stopping=scvi_early_stopping,
        n_neighbors=n_neighbors,
    )
    return CacheKey(
        data_hash=compute_data_hash(baseline_embedding),
        batch_key=batch_key,
        scvi_config=scvi_config,
        code_version=code_version,
        gate_params=gate_params,
        batch_state_hash=(
            compute_batch_state_hash(batch_labels, obs_names)
            if batch_labels is not None else ""
        ),
    )


def _encode_gate_params(
    *,
    margin_mix: float | None,
    harmony_theta: float | None,
    harmony_sigma: float | None,
    harmony_max_iter: int | None,
    scvi_max_epochs: int | None,
    scvi_n_latent: int | None,
    scvi_early_stopping: bool | None,
    n_neighbors: int | None,
) -> str:
    """Encode the gate-affecting params as a stable, order-fixed string.

    Only non-``None`` params are emitted (so the legacy all-``None`` call keeps
    the pre-W12 digest). Floats are rounded to 6 decimals so float-repr jitter
    never produces spurious cache misses.
    """
    fields: list[str] = []
    if margin_mix is not None:
        fields.append(f"mm={round(float(margin_mix), 6)}")
    if harmony_theta is not None:
        fields.append(f"ht={round(float(harmony_theta), 6)}")
    if harmony_sigma is not None:
        fields.append(f"hs={round(float(harmony_sigma), 6)}")
    if harmony_max_iter is not None:
        fields.append(f"hmi={int(harmony_max_iter)}")
    if scvi_max_epochs is not None:
        fields.append(f"sme={int(scvi_max_epochs)}")
    if scvi_n_latent is not None:
        fields.append(f"snl={int(scvi_n_latent)}")
    if scvi_early_stopping is not None:
        fields.append(f"ses={int(bool(scvi_early_stopping))}")
    if n_neighbors is not None:
        fields.append(f"nn={int(n_neighbors)}")
    return ";".join(fields)


@dataclass
class CacheResult:
    """Outcome of a cache lookup."""

    hit: bool
    key: CacheKey
    payload: dict[str, Any] | None = None
    path: Path | None = None
    notes: list[str] = field(default_factory=list)


class IntegrationCache:
    """Filesystem cache of prior recommendations, default-ON (§D7).

    Entries live under ``<cache_dir>/<digest>.json``. The cache dir is resolved
    UNDER the project root (never inside the factory tree).
    """

    def __init__(self, cache_dir: Path):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _path_for(self, key: CacheKey) -> Path:
        return self.cache_dir / f"{key.digest()}.json"

    def lookup(self, key: CacheKey, *, force_refresh: bool = False) -> CacheResult:
        """Return a CacheResult. ``force_refresh`` forces a MISS (loud)."""
        path = self._path_for(key)
        if force_refresh:
            logger.warning(
                "--force-refresh set: ignoring cache entry %s (recomputing).",
                path.name)
            return CacheResult(hit=False, key=key, path=path,
                               notes=["force_refresh: cache bypassed (loud)."])
        if path.is_file():
            try:
                payload = json.loads(path.read_text())
                logger.info("cache HIT %s — returning prior recommendation, "
                            "NO scVI retrain.", path.name)
                return CacheResult(hit=True, key=key, payload=payload, path=path,
                                   notes=["cache hit: no scVI retrain."])
            except Exception as exc:  # noqa: BLE001
                logger.warning("cache entry %s unreadable (%r); treating as MISS.",
                               path.name, exc)
        return CacheResult(hit=False, key=key, path=path,
                           notes=["cache miss: fresh evaluation required."])

    def store(self, key: CacheKey, payload: dict[str, Any]) -> Path:
        """Persist a recommendation payload under the key digest."""
        path = self._path_for(key)
        enriched = dict(payload)
        enriched["cache_key"] = key.as_dict()
        path.write_text(json.dumps(enriched, indent=2) + "\n")
        logger.info("cache STORE %s", path.name)
        return path


def loud_mode_banner(*, fast: bool, skip_scvi: bool, force_refresh: bool) -> list[str]:
    """Return loud-mode warning lines (also logged at WARNING)."""
    notes: list[str] = []
    if fast:
        msg = ("FAST mode: single scVI seed, NO variance band — NOT FOR FINAL "
               "CLAIM. Re-run without --fast for a reportable recommendation.")
        logger.warning(msg)
        notes.append(msg)
    if skip_scvi:
        msg = ("--skip-scvi: scVI EXCLUDED (GPU unavailable); recommendation is "
               "Harmony-vs-baseline ONLY. Loud: a scVI-eligible dataset was not "
               "fully evaluated.")
        logger.warning(msg)
        notes.append(msg)
    if force_refresh:
        msg = "--force-refresh: cache ignored; full recompute."
        logger.warning(msg)
        notes.append(msg)
    return notes


def resolve_cache_dir(project_root: Path, subdir: str = "discovery_gate/cache") -> Path:
    """Resolve + create the cache dir UNDER project_root (never the factory tree)."""
    factory_root = Path(__file__).resolve().parents[3]
    pr = Path(project_root).resolve()
    if factory_root in pr.parents or pr == factory_root:
        raise ValueError(
            f"--project-root {pr} is inside the factory tree {factory_root}; "
            "the cache must never live inside the factory (CLAUDE.md).")
    cache_dir = pr / subdir
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir
