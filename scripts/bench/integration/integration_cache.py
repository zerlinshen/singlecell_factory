"""Default-ON cache for the discovery integration-selection gate (plan §D7).

"Per-run" means "per dataset VERSION". A literal per-run scVI retrain (3 seeds x
3600s on large LUSC cohorts) is infeasible, so the cache IS the per-run
semantics: a recommendation is keyed on
``(data_hash, batch_key, scVI_config, code_version)``. A cache HIT returns the
prior recommendation with NO scVI retrain; a MISS triggers a fresh evaluation.

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
GATE_CODE_VERSION = "discovery-gate-v1"

DEFAULT_SCVI_SEEDS: tuple[int, ...] = (0, 1, 2)


def compute_data_hash(X: np.ndarray) -> str:
    """Stable sha256 (truncated) of the baseline embedding bytes."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(X, dtype=np.float32).tobytes())
    return h.hexdigest()


@dataclass
class CacheKey:
    """The (data_hash, batch_key, scVI_config, code_version) cache key (§D7)."""

    data_hash: str
    batch_key: str
    scvi_config: str
    code_version: str = GATE_CODE_VERSION

    def digest(self) -> str:
        raw = "|".join([
            self.data_hash, self.batch_key, self.scvi_config, self.code_version])
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]

    def as_dict(self) -> dict[str, str]:
        return {
            "data_hash": self.data_hash,
            "batch_key": self.batch_key,
            "scvi_config": self.scvi_config,
            "code_version": self.code_version,
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
) -> CacheKey:
    """Build a CacheKey. scVI config encodes seeds + FAST/skip flags so a FAST
    or skip-scvi recommendation never aliases a full recommendation."""
    if skip_scvi:
        scvi_config = "skip_scvi"
    elif fast:
        scvi_config = f"fast_seed{scvi_seeds[0]}"
    else:
        scvi_config = "seeds=" + ",".join(str(s) for s in scvi_seeds)
    return CacheKey(
        data_hash=compute_data_hash(baseline_embedding),
        batch_key=batch_key,
        scvi_config=scvi_config,
        code_version=code_version,
    )


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
