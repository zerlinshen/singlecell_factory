"""Label-agnostic discovery metric set for the integration-selection gate (plan §D1).

This module is a NEW layer over the reviewed scoring engine
(``score_integration.py`` + ``methods.py``); it does NOT modify either. It
implements the per-batch-cluster reference and the three label-free discovery
metrics that drive the per-run integration recommendation:

  (a) GOOD-MERGE cross-batch correspondence (COND #1, COND #5)
      Per-batch Leiden clusters are matched across batches by
      MUTUAL-NEAREST-CENTROID in the SHARED full-cohort pre-integration
      ``X_pca`` (per-batch PCA bases are incomparable, so the common baseline
      embedding is the only valid comparison space). A per-batch cluster with no
      mutual match is EXCLUDED from the denominator (never scored as a failed
      merge). The metric is the fraction of corresponded cross-batch pairs whose
      cells co-localize into ONE cluster of the embedding ``E`` under test.

  (b) WITHIN-BATCH DISTINCTNESS (PRIMARY discovery metric, COND #2)
      ``E`` restricted to each batch's cells, scored AGAINST that batch's FROZEN
      pre-integration per-batch Leiden labels. ``E`` is NEVER re-clustered to
      define its own reference (that would be circular).

  (c) WITHIN-BATCH ISOLATION (rare/novel anchor)
      Per-batch silhouette / isolated_labels on the FROZEN per-batch labels
      evaluated on ``E``'s geometry, aggregated across batches.

Determinism (COND #4): per-batch Leiden uses FIXED ``resolution=1.0``,
``random_state=0``, ``n_neighbors=15``, ``flavor="leidenalg"``. Per-batch labels
are computed ONCE and FROZEN. A version-assertion EXTENDS the engine's pinned
check to also assert ``leidenalg`` and ``igraph`` versions, fail-loud with the
``SC_ALLOW_SCIB_VERSION_DRIFT`` escape.

The old frozen-global ``L_ref`` ARI/iso quantity is retained as a LOGGED
SECONDARY cross-check only (``compute_global_lref_crosscheck``); it never drives
the gate.
"""
from __future__ import annotations

import importlib.metadata as _im
import logging
import os
from dataclasses import dataclass, field
from typing import Any

import numpy as np

# Reuse the reviewed engine (import only — never edit it).
from scripts.bench.integration.score_integration import (
    SC_ALLOW_SCIB_VERSION_DRIFT,
    assert_pinned_metric_engine,
    build_neighbors,
)

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# COND #4 — pinned Leiden determinism parameters. These are reproducibility
# invariants recorded in the audit. Do NOT vary per run.
# --------------------------------------------------------------------------
LEIDEN_RESOLUTION = 1.0
LEIDEN_RANDOM_STATE = 0
LEIDEN_N_NEIGHBORS = 15
LEIDEN_FLAVOR = "leidenalg"

# Resolution stability sweep (logged-secondary only, never gate-driving).
RESOLUTION_SWEEP = (0.5, 1.0, 2.0)

# COND #4 — Leiden reproducibility additionally depends on leidenalg + igraph.
# These pins EXTEND GATE G1 (which pins scib-metrics/chex/plottable). The dist
# name for the C-extension is ``igraph`` in this env; ``python-igraph`` is the
# legacy alias and is checked as a fallback so the assert is robust across envs.
PINNED_LEIDEN_VERSIONS: dict[str, str] = {
    "leidenalg": "0.11.0",
    "igraph": "1.0.0",
}
_IGRAPH_DIST_ALIASES = ("igraph", "python-igraph")

# A batch needs at least this many per-batch clusters for within-batch
# distinctness (D1b) to be defined (>=2 populations). Fewer -> loud note.
MIN_CLUSTERS_FOR_DISTINCTNESS = 2
# A per-batch cluster smaller than this is too small for a stable centroid /
# silhouette; flagged loud, never silently scored away.
MIN_CELLS_PER_CLUSTER = 10


# ==========================================================================
# Version gate extension (COND #4) — EXTENDS GATE G1
# ==========================================================================

def _installed_version_any(dist_aliases: tuple[str, ...]) -> str | None:
    for dist in dist_aliases:
        try:
            return _im.version(dist)
        except _im.PackageNotFoundError:
            continue
    return None


def assert_pinned_discovery_engine(*, allow_drift_env: bool = True) -> dict[str, str]:
    """Assert scib-metrics pins (GATE G1) PLUS leidenalg/igraph pins (COND #4).

    Returns the resolved version map (scib pins + leiden pins). Raises
    ``RuntimeError`` on any mismatch unless ``SC_ALLOW_SCIB_VERSION_DRIFT=1`` is
    set (opt-in, loud) — the SAME escape hatch the engine uses, so a single env
    var governs all reproducibility pins.
    """
    # First run the reviewed engine's scib/chex/plottable gate (fail-loud).
    resolved = dict(assert_pinned_metric_engine(allow_drift_env=allow_drift_env))

    problems: list[str] = []
    for dist, pinned in PINNED_LEIDEN_VERSIONS.items():
        aliases = _IGRAPH_DIST_ALIASES if dist == "igraph" else (dist,)
        got = _installed_version_any(aliases)
        resolved[dist] = got if got is not None else "MISSING"
        if got is None:
            problems.append(f"{dist}: NOT INSTALLED (pinned {pinned})")
        elif got != pinned:
            problems.append(f"{dist}: installed {got} != pinned {pinned}")

    if problems:
        opt_in = allow_drift_env and os.environ.get(
            SC_ALLOW_SCIB_VERSION_DRIFT, "").strip() == "1"
        msg = (
            "discovery gate: Leiden-engine version gate FAILED (COND #4) — "
            + "; ".join(problems)
            + ". Per-batch Leiden determinism depends on leidenalg + igraph "
            f"versions (pinned {PINNED_LEIDEN_VERSIONS}). Reinstall the pins "
            f"into sc_gpu, or set {SC_ALLOW_SCIB_VERSION_DRIFT}=1 to opt in to "
            "drift (NOT permitted for any reported run)."
        )
        if not opt_in:
            raise RuntimeError(msg)
        logger.error(
            "%s; %s=1 acknowledged — proceeding with drifted Leiden engine.",
            msg, SC_ALLOW_SCIB_VERSION_DRIFT,
        )
    return resolved


# ==========================================================================
# Per-batch Leiden reference (COND #2, COND #4) — computed ONCE and FROZEN
# ==========================================================================

@dataclass
class PerBatchReference:
    """Frozen per-batch Leiden partition + cross-batch correspondence.

    Attributes
    ----------
    batches:
        ordered batch ids (sorted for determinism).
    labels:
        ``{batch -> (n_cells_in_batch,) int array}`` of FROZEN per-batch Leiden
        labels (0-indexed within each batch). Computed once; never recomputed
        per embedding under test.
    cluster_keys:
        ``{batch -> [global_cluster_key, ...]}`` aligned to the sorted unique
        labels of that batch. A global key is ``f"{batch}::{local_label}"``.
    centroids:
        ``{global_cluster_key -> (n_pca_dims,) centroid}`` in the SHARED
        full-cohort pre-integration ``X_pca`` (COND #1 common space).
    sizes:
        ``{global_cluster_key -> n_cells}``.
    correspondence:
        list of ``(key_a, key_b)`` mutual-nearest-centroid matched cross-batch
        cluster pairs (deterministic; sorted). COND #5: clusters with no mutual
        match are absent here.
    matched_keys:
        set of global keys that participate in >=1 mutual match.
    notes:
        loud notes (e.g. a batch too small for distinctness).
    """

    batches: list[str]
    labels: dict[str, np.ndarray]
    cluster_keys: dict[str, list[str]]
    centroids: dict[str, np.ndarray]
    sizes: dict[str, int]
    correspondence: list[tuple[str, str]]
    matched_keys: set[str]
    notes: list[str] = field(default_factory=list)


def _leiden_labels_for_subset(
    X_sub: np.ndarray,
    *,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
    resolution: float = LEIDEN_RESOLUTION,
    random_state: int = LEIDEN_RANDOM_STATE,
) -> np.ndarray:
    """Deterministic Leiden labels for one batch's cells (within-batch space).

    Clusters the within-batch slice of the SHARED pre-integration embedding
    DIRECTLY (``use_rep`` on the slice itself — no PCA recompute). Recomputing a
    fresh PCA on a within-batch subset manufactures spurious noise axes that
    over-fragment isotropic populations; clustering the slice as-is reflects the
    genuine within-batch structure (the per-batch reference is method-independent
    of the integration under test). Builds the kNN graph with scanpy + runs
    ``sc.tl.leiden(flavor="leidenalg")`` with the pinned seed/resolution. Returns
    0-indexed integer labels. Fully deterministic for COND #4 / TM-16.
    """
    import anndata as ad
    import scanpy as sc

    n = X_sub.shape[0]
    if n < 2:
        return np.zeros(n, dtype=np.int64)

    sub = ad.AnnData(np.asarray(X_sub, dtype=np.float32).copy())
    sub.obsm["X_sub"] = np.asarray(X_sub, dtype=np.float32)
    k = int(min(n_neighbors, n - 1))
    sc.pp.neighbors(sub, n_neighbors=k, use_rep="X_sub", random_state=random_state)
    sc.tl.leiden(
        sub,
        resolution=resolution,
        random_state=random_state,
        flavor=LEIDEN_FLAVOR,
        directed=False,
        n_iterations=2,
    )
    return sub.obs["leiden"].astype(int).to_numpy().astype(np.int64)


def compute_per_batch_reference(
    X_pca: np.ndarray,
    batch: np.ndarray,
    *,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
    resolution: float = LEIDEN_RESOLUTION,
    random_state: int = LEIDEN_RANDOM_STATE,
) -> PerBatchReference:
    """Build the FROZEN per-batch Leiden reference + cross-batch correspondence.

    ``X_pca`` is the SHARED full-cohort PRE-INTEGRATION embedding (the baseline).
    Steps (COND #1):
      1. Per-batch Leiden on each batch's within-batch PCA (frozen labels).
      2. Each per-batch cluster -> centroid in the COMMON full-cohort ``X_pca``.
      3. Cross-batch match by MUTUAL-NEAREST-CENTROID; ties broken by lowest
         cluster index (deterministic). COND #5: unmatched clusters excluded.
    """
    X_pca = np.asarray(X_pca, dtype=np.float64)
    batch = np.asarray(batch)
    batches = sorted({str(b) for b in batch})

    labels: dict[str, np.ndarray] = {}
    cluster_keys: dict[str, list[str]] = {}
    centroids: dict[str, np.ndarray] = {}
    sizes: dict[str, int] = {}
    notes: list[str] = []

    batch_str = np.asarray([str(b) for b in batch])

    for b in batches:
        mask = batch_str == b
        idx = np.where(mask)[0]
        X_sub = X_pca[idx]
        lab = _leiden_labels_for_subset(
            X_sub, n_neighbors=n_neighbors, resolution=resolution,
            random_state=random_state)
        labels[b] = lab
        uniq = sorted(set(int(x) for x in lab))
        keys: list[str] = []
        for u in uniq:
            key = f"{b}::{u}"
            keys.append(key)
            member = idx[lab == u]
            sizes[key] = int(member.size)
            # COND #1: centroid in the COMMON full-cohort X_pca.
            centroids[key] = X_pca[member].mean(axis=0)
            if member.size < MIN_CELLS_PER_CLUSTER:
                notes.append(
                    f"per-batch cluster {key} has only {member.size} cells "
                    f"(< {MIN_CELLS_PER_CLUSTER}); centroid/silhouette low-power "
                    "(loud note, not silently dropped).")
        cluster_keys[b] = keys
        if len(uniq) < MIN_CLUSTERS_FOR_DISTINCTNESS:
            notes.append(
                f"batch {b} has only {len(uniq)} per-batch cluster(s); "
                "within-batch distinctness (D1b) is undefined for this batch "
                "(loud note).")

    correspondence, matched_keys = _mutual_nearest_centroid(
        batches, cluster_keys, centroids)

    return PerBatchReference(
        batches=batches, labels=labels, cluster_keys=cluster_keys,
        centroids=centroids, sizes=sizes, correspondence=correspondence,
        matched_keys=matched_keys, notes=notes,
    )


def _local_index(key: str) -> int:
    """Local cluster index from a global key ``batch::idx`` (tie-break key)."""
    return int(key.rsplit("::", 1)[1])


def _mutual_nearest_centroid(
    batches: list[str],
    cluster_keys: dict[str, list[str]],
    centroids: dict[str, np.ndarray],
) -> tuple[list[tuple[str, str]], set[str]]:
    """Mutual-nearest-centroid cross-batch matching (COND #1, deterministic).

    For every ordered pair of distinct batches, compute each cluster's nearest
    centroid in the other batch (Euclidean in the common ``X_pca``). A pair
    ``(i in A, j in B)`` is a match iff each is the other's nearest. Ties
    (equidistant centroids) are broken by LOWEST cluster index. Returns a sorted
    list of matched ``(key_a, key_b)`` (a<b lexicographically per pair) and the
    set of all matched keys. COND #5: a cluster with no mutual match is simply
    absent from the result.
    """
    matches: set[tuple[str, str]] = set()

    def nearest(src_key: str, dst_keys: list[str]) -> str | None:
        if not dst_keys:
            return None
        c = centroids[src_key]
        best_key = None
        best_d = None
        # Iterate in index order so ties resolve to the LOWEST cluster index.
        for dk in sorted(dst_keys, key=_local_index):
            d = float(np.sum((c - centroids[dk]) ** 2))
            if best_d is None or d < best_d:
                best_d = d
                best_key = dk
        return best_key

    for ai in range(len(batches)):
        for bi in range(ai + 1, len(batches)):
            ba, bb = batches[ai], batches[bi]
            keys_a, keys_b = cluster_keys[ba], cluster_keys[bb]
            for ka in keys_a:
                nb = nearest(ka, keys_b)
                if nb is None:
                    continue
                na = nearest(nb, keys_a)
                if na == ka:
                    # mutual match. order pair deterministically.
                    pair = tuple(sorted((ka, nb)))
                    matches.add(pair)  # type: ignore[arg-type]

    matched_keys: set[str] = set()
    for a, b in matches:
        matched_keys.add(a)
        matched_keys.add(b)
    return sorted(matches), matched_keys


# ==========================================================================
# D1a — GOOD-MERGE cross-batch correspondence on embedding E
# ==========================================================================

def _embedding_cluster_assignment(
    E: np.ndarray,
    *,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
    resolution: float = LEIDEN_RESOLUTION,
    random_state: int = LEIDEN_RANDOM_STATE,
) -> np.ndarray:
    """Cluster the WHOLE embedding E (all batches together) for co-localization.

    Used ONLY by D1a to ask "do corresponded cross-batch clusters co-localize
    into ONE E-cluster?" This is NOT a reference for D1b/D1c (those use frozen
    per-batch labels). Deterministic Leiden with the pinned seed/resolution.
    """
    import anndata as ad
    import scanpy as sc

    n = E.shape[0]
    if n < 2:
        return np.zeros(n, dtype=np.int64)
    a = ad.AnnData(np.asarray(E, dtype=np.float32).copy())
    # Use E directly as the representation (it is already an embedding).
    a.obsm["X_emb"] = np.asarray(E, dtype=np.float32)
    k = int(min(n_neighbors, n - 1))
    sc.pp.neighbors(a, n_neighbors=k, use_rep="X_emb", random_state=random_state)
    sc.tl.leiden(
        a, resolution=resolution, random_state=random_state,
        flavor=LEIDEN_FLAVOR, directed=False, n_iterations=2)
    return a.obs["leiden"].astype(int).to_numpy().astype(np.int64)


# A corresponded pair co-localizes iff their E-cluster distributions overlap by
# at least this much (overlap coefficient). Robust to Leiden over-fragmenting a
# single true population into several E-clusters: a GOOD merge places both
# per-batch copies into the SAME set of E-clusters even if that set has >1
# member, so the overlap is high; a BAD merge / non-merge sends them to disjoint
# E-clusters -> low overlap. Faithful to "co-localize into one cluster of E"
# while not brittle to even fragment splits.
CO_LOCALIZE_OVERLAP_THRESHOLD = 0.5


def good_merge_correspondence(
    E: np.ndarray,
    ref: PerBatchReference,
    batch: np.ndarray,
    *,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
    resolution: float = LEIDEN_RESOLUTION,
    random_state: int = LEIDEN_RANDOM_STATE,
    overlap_threshold: float = CO_LOCALIZE_OVERLAP_THRESHOLD,
) -> dict[str, Any]:
    """D1a: fraction of corresponded cross-batch pairs that co-localize in E.

    For each mutually-matched cross-batch cluster pair ``(key_a, key_b)``, the
    pair "co-localizes" iff their E-cluster occupancy distributions OVERLAP by at
    least ``overlap_threshold`` (overlap coefficient = sum of per-E-cluster
    minima of the two normalized distributions). This is robust to Leiden
    over-fragmenting a single true population into several E-clusters: a GOOD
    merge places both per-batch copies into the SAME E-clusters (high overlap)
    while a non-merge / bad merge sends them to disjoint E-clusters (low
    overlap). COND #5: only matched pairs are in the denominator; single-batch
    clusters are excluded. Returns the fraction + per-pair detail.
    """
    batch_str = np.asarray([str(b) for b in batch])
    E = np.asarray(E, dtype=np.float64)
    e_labels = _embedding_cluster_assignment(
        E, n_neighbors=n_neighbors, resolution=resolution,
        random_state=random_state)
    e_cluster_ids = sorted(set(int(x) for x in e_labels))

    # Map each global cluster key -> the global row indices of its cells.
    key_to_rows: dict[str, np.ndarray] = {}
    for b in ref.batches:
        mask = batch_str == b
        idx = np.where(mask)[0]
        lab = ref.labels[b]
        for u in sorted(set(int(x) for x in lab)):
            key_to_rows[f"{b}::{u}"] = idx[lab == u]

    def e_distribution(rows: np.ndarray) -> dict[int, float] | None:
        if rows.size == 0:
            return None
        vals, counts = np.unique(e_labels[rows], return_counts=True)
        total = float(counts.sum())
        return {int(v): float(c) / total for v, c in zip(vals, counts)}

    def overlap(da: dict[int, float], db: dict[int, float]) -> float:
        return float(sum(min(da.get(c, 0.0), db.get(c, 0.0)) for c in e_cluster_ids))

    pairs_detail: list[dict[str, Any]] = []
    n_co = 0
    for ka, kb in ref.correspondence:
        da = e_distribution(key_to_rows.get(ka, np.array([], dtype=int)))
        db = e_distribution(key_to_rows.get(kb, np.array([], dtype=int)))
        if da is None or db is None:
            ov = 0.0
        else:
            ov = overlap(da, db)
        co = ov >= overlap_threshold
        if co:
            n_co += 1
        pairs_detail.append({
            "pair": [ka, kb], "e_overlap": round(ov, 4),
            "co_localized": bool(co)})

    denom = len(ref.correspondence)
    frac = float(n_co) / denom if denom > 0 else None
    return {
        "good_merge_fraction": frac,
        "n_corresponded_pairs": denom,
        "n_co_localized": n_co,
        "overlap_threshold": overlap_threshold,
        "pairs": pairs_detail,
    }


# ==========================================================================
# D1b — WITHIN-BATCH DISTINCTNESS (primary) on frozen per-batch labels
# ==========================================================================

def within_batch_distinctness(
    E: np.ndarray,
    ref: PerBatchReference,
    batch: np.ndarray,
) -> dict[str, Any]:
    """D1b: per-batch silhouette of E vs the FROZEN per-batch labels (COND #2).

    E restricted to each batch's cells is scored against that batch's frozen
    pre-integration Leiden labels. A batch with <2 populations is undefined
    (loud note, excluded from the aggregate). Returns the cross-batch mean +
    per-batch detail. Higher = distinct populations stay distinct.
    """
    import scib_metrics as sm

    batch_str = np.asarray([str(b) for b in batch])
    E = np.asarray(E, dtype=np.float64)
    per_batch: dict[str, float | None] = {}
    notes: list[str] = []
    present: list[float] = []

    for b in ref.batches:
        mask = batch_str == b
        idx = np.where(mask)[0]
        lab = ref.labels[b]
        n_clusters = len(set(int(x) for x in lab))
        if n_clusters < MIN_CLUSTERS_FOR_DISTINCTNESS:
            per_batch[b] = None
            notes.append(
                f"D1b batch {b}: only {n_clusters} cluster(s); distinctness "
                "undefined (excluded from aggregate).")
            continue
        try:
            val = float(sm.silhouette_label(
                E[idx], lab.astype(str), rescale=True))
            per_batch[b] = val
            present.append(val)
        except Exception as exc:  # noqa: BLE001 — recorded loud
            per_batch[b] = None
            notes.append(f"D1b batch {b}: FAILED ({exc!r})")

    agg = float(np.mean(present)) if present else None
    return {"within_batch_distinctness": agg, "per_batch": per_batch,
            "notes": notes}


# ==========================================================================
# D1c — WITHIN-BATCH ISOLATION (rare/novel anchor) on frozen per-batch labels
# ==========================================================================

def within_batch_isolation(
    E: np.ndarray,
    ref: PerBatchReference,
    batch: np.ndarray,
) -> dict[str, Any]:
    """D1c: WORST-isolated frozen per-batch population on E's geometry.

    The rare/novel anchor. For each batch we compute the per-LABEL silhouette
    (rescaled to [0,1]) of every frozen per-batch population on E's geometry and
    take the MINIMUM — the WORST-isolated population in that batch. A mean over
    all populations would DROWN a rare type (plan §0 / A-PRIMARY): erasing one
    20-cell population barely moves a mean over ~440 cells, but it dominates the
    minimum. So a rare/novel type present within ONE batch that gets ERASED in E
    (its cells absorbed by a neighbor) collapses its own silhouette -> becomes
    the minimum -> D1c drops. D1c is the cross-batch MEAN of the per-batch
    worst-isolation, kept as its own gated quantity (COND #3).
    """
    from sklearn.metrics import silhouette_samples

    batch_str = np.asarray([str(b) for b in batch])
    E = np.asarray(E, dtype=np.float64)
    per_batch: dict[str, float | None] = {}
    per_batch_worst_label: dict[str, Any] = {}
    notes: list[str] = []
    present: list[float] = []

    for b in ref.batches:
        mask = batch_str == b
        idx = np.where(mask)[0]
        lab = ref.labels[b]
        uniq = sorted(set(int(x) for x in lab))
        if len(uniq) < 2:
            # Single-population batch: isolation is trivially undefined as a
            # contrast; record loud, exclude from aggregate.
            per_batch[b] = None
            notes.append(
                f"D1c batch {b}: only {len(uniq)} cluster(s); isolation "
                "undefined (excluded from aggregate).")
            continue
        try:
            s = silhouette_samples(E[idx], lab, metric="euclidean")
            per_label = {}
            for u in uniq:
                # rescale to [0,1] to match scib silhouette_label convention.
                per_label[u] = float((np.mean(s[lab == u]) + 1.0) / 2.0)
            worst_label = min(per_label, key=lambda k: per_label[k])
            worst = per_label[worst_label]
            per_batch[b] = worst
            per_batch_worst_label[b] = {
                "worst_label": int(worst_label),
                "worst_isolation": worst,
                "worst_label_size": int((lab == worst_label).sum()),
            }
            present.append(worst)
        except Exception as exc:  # noqa: BLE001 — recorded loud
            per_batch[b] = None
            notes.append(f"D1c batch {b}: FAILED ({exc!r})")

    agg = float(np.mean(present)) if present else None
    return {"within_batch_isolation": agg, "per_batch": per_batch,
            "per_batch_worst": per_batch_worst_label, "notes": notes}


# ==========================================================================
# Logged-secondary cross-check: frozen GLOBAL L_ref (never gate-driving)
# ==========================================================================

def compute_global_lref_crosscheck(
    E: np.ndarray,
    L_ref: np.ndarray,
    *,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
) -> dict[str, Any]:
    """LOGGED SECONDARY ONLY (never gate-driving): the v1 frozen-global L_ref.

    Retained purely so the audit shows the v1 quantity and its known
    batch-contamination bias for traceability. ARI/NMI of E (kmeans-clustered)
    vs a single global pre-integration partition L_ref.
    """
    import scib_metrics as sm

    out: dict[str, Any] = {"note": "LOGGED SECONDARY — never drives the gate "
                           "(known batch-contamination bias, plan §0)."}
    E = np.asarray(E, dtype=np.float64)
    L_ref = np.asarray(L_ref)
    try:
        d = sm.nmi_ari_cluster_labels_kmeans(E, L_ref)
        out["ari_vs_global_lref"] = float(d.get("ari")) if d.get("ari") is not None else None
        out["nmi_vs_global_lref"] = float(d.get("nmi")) if d.get("nmi") is not None else None
    except Exception as exc:  # noqa: BLE001
        out["error"] = repr(exc)
    return out


# ==========================================================================
# Per-embedding discovery bundle
# ==========================================================================

@dataclass
class DiscoveryMetrics:
    """All label-agnostic discovery metrics for ONE embedding E."""

    method: str
    # batch mixing (reuse engine's scib mixing aggregate; higher=better mixed)
    batch_mixing_score: float | None = None
    ilisi: float | None = None
    kbet_acceptance: float | None = None
    batch_asw: float | None = None
    # discovery metrics
    good_merge_fraction: float | None = None        # D1a (bonus / ranking)
    within_batch_distinctness: float | None = None   # D1b (gate)
    within_batch_isolation: float | None = None      # D1c (gate, COND #3)
    # detail
    d1a_detail: dict[str, Any] = field(default_factory=dict)
    d1b_detail: dict[str, Any] = field(default_factory=dict)
    d1c_detail: dict[str, Any] = field(default_factory=dict)
    lref_crosscheck: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def compute_batch_mixing(
    E: np.ndarray,
    batch: np.ndarray,
    *,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
) -> dict[str, float | None]:
    """Batch-mixing aggregate (iLISI + kBET) reusing the engine's kNN builder.

    Note: batch-ASW from the engine requires a cell-type label; for the
    label-agnostic gate we use only iLISI + kBET (both label-free) so mixing is
    never contaminated by labels. Higher = better mixed.
    """
    import scib_metrics as sm

    E = np.asarray(E, dtype=np.float64)
    batch_arr = np.asarray(batch)
    out: dict[str, float | None] = {"ilisi": None, "kbet_acceptance": None}
    neigh = build_neighbors(E, n_neighbors=n_neighbors)
    try:
        out["ilisi"] = float(sm.ilisi_knn(neigh, batch_arr))
    except Exception:  # noqa: BLE001
        out["ilisi"] = None
    try:
        out["kbet_acceptance"] = float(sm.kbet(neigh, batch_arr)[0])
    except Exception:  # noqa: BLE001
        out["kbet_acceptance"] = None
    present = [v for v in (out["ilisi"], out["kbet_acceptance"]) if v is not None]
    out["batch_mixing_score"] = float(np.mean(present)) if present else None
    return out


def score_embedding_discovery(
    method: str,
    E: np.ndarray,
    ref: PerBatchReference,
    batch: np.ndarray,
    *,
    global_lref: np.ndarray | None = None,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
) -> DiscoveryMetrics:
    """Score one embedding against the FROZEN per-batch reference (D1a/b/c)."""
    dm = DiscoveryMetrics(method=method)

    mix = compute_batch_mixing(E, batch, n_neighbors=n_neighbors)
    dm.ilisi = mix["ilisi"]
    dm.kbet_acceptance = mix["kbet_acceptance"]
    dm.batch_mixing_score = mix["batch_mixing_score"]

    d1a = good_merge_correspondence(E, ref, batch, n_neighbors=n_neighbors)
    dm.good_merge_fraction = d1a["good_merge_fraction"]
    dm.d1a_detail = d1a

    d1b = within_batch_distinctness(E, ref, batch)
    dm.within_batch_distinctness = d1b["within_batch_distinctness"]
    dm.d1b_detail = d1b
    dm.notes.extend(d1b.get("notes", []))

    d1c = within_batch_isolation(E, ref, batch)
    dm.within_batch_isolation = d1c["within_batch_isolation"]
    dm.d1c_detail = d1c
    dm.notes.extend(d1c.get("notes", []))

    if global_lref is not None:
        dm.lref_crosscheck = compute_global_lref_crosscheck(
            E, global_lref, n_neighbors=n_neighbors)

    return dm


def score_all_embeddings_discovery(
    embeddings: dict[str, np.ndarray],
    X_pca_baseline: np.ndarray,
    batch: np.ndarray,
    *,
    global_lref: np.ndarray | None = None,
    n_neighbors: int = LEIDEN_N_NEIGHBORS,
) -> tuple[list[DiscoveryMetrics], PerBatchReference]:
    """Build the frozen per-batch reference ONCE, score every embedding on it.

    ``X_pca_baseline`` is the shared full-cohort pre-integration PCA (COND #1
    common space) used to build the reference. The SAME reference object is used
    for every embedding -> swapping E never changes the reference (COND #2 /
    TM-14).
    """
    ref = compute_per_batch_reference(
        X_pca_baseline, batch, n_neighbors=n_neighbors)
    results: list[DiscoveryMetrics] = []
    for method, E in embeddings.items():
        results.append(score_embedding_discovery(
            method, E, ref, batch, global_lref=global_lref,
            n_neighbors=n_neighbors))
    return results, ref
