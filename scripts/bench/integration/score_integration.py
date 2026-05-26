"""Integration calibration harness — scib-metrics scoreboard + audit.

CALIBRATION / RELATIVE lane ONLY. This harness scores batch-integration
embeddings (e.g. baseline PCA vs Harmony vs scVI) on a fixed dataset and emits
a *relative* scoreboard. It MUST NEVER emit a public "method A beats method B"
concordance verdict — concordance is gated behind a future curated multi-donor
PBMC dataset with author cell-type labels (GATE G6, see
.omc/plans/open-questions.md). The verdict_tier mechanism (T3) enforces this in
code.

Metric engine (GATE G1)
-----------------------
scib-metrics is the DEFAULT and authoritative engine. At import/run time we
ASSERT the installed scib-metrics version equals the pinned version and FAIL
LOUD on mismatch (mirroring the SC_ALLOW_* fail-loud idiom at
workflow/modular/modules/batch_correction.py:190-220). chex and plottable are
pinned transitive deps and are likewise asserted. Native sklearn metrics are
computed ONLY as a cross-check and are never the default source of truth.
Resolved jax / numpy / scikit-learn versions are RECORDED (not asserted) into
the audit JSON so environment drift is observable without being fatal.

Output shape (mirrors the hgmm benchmark convention)
---------------------------------------------------
  - ``<out_dir>/integration_scoreboard.csv`` : one row per scored method.
  - ``<out_dir>/integration_audit.json``     : machine-readable audit (env
    provenance, per-method raw metrics, native-vs-scib agreement, verdict
    tier).
  - ``<out_dir>/integration_audit.md``       : human-readable summary.

Run the built-in synthetic self-test (no external data required)::

    SC_ALLOW... not needed; just:
    python -m scripts.bench.integration.score_integration --self-test --out-dir /tmp/intg
"""
from __future__ import annotations

import argparse
import importlib.metadata as _im
import json
import platform
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

# --------------------------------------------------------------------------
# GATE G1: pinned metric-engine versions. These are the single source of truth
# for the version assertion below. Bump deliberately + re-run the fixture
# agreement check when changing.
# --------------------------------------------------------------------------
PINNED_VERSIONS: dict[str, str] = {
    "scib-metrics": "0.5.9",
    "chex": "0.1.91",
    "plottable": "0.1.5",
}

# Opt-in escape hatch, mirroring the SC_ALLOW_BATCH_BACKEND_SKIP idiom at
# workflow/modular/modules/batch_correction.py:190-220. Banned by default;
# only set with explicit human approval (e.g. a deliberate version bump in
# progress). When set, a version mismatch logs LOUDLY + is recorded in the
# audit instead of raising.
SC_ALLOW_SCIB_VERSION_DRIFT = "SC_ALLOW_SCIB_VERSION_DRIFT"

# Number of neighbors used to build the shared kNN graph that LISI/kBET
# consume. Fixed for determinism.
DEFAULT_N_NEIGHBORS = 15


# ==========================================================================
# Version gate (GATE G1)
# ==========================================================================

def _installed_version(dist_name: str) -> str | None:
    try:
        return _im.version(dist_name)
    except _im.PackageNotFoundError:
        return None


def assert_pinned_metric_engine(*, allow_drift_env: bool = True) -> dict[str, str]:
    """Assert installed metric-engine versions match the pins. FAIL LOUD.

    Returns the resolved (installed) version map. Raises ``RuntimeError`` on
    any mismatch or missing package unless ``SC_ALLOW_SCIB_VERSION_DRIFT=1`` is
    set (opt-in, loud).
    """
    import os

    resolved: dict[str, str] = {}
    problems: list[str] = []
    for dist, pinned in PINNED_VERSIONS.items():
        got = _installed_version(dist)
        resolved[dist] = got if got is not None else "MISSING"
        if got is None:
            problems.append(f"{dist}: NOT INSTALLED (pinned {pinned})")
        elif got != pinned:
            problems.append(f"{dist}: installed {got} != pinned {pinned}")

    if problems:
        opt_in = allow_drift_env and os.environ.get(
            SC_ALLOW_SCIB_VERSION_DRIFT, "").strip() == "1"
        msg = (
            "integration scoreboard: metric-engine version gate FAILED — "
            + "; ".join(problems)
            + f". scib-metrics is the DEFAULT engine and its version is a "
            f"reproducibility invariant. Reinstall the pins "
            f"(scib-metrics=={PINNED_VERSIONS['scib-metrics']}, "
            f"chex=={PINNED_VERSIONS['chex']}, "
            f"plottable=={PINNED_VERSIONS['plottable']}) into the sc_gpu env, "
            f"or set {SC_ALLOW_SCIB_VERSION_DRIFT}=1 to opt in to drift "
            "(NOT permitted for any reported run)."
        )
        if not opt_in:
            raise RuntimeError(msg)
        import logging
        logging.getLogger(__name__).error(
            "%s; %s=1 acknowledged — proceeding with drifted engine.",
            msg, SC_ALLOW_SCIB_VERSION_DRIFT,
        )
    return resolved


def record_resolved_runtime_versions() -> dict[str, str]:
    """RECORD (do NOT assert) jax / numpy / scikit-learn versions.

    These can shift with the environment without invalidating the run; they are
    captured for provenance only (GATE G1 distinguishes asserted pins from
    recorded runtime deps).
    """
    out: dict[str, str] = {}
    for dist in ("jax", "numpy", "scikit-learn"):
        out[dist] = _installed_version(dist) or "MISSING"
    out["python"] = platform.python_version()
    return out


# ==========================================================================
# Shared kNN graph (deterministic; feeds BOTH scib and native cross-checks)
# ==========================================================================

def build_neighbors(X: np.ndarray, n_neighbors: int = DEFAULT_N_NEIGHBORS):
    """Build a deterministic kNN graph as a scib-metrics ``NeighborsResults``.

    Uses sklearn's exact NearestNeighbors so the graph is reproducible and is
    the *same* graph both metric engines consume (apples-to-apples agreement).
    scib-metrics expects indices INCLUDING self as the first neighbor.
    """
    from sklearn.neighbors import NearestNeighbors
    from scib_metrics.nearest_neighbors import NeighborsResults

    n = X.shape[0]
    k = min(n_neighbors, n - 1)
    nn = NearestNeighbors(n_neighbors=k, algorithm="brute", metric="euclidean")
    nn.fit(X)
    distances, indices = nn.kneighbors(X, return_distance=True)
    return NeighborsResults(indices=indices, distances=distances)


# ==========================================================================
# Metric layer
# ==========================================================================

@dataclass
class MetricResult:
    """Per-method metric bundle for one embedding."""

    method: str
    # batch mixing (higher = better mixed)
    ilisi: float | None = None
    kbet_acceptance: float | None = None
    batch_asw: float | None = None
    # bio conservation (higher = better preserved)
    clisi: float | None = None
    celltype_asw: float | None = None
    ari: float | None = None
    nmi: float | None = None
    # aggregates
    batch_mixing_score: float | None = None
    bio_conservation_score: float | None = None
    # native cross-check raw values
    native: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def _safe(fn: Callable[[], float], notes: list[str], label: str) -> float | None:
    try:
        val = float(fn())
        if np.isnan(val):
            notes.append(f"{label}: NaN")
            return None
        return val
    except Exception as exc:  # noqa: BLE001 - metric robustness, recorded loudly
        notes.append(f"{label}: FAILED ({exc!r})")
        return None


def compute_scib_metrics(
    X: np.ndarray,
    batch: np.ndarray,
    cell_type: np.ndarray,
    *,
    method: str,
    n_neighbors: int = DEFAULT_N_NEIGHBORS,
) -> MetricResult:
    """Compute the default (scib-metrics) batch-mixing + bio-conservation set.

    Batch mixing : iLISI, kBET acceptance, batch-ASW.
    Bio conservation : cLISI, celltype-ASW, ARI, NMI (kmeans-clustered).
    """
    import scib_metrics as sm

    res = MetricResult(method=method)
    notes = res.notes

    neigh = build_neighbors(X, n_neighbors=n_neighbors)
    batch_arr = np.asarray(batch)
    ct_arr = np.asarray(cell_type)

    # --- batch mixing ---
    res.ilisi = _safe(lambda: sm.ilisi_knn(neigh, batch_arr), notes, "ilisi")
    # kbet returns (acceptance_rate, stats, pvalues); element 0 is the scalar
    # acceptance rate in [0,1] (higher = better mixed).
    res.kbet_acceptance = _safe(
        lambda: sm.kbet(neigh, batch_arr)[0], notes, "kbet")
    res.batch_asw = _safe(
        lambda: sm.silhouette_batch(X, ct_arr, batch_arr), notes, "batch_asw")

    # --- bio conservation ---
    res.clisi = _safe(lambda: sm.clisi_knn(neigh, ct_arr), notes, "clisi")
    res.celltype_asw = _safe(
        lambda: sm.silhouette_label(X, ct_arr), notes, "celltype_asw")
    nmi_ari = _safe_dict(
        lambda: sm.nmi_ari_cluster_labels_kmeans(X, ct_arr), notes, "nmi_ari")
    if nmi_ari is not None:
        res.ari = float(nmi_ari.get("ari")) if nmi_ari.get("ari") is not None else None
        res.nmi = float(nmi_ari.get("nmi")) if nmi_ari.get("nmi") is not None else None

    res.batch_mixing_score = _mean_present([res.ilisi, res.kbet_acceptance, res.batch_asw])
    res.bio_conservation_score = _mean_present(
        [res.clisi, res.celltype_asw, res.ari, res.nmi])
    return res


def _safe_dict(fn: Callable[[], dict], notes: list[str], label: str) -> dict | None:
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001
        notes.append(f"{label}: FAILED ({exc!r})")
        return None


def _mean_present(vals: list[float | None]) -> float | None:
    present = [v for v in vals if v is not None]
    if not present:
        return None
    return float(np.mean(present))


def compute_native_crosscheck(
    X: np.ndarray,
    batch: np.ndarray,
    cell_type: np.ndarray,
) -> dict[str, float]:
    """Native sklearn cross-check (NEVER the default source of truth).

    Provides independent estimates for batch-ASW, celltype-ASW, ARI, NMI so we
    can report agreement vs the scib-metrics engine on the fixture.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import (
        adjusted_rand_score,
        normalized_mutual_info_score,
        silhouette_samples,
    )

    ct_arr = np.asarray(cell_type)
    batch_arr = np.asarray(batch)
    out: dict[str, float] = {}

    # Native ASW: silhouette over celltype labels, rescaled to [0,1] to match
    # scib's silhouette_label convention ((s+1)/2).
    try:
        sil_ct = silhouette_samples(X, ct_arr, metric="euclidean")
        out["celltype_asw_native"] = float((np.mean(sil_ct) + 1.0) / 2.0)
    except Exception:
        pass
    # Native batch-ASW: scib computes per-celltype, 1-|silhouette| over batch,
    # rescaled. Mirror that: for each celltype group, 1-|s(batch)| then average.
    try:
        accum = []
        for ct in np.unique(ct_arr):
            mask = ct_arr == ct
            if mask.sum() < 3 or len(np.unique(batch_arr[mask])) < 2:
                continue
            s = silhouette_samples(X[mask], batch_arr[mask], metric="euclidean")
            accum.append(float(np.mean(1.0 - np.abs(s))))
        if accum:
            out["batch_asw_native"] = float(np.mean(accum))
    except Exception:
        pass
    # Native ARI/NMI from KMeans(k=n_celltypes).
    try:
        k = len(np.unique(ct_arr))
        km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
        out["ari_native"] = float(adjusted_rand_score(ct_arr, km.labels_))
        out["nmi_native"] = float(
            normalized_mutual_info_score(ct_arr, km.labels_))
    except Exception:
        pass
    return out


def native_scib_agreement(res: MetricResult) -> dict[str, float]:
    """Absolute differences between scib defaults and native cross-check."""
    n = res.native
    agree: dict[str, float] = {}
    pairs = [
        ("celltype_asw", res.celltype_asw, n.get("celltype_asw_native")),
        ("batch_asw", res.batch_asw, n.get("batch_asw_native")),
        ("ari", res.ari, n.get("ari_native")),
        ("nmi", res.nmi, n.get("nmi_native")),
    ]
    for label, scib_v, native_v in pairs:
        if scib_v is not None and native_v is not None:
            agree[f"{label}_abs_diff"] = abs(float(scib_v) - float(native_v))
    return agree


# ==========================================================================
# GATE G3: verdict_tier — enforced contrast classification
# ==========================================================================
# Pinned mapping from a comparison's *contrast type* to its allowed verdict
# tier. This is a reproducibility invariant: a global / sample-level contrast
# (or one whose labels are derived from the pipeline itself) can only support a
# RELATIVE RANKING, never a concordance claim, because there is no independent
# ground truth to be concordant *with*. A within-age replicate-pair contrast
# (matched age band, author/ground-truth cell-type labels) supports the
# stronger CONDITIONED CONCORDANCE tier — concordance conditioned on the age
# band. No contrast in this calibration lane may be elevated above
# "conditioned_concordance" (public Harmony-vs-scVI concordance is GATE G6,
# deferred to a future curated PBMC dataset).
VERDICT_TIER_MAP: dict[str, str] = {
    "global": "relative_ranking_only",
    "sample_level": "relative_ranking_only",
    "pipeline_derived_labels": "relative_ranking_only",
    "within_age_replicate_pair": "conditioned_concordance",
}

# Pre-registered OVER-CORRECTING threshold T. A corrected method whose bio
# aggregate falls more than T BELOW the uncorrected baseline is flagged as
# over-correcting, REGARDLESS of how well it mixes batches — a method that
# shreds biology to win on mixing is failing, not winning.
#
# How T is set (pre-registered, not tuned post-hoc): T is the bio-conservation
# spread between the uncorrected baseline and the negative control (a
# label-permuted / structure-destroying control). That spread is the largest
# bio drop attributable to NON-biological causes on this dataset; any real
# method dropping more than that is destroying genuine signal. T defaults to a
# conservative floor (0.05) when no negative control is supplied, and is
# recomputed as ``baseline_bio - neg_control_bio`` (clipped to >= floor) when a
# neg-control method is present. Documented in the audit JSON under
# ``over_correcting_threshold`` with its provenance.
OVER_CORRECTING_THRESHOLD_FLOOR = 0.05


class BioConservationLabelError(RuntimeError):
    """GATE G2: bio-conservation label column missing or join is invalid."""


@dataclass(frozen=True)
class ContrastSpec:
    """Describes ONE comparison so its verdict tier can be classified.

    Attributes
    ----------
    contrast_type:
        One of the keys in VERDICT_TIER_MAP. Unknown values fail loud.
    labels_source:
        "author"/"ground_truth" or "pipeline" — recorded for provenance and
        used as an extra guard: a contrast that claims concordance
        ("within_age_replicate_pair") while using pipeline-derived labels is
        downgraded to relative_ranking_only (you cannot be concordant with
        your own clustering output).
    """

    contrast_type: str
    labels_source: str = "author"


def classify_verdict_tier(spec: ContrastSpec) -> str:
    """Return the enforced verdict tier for a contrast (GATE G3). Fail loud.

    Raises ``ValueError`` for an unknown contrast_type so a typo can never
    silently produce an over-strong tier.
    """
    if spec.contrast_type not in VERDICT_TIER_MAP:
        raise ValueError(
            f"unknown contrast_type {spec.contrast_type!r}; allowed: "
            f"{sorted(VERDICT_TIER_MAP)}"
        )
    tier = VERDICT_TIER_MAP[spec.contrast_type]
    # Concordance requires INDEPENDENT labels. If the labels came from the
    # pipeline itself, force the weaker tier regardless of contrast_type.
    if tier == "conditioned_concordance" and spec.labels_source not in {
        "author", "ground_truth",
    }:
        return "relative_ranking_only"
    return tier


def require_bio_conservation_labels(
    obs_columns,
    label_column: str,
    *,
    cluster_names_join: dict | None = None,
    cluster_names_path: "Path | None" = None,
) -> None:
    """GATE G2 fail-loud guard for the bio-conservation label join.

    Raises ``BioConservationLabelError`` if:
      - ``label_column`` is absent from ``obs_columns`` (no ground-truth
        cell-type labels => bio conservation is unmeasurable, must not be
        silently scored as 0/None), OR
      - a ``cluster_names -> celltype`` join is attempted through a path that
        is a DEAD SYMLINK (or otherwise non-resolvable), OR an empty/None join
        map — a broken external annotation must fail loud, never join to NaN.
    """
    cols = set(obs_columns)
    if label_column not in cols:
        raise BioConservationLabelError(
            f"GATE G2: bio-conservation label column {label_column!r} absent "
            f"from obs (have: {sorted(cols)}). Without ground-truth cell-type "
            "labels the bio-conservation metrics are unmeasurable and MUST NOT "
            "be silently scored — supply author/ground-truth labels or use a "
            "dataset that carries them."
        )
    # If a cluster_names->celltype join is requested, the annotation source
    # must resolve. A dead symlink is the classic silent-NaN trap.
    if cluster_names_path is not None:
        p = Path(cluster_names_path)
        if p.is_symlink() and not p.exists():
            raise BioConservationLabelError(
                f"GATE G2: cluster_names annotation {p} is a DEAD SYMLINK; "
                "the cluster_names->celltype join would map every cluster to "
                "NaN. Fix the symlink target before scoring bio conservation."
            )
        if not p.exists():
            raise BioConservationLabelError(
                f"GATE G2: cluster_names annotation {p} does not exist; "
                "the cluster_names->celltype join cannot be performed."
            )
    if cluster_names_join is not None and len(cluster_names_join) == 0:
        raise BioConservationLabelError(
            "GATE G2: cluster_names->celltype join map is empty; the join "
            "would map every cluster to NaN. Refuse to score against an empty "
            "annotation."
        )


def pre_register_over_correcting_threshold(
    baseline_bio: float | None,
    neg_control_bio: float | None,
    *,
    floor: float = OVER_CORRECTING_THRESHOLD_FLOOR,
) -> tuple[float, str]:
    """Return (T, provenance) for the OVER-CORRECTING flag.

    T = max(floor, baseline_bio - neg_control_bio) when a neg control is
    present; otherwise the conservative floor. Pre-registered: it is a function
    of the baseline + neg-control spread, NOT tuned against the methods under
    test.
    """
    if baseline_bio is not None and neg_control_bio is not None:
        spread = float(baseline_bio) - float(neg_control_bio)
        T = max(floor, spread)
        prov = (
            f"baseline_bio({baseline_bio:.4f}) - neg_control_bio"
            f"({neg_control_bio:.4f}) = {spread:.4f}, clipped to floor {floor}"
            f" => T={T:.4f}"
        )
        return T, prov
    return floor, (
        f"no negative control supplied; conservative floor T={floor} "
        "(pre-registered default)"
    )


def flag_over_correcting(
    results: list[MetricResult],
    *,
    baseline_method: str = "baseline",
    neg_control_method: str | None = None,
) -> tuple[list[str], dict[str, Any]]:
    """Flag corrected methods that over-correct (bio drop > T below baseline).

    A method is flagged if its bio-conservation aggregate is more than T below
    the uncorrected baseline's, REGARDLESS of its batch-mixing score. Returns
    (flag_messages, threshold_audit).
    """
    by_method = {r.method: r for r in results}
    baseline = by_method.get(baseline_method)
    baseline_bio = baseline.bio_conservation_score if baseline else None

    neg_bio = None
    if neg_control_method and neg_control_method in by_method:
        neg_bio = by_method[neg_control_method].bio_conservation_score

    T, provenance = pre_register_over_correcting_threshold(baseline_bio, neg_bio)

    flags: list[str] = []
    if baseline_bio is None:
        return flags, {
            "over_correcting_threshold": T,
            "threshold_provenance": provenance,
            "note": "no baseline method present; over-correcting check skipped",
        }

    for r in results:
        if r.method in {baseline_method, neg_control_method}:
            continue
        if r.bio_conservation_score is None:
            continue
        drop = float(baseline_bio) - float(r.bio_conservation_score)
        if drop > T:
            flags.append(
                f"OVER-CORRECTING: method '{r.method}' bio-conservation "
                f"{r.bio_conservation_score:.4f} is {drop:.4f} below baseline "
                f"{baseline_bio:.4f} (> T={T:.4f}); flagged regardless of "
                f"batch-mixing ({r.batch_mixing_score})."
            )
    return flags, {
        "over_correcting_threshold": T,
        "threshold_provenance": provenance,
        "baseline_bio_conservation": baseline_bio,
        "neg_control_bio_conservation": neg_bio,
    }


# ==========================================================================
# Output shape: scoreboard CSV + audit JSON + audit MD
# ==========================================================================

_SCOREBOARD_COLUMNS = [
    "method",
    "ilisi", "kbet_acceptance", "batch_asw", "batch_mixing_score",
    "clisi", "celltype_asw", "ari", "nmi", "bio_conservation_score",
]


def write_scoreboard_csv(results: list[MetricResult], path: Path) -> None:
    import csv

    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(_SCOREBOARD_COLUMNS)
        for r in results:
            writer.writerow([
                r.method,
                _fmt(r.ilisi), _fmt(r.kbet_acceptance), _fmt(r.batch_asw),
                _fmt(r.batch_mixing_score),
                _fmt(r.clisi), _fmt(r.celltype_asw), _fmt(r.ari), _fmt(r.nmi),
                _fmt(r.bio_conservation_score),
            ])


def _fmt(v: float | None) -> str:
    return "" if v is None else f"{v:.6f}"


def build_audit(
    results: list[MetricResult],
    resolved_pins: dict[str, str],
    runtime_versions: dict[str, str],
    *,
    contrast: ContrastSpec | None = None,
    baseline_method: str = "baseline",
    neg_control_method: str | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    audit: dict[str, Any] = {
        "schema_version": "integration-bench-v1",
        "computed_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "lane": "calibration_relative_only",
        "engine": {
            "default": "scib-metrics",
            "pinned_versions": PINNED_VERSIONS,
            "resolved_pins": resolved_pins,
            "cross_check": "sklearn (non-default)",
        },
        "recorded_runtime_versions": runtime_versions,
        "n_neighbors": DEFAULT_N_NEIGHBORS,
        "methods": [],
        "native_scib_agreement": {},
    }
    # GATE G3: verdict tier. Default to the most conservative contrast (global)
    # when none is supplied so the audit can NEVER silently omit the tier.
    spec = contrast or ContrastSpec(contrast_type="global", labels_source="pipeline")
    audit["contrast"] = {
        "contrast_type": spec.contrast_type,
        "labels_source": spec.labels_source,
    }
    audit["verdict_tier"] = classify_verdict_tier(spec)

    # OVER-CORRECTING flags + pre-registered threshold provenance.
    flags, threshold_audit = flag_over_correcting(
        results,
        baseline_method=baseline_method,
        neg_control_method=neg_control_method,
    )
    audit["over_correcting_flags"] = flags
    audit["over_correcting"] = threshold_audit

    for r in results:
        audit["methods"].append({
            "method": r.method,
            "batch_mixing": {
                "ilisi": r.ilisi,
                "kbet_acceptance": r.kbet_acceptance,
                "batch_asw": r.batch_asw,
                "aggregate": r.batch_mixing_score,
            },
            "bio_conservation": {
                "clisi": r.clisi,
                "celltype_asw": r.celltype_asw,
                "ari": r.ari,
                "nmi": r.nmi,
                "aggregate": r.bio_conservation_score,
            },
            "native_crosscheck": r.native,
            "notes": r.notes,
        })
        audit["native_scib_agreement"][r.method] = native_scib_agreement(r)
    if extra:
        audit.update(extra)
    return audit


def write_audit_json(audit: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(audit, indent=2) + "\n")


def write_audit_md(audit: dict[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Integration Calibration Audit")
    lines.append("")
    lines.append(f"- computed: `{audit['computed_iso']}`")
    lines.append(f"- lane: **{audit['lane']}** (relative ranking; NOT a "
                 "Harmony-vs-scVI concordance verdict)")
    lines.append(f"- engine: **{audit['engine']['default']}** "
                 f"(pinned {audit['engine']['pinned_versions']})")
    lines.append(f"- recorded runtime deps: {audit['recorded_runtime_versions']}")
    if "verdict_tier" in audit:
        lines.append(f"- verdict_tier: **{audit['verdict_tier']}**")
    lines.append("")
    lines.append("## Scoreboard")
    lines.append("")
    lines.append("| method | iLISI | kBET | batch-ASW | mixing | "
                 "cLISI | ct-ASW | ARI | NMI | bio |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for m in audit["methods"]:
        bm, bc = m["batch_mixing"], m["bio_conservation"]
        lines.append(
            f"| {m['method']} | {_md(bm['ilisi'])} | {_md(bm['kbet_acceptance'])} "
            f"| {_md(bm['batch_asw'])} | {_md(bm['aggregate'])} "
            f"| {_md(bc['clisi'])} | {_md(bc['celltype_asw'])} | {_md(bc['ari'])} "
            f"| {_md(bc['nmi'])} | {_md(bc['aggregate'])} |")
    lines.append("")
    lines.append("## Native (sklearn) vs scib-metrics agreement (abs diff)")
    lines.append("")
    for method, agree in audit["native_scib_agreement"].items():
        if agree:
            parts = ", ".join(f"{k}={v:.4f}" for k, v in agree.items())
            lines.append(f"- **{method}**: {parts}")
    if audit.get("over_correcting_flags"):
        lines.append("")
        lines.append("## OVER-CORRECTING flags")
        lines.append("")
        for f in audit["over_correcting_flags"]:
            lines.append(f"- {f}")
    lines.append("")
    path.write_text("\n".join(lines))


def _md(v: float | None) -> str:
    return "—" if v is None else f"{v:.4f}"


# ==========================================================================
# Self-test on the synthetic fixture
# ==========================================================================

def score_embeddings(
    embeddings: dict[str, np.ndarray],
    batch: np.ndarray,
    cell_type: np.ndarray,
    *,
    n_neighbors: int = DEFAULT_N_NEIGHBORS,
    with_native: bool = True,
) -> list[MetricResult]:
    """Score a set of named embeddings against shared batch/celltype labels."""
    results: list[MetricResult] = []
    for method, X in embeddings.items():
        X = np.asarray(X, dtype=np.float64)
        res = compute_scib_metrics(
            X, batch, cell_type, method=method, n_neighbors=n_neighbors)
        if with_native:
            res.native = compute_native_crosscheck(X, batch, cell_type)
        results.append(res)
    return results


def load_synthetic_fixture(*, seed: int = 0):
    """Import + build the synthetic 2-batch fixture by file path.

    ``tests`` is a pytest namespace dir (no top-level ``__init__.py``), so a
    plain ``import tests.fixtures...`` is unreliable when this harness runs as a
    standalone module. Load the generator directly from its file location.
    """
    import importlib.util
    import sys

    repo_root = Path(__file__).resolve().parents[3]
    fixture_path = repo_root / "tests" / "fixtures" / "integration_synthetic.py"
    mod_name = "integration_synthetic_fixture"
    spec = importlib.util.spec_from_file_location(mod_name, fixture_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load fixture at {fixture_path}")
    mod = importlib.util.module_from_spec(spec)
    # Register before exec so the frozen dataclass can resolve its own module
    # under `from __future__ import annotations`.
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod.make_integration_fixture(seed=seed)


def run_self_test(out_dir: Path, *, seed: int = 0) -> dict[str, Any]:
    """Score the synthetic fixture and write the scoreboard + audit.

    The fixture's ``X_integrated`` should out-mix ``X_unintegrated`` while
    preserving bio conservation. Returns the audit dict.
    """
    fx = load_synthetic_fixture(seed=seed)
    embeddings = {
        "unintegrated": fx.X_unintegrated,
        "integrated": fx.X_integrated,
    }
    results = score_embeddings(embeddings, fx.batch, fx.cell_type)

    resolved_pins = assert_pinned_metric_engine()
    runtime_versions = record_resolved_runtime_versions()
    audit = build_audit(results, resolved_pins, runtime_versions,
                        extra={"fixture_seed": seed,
                               "n_cell_types": fx.n_cell_types})

    out_dir.mkdir(parents=True, exist_ok=True)
    write_scoreboard_csv(results, out_dir / "integration_scoreboard.csv")
    write_audit_json(audit, out_dir / "integration_audit.json")
    write_audit_md(audit, out_dir / "integration_audit.md")
    return audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true",
                        help="Run the synthetic-fixture self-test.")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Directory for scoreboard CSV + audit md/json.")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    # GATE G1 always runs first — fail loud before any compute.
    assert_pinned_metric_engine()

    if args.self_test:
        audit = run_self_test(args.out_dir, seed=args.seed)
        print(f"[integration] self-test scoreboard -> "
              f"{args.out_dir / 'integration_scoreboard.csv'}")
        for m in audit["methods"]:
            print(f"  {m['method']}: mixing="
                  f"{m['batch_mixing']['aggregate']}, "
                  f"bio={m['bio_conservation']['aggregate']}")
        return 0

    parser.error("no action requested (use --self-test for now)")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
