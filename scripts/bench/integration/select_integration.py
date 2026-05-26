"""Per-run discovery integration-selection gate — recommendation logic + CLI (§D4/§D5).

This is a NEW layer over the reviewed engine (``score_integration.py`` /
``methods.py``) and the discovery metrics (``discovery_metrics.py``). It
implements the pre-registered DISQUALIFYING-GATE recommendation rule and a CLI
that selects ``cfg.batch.method`` AFTER clustering and BEFORE batch_correction.

Recommendation rule (plan §D4; floors are GATES, bonus only ranks survivors):
  1. MIXING FLOOR (gate): candidate batch-mixing > baseline mixing + MARGIN_MIX.
  2. WITHIN-BATCH DISTINCTNESS FLOOR (gate, D1b): >= pre-registered floor.
  3. WITHIN-BATCH ISOLATION FLOOR (gate, D1c, COND #3): >= pre-registered floor.
     A candidate below D1c floor is REJECTED EVEN IF its D1a is the highest.
  4. RANK survivors by D1a good-merge correspondence (bonus only).
  5. Band-aware ties via band_aware_rank / ScviBand.
  6. If no candidate survives gates 1-3 -> recommend baseline (loudly).
  7. verdict_tier = relative_ranking_only; scope = "for THIS dataset".

Over-correction anchor policy (plan §D3): the SHUFFLE control is the REQUIRED
falsifiability anchor (must collapse); the extreme-theta Harmony control is
REPORTED-BUT-NEVER-GATING (known-weak on real data). Gating keys on the shuffle
control.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from scripts.bench.integration import discovery_metrics as dm
from scripts.bench.integration.methods import (
    BASELINE_METHOD,
    HARMONY_METHOD,
    NEG_CONTROL_METHOD,
    SCVI_METHOD,
    band_aware_rank,
    summarize_scvi_band,
)
from scripts.bench.integration.score_integration import (
    OVER_CORRECTING_THRESHOLD_FLOOR,
    pre_register_over_correcting_threshold,
)

logger = logging.getLogger(__name__)

# Pre-registered margin the candidate batch-mixing must beat baseline by
# (gate 1). Conservative default; recorded in audit. A larger margin demands a
# more decisive mixing gain before recommending integration.
MARGIN_MIX = 0.05

# Pre-registered floor position along the baseline->neg-control spread (gates 2
# & 3). The disqualifying floor is ``baseline - FLOOR_FRACTION * spread`` where
# ``spread`` is the engine's pre-registered baseline-minus-neg-control quantity
# (pre_register_over_correcting_threshold). FLOOR_FRACTION=0.5 places the line at
# the MIDPOINT between a structure-preserving baseline and the
# structure-destroying control: a candidate that has degraded HALFWAY toward the
# control's collapse is fusing distinct populations / eroding rare isolation and
# is disqualified. The engine's own over-correcting FLAG uses the full spread
# (FLOOR_FRACTION=1.0, i.e. "shredded to control level"); the gate uses the
# stricter midpoint so over-correction milder than total collapse is still
# caught (COND #3 / TM-3 / TM-15). Recorded in the audit.
FLOOR_FRACTION = 0.5

# Minimum collapse (baseline - control) on a discovery metric for the control to
# count as "fired". The shuffle control destroys cluster structure, so its
# distinctness/isolation drop far exceeds this; the bar is intentionally low (a
# fired anchor is unambiguous, and a NON-fired shuffle would itself be alarming).
# Named explicitly (not an inline literal) so every gate constant is auditable.
CONTROL_FIRING_MIN_COLLAPSE = 0.05

# Shuffle control method name (the REQUIRED falsifiability anchor, §D3).
SHUFFLE_CONTROL_METHOD = "neg_control_shuffle_label"

# Control methods are never recommended; they exist only as anchors.
CONTROL_METHODS = {SHUFFLE_CONTROL_METHOD, NEG_CONTROL_METHOD}


# ==========================================================================
# Floor pre-registration (reuses pre_register_over_correcting_threshold)
# ==========================================================================

def _midpoint_clip_note(baseline: float | None, neg_control: float | None) -> str:
    """Flag when T was clipped, so the floor is NOT the true baseline<->control midpoint.

    ``pre_register_over_correcting_threshold`` returns ``T = max(floor, baseline -
    neg_control)``. When the raw baseline->control spread is below
    ``OVER_CORRECTING_THRESHOLD_FLOOR`` the clip is ACTIVE and ``baseline -
    0.5*T`` is a conservative FIXED offset, not the literal midpoint to the
    control. Surfacing this keeps the audit honest (review condition a)."""
    if baseline is None or neg_control is None:
        return ""
    spread = float(baseline) - float(neg_control)
    if spread < OVER_CORRECTING_THRESHOLD_FLOOR:
        return (f" [NOTE: raw spread {spread:.4f} < floor "
                f"{OVER_CORRECTING_THRESHOLD_FLOOR}; T clipped -> floor is a "
                "conservative fixed offset, NOT the baseline<->control midpoint]")
    return ""


def pre_register_distinctness_floor(
    baseline_distinctness: float | None,
    neg_control_distinctness: float | None,
    *,
    floor_fraction: float = FLOOR_FRACTION,
) -> tuple[float, str]:
    """D1b floor = baseline - floor_fraction * spread. Pre-registered (§D3).

    ``spread`` is the baseline-minus-neg(shuffle)-control quantity from the
    engine's ``pre_register_over_correcting_threshold`` helper. The floor sits
    ``floor_fraction`` of the way down that spread (default 0.5 = midpoint): a
    candidate that has eroded HALFWAY toward the structure-destroying control is
    fusing distinct populations (over-correction) -> disqualified. Stricter than
    the engine's full-spread over-correcting flag so partial over-merging is
    still caught (TM-3 / TM-15).
    """
    T, prov = pre_register_over_correcting_threshold(
        baseline_distinctness, neg_control_distinctness)
    if baseline_distinctness is None:
        return 0.0, "no baseline distinctness; floor=0.0 (degenerate)"
    floor = max(0.0, float(baseline_distinctness) - floor_fraction * T)
    clip_note = _midpoint_clip_note(baseline_distinctness, neg_control_distinctness)
    return floor, (f"baseline_distinctness({baseline_distinctness:.4f}) - "
                   f"{floor_fraction}*T[{prov}] => floor={floor:.4f}{clip_note}")


def pre_register_isolation_floor(
    baseline_isolation: float | None,
    neg_control_isolation: float | None,
    *,
    floor_fraction: float = FLOOR_FRACTION,
) -> tuple[float, str]:
    """D1c floor = baseline - floor_fraction * spread. Pre-registered (COND #3).

    Same midpoint pre-registration as the distinctness floor, applied to the
    WORST-isolated per-batch population (D1c). A rare/novel type eroded halfway
    toward the control's collapse is disqualified EVEN IF D1a is highest
    (COND #3 / TM-15).
    """
    T, prov = pre_register_over_correcting_threshold(
        baseline_isolation, neg_control_isolation)
    if baseline_isolation is None:
        return 0.0, "no baseline isolation; floor=0.0 (degenerate)"
    floor = max(0.0, float(baseline_isolation) - floor_fraction * T)
    clip_note = _midpoint_clip_note(baseline_isolation, neg_control_isolation)
    return floor, (f"baseline_isolation({baseline_isolation:.4f}) - "
                   f"{floor_fraction}*T[{prov}] => floor={floor:.4f}{clip_note}")


# ==========================================================================
# Control firing check (§D3)
# ==========================================================================

def check_controls_fire(
    results: list[dm.DiscoveryMetrics],
    *,
    baseline_method: str = BASELINE_METHOD,
) -> dict[str, Any]:
    """Verify the REQUIRED shuffle anchor FIRES (collapses) on discovery metrics.

    The shuffle control destroys all structure -> its within-batch distinctness
    AND isolation must collapse well below baseline. The extreme-theta control
    (if present) is REPORTED verbatim but never gates. Returns the audit dict;
    ``at_least_one_control_fires`` is driven by the shuffle control ONLY.
    """
    by_method = {r.method: r for r in results}
    baseline = by_method.get(baseline_method)
    base_dist = baseline.within_batch_distinctness if baseline else None
    base_iso = baseline.within_batch_isolation if baseline else None

    out: dict[str, Any] = {
        "shuffle_control_present": SHUFFLE_CONTROL_METHOD in by_method,
        "extreme_theta_control_present": NEG_CONTROL_METHOD in by_method,
        "shuffle_fired": False,
        "extreme_theta_fired": None,  # logged-only
        "at_least_one_control_fires": False,
    }

    shuf = by_method.get(SHUFFLE_CONTROL_METHOD)
    if shuf is not None and base_dist is not None:
        # Fired = shuffle distinctness OR isolation collapses substantially
        # below baseline (structure destroyed). Use the engine's spread idea: a
        # collapse of more than the conservative floor.
        sd = shuf.within_batch_distinctness
        si = shuf.within_batch_isolation
        dist_collapse = (sd is not None and base_dist is not None
                         and (base_dist - sd) > CONTROL_FIRING_MIN_COLLAPSE)
        iso_collapse = (si is not None and base_iso is not None
                        and (base_iso - si) > CONTROL_FIRING_MIN_COLLAPSE)
        out["shuffle_distinctness"] = sd
        out["shuffle_isolation"] = si
        out["baseline_distinctness"] = base_dist
        out["baseline_isolation"] = base_iso
        out["shuffle_fired"] = bool(dist_collapse or iso_collapse)

    theta = by_method.get(NEG_CONTROL_METHOD)
    if theta is not None:
        td = theta.within_batch_distinctness
        ti = theta.within_batch_isolation
        out["extreme_theta_distinctness"] = td
        out["extreme_theta_isolation"] = ti
        # Logged verbatim; non-firing is acceptable (harmonypy robust).
        out["extreme_theta_fired"] = bool(
            (td is not None and base_dist is not None
             and (base_dist - td) > CONTROL_FIRING_MIN_COLLAPSE)
            or (ti is not None and base_iso is not None
                and (base_iso - ti) > CONTROL_FIRING_MIN_COLLAPSE))
        out["extreme_theta_note"] = (
            "REPORTED-NOT-GATING: harmonypy robust; theta=100 may not collapse "
            "structure on real data (known-weak control, §D3).")

    out["at_least_one_control_fires"] = bool(out["shuffle_fired"])
    return out


# ==========================================================================
# Recommendation (§D4)
# ==========================================================================

@dataclass
class Recommendation:
    chosen_method: str           # cfg.batch.method value
    reason: str
    rule_constants: dict[str, Any]
    gate_results: dict[str, dict[str, Any]]  # per-candidate gate pass/fail
    ranking: list[dict[str, Any]]
    verdict_tier: str = "relative_ranking_only"
    scope: str = "recommended FOR THIS DATASET ONLY (per-run, not a public claim)"
    notes: list[str] = field(default_factory=list)


# Map a discovery-method name to the cfg.batch.method value the runner consumes.
_CFG_METHOD = {
    BASELINE_METHOD: "none",
    HARMONY_METHOD: "harmony",
    SCVI_METHOD: "scvi",
}


def _cfg_method_for(method: str) -> str:
    if method.startswith(SCVI_METHOD):
        return "scvi"
    return _CFG_METHOD.get(method, method)


def recommend(
    results: list[dm.DiscoveryMetrics],
    *,
    baseline_method: str = BASELINE_METHOD,
    margin_mix: float = MARGIN_MIX,
) -> Recommendation:
    """Apply the §D4 DISQUALIFYING-GATE rule and return a Recommendation.

    Floors are gates (steps 1-3). D1a only ranks survivors (step 4). scVI seeds
    are collapsed into a band; the band is treated as a single candidate using
    its MEAN for gating and ranking, with band-aware tie annotation.
    """
    by_method = {r.method: r for r in results}
    baseline = by_method.get(baseline_method)
    if baseline is None:
        raise ValueError(
            f"baseline method {baseline_method!r} absent from results; cannot "
            "compute mixing floor (fail loud).")

    base_mix = baseline.batch_mixing_score or 0.0
    base_dist = baseline.within_batch_distinctness
    base_iso = baseline.within_batch_isolation

    # Pre-register floors from baseline vs the shuffle (neg) control.
    shuf = by_method.get(SHUFFLE_CONTROL_METHOD)
    neg_dist = shuf.within_batch_distinctness if shuf else None
    neg_iso = shuf.within_batch_isolation if shuf else None
    dist_floor, dist_floor_prov = pre_register_distinctness_floor(base_dist, neg_dist)
    iso_floor, iso_floor_prov = pre_register_isolation_floor(base_iso, neg_iso)

    rule_constants = {
        "MARGIN_MIX": margin_mix,
        "baseline_mixing": base_mix,
        "mixing_floor": base_mix + margin_mix,
        "distinctness_floor": dist_floor,
        "distinctness_floor_provenance": dist_floor_prov,
        "isolation_floor": iso_floor,
        "isolation_floor_provenance": iso_floor_prov,
        "baseline_distinctness": base_dist,
        "baseline_isolation": base_iso,
    }

    # Candidate set: everything that is NOT the baseline and NOT a control.
    # scVI seeds are collapsed to a single band candidate.
    scvi_results = [r for r in results if r.method.startswith(SCVI_METHOD)]
    non_scvi_candidates = [
        r for r in results
        if r.method != baseline_method
        and r.method not in CONTROL_METHODS
        and not r.method.startswith(SCVI_METHOD)
    ]

    gate_results: dict[str, dict[str, Any]] = {}
    survivors: dict[str, float] = {}  # method -> D1a (ranking score)

    def evaluate(method: str, mix: float | None, dist: float | None,
                 iso: float | None, d1a: float | None) -> None:
        passed = {}
        mix_ok = mix is not None and mix > base_mix + margin_mix
        dist_ok = dist is not None and dist >= dist_floor
        iso_ok = iso is not None and iso >= iso_floor
        passed["mixing_floor"] = bool(mix_ok)
        passed["distinctness_floor"] = bool(dist_ok)
        passed["isolation_floor"] = bool(iso_ok)
        passed["batch_mixing"] = mix
        passed["within_batch_distinctness"] = dist
        passed["within_batch_isolation"] = iso
        passed["good_merge_fraction"] = d1a
        survived = bool(mix_ok and dist_ok and iso_ok)
        passed["survived_gates"] = survived
        gate_results[method] = passed
        if survived:
            survivors[method] = float(d1a) if d1a is not None else 0.0

    for r in non_scvi_candidates:
        evaluate(r.method, r.batch_mixing_score, r.within_batch_distinctness,
                 r.within_batch_isolation, r.good_merge_fraction)

    # scVI band candidate (mean across seeds for each metric).
    scvi_band = None
    if scvi_results:
        mix_vals = [r.batch_mixing_score for r in scvi_results if r.batch_mixing_score is not None]
        dist_vals = [r.within_batch_distinctness for r in scvi_results if r.within_batch_distinctness is not None]
        iso_vals = [r.within_batch_isolation for r in scvi_results if r.within_batch_isolation is not None]
        d1a_vals = [r.good_merge_fraction for r in scvi_results if r.good_merge_fraction is not None]
        scvi_mix = float(np.mean(mix_vals)) if mix_vals else None
        scvi_dist = float(np.mean(dist_vals)) if dist_vals else None
        scvi_iso = float(np.mean(iso_vals)) if iso_vals else None
        scvi_d1a = float(np.mean(d1a_vals)) if d1a_vals else None
        scvi_band = summarize_scvi_band(d1a_vals, metric="good_merge_fraction") if d1a_vals else None
        evaluate(SCVI_METHOD, scvi_mix, scvi_dist, scvi_iso, scvi_d1a)

    # Step 4: rank survivors by D1a (bonus). Step 5: band-aware ties.
    ranking = band_aware_rank(survivors, scvi_band=scvi_band, scvi_label=SCVI_METHOD)

    notes: list[str] = []
    if not survivors:
        # Step 6: baseline fallback (loud).
        notes.append(
            "NO candidate survived the disqualifying gates (mixing/distinctness/"
            "isolation floors). Recommending BASELINE (no integration) — "
            "legitimate when there is no real batch effect, OR every candidate "
            "over-corrected. This is the discovery-safe choice.")
        return Recommendation(
            chosen_method="none",
            reason="baseline_fallback_no_survivor",
            rule_constants=rule_constants,
            gate_results=gate_results,
            ranking=ranking,
            notes=notes,
        )

    # Top survivor by D1a; band-aware tie note if a deterministic method ties scVI.
    top = ranking[0]
    chosen = top["method"]
    notes.append(
        f"Chosen '{chosen}' = highest good-merge correspondence (D1a) AMONG "
        f"candidates that PASSED all three disqualifying gates. D1a is a "
        f"ranking bonus, NOT a rescue: gates were applied first (§D4).")
    if top.get("tie_with_scvi_band"):
        notes.append(
            f"'{chosen}' falls inside the scVI seed band -> band-aware TIE with "
            "scVI; no forced winner across the band (seed variance).")
    return Recommendation(
        chosen_method=_cfg_method_for(chosen),
        reason="highest_d1a_among_gate_survivors",
        rule_constants=rule_constants,
        gate_results=gate_results,
        ranking=ranking,
        notes=notes,
    )


# ==========================================================================
# Output writers (§D5) — NEVER inside the factory tree
# ==========================================================================

_SCOREBOARD_COLUMNS = [
    "method",
    "batch_mixing_score", "ilisi", "kbet_acceptance",
    "good_merge_fraction", "within_batch_distinctness", "within_batch_isolation",
    "n_corresponded_pairs", "n_co_localized",
]


def write_discovery_scoreboard(
    results: list[dm.DiscoveryMetrics], path: Path) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(_SCOREBOARD_COLUMNS)
        for r in results:
            d1a = r.d1a_detail or {}
            w.writerow([
                r.method,
                _fmt(r.batch_mixing_score), _fmt(r.ilisi), _fmt(r.kbet_acceptance),
                _fmt(r.good_merge_fraction), _fmt(r.within_batch_distinctness),
                _fmt(r.within_batch_isolation),
                d1a.get("n_corresponded_pairs", ""),
                d1a.get("n_co_localized", ""),
            ])


def _fmt(v: float | None) -> str:
    return "" if v is None else f"{v:.6f}"


def build_recommendation_json(
    rec: Recommendation,
    ref: dm.PerBatchReference,
    resolved_versions: dict[str, str],
    *,
    batch_key: str,
    cache_key: str | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "schema_version": "discovery-integration-gate-v1",
        "computed_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "verdict_tier": rec.verdict_tier,
        "scope": rec.scope,
        "batch_key": batch_key,
        "cfg.batch.method": rec.chosen_method,
        "chosen_method": rec.chosen_method,
        "reason": rec.reason,
        "rule_constants": rec.rule_constants,
        "gate_results": rec.gate_results,
        "ranking": rec.ranking,
        "notes": rec.notes,
        "leiden_params": {
            "resolution": dm.LEIDEN_RESOLUTION,
            "random_state": dm.LEIDEN_RANDOM_STATE,
            "n_neighbors": dm.LEIDEN_N_NEIGHBORS,
            "flavor": dm.LEIDEN_FLAVOR,
        },
        "pinned_versions": resolved_versions,
        "correspondence_map": [list(p) for p in ref.correspondence],
        "per_batch_cluster_counts": {
            b: len(ref.cluster_keys[b]) for b in ref.batches},
        "reference_notes": ref.notes,
        "harmonypy_direct_caveat": (
            "If Harmony is chosen, evaluation used harmonypy.run_harmony DIRECT "
            "(production scanpy/rapids Harmony backend is broken in sc_gpu; "
            "Q15/Q19). The runner must route Harmony through the direct path "
            "until the production backend fix lands."),
        "cache_key": cache_key,
    }
    if extra:
        out.update(extra)
    return out


def write_recommendation_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n")


def write_audit_md(payload: dict[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Discovery Integration-Selection Gate — Audit")
    lines.append("")
    lines.append(f"- computed: `{payload['computed_iso']}`")
    lines.append(f"- verdict_tier: **{payload['verdict_tier']}**")
    lines.append(f"- scope: {payload['scope']}")
    lines.append(f"- batch_key: `{payload['batch_key']}`")
    lines.append(f"- **chosen cfg.batch.method: `{payload['cfg.batch.method']}`** "
                 f"(reason: {payload['reason']})")
    lines.append("")
    lines.append("## Rule constants (pre-registered)")
    lines.append("")
    for k, v in payload["rule_constants"].items():
        lines.append(f"- {k}: `{v}`")
    lines.append("")
    lines.append("## Per-candidate gate pass/fail")
    lines.append("")
    lines.append("| method | mixing | dist | iso | mix_floor | dist_floor | iso_floor | survived |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for m, g in payload["gate_results"].items():
        lines.append(
            f"| {m} | {_md(g.get('batch_mixing'))} | "
            f"{_md(g.get('within_batch_distinctness'))} | "
            f"{_md(g.get('within_batch_isolation'))} | "
            f"{g.get('mixing_floor')} | {g.get('distinctness_floor')} | "
            f"{g.get('isolation_floor')} | {g.get('survived_gates')} |")
    lines.append("")
    lines.append("## Ranking (survivors only; D1a bonus)")
    lines.append("")
    for r in payload["ranking"]:
        lines.append(f"- {r}")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    for n in payload["notes"]:
        lines.append(f"- {n}")
    lines.append("")
    lines.append(f"- harmonypy-direct caveat: {payload['harmonypy_direct_caveat']}")
    lines.append(f"- correspondence map: `{payload['correspondence_map']}`")
    lines.append(f"- leiden params: `{payload['leiden_params']}`")
    lines.append(f"- pinned versions: `{payload['pinned_versions']}`")
    if payload.get("controls"):
        lines.append("")
        lines.append("## Controls (shuffle = REQUIRED anchor; extreme-theta = logged-only)")
        lines.append("")
        for k, v in payload["controls"].items():
            lines.append(f"- {k}: `{v}`")
    if payload.get("lref_crosscheck"):
        lines.append("")
        lines.append("## LOGGED SECONDARY: global L_ref cross-check (never gate-driving)")
        lines.append("")
        for m, v in payload["lref_crosscheck"].items():
            lines.append(f"- {m}: `{v}`")
    lines.append("")
    path.write_text("\n".join(lines))


def _md(v: float | None) -> str:
    return "—" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v))


# ==========================================================================
# Project-root resolution (CLAUDE.md agent rule: never write in factory tree)
# ==========================================================================

def resolve_out_dir(project_root: Path, subdir: str = "discovery_gate") -> Path:
    """Resolve + create the output dir UNDER project_root (never the factory)."""
    factory_root = Path(__file__).resolve().parents[3]
    pr = Path(project_root).resolve()
    if factory_root in pr.parents or pr == factory_root:
        raise ValueError(
            f"--project-root {pr} is inside the factory tree {factory_root}; "
            "outputs must never be written inside the factory (CLAUDE.md). "
            "Pass a project directory under /home/zerlinshen/projects/.")
    out = pr / subdir
    out.mkdir(parents=True, exist_ok=True)
    return out


# ==========================================================================
# CLI
# ==========================================================================

def _compute_data_hash(X: np.ndarray) -> str:
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(X, dtype=np.float32).tobytes())
    return h.hexdigest()[:16]


def run_gate(
    embeddings: dict[str, np.ndarray],
    batch: np.ndarray,
    *,
    baseline_method: str = BASELINE_METHOD,
    global_lref: np.ndarray | None = None,
    margin_mix: float = MARGIN_MIX,
) -> tuple[Recommendation, list[dm.DiscoveryMetrics], dm.PerBatchReference, dict[str, Any]]:
    """Score embeddings, check controls, produce the recommendation.

    ``embeddings`` MUST contain ``baseline_method`` (the shared pre-integration
    X_pca, used both as the baseline candidate AND to build the per-batch
    reference). Fail loud if absent.
    """
    if baseline_method not in embeddings:
        raise ValueError(
            f"baseline embedding {baseline_method!r} absent; required to build "
            "the per-batch reference (shared X_pca) and the mixing floor.")
    X_pca_baseline = embeddings[baseline_method]
    results, ref = dm.score_all_embeddings_discovery(
        embeddings, X_pca_baseline, batch, global_lref=global_lref)
    controls = check_controls_fire(results, baseline_method=baseline_method)
    rec = recommend(results, baseline_method=baseline_method, margin_mix=margin_mix)
    return rec, results, ref, controls


def _load_embeddings_from_h5ad(
    h5ad_path: Path, batch_key: str
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray | None]:
    """Load obsm embeddings + batch labels from a persisted h5ad.

    obsm keys are normalized: ``X_baseline`` -> ``baseline``, ``X_harmony`` ->
    ``harmony``, ``X_scvi_seed{n}`` -> ``scvi_seed{n}``, the shuffle/neg keys map
    to their method names. Fail loud if the batch_key column is missing.
    """
    import anndata as ad

    a = ad.read_h5ad(h5ad_path)
    if batch_key not in a.obs.columns:
        raise ValueError(
            f"batch_key {batch_key!r} absent from obs (have "
            f"{list(a.obs.columns)}); fail loud.")
    batch = a.obs[batch_key].astype(str).to_numpy()

    keymap = {
        "X_baseline": BASELINE_METHOD,
        "X_pca": BASELINE_METHOD,
        "X_harmony": HARMONY_METHOD,
        "X_pca_harmony": HARMONY_METHOD,
        "X_neg_control_shuffle_label": SHUFFLE_CONTROL_METHOD,
        "X_shuffle": SHUFFLE_CONTROL_METHOD,
        "X_neg_control": NEG_CONTROL_METHOD,
        "X_neg_control_harmony_extreme_theta": NEG_CONTROL_METHOD,
    }
    embeddings: dict[str, np.ndarray] = {}
    for k in a.obsm.keys():
        if k in keymap:
            embeddings[keymap[k]] = np.asarray(a.obsm[k], dtype=np.float64)
        elif k.startswith("X_scvi_seed"):
            embeddings[k.replace("X_", "")] = np.asarray(a.obsm[k], dtype=np.float64)
    if BASELINE_METHOD not in embeddings:
        raise ValueError(
            f"no baseline embedding found in obsm (have {list(a.obsm.keys())}); "
            "fail loud (the gate requires the shared pre-integration X_pca).")

    # Optional global L_ref for the logged-secondary cross-check.
    global_lref = None
    for col in ("seurat_clusters", "leiden"):
        if col in a.obs.columns:
            global_lref = a.obs[col].astype(str).to_numpy()
            break
    return embeddings, batch, global_lref


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True,
                        help="Project dir (outputs land under <root>/discovery_gate/).")
    parser.add_argument("--batch-key", type=str, required=True,
                        help="obs column to use as the batch key (fail loud if absent).")
    parser.add_argument("--embeddings-h5ad", type=Path, default=None,
                        help="Persisted h5ad with obsm per method (re-score, no scVI retrain).")
    parser.add_argument("--margin-mix", type=float, default=MARGIN_MIX)
    parser.add_argument("--out-subdir", type=str, default="discovery_gate")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO)
    # GATE G1 + COND #4 version assert first — fail loud before compute.
    resolved = dm.assert_pinned_discovery_engine()

    out_dir = resolve_out_dir(args.project_root, subdir=args.out_subdir)

    if args.embeddings_h5ad is None:
        parser.error(
            "--embeddings-h5ad is required (re-score persisted embeddings; "
            "literal per-run scVI retrain is gated through the cache layer).")
        return 2

    embeddings, batch, global_lref = _load_embeddings_from_h5ad(
        args.embeddings_h5ad, args.batch_key)
    data_hash = _compute_data_hash(embeddings[BASELINE_METHOD])
    cache_key = f"{data_hash}:{args.batch_key}"

    rec, results, ref, controls = run_gate(
        embeddings, batch, global_lref=global_lref, margin_mix=args.margin_mix)

    payload = build_recommendation_json(
        rec, ref, resolved, batch_key=args.batch_key, cache_key=cache_key,
        extra={
            "controls": controls,
            "lref_crosscheck": {r.method: r.lref_crosscheck for r in results
                                if r.lref_crosscheck},
        })
    write_recommendation_json(payload, out_dir / "integration_recommendation.json")
    write_discovery_scoreboard(results, out_dir / "integration_scoreboard.csv")
    write_audit_md(payload, out_dir / "integration_audit.md")
    (out_dir / "integration_audit.json").write_text(json.dumps(payload, indent=2) + "\n")

    print(f"[discovery-gate] chosen cfg.batch.method = {rec.chosen_method}")
    print(f"[discovery-gate] reason = {rec.reason}")
    print(f"[discovery-gate] shuffle control fired = {controls['shuffle_fired']}")
    print(f"[discovery-gate] outputs -> {out_dir}")
    for m, g in rec.gate_results.items():
        print(f"  {m}: mixing_floor={g['mixing_floor']} "
              f"distinctness_floor={g['distinctness_floor']} "
              f"isolation_floor={g['isolation_floor']} survived={g['survived_gates']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
