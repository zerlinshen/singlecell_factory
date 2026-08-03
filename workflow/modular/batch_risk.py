"""Plan-time multi-batch detection and clustering claim-status resolution.

Deliberately dependency-light for the same reason as module_catalog.py: every
surface that can start an analysis — the canonical CLI, ``scfactory run``, and
``scfactory run --dry-run`` — must be able to detect and announce batch
structure without importing scanpy/anndata. Nothing beyond the standard library
is imported at module scope (``module_catalog`` imports the strategy vocabulary
from here, so this file must stay importable on a bare help path);
``numpy``/``h5py``/``zarr`` are imported lazily inside the reader.

Why this exists
---------------
The default optional-module set (``module_catalog.DEFAULT_OPTIONAL_MODULES``)
contains no integration step, so a default run on multi-batch input produces
clusters that track the batch rather than the biology — and it produced them
silently. On the LuCA LUSC cohort (92,430 cells, 87 samples, 9 datasets) the
unintegrated default scored ARI 0.216 / matched accuracy 0.381 / kBET 0.117 and
split the 24 published classes into 47 clusters; the same input with Harmony
scored ARI 0.471 / accuracy 0.581 / kBET 0.554 with 23 clusters
(governance/realrun_gt_concordance_lusc_2026-07-28.md).

The fix is NOT to force Harmony on by default — silently changing the analysis
is the same class of error in the other direction. Instead this module detects
the structure before any compute happens, warns loudly and machine-readably,
and downgrades the affected claim to ``exploratory`` unless the operator makes
an affirmative declaration. Silence never clears the warning.

``clustering.py::_record_batch_confounding_risk`` records the same condition at
*runtime*, after clustering has already run, and shares the candidate-column
tuple defined here so the two readings cannot drift apart.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

logger = logging.getLogger(__name__)


# obs columns that, with more than one level, indicate the run spans multiple
# technical units and is therefore exposed to batch effects. This tuple is the
# single source of truth: clustering.py imports it for its runtime reading.
BATCH_CANDIDATE_COLUMNS: tuple[str, ...] = (
    "batch", "sample", "sample_id", "donor_id", "donor",
    "dataset", "study", "patient", "patient_id", "platform", "assay",
    # Seurat's standard per-object sample column, and this suite's own R
    # default: r_multiomics_factory/R/pipeline_steps.R:18-19 resolves it as THE
    # batch column when present, and plotting_factory/r/general/qc_plots.R:62
    # uses it as the default group.by. The Python detector was blind to it, so
    # any Seurat-derived object read as single-batch.
    "orig.ident",
)

# Case-insensitive lookup: `Sample`, `BATCH` and `Orig.Ident` are the same axis
# as their lowercase spellings, and missing them used to be reported as an
# affirmative single-batch answer rather than as "not found".
_CANDIDATE_ORDER: dict[str, int] = {
    name.lower(): index for index, name in enumerate(BATCH_CANDIDATE_COLUMNS)
}

# Modules that actually correct batch structure. ``integration_select`` only
# SCORES candidates and sets cfg.batch.method — on its own it corrects nothing,
# so it does not satisfy an "integrate" declaration.
INTEGRATION_MODULES: tuple[str, ...] = ("batch_correction",)

# Modules whose primary output is derived from the latent structure, or from the
# leiden labels computed on it. When the embedding is batch-driven, so are these,
# so they carry the downgraded claim status. Deliberately narrow: modules that
# consume an explicit, separately-declared design (pseudobulk_de, composition
# with a condition column) carry their own inference contracts.
BATCH_SENSITIVE_MODULES: tuple[str, ...] = (
    "clustering",
    "differential_expression",
    "annotation",
    "context_aware_annotation",
    "composition",
)

# Declared batch strategies. ``auto`` is the absence of a declaration, not a
# declaration of "no batches" — it is the only value that can be reached by
# silence, and it is the only one that leaves the claim downgraded.
BATCH_STRATEGY_AUTO = "auto"
BATCH_STRATEGY_SINGLE_BATCH = "single-batch"
BATCH_STRATEGY_INTEGRATE = "integrate"
BATCH_STRATEGY_ACCEPT_UNCORRECTED = "accept-uncorrected"
BATCH_STRATEGY_CHOICES: tuple[str, ...] = (
    BATCH_STRATEGY_AUTO,
    BATCH_STRATEGY_SINGLE_BATCH,
    BATCH_STRATEGY_INTEGRATE,
    BATCH_STRATEGY_ACCEPT_UNCORRECTED,
)

# Claim vocabulary. ``exploratory`` matches velocity_render_contract.CLAIM_CLASS;
# the long-form reason strings match the ``exploratory_nonclaimable_*`` idiom
# already used by pseudobulk_de / composition / cell_communication.
CLAIM_STATUS_CLAIMABLE = "claimable"
CLAIM_STATUS_EXPLORATORY = "exploratory"
CLAIM_STATUS_UNDETERMINED = "undetermined"

DETECTION_MULTI_BATCH = "multi_batch"
DETECTION_SINGLE_BATCH = "single_batch"
DETECTION_UNAVAILABLE = "unavailable"

# Runtime (manifest-time) reading of the object the run actually produced. This
# is what makes the claim trustworthy: plan-time detection is UNAVAILABLE for
# any run that starts from a raw sample root, so a declaration made there is
# accepted on the operator's word alone until this reading checks it.
OBSERVATION_MULTI_BATCH = "multi_batch"
OBSERVATION_SINGLE_BATCH = "single_batch"
OBSERVATION_UNOBSERVED = "unobserved"

# Read obs columns in bounded blocks so a 10M-cell cohort cannot materialise a
# whole column at once, and stop tracking distinct levels past a sane ceiling
# (the exact count stops mattering long before this).
_READ_BLOCK = 1_000_000
_MAX_TRACKED_LEVELS = 10_000

# Per-column read outcomes. "absent" and "unreadable" are distinct because only
# the second means the reader positively KNOWS it failed.
_COLUMN_READ = "read"
_COLUMN_ABSENT = "absent"
_COLUMN_UNREADABLE = "unreadable"


class BatchStrategyConflict(ValueError):
    """A declared batch strategy contradicts the data or the planned modules."""


# ---------------------------------------------------------------------------
# obs-only reading
# ---------------------------------------------------------------------------


def resolve_obs_source(
    *,
    input_h5ad: Optional[Path] = None,
    sample_root: Optional[Path] = None,
) -> Optional[Path]:
    """Return the prepared-input path the run will actually ingest, if any.

    Mirrors ``modules/cellranger.py``: an explicit ``--input-h5ad`` wins,
    otherwise ``<sample-root>/prepared_input.h5ad`` then ``prepared_input.zarr``.
    A Cell Ranger sample-root with neither has no obs to read before the run, so
    the caller gets ``None`` (detection status ``unavailable``) rather than a
    guess.
    """
    if input_h5ad is not None:
        candidate = Path(input_h5ad)
        return candidate if candidate.is_file() else None
    if sample_root is None:
        return None
    root = Path(sample_root)
    prepared_h5ad = root / "prepared_input.h5ad"
    if prepared_h5ad.is_file():
        return prepared_h5ad
    prepared_zarr = root / "prepared_input.zarr"
    if prepared_zarr.is_dir():
        return prepared_zarr
    return None


def _open_obs_group(path: Path):
    """Open ``<path>/obs`` for h5ad or zarr. Returns (group, closer) or None.

    h5py Groups/Datasets and zarr Groups/Arrays expose the same subset of the
    mapping + array protocol used below, so one traversal serves both.
    """
    if path.is_dir():
        try:
            import zarr  # type: ignore
        except ImportError:
            logger.warning(
                "batch detection: %s is a zarr store but zarr is not importable; "
                "batch structure cannot be read before the run.",
                path,
            )
            return None
        root = zarr.open(str(path), mode="r")
        if "obs" not in root:
            return None
        return root["obs"], (lambda: None)

    import h5py  # type: ignore

    handle = h5py.File(str(path), "r")
    if "obs" not in handle:
        handle.close()
        return None
    return handle["obs"], handle.close


def _is_group(node: Any) -> bool:
    """True for a container node (h5py.Group / zarr.Group), False for an array."""
    return not hasattr(node, "dtype")


def _iter_blocks(dataset: Any) -> Iterable[Any]:
    import numpy as np

    shape = getattr(dataset, "shape", None)
    if not shape:
        return
    total = int(shape[0])
    for start in range(0, total, _READ_BLOCK):
        yield np.asarray(dataset[start:start + _READ_BLOCK])


def _count_levels(dataset: Any, *, drop_negative: bool) -> int:
    """Count distinct non-null values, matching ``nunique(dropna=True)``."""
    import numpy as np

    levels: set[Any] = set()
    for block in _iter_blocks(dataset):
        for value in np.unique(block):
            item = value.item() if hasattr(value, "item") else value
            if drop_negative and isinstance(item, int) and item < 0:
                continue
            if item is None or item == b"" or item == "":
                continue
            if isinstance(item, float) and np.isnan(item):
                continue
            levels.add(item)
            if len(levels) >= _MAX_TRACKED_LEVELS:
                return _MAX_TRACKED_LEVELS
    return len(levels)


def _candidate_rank(name: Any) -> int:
    """Tie-break position of ``name`` in the canonical candidate order."""
    return _CANDIDATE_ORDER.get(str(name).lower(), len(_CANDIDATE_ORDER))


def _resolve_candidate_columns(
    available: Iterable[Any], preferred_key: str = ""
) -> list[str]:
    """Candidate columns actually present, matched case-insensitively.

    Driven by what the object HAS rather than by the canonical spellings, so
    ``Sample`` / ``BATCH`` / ``Orig.Ident`` resolve to the same axes as their
    lowercase forms instead of silently going unfound.
    """
    names = [str(name) for name in available]
    resolved = [name for name in names if name.lower() in _CANDIDATE_ORDER]
    if preferred_key and preferred_key in names and preferred_key not in resolved:
        resolved.insert(0, preferred_key)
    return sorted(resolved, key=_candidate_rank)


def _column_levels(obs_group: Any, column: str) -> tuple[str, Optional[int]]:
    """Distinct level count for one obs column, as an explicit tri-state.

    Returns ``("read", n)``, ``("absent", None)`` or ``("unreadable", None)``.
    These are deliberately NOT collapsed: a single ``None`` for both "the column
    is not there" and "the read failed" is what let a known read failure be
    reported as an affirmative single-batch answer.
    """
    if column not in obs_group:
        return _COLUMN_ABSENT, None
    node = obs_group[column]
    try:
        if _is_group(node):
            # AnnData categorical: {categories, codes}. Count used codes rather
            # than declared categories — a subset cohort routinely keeps unused
            # categories, and those are not batches present in this run.
            if "codes" not in node:
                logger.warning(
                    "batch detection: obs[%r] is a group without 'codes'; "
                    "the column is present but cannot be read.",
                    column,
                )
                return _COLUMN_UNREADABLE, None
            return _COLUMN_READ, _count_levels(node["codes"], drop_negative=True)
        return _COLUMN_READ, _count_levels(node, drop_negative=False)
    except Exception as exc:
        logger.warning("batch detection: obs[%r] unreadable: %s", column, exc)
        return _COLUMN_UNREADABLE, None


def detect_batch_structure(
    obs_source: Optional[Path],
    *,
    preferred_key: str = "",
) -> dict:
    """Detect multi-batch structure from obs alone. Never raises.

    Reads only the candidate obs columns (categorical codes where available), so
    cost is bounded by the number of cells times a handful of narrow columns —
    the count matrix is never touched, and the lazy/backed ingest path is
    unaffected because nothing is loaded into the run's AnnData here.
    """
    unavailable = {
        "detection_status": DETECTION_UNAVAILABLE,
        "detection_source": str(obs_source) if obs_source is not None else "",
        "detected_batches": 0,
        "key": "",
        "candidate_batch_columns": {},
        "candidate_columns_present": [],
        "candidate_columns_read": [],
        "unreadable_columns": [],
    }
    if obs_source is None:
        unavailable["detection_reason"] = "no prepared obs available before the run"
        return unavailable

    try:
        opened = _open_obs_group(Path(obs_source))
    except Exception as exc:
        unavailable["detection_reason"] = f"obs could not be opened: {exc}"
        logger.warning("batch detection: could not open %s: %s", obs_source, exc)
        return unavailable
    if opened is None:
        unavailable["detection_reason"] = "input carries no obs group"
        return unavailable

    obs_group, closer = opened
    levels: dict[str, int] = {}
    unreadable: list[str] = []
    present: list[str] = []
    try:
        present = _resolve_candidate_columns(obs_group.keys(), preferred_key)
        for column in present:
            state, count = _column_levels(obs_group, column)
            if state == _COLUMN_READ and count is not None:
                levels[column] = count
            elif state == _COLUMN_UNREADABLE:
                unreadable.append(column)
    finally:
        try:
            closer()
        except Exception:  # pragma: no cover - close failures are not run errors
            pass

    return _classify_levels(
        levels,
        unreadable=unreadable,
        present=present,
        preferred_key=preferred_key,
        source=str(obs_source),
    )


def _classify_levels(
    levels: dict[str, int],
    *,
    unreadable: list[str],
    present: list[str],
    preferred_key: str,
    source: str,
) -> dict:
    """Turn per-column level counts into a detection verdict.

    The ordering of the three outcomes is the whole point of this function:

    * any readable column with >1 level -> multi_batch. An unreadable column
      cannot undo that, so this decision is safe to make first.
    * nothing readable, or ANY candidate that failed to read -> unavailable.
      Not single_batch: the unread column could have been the multi-level axis,
      and a detector that knows it failed must not answer affirmatively.
    * at least one column read, all single-level, none unreadable -> the only
      genuine single_batch.
    """
    multi = {name: count for name, count in levels.items() if count > 1}
    read = sorted(levels, key=_candidate_rank)
    base = {
        "detection_source": source,
        "candidate_columns_present": list(present),
        "candidate_columns_read": read,
        "unreadable_columns": sorted(unreadable),
    }

    if multi:
        if preferred_key and preferred_key in multi:
            key = preferred_key
        else:
            key = min(multi, key=lambda name: (-multi[name], _candidate_rank(name)))
        return {
            **base,
            "detection_status": DETECTION_MULTI_BATCH,
            "detected_batches": int(multi[key]),
            "key": key,
            "candidate_batch_columns": dict(sorted(multi.items())),
        }

    if unreadable or not levels:
        if unreadable:
            reason = (
                "batch-like obs column(s) present but unreadable: "
                + ", ".join(sorted(unreadable))
            )
        else:
            reason = (
                "no batch-like obs column found (searched "
                f"{', '.join(BATCH_CANDIDATE_COLUMNS)}, case-insensitively)"
            )
        logger.warning("batch detection: %s; batch structure is UNDETERMINED.", reason)
        return {
            **base,
            "detection_status": DETECTION_UNAVAILABLE,
            "detected_batches": 0,
            "key": "",
            "candidate_batch_columns": {},
            "detection_reason": reason,
        }

    return {
        **base,
        "detection_status": DETECTION_SINGLE_BATCH,
        "detected_batches": 1,
        "key": read[0] if read else "",
        "candidate_batch_columns": {},
    }


# ---------------------------------------------------------------------------
# strategy resolution
# ---------------------------------------------------------------------------


def resolve_batch_risk(
    *,
    detection: dict,
    declared_strategy: str,
    planned_modules: Sequence[str],
) -> dict:
    """Combine detection, declaration, and plan into the manifest envelope.

    Raises :class:`BatchStrategyConflict` when the declaration contradicts what
    is on disk or what the plan can deliver. A false declaration is a scientific
    error, not a preference, so it fails before any compute happens.
    """
    if declared_strategy not in BATCH_STRATEGY_CHOICES:
        raise BatchStrategyConflict(
            f"unknown batch strategy {declared_strategy!r}; "
            f"expected one of: {', '.join(BATCH_STRATEGY_CHOICES)}"
        )

    planned = list(planned_modules)
    integration_planned = [name for name in INTEGRATION_MODULES if name in planned]
    status = detection.get("detection_status", DETECTION_UNAVAILABLE)
    detected = int(detection.get("detected_batches", 0))
    key = str(detection.get("key", ""))

    if declared_strategy == BATCH_STRATEGY_SINGLE_BATCH and status == DETECTION_MULTI_BATCH:
        raise BatchStrategyConflict(
            f"--batch-strategy single-batch was declared, but obs[{key!r}] carries "
            f"{detected} levels "
            f"({_format_candidates(detection.get('candidate_batch_columns', {}))}). "
            "Declare --batch-strategy integrate with batch_correction in "
            "--optional-modules, or --batch-strategy accept-uncorrected to record "
            "an explicitly uncorrected exploratory run."
        )
    if declared_strategy == BATCH_STRATEGY_INTEGRATE and not integration_planned:
        raise BatchStrategyConflict(
            "--batch-strategy integrate was declared, but the plan contains no "
            f"integration module ({', '.join(INTEGRATION_MODULES)}). Add "
            "batch_correction to --optional-modules or drop the declaration."
        )

    if declared_strategy == BATCH_STRATEGY_AUTO and integration_planned:
        # Planning an integration module IS an affirmative strategy: it is
        # explicit, it is in the manifest, and it changes the analysis.
        strategy = BATCH_STRATEGY_INTEGRATE
        strategy_declared = True
        strategy_source = "planned_modules"
    elif declared_strategy == BATCH_STRATEGY_AUTO:
        strategy = BATCH_STRATEGY_AUTO
        strategy_declared = False
        strategy_source = "undeclared"
    else:
        strategy = declared_strategy
        strategy_declared = True
        strategy_source = "operator_declaration"

    claim_status, claim_reason = _resolve_claim(
        status, strategy, read_failed=bool(detection.get("unreadable_columns"))
    )
    sensitive_planned = [name for name in BATCH_SENSITIVE_MODULES if name in planned]

    return {
        "detected_batches": detected,
        "key": key,
        "candidate_batch_columns": dict(detection.get("candidate_batch_columns", {})),
        "detection_status": status,
        "detection_source": str(detection.get("detection_source", "")),
        "detection_reason": str(detection.get("detection_reason", "")),
        "strategy": strategy,
        "declared_strategy": declared_strategy,
        "strategy_declared": strategy_declared,
        "strategy_source": strategy_source,
        "integration_modules_planned": integration_planned,
        # Retained so the manifest-time reconciliation can recompute the
        # downgrade set without being handed the plan again.
        "batch_sensitive_modules_planned": sensitive_planned,
        "clustering_claim_status": claim_status,
        "clustering_claim_reason": claim_reason,
        "downgraded_modules": (
            sensitive_planned if claim_status != CLAIM_STATUS_CLAIMABLE else []
        ),
        "claim_basis": "plan_time_detection",
        "observation_status": OBSERVATION_UNOBSERVED,
    }


def _resolve_claim(
    detection_status: str, strategy: str, *, read_failed: bool = False
) -> tuple[str, str]:
    if detection_status == DETECTION_SINGLE_BATCH:
        return CLAIM_STATUS_CLAIMABLE, "single_batch_input"
    if detection_status == DETECTION_UNAVAILABLE:
        if strategy == BATCH_STRATEGY_INTEGRATE:
            return CLAIM_STATUS_CLAIMABLE, "integration_planned"
        if strategy == BATCH_STRATEGY_SINGLE_BATCH and not read_failed:
            # An affirmative declaration is the operator's to make when there is
            # simply no batch annotation to check it against. It is NOT theirs to
            # make when the detector positively failed to read a column that was
            # there — that column could be the one contradicting them.
            return CLAIM_STATUS_CLAIMABLE, "single_batch_declared"
        if strategy == BATCH_STRATEGY_ACCEPT_UNCORRECTED:
            return (
                CLAIM_STATUS_EXPLORATORY,
                "exploratory_nonclaimable_uncorrected_declared",
            )
        return (
            CLAIM_STATUS_UNDETERMINED,
            "batch_structure_unreadable_before_run",
        )
    # detection_status == DETECTION_MULTI_BATCH. A single-batch declaration
    # cannot reach here: resolve_batch_risk raises on that contradiction.
    if strategy == BATCH_STRATEGY_INTEGRATE:
        return CLAIM_STATUS_CLAIMABLE, "multi_batch_integrated"
    if strategy == BATCH_STRATEGY_ACCEPT_UNCORRECTED:
        return (
            CLAIM_STATUS_EXPLORATORY,
            "exploratory_nonclaimable_unintegrated_multi_batch_declared",
        )
    return (
        CLAIM_STATUS_EXPLORATORY,
        "exploratory_nonclaimable_unintegrated_multi_batch",
    )


def observe_batch_structure(obs: Any, *, preferred_key: str = "") -> dict:
    """Read batch structure from a loaded ``adata.obs``. Never raises.

    The in-memory counterpart of :func:`detect_batch_structure`, using the same
    candidate columns and the same ``nunique(dropna=True)`` semantics so the
    plan-time and runtime readings are directly comparable. Shared with
    ``clustering.py::_record_batch_confounding_risk``.
    """
    unobserved = {
        "observation_status": OBSERVATION_UNOBSERVED,
        "observed_batches": 0,
        "observed_key": "",
        "observed_candidate_batch_columns": {},
        "observed_columns_present": [],
        "observed_columns_read": [],
        "observed_unreadable_columns": [],
    }
    if obs is None:
        return unobserved

    try:
        present = _resolve_candidate_columns(obs.columns, preferred_key)
    except Exception as exc:  # pragma: no cover - never block a manifest write
        logger.warning("batch observation failed: %s", exc)
        return {**unobserved, "observation_error": repr(exc)}

    levels: dict[str, int] = {}
    unreadable: list[str] = []
    for column in present:
        try:
            levels[column] = int(obs[column].nunique(dropna=True))
        except Exception as exc:
            logger.warning("batch observation: obs[%r] unreadable: %s", column, exc)
            unreadable.append(column)

    verdict = _classify_levels(
        levels,
        unreadable=unreadable,
        present=present,
        preferred_key=preferred_key,
        source="runtime_obs",
    )
    status = {
        DETECTION_MULTI_BATCH: OBSERVATION_MULTI_BATCH,
        DETECTION_SINGLE_BATCH: OBSERVATION_SINGLE_BATCH,
        # A runtime reading that found nothing readable is NOT an observed
        # single batch; it is an absence of observation, and reconcile keeps the
        # plan-time claim rather than upgrading on it.
        DETECTION_UNAVAILABLE: OBSERVATION_UNOBSERVED,
    }[verdict["detection_status"]]
    return {
        "observation_status": status,
        "observed_batches": verdict["detected_batches"],
        "observed_key": verdict["key"],
        "observed_candidate_batch_columns": verdict["candidate_batch_columns"],
        "observed_columns_present": verdict["candidate_columns_present"],
        "observed_columns_read": verdict["candidate_columns_read"],
        "observed_unreadable_columns": verdict["unreadable_columns"],
    }


def reconcile_batch_risk(
    risk: dict,
    *,
    observation: dict,
    integration_completed: bool,
) -> dict:
    """Re-resolve the authoritative claim against the object the run produced.

    There is exactly ONE authoritative claim field,
    ``clustering_claim_status``. Plan-time detection can be UNAVAILABLE (every
    run starting from a raw sample root), in which case a declared strategy is
    granted on the operator's word; this pass checks that word against the data
    and against whether integration actually completed, and downgrades when it
    does not hold. The plan-time reading is preserved under ``plan_time_claim``
    rather than overwritten — a disagreement between the two timepoints is
    itself evidence, the same contract ``_restate_raw_axis_at_manifest_time``
    uses for the ``.raw`` gene axis.

    Deliberately never upgrades a claim to ``claimable`` on the strength of a
    declaration alone, and never raises: the compute has already happened, so
    destroying the artifacts would cost more than recording the truth about
    them.
    """
    reconciled = dict(risk)
    reconciled.update(
        {key: value for key, value in observation.items() if key != "observation_status"}
    )
    status = observation.get("observation_status", OBSERVATION_UNOBSERVED)
    reconciled["observation_status"] = status
    reconciled["integration_completed"] = bool(integration_completed)
    reconciled["plan_time_claim"] = {
        "clustering_claim_status": risk.get("clustering_claim_status"),
        "clustering_claim_reason": risk.get("clustering_claim_reason"),
    }

    if status == OBSERVATION_UNOBSERVED:
        # Nothing was loaded (crash before ingest, or a metadata-only run).
        # The plan-time claim stands, and says so.
        reconciled["claim_basis"] = "plan_time_detection"
        return reconciled

    reconciled["claim_basis"] = "runtime_observation"
    if status == OBSERVATION_SINGLE_BATCH:
        claim, reason = CLAIM_STATUS_CLAIMABLE, "single_batch_input_observed"
    elif integration_completed:
        claim, reason = CLAIM_STATUS_CLAIMABLE, "multi_batch_integrated"
    else:
        claim = CLAIM_STATUS_EXPLORATORY
        reason = _falsified_reason(risk)

    reconciled["clustering_claim_status"] = claim
    reconciled["clustering_claim_reason"] = reason
    reconciled["downgraded_modules"] = (
        list(risk.get("batch_sensitive_modules_planned", []))
        if claim != CLAIM_STATUS_CLAIMABLE
        else []
    )
    return reconciled


def _falsified_reason(risk: dict) -> str:
    """Name WHY an observed multi-batch run is not claimable, precisely."""
    plan_reason = str(risk.get("clustering_claim_reason", ""))
    strategy = risk.get("strategy")
    if strategy == BATCH_STRATEGY_ACCEPT_UNCORRECTED:
        return "exploratory_nonclaimable_unintegrated_multi_batch_declared"
    if plan_reason == "single_batch_declared":
        return "exploratory_nonclaimable_single_batch_declaration_falsified"
    if strategy == BATCH_STRATEGY_INTEGRATE:
        return "exploratory_nonclaimable_integration_planned_but_not_completed"
    return "exploratory_nonclaimable_unintegrated_multi_batch"


def reconciliation_warning(risk: dict) -> str:
    """Warning for a claim the runtime reading made WORSE. "" otherwise."""
    before = (risk.get("plan_time_claim") or {}).get("clustering_claim_status")
    after = risk.get("clustering_claim_status")
    if before != CLAIM_STATUS_CLAIMABLE or after == CLAIM_STATUS_CLAIMABLE:
        return ""
    reason = risk.get("clustering_claim_reason", "")
    if reason == "exploratory_nonclaimable_single_batch_declaration_falsified":
        head = (
            "BATCH_RISK_FALSIFIED: this run was launched with "
            "--batch-strategy single-batch, which plan-time detection could not "
            "check (no prepared obs before the run). The object the run produced "
            f"carries {risk.get('observed_batches')} levels in "
            f"obs[{risk.get('observed_key')!r}] "
            f"({_format_candidates(risk.get('observed_candidate_batch_columns', {}))})."
        )
    elif reason == (
        "exploratory_nonclaimable_integration_planned_but_not_completed"
    ):
        head = (
            "BATCH_RISK_FALSIFIED: integration was declared and planned, but "
            "batch_correction did not complete on the corrected representation, "
            f"and the object carries {risk.get('observed_batches')} batches in "
            f"obs[{risk.get('observed_key')!r}]."
        )
    else:
        head = (
            "BATCH_RISK_FALSIFIED: the run produced an object carrying "
            f"{risk.get('observed_batches')} batches in "
            f"obs[{risk.get('observed_key')!r}] with no completed integration."
        )
    return head + (
        " The clustering claim recorded at plan time was 'claimable'; it is "
        "DOWNGRADED to 'exploratory' "
        f"({reason}). Output of "
        f"{', '.join(risk.get('downgraded_modules') or []) or 'the affected modules'}"
        " must not be reported as integrated. The run's artifacts are kept: the "
        "compute already happened, and the honest record of it is more useful "
        "than deleting it."
    )


def plan_batch_risk(cfg: Any, planned_modules: Sequence[str]) -> dict:
    """Resolve the batch-risk envelope for a fully-built ``PipelineConfig``."""
    cellranger = getattr(cfg, "cellranger", None)
    source = resolve_obs_source(
        input_h5ad=getattr(cellranger, "input_h5ad", None),
        sample_root=getattr(cellranger, "sample_root", None),
    )
    detection = detect_batch_structure(
        source, preferred_key=getattr(getattr(cfg, "batch", None), "batch_key", "")
    )
    return resolve_batch_risk(
        detection=detection,
        declared_strategy=getattr(cfg, "batch_strategy", BATCH_STRATEGY_AUTO),
        planned_modules=planned_modules,
    )


# ---------------------------------------------------------------------------
# announcement
# ---------------------------------------------------------------------------


def _format_candidates(candidates: dict) -> str:
    if not candidates:
        return "no batch-like obs columns"
    return ", ".join(f"{name}={count}" for name, count in sorted(candidates.items()))


def batch_risk_warning(risk: dict) -> str:
    """Return the operator-facing warning, or "" when nothing needs saying."""
    status = risk.get("clustering_claim_status")
    if status == CLAIM_STATUS_CLAIMABLE:
        return ""

    if status == CLAIM_STATUS_UNDETERMINED:
        return (
            "BATCH_RISK: batch structure could not be read before the run "
            f"({risk.get('detection_reason') or 'no prepared obs input'}), and no "
            "batch strategy was declared. Clustering claim status is recorded as "
            "'undetermined' in the manifest (batch_risk.clustering_claim_status); "
            "clustering.py re-checks the loaded object at runtime and records "
            "batch_confounding_risk if the input turns out to span batches. "
            "Declare --batch-strategy single-batch / accept-uncorrected, or add "
            "batch_correction to --optional-modules, to resolve this up front."
        )

    downgraded = risk.get("downgraded_modules") or []
    detected = risk.get("detected_batches", 0)
    key = risk.get("key") or "<unknown>"
    declared = risk.get("strategy_declared")
    head = (
        f"BATCH_RISK: input spans {detected} batches (obs column {key!r}; "
        f"{_format_candidates(risk.get('candidate_batch_columns', {}))})"
    )
    if declared:
        head += " and the run was explicitly declared uncorrected."
    else:
        head += " and no batch strategy was declared."
    body = (
        " Clustering on unintegrated multi-batch input tracks the batch, not the "
        "biology: on the LuCA LUSC cohort (92,430 cells) the unintegrated default "
        "scored ARI 0.216 / kBET 0.117 with 47 clusters against 24 published "
        "classes, while Harmony scored ARI 0.471 / kBET 0.554 with 23 clusters "
        "(governance/realrun_gt_concordance_lusc_2026-07-28.md). Output of "
        f"{', '.join(downgraded) or 'the affected modules'} is recorded as "
        "EXPLORATORY in the manifest (batch_risk.clustering_claim_status)."
    )
    if declared:
        return head + body
    return head + body + (
        " To change that, declare a strategy:"
        f"  --optional-modules ...,batch_correction --batch-key {key}  (integrate);"
        "  --batch-strategy accept-uncorrected  (keep uncorrected; stays exploratory);"
        "  --batch-strategy single-batch  (assert single batch; fails if contradicted)."
    )


def _logging_is_configured() -> bool:
    """True when some ancestor logger owns a handler.

    When nothing is configured, ``logging.lastResort`` already writes WARNING
    records to stderr, so logging AND printing would emit the banner twice.
    """
    current: Optional[logging.Logger] = logger
    while current is not None:
        if current.handlers:
            return True
        if not current.propagate:
            return False
        current = current.parent
    return False


def _emit_once(message: str, stream, level: int) -> str:
    """Emit ``message`` to stderr exactly once, logging it too when useful.

    stderr is the contract — the operator must see this even when the launcher
    never configured logging. The logger call is additive, for hosts that do
    configure one, and is skipped when it would only duplicate the print.
    """
    if not message:
        return ""
    if _logging_is_configured():
        logger.log(level, "%s", message)
    print(message, file=stream if stream is not None else sys.stderr)
    return message


def announce_batch_risk(risk: dict, stream=None) -> str:
    """Emit the plan-time warning for ``risk`` exactly once."""
    return _emit_once(batch_risk_warning(risk), stream, logging.WARNING)


def announce_reconciliation(risk: dict, stream=None) -> str:
    """Emit the manifest-time downgrade notice for ``risk`` exactly once.

    ERROR level, not WARNING: a claim the data falsified is a stronger event
    than the plan-time "you have not declared anything yet" notice.
    """
    return _emit_once(reconciliation_warning(risk), stream, logging.ERROR)
