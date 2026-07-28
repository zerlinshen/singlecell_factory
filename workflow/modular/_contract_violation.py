"""Wave 3 / US-W3-2 — Pipeline contract violation poison.

When a GPU module fails post-host-mutation and the configured
SC_GPU_FAILURE_POLICY is `raise`, the implementation MUST:

  1. Poison the adata: ``adata.X = None`` AND ``adata.uns["__corrupted__"] = True``.
  2. Raise ``ClusteringContractViolation`` (an ``Exception`` subclass — not
     ``BaseException``, so legitimate ``finally`` cleanup still runs).

The pipeline orchestrator (`pipeline._run_module`) checks
``ctx.adata.uns.get("__corrupted__")`` at every module entry and aborts
immediately with a fresh ``ClusteringContractViolation`` if poisoned. This
prevents any caller that accidentally swallows the original exception from
sneaking a corrupted adata to a downstream module.

Wave 5 / US-W5-FOUNDATION adds ``ModuleContractError`` as a sibling error
class for precondition failures across multiple optional modules. Wave 6+ may
consolidate ``ClusteringContractViolation`` and ``ModuleContractError`` into a
shared ``ContractRegistry``.
"""
from __future__ import annotations


class ClusteringContractViolation(Exception):
    """Raised when a GPU clustering failure leaves adata in an unsafe state."""


class ModuleContractError(Exception):
    """Raised when a module's required input contract is violated at runtime.

    Canonical violation sites (Wave 5):

    - ``atac_lsi``: ``adata.layers["atac_peaks"]`` is absent.
    - ``multimodal_integration``: explicit-config mode requires a second obsm
      key (SC_MULTIMODAL_SECOND_OBSM) that is absent from ``adata.obsm``.
    - ``peak_to_gene``: cell axis of the linkage matrix is dense (sparse-axis
      discipline violated).
    - ``trajectory``: neither ``adata.obsm["X_wnn"]`` nor
      ``adata.obsm["X_pca"]`` is present, so no joint embedding is available
      for pseudotime computation.

    Relationship to ``ClusteringContractViolation``: sibling errors. Both are
    plain ``Exception`` subclasses raised when pipeline preconditions are
    unmet. Wave 6+ may consolidate these into a ``ContractRegistry``.
    """


CORRUPTED_FLAG = "__corrupted__"


def poison_adata(adata, reason: str) -> None:
    """Mark an adata as corrupted: null X and set the orchestrator poison flag."""
    if adata is None:
        return
    try:
        adata.X = None
    except Exception:
        pass
    try:
        adata.uns[CORRUPTED_FLAG] = True
        adata.uns[CORRUPTED_FLAG + "_reason"] = reason
    except Exception:
        pass


def assert_not_corrupted(adata, module_name: str) -> None:
    """Orchestrator guard: raise if adata carries the poison flag.

    Called by `pipeline._run_module` immediately before delegating to a
    module's `.run(ctx)`. Layered defense: if a caller catches the original
    ClusteringContractViolation and tries to continue, the next module
    invocation aborts here.
    """
    if adata is None:
        return
    try:
        if adata.uns.get(CORRUPTED_FLAG):
            reason = adata.uns.get(CORRUPTED_FLAG + "_reason", "unknown")
            raise ClusteringContractViolation(
                f"Cannot run module '{module_name}': adata is poisoned "
                f"({CORRUPTED_FLAG}=True; reason: {reason}). "
                f"A prior module raised ClusteringContractViolation; the pipeline "
                f"MUST NOT continue. Restart with --checkpoint + --resume-from to recover."
            )
    except ClusteringContractViolation:
        raise
    except Exception:
        # Defensive: if uns is not subscriptable for any reason, do not crash
        # the guard itself — best-effort enforcement.
        pass


def resolve_gpu_failure_policy(cfg) -> str:
    """Resolve the active GPU-failure policy.

    Precedence: SC_GPU_FAILURE_POLICY env var > cfg.gpu_failure_policy >
    a default derived from ``cfg.gpu_mode``.

    The derived default matters because ``--gpu-mode auto`` is an *opportunistic*
    choice: the pipeline elects the GPU on the operator's behalf, so it must also
    be able to retract that choice on the operator's behalf. Pairing an
    opportunistic selection with an unconditional ``raise`` meant a GPU-only
    fault destroyed the whole scientific payload — on this workstation
    ``rsc.pp.pca`` is a documented cuSOLVER incompatibility
    (docs/HOTSPOT1_DIAGNOSIS.md), so a default run lost clustering, annotation
    and differential expression to a failure the CPU path handles fine.

    ``auto``  -> ``restore-cpu``  (retract the automatic choice; M2 preserves
                                   ``adata.raw`` precisely so this is safe)
    ``force`` -> ``raise``        (the operator demanded GPU; do not silently
                                   deliver something else)

    Valid values: raise, restore-cpu, reload-checkpoint.
    """
    import os
    raw = os.environ.get("SC_GPU_FAILURE_POLICY", "").strip().lower()
    if not raw:
        raw = (getattr(cfg, "gpu_failure_policy", "") or "").strip().lower()
    if not raw:
        raw = "raise" if (getattr(cfg, "gpu_mode", "auto") or "auto").lower() == "force" else "restore-cpu"
    if raw not in {"raise", "restore-cpu", "reload-checkpoint"}:
        raise ValueError(
            f"Invalid SC_GPU_FAILURE_POLICY={raw!r}. "
            f"Expected one of: raise, restore-cpu, reload-checkpoint."
        )
    return raw
