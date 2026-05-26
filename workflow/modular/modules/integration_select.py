"""Per-run discovery integration-selection gate as a pipeline module (plan §D5/§D6, TM-10).

This OPT-IN module runs AFTER clustering and BEFORE batch_correction. It scores
the candidate integration embeddings label-free (baseline X_pca, Harmony, a scVI
seed sweep, plus a shuffle falsifiability control) on the discovery metric set
and applies the pre-registered disqualifying-gate rule to RECOMMEND a
``cfg.batch.method`` for THIS dataset. The recommendation is APPLIED to
``ctx.cfg.batch.method`` so the downstream batch_correction module honors it.

When the gate recommends Harmony it ALSO sets
``ctx.cfg.batch.harmony_backend = "direct"`` so batch_correction uses the
proven-working ``harmonypy.run_harmony`` path (the production scanpy/rapids
Harmony backends are broken in the sc_gpu env — see batch_correction
``_run_harmony_direct``).

Cost control (plan §D7): the scVI seed sweep is expensive, so the per-run cost
is bounded by a default-ON filesystem cache keyed on
``(data_hash, batch_key, scVI_config, code_version)``. A cache HIT replays the
prior recommendation with NO scVI retrain.

Reuses (never reimplements) the reviewed engine in scripts/bench/integration:
``methods`` (production-backed embedding helpers + _CtxShim),
``discovery_metrics`` (label-free scoring), ``select_integration`` (the gate
rule + output writers), ``integration_cache`` (the per-run cache).

Determinism: scVI sweep uses the fixed seeds in cfg.batch.integration_scvi_seeds;
the shuffle control uses seed=0.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ..context import PipelineContext

logger = logging.getLogger(__name__)


# Factory tree root (workflow/modular/modules/integration_select.py -> repo root).
_FACTORY_ROOT = Path(__file__).resolve().parents[3]


class IntegrationSelectModule:
    """Optional module: discovery-driven integration method selection gate.

    Runs after clustering and (ordering-only) before batch_correction. Sets
    ``cfg.batch.method`` from the per-run recommendation; opt-in via
    ``--select-integration`` / ``cfg.batch.select_integration``.
    """

    name = "integration_select"
    # NOTE: intentionally NO framework-level `requires_keys` for X_pca. Declaring
    # it would make _run_module SKIP this module (logger.warning + _SkipModule)
    # when X_pca is absent, silently falling through to the default
    # cfg.batch.method (Harmony on the broken non-direct backend). The module
    # validates X_pca itself in run() and FAILS LOUD instead (review MEDIUM-1).

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("integration_select requires AnnData.")

        cfg_batch = ctx.cfg.batch
        batch_key = cfg_batch.batch_key

        # --- Skip gracefully when there is no multi-batch structure to gate on.
        if batch_key not in adata.obs.columns:
            ctx.metadata["integration_select_status"] = "skipped_missing_batch_key"
            ctx.status(
                self.name,
                "skipped",
                f"Batch key '{batch_key}' not found in adata.obs; leaving "
                f"cfg.batch.method='{cfg_batch.method}' unchanged.",
            )
            return
        n_batches = int(adata.obs[batch_key].nunique())
        if n_batches < 2:
            ctx.metadata["integration_select_status"] = "skipped_single_batch"
            ctx.status(
                self.name,
                "skipped",
                f"Only one batch in '{batch_key}'; no integration needed "
                f"(leaving cfg.batch.method='{cfg_batch.method}').",
            )
            return

        # --- Fail loud: clustering must have produced PCA first.
        if "X_pca" not in adata.obsm:
            raise ValueError(
                "integration_select requires PCA in adata.obsm['X_pca'] "
                "(run clustering before integration_select)."
            )

        # Reuse the reviewed engine (import here to keep package import light).
        from scripts.bench.integration import discovery_metrics as dm
        from scripts.bench.integration import integration_cache as ic
        from scripts.bench.integration import select_integration as sel

        out_dir = self._resolve_out_dir(ctx)
        cache_dir = out_dir / "cache"
        cache_dir.mkdir(parents=True, exist_ok=True)

        batch = adata.obs[batch_key].astype(str).to_numpy()
        baseline = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
        scvi_seeds = tuple(int(s) for s in cfg_batch.integration_scvi_seeds)
        margin_mix = float(cfg_batch.integration_margin_mix)

        cache = ic.IntegrationCache(cache_dir)
        cache_key = ic.make_cache_key(baseline, batch_key, scvi_seeds=scvi_seeds)
        cache_res = cache.lookup(cache_key)

        if cache_res.hit and cache_res.payload is not None:
            chosen = str(cache_res.payload.get("chosen_method", "none"))
            self._apply_choice(ctx, chosen)
            ctx.metadata["integration_select"] = {
                "chosen_method": chosen,
                "source": "cache",
                "cache_key": cache_key.digest(),
                "rule_constants": cache_res.payload.get("rule_constants"),
                "reason": cache_res.payload.get("reason"),
            }
            ctx.metadata["integration_select_status"] = "completed_cache_hit"
            ctx.status(
                self.name,
                "ok",
                f"cache HIT: applied cfg.batch.method='{chosen}' "
                f"(no scVI retrain).",
            )
            return

        # --- Cache MISS: compute candidate embeddings on the production backends.
        embeddings = self._compute_embeddings(adata, batch_key, scvi_seeds)

        # Score + recommend via the reviewed gate.
        results, ref = dm.score_all_embeddings_discovery(
            embeddings, embeddings[sel.BASELINE_METHOD], batch
        )
        controls = sel.check_controls_fire(results)
        rec = sel.recommend(results, margin_mix=margin_mix)

        try:
            resolved_versions = dm.assert_pinned_discovery_engine()
        except Exception as exc:  # noqa: BLE001 — record drift, do not crash the gate write
            resolved_versions = {"version_assert_error": str(exc)}

        payload = sel.build_recommendation_json(
            rec, ref, resolved_versions, batch_key=batch_key,
            cache_key=cache_key.digest(),
            extra={
                "controls": controls,
                "n_batches": n_batches,
                "scvi_seeds": list(scvi_seeds),
            },
        )
        sel.write_recommendation_json(payload, out_dir / "integration_recommendation.json")
        sel.write_discovery_scoreboard(results, out_dir / "integration_scoreboard.csv")
        sel.write_audit_md(payload, out_dir / "integration_audit.md")
        (out_dir / "integration_audit.json").write_text(
            __import__("json").dumps(payload, indent=2) + "\n"
        )
        cache.store(cache_key, payload)

        self._apply_choice(ctx, rec.chosen_method)
        ctx.metadata["integration_select"] = {
            "chosen_method": rec.chosen_method,
            "source": "computed",
            "reason": rec.reason,
            "cache_key": cache_key.digest(),
            "rule_constants": rec.rule_constants,
            "shuffle_control_fired": bool(controls.get("shuffle_fired")),
        }
        ctx.metadata["integration_select_status"] = "completed"
        ctx.status(
            self.name,
            "ok",
            f"applied cfg.batch.method='{rec.chosen_method}' "
            f"(reason={rec.reason}); outputs -> {out_dir}",
        )

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _resolve_out_dir(self, ctx: PipelineContext) -> Path:
        """Resolve the module output dir; fail loud if it lands in the factory tree.

        Uses the codebase's per-module dir convention (ctx.set_module_dir), which
        lives under ``run_dir`` (always resolved under the project root, never
        the factory). The factory-tree guard makes the CLAUDE.md "never write in
        the factory tree" invariant explicit and fail-loud.
        """
        ctx.set_module_dir(self.name)
        out_dir = Path(ctx.table_dir).resolve()
        if _FACTORY_ROOT == out_dir or _FACTORY_ROOT in out_dir.parents:
            raise ValueError(
                f"integration_select output dir {out_dir} is inside the factory "
                f"tree {_FACTORY_ROOT}; outputs must never be written inside the "
                "factory (CLAUDE.md). Pass --project-root pointing under "
                "/home/zerlinshen/projects/."
            )
        return out_dir

    @staticmethod
    def _compute_embeddings(adata, batch_key: str, scvi_seeds: tuple[int, ...]) -> dict:
        """Build the candidate embeddings dict via the production-backed helpers.

        baseline = X_pca (already present); harmony via harmonypy-direct; scVI
        seed sweep; shuffle control (required falsifiability anchor). Reuses
        scripts/bench/integration/methods.py call signatures verbatim.
        """
        from scripts.bench.integration import methods as M

        embeddings: dict[str, np.ndarray] = {
            M.BASELINE_METHOD: np.asarray(adata.obsm["X_pca"], dtype=np.float64),
        }

        harmony = M.compute_harmonypy_direct(adata, batch_key)
        if harmony.X is not None:
            embeddings[M.HARMONY_METHOD] = np.asarray(harmony.X, dtype=np.float64)
        else:
            logger.warning(
                "integration_select: Harmony embedding failed (%s); harmony "
                "candidate excluded from the gate.", harmony.error,
            )

        scvi_runs = M.compute_scvi_seed_sweep(adata, batch_key, seeds=scvi_seeds)
        for run in scvi_runs:
            if run.X is not None:
                embeddings[run.method] = np.asarray(run.X, dtype=np.float64)
            else:
                logger.warning(
                    "integration_select: scVI %s failed (%s); excluded.",
                    run.method, run.error,
                )

        shuffle = M.compute_shuffle_label_control(adata, batch_key, seed=0)
        if shuffle.X is not None:
            from scripts.bench.integration import select_integration as sel
            embeddings[sel.SHUFFLE_CONTROL_METHOD] = np.asarray(
                shuffle.X, dtype=np.float64
            )

        return embeddings

    @staticmethod
    def _apply_choice(ctx: PipelineContext, chosen_method: str) -> None:
        """Apply the recommendation to cfg.batch (and route Harmony to direct)."""
        ctx.cfg.batch.method = chosen_method
        if chosen_method == "harmony":
            # Use the proven-working direct backend (W1); the auto/gpu/cpu paths
            # are broken in sc_gpu on this hardware.
            ctx.cfg.batch.harmony_backend = "direct"
