"""Ambient RNA correction via DecontX (celda) -- conditional QC remediation.

Lives in singlecell_factory (Python pipeline) and dispatches to an R driver
(``r_multiomics_factory/scripts/decontx_run.R``) in the ``r_multiomics`` conda
env via subprocess. Mirrors the cross-language pattern established by
``multimodal_integration.py`` + ``run_wnn.R``.

Policy: ``ops/policy/ambient_correction_policy.md`` (Phase 2 decided 2026-05-20).

Activation:
  * Conditional per-sample evaluation against 4 trigger rules (T1, T2, T4, T5).
  * If no trigger fires -> module records ``decision="skipped_no_trigger"``
    in ``adata.uns["ambient_correction"]`` and returns without invoking R.
  * If any trigger fires -> R driver runs DecontX; corrected counts replace
    ``.X``; downstream modules (doublet, clustering) operate on corrected data.

Trigger T3 (doublet-rate excess) was REMOVED 2026-05-20: the DAG ordering
places ``ambient_correction`` BEFORE ``doublet_detection`` (per
``module_catalog.py``), so ``obs["predicted_doublet"]`` is never populated
when triggers are evaluated. A doublet-driven ambient re-trigger is left to
future Phase 3 work if needed.

Cross-language debugging design:
  * **Env isolation**: R driver runs in ``r_multiomics`` conda env via
    ``conda run -n r_multiomics Rscript ...`` -- the Python parent env stays
    clean; no rpy2 bridge.
  * **Cross-factory R script location**: the R driver lives in
    ``r_multiomics_factory/scripts/decontx_run.R`` (sibling factory).
    Resolved via env var ``R_MULTIOMICS_SCRIPT_DIR`` when set, else via
    repo-relative path computation from this file.
  * **Stderr pass-through**: R stderr is captured and logged at WARN/ERROR
    levels, not swallowed.
  * **Keep-temp-on-failure**: when the R subprocess exits non-zero, the
    temp h5ad input + out-dir are preserved (not cleaned up) and the path
    is logged so the user can reproduce manually with the printed command.
  * **Repro command logging**: the exact ``conda run -n r_multiomics
    Rscript decontx_run.R --in-h5ad ... --out-dir ...`` command line is
    logged at INFO before invocation -- copy/paste reproducible.
  * **Trigger dry-run**: setting ``cfg.ambient.dry_run_triggers_only=True``
    evaluates triggers in pure Python and writes decision-only provenance,
    skipping R entirely (useful for "would this fire?" introspection).
  * **No silent skip on R-missing**: unlike multimodal_integration (which
    skips cleanly when Rscript is missing because WNN is optional/visual),
    ambient_correction RAISES when triggers fire and R is unavailable --
    silent skip on a scientific correction step would be a Principle 9
    -style silent compromise.
  * **No silent skip on missing QC**: when required QC columns
    (``pct_counts_in_top_50`` AND ``pct_counts_mt``) are absent, the module
    raises ``RuntimeError`` with ``decision="qc_metrics_unavailable"``
    recorded in ``adata.uns``. Principle 9 forbids "no triggers fired" as
    a label for "no triggers evaluable".

Status: module wired into pipeline.py / module_catalog.py as MANDATORY
(self-skips internally per trigger evaluation).
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import signal
import subprocess
import tempfile
from pathlib import Path
from typing import Any, NoReturn

import numpy as np
import pandas as pd

from ..config import AmbientCorrectionConfig
from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "yang_decontx_2020": {
        "title": "Decontamination of ambient RNA in single-cell RNA-seq with DecontX",
        "authors": "Yang S, Corbett SE, Koga Y, Wang Z, Johnson WE, Yajima M, Campbell JD",
        "journal": "Genome Biology",
        "year": "2020",
        "doi": "10.1186/s13059-020-1950-6",
        "description": "DecontX: EM-based per-cell ambient RNA contamination estimation; chosen tool per ambient_correction_policy.md Phase 2.",
    },
}


# Resolve R driver path:
#   __file__              = .../singlecell_factory/workflow/modular/modules/ambient_correction.py
#   parents[0]            = .../singlecell_factory/workflow/modular/modules
#   parents[1]            = .../singlecell_factory/workflow/modular
#   parents[2]            = .../singlecell_factory/workflow
#   parents[3]            = .../singlecell_factory
#   parents[4]            = .../Bioinformatics Research Pipeline (suite root)
_SUITE_ROOT = Path(__file__).resolve().parents[4]
_DEFAULT_R_SCRIPT_DIR = _SUITE_ROOT / "r_multiomics_factory" / "scripts"


def _resolve_r_driver_path() -> Path:
    """Resolve the DecontX R driver script path.

    Priority order:
      1. ``R_MULTIOMICS_SCRIPT_DIR`` env var (explicit override).
      2. Suite-root sibling: ``<suite>/r_multiomics_factory/scripts/decontx_run.R``.
    """
    env_dir = os.environ.get("R_MULTIOMICS_SCRIPT_DIR")
    if env_dir:
        return Path(env_dir) / "decontx_run.R"
    return _DEFAULT_R_SCRIPT_DIR / "decontx_run.R"


R_SCRIPT_PATH = _resolve_r_driver_path()

_TRIGGER_IDS = ("T1_top50_high", "T2_mt_excess",
                "T4_low_count_correlation", "T5_cross_sample_variance")


class AmbientCorrectionModule:
    """Conditional DecontX ambient correction (per-sample trigger evaluation).

    DAG position: between ``qc`` and ``doublet_detection``.
    """

    name = "ambient_correction"
    # MANDATORY-listed in module_catalog.MANDATORY_MODULES (always invoked).
    # The conditional policy decides per-sample whether DecontX actually runs.
    required = True
    mutates_structure = True  # replaces .X with decontXcounts when triggered
    # Soft contract: QC metrics are EXPECTED (produced by qc.py). When BOTH
    # pct_counts_in_top_50 AND pct_counts_mt are absent the module FAIL-LOUDS
    # (Principle 9: no silent algorithmic compromise) -- it cannot honestly
    # claim "no trigger fired" if no trigger was evaluable.
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "uns": ["ambient_correction"],
        # When triggered: also writes .obs["decontx_contamination"]
    }

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        cfg = self._resolve_config(ctx)
        logger.info("Running ambient_correction (DecontX, conditional)...")

        # ---- Step 1: hard-disable check ----
        if cfg.disable_triggers or os.environ.get("SC_AMBIENT_TRIGGERS_DISABLE") == "1":
            self._record_skip(adata, cfg, reason="disabled_by_user",
                              triggers_fired=[], qc_pre=None, ctx=ctx)
            logger.info("ambient_correction: disabled (SC_AMBIENT_TRIGGERS_DISABLE / cfg).")
            return

        # ---- Step 2: evaluate triggers ----
        qc_pre, triggers_fired = self._evaluate_triggers(adata, cfg)

        if not triggers_fired:
            self._record_skip(adata, cfg, reason="skipped_no_trigger",
                              triggers_fired=[], qc_pre=qc_pre, ctx=ctx)
            logger.info(
                "ambient_correction: no triggers fired -> SKIP (qc_pre=%s)",
                {k: round(v, 3) if isinstance(v, float) else v
                 for k, v in qc_pre.items()},
            )
            return

        logger.info("ambient_correction: triggers fired: %s -> running DecontX",
                    triggers_fired)

        # ---- Step 3: dry-run mode (skip R, record-only) ----
        if cfg.dry_run_triggers_only:
            self._record_skip(adata, cfg, reason="dry_run_triggers_only",
                              triggers_fired=triggers_fired, qc_pre=qc_pre, ctx=ctx)
            logger.warning("ambient_correction: dry-run mode -- triggers fired "
                           "but R not invoked.")
            return

        # ---- Step 4: invoke R driver (DecontX) ----
        contamination, summary, qc_post = self._invoke_decontx(adata, cfg, ctx)

        # ---- Step 5: record provenance ----
        adata.uns["ambient_correction"] = {
            "engine": "decontx",
            "tool_version": summary.get("celda_version", "unknown"),
            "decision": "triggered",
            "triggers_fired": triggers_fired,
            "trigger_thresholds": self._snapshot_thresholds(cfg),
            "qc_pre": qc_pre,
            "qc_post": qc_post,
            "contamination_summary": {
                "median": summary["median"],
                "q25": summary["q25"],
                "q75": summary["q75"],
                "max": summary["max"],
            },
            "decontx_runtime_seconds": summary.get("decontx_runtime_seconds"),
            "decontx_max_iter": cfg.decontx_max_iter,
            "decontx_seed": cfg.decontx_seed,
        }
        # Surface key numbers in manifest
        ctx.metadata["ambient_correction_engine"] = "decontx"
        ctx.metadata["ambient_correction_decision"] = "triggered"
        ctx.metadata["ambient_correction_triggers"] = ",".join(triggers_fired)
        ctx.metadata["ambient_contamination_median"] = float(summary["median"])
        ctx.metadata["ambient_contamination_max"] = float(summary["max"])

        logger.info(
            "ambient_correction: DecontX done, median contamination=%.3f, "
            "max=%.3f, runtime=%.1fs",
            summary["median"], summary["max"],
            summary.get("decontx_runtime_seconds", -1),
        )

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _resolve_config(self, ctx: PipelineContext) -> AmbientCorrectionConfig:
        """Read AmbientCorrectionConfig from ctx.cfg. Env-var override is
        applied inline at the disable check (no mutation of shared cfg)."""
        return getattr(ctx.cfg, "ambient", None) or AmbientCorrectionConfig()

    def _evaluate_triggers(
        self, adata: Any, cfg: AmbientCorrectionConfig
    ) -> tuple[dict[str, float], list[str]]:
        """Compute QC metrics summary and return list of triggers that fired.

        Fail-loud on missing QC: if BOTH ``pct_counts_in_top_50`` and
        ``pct_counts_mt`` are absent from ``adata.obs``, raise RuntimeError.
        Principle 9 ("no silent algorithmic compromise") forbids recording
        "no trigger fired" as a proxy for "no trigger was evaluable" -- the
        downstream reader of ``adata.uns["ambient_correction"]`` would be
        actively misled.
        """
        from scipy.stats import spearmanr

        obs = adata.obs
        qc_pre: dict[str, float] = {}
        fired: list[str] = []

        has_top50 = "pct_counts_in_top_50" in obs.columns
        has_mt = "pct_counts_mt" in obs.columns
        if not has_top50 and not has_mt:
            # Record the fail-loud decision in adata.uns BEFORE raising so
            # downstream forensics can read the provenance from a partial
            # checkpoint.
            adata.uns["ambient_correction"] = {
                "engine": "decontx",
                "tool_version": "celda 1.26.0 (not invoked)",
                "decision": "qc_metrics_unavailable",
                "triggers_fired": [],
                "trigger_thresholds": self._snapshot_thresholds(cfg),
                "qc_pre": {},
                "qc_post": None,
                "contamination_summary": None,
                "decontx_runtime_seconds": None,
            }
            raise RuntimeError(
                "ambient_correction: required QC metrics are unavailable -- "
                "neither 'pct_counts_in_top_50' nor 'pct_counts_mt' was found "
                "in adata.obs. Run the qc module before ambient_correction, "
                "or disable triggers explicitly via "
                "SC_AMBIENT_TRIGGERS_DISABLE=1 / --ambient-disable-triggers. "
                "Silent-skip is forbidden (Principle 9: no silent "
                "algorithmic compromise)."
            )

        if has_top50:
            top50 = float(np.median(obs["pct_counts_in_top_50"]))
            qc_pre["top50_median"] = top50
            if top50 > cfg.trigger_top50:
                fired.append("T1_top50_high")
        if has_mt:
            mt = float(np.median(obs["pct_counts_mt"]))
            qc_pre["mt_median"] = mt
            if mt > cfg.trigger_mt:
                fired.append("T2_mt_excess")
        if "total_counts" in obs.columns and "n_genes_by_counts" in obs.columns:
            rho_total_genes, _ = spearmanr(
                obs["total_counts"], obs["n_genes_by_counts"]
            )
            if np.isnan(rho_total_genes):
                rho_total_genes = 0.0
            qc_pre["count_gene_spearman"] = float(rho_total_genes)
            if rho_total_genes < cfg.trigger_count_correlation:
                fired.append("T4_low_count_correlation")

        # T3 removed 2026-05-20 -- DAG ordering (qc -> ambient_correction ->
        # doublet_detection) makes obs["predicted_doublet"] always absent when
        # this evaluator runs. See ops/policy/ambient_correction_policy.md
        # Phase 2 notes.

        # T5: cross-sample CV (cohort-level). Only evaluable if a sample obs
        # column exists. Skip if not present -- the per-cohort orchestrator
        # owns this trigger in multi-sample runs.
        if "sample_id" in obs.columns or "batch" in obs.columns:
            sample_col = "sample_id" if "sample_id" in obs.columns else "batch"
            housekeeping = [g for g in ("ACTB", "GAPDH", "B2M", "TUBB", "PPIA")
                            if g in adata.var_names]
            if housekeeping:
                expr = adata[:, housekeeping].X
                if hasattr(expr, "toarray"):
                    expr = expr.toarray()  # densify-allowed: bounded to <=5 housekeeping gene columns
                means_per_sample = (
                    np.asarray([
                        expr[obs[sample_col].values == s].mean()
                        for s in obs[sample_col].unique()
                    ])
                )
                if means_per_sample.size >= 2 and means_per_sample.mean() > 0:
                    cv = float(means_per_sample.std() / means_per_sample.mean())
                    qc_pre["housekeeping_cv"] = cv
                    if cv > cfg.trigger_cross_sample_cv:
                        fired.append("T5_cross_sample_variance")

        return qc_pre, fired

    def _invoke_decontx(
        self, adata: Any, cfg: AmbientCorrectionConfig, ctx: PipelineContext
    ) -> tuple[np.ndarray, dict[str, Any], dict[str, float]]:
        """Write h5ad -> Rscript decontx -> read back. Raises on R failure.

        Uses Popen + communicate(timeout=...) + os.killpg() so that an R
        subprocess wedged in celda::decontX is reliably terminated together
        with its grandchildren (Rscript launched via ``conda run``).
        """
        import anndata as ad

        # Preflight: conda must be on PATH; otherwise FileNotFoundError from
        # subprocess surfaces as a cryptic message at module mid-run.
        if shutil.which("conda") is None:
            raise RuntimeError(
                "ambient_correction: 'conda' executable not found on PATH. "
                "DecontX runs in the r_multiomics conda env via subprocess; "
                "either install conda + create the env via "
                "'conda env create -f r_multiomics_factory/envs/r_multiomics.yml', "
                "or set SC_AMBIENT_TRIGGERS_DISABLE=1 to bypass."
            )

        r_driver = _resolve_r_driver_path()
        if not r_driver.exists():
            raise RuntimeError(
                f"ambient_correction: DecontX R driver not found at {r_driver}. "
                "Expected location: <suite>/r_multiomics_factory/scripts/decontx_run.R. "
                "Override via env var R_MULTIOMICS_SCRIPT_DIR=<dir-containing-decontx_run.R>."
            )

        keep_temp = cfg.keep_temp_on_failure
        tmp_root = Path(tempfile.mkdtemp(prefix="ambient_decontx_"))
        mtx_dir = tmp_root / "mtx"
        mtx_dir.mkdir(parents=True, exist_ok=True)
        out_dir = tmp_root / "out"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Write filtered counts as 10X mtx triplet. We deliberately avoid
        # h5ad here: the previous contract round-tripped through
        # zellkonverter::readH5AD in R, which relies on basilisk and was
        # observed attempting to compile Python 3.14 from source via
        # pyenv during the 2026-05-21 hgmm_10k_v3 benchmark, hanging the
        # subprocess indefinitely. mtx is native to Seurat::ReadMtx with
        # no Python dependency on the R side.
        from scipy import io as _sio
        from scipy import sparse as _sp
        X = adata.X
        if not _sp.issparse(X):
            X = _sp.csr_matrix(X)
        # 10X convention: rows = genes, cols = cells.
        _sio.mmwrite(str(mtx_dir / "matrix.mtx"), X.T.tocoo())
        with open(mtx_dir / "features.tsv", "w") as f:
            for g in adata.var_names:
                f.write(f"{g}\t{g}\tGene Expression\n")
        with open(mtx_dir / "barcodes.tsv", "w") as f:
            for b in adata.obs_names:
                f.write(f"{b}\n")

        cmd = [
            "conda", "run", "-n", cfg.r_conda_env,
            "Rscript", str(r_driver),
            "--mtx-dir", str(mtx_dir),
            "--out-dir", str(out_dir),
            "--max-iter", str(cfg.decontx_max_iter),
            "--seed", str(cfg.decontx_seed),
        ]
        # Optional per-sample DecontX via a barcode\tbatch TSV.
        if cfg.batch_obs_column and cfg.batch_obs_column in adata.obs.columns:
            batch_tsv = tmp_root / "batch.tsv"
            with open(batch_tsv, "w") as f:
                for bc, b in zip(adata.obs_names,
                                 adata.obs[cfg.batch_obs_column].astype(str)):
                    f.write(f"{bc}\t{b}\n")
            cmd += ["--batch-tsv", str(batch_tsv)]

        repro_cmd_str = " ".join(cmd)
        logger.info("ambient_correction: repro command -- %s", repro_cmd_str)

        # Popen + communicate(timeout=...) so we can kill the entire process
        # group (conda + Rscript + R + grandchildren) on timeout. plain
        # subprocess.run(...).timeout only kills the immediate child.
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,  # detaches into its own process group
        )
        stdout = ""
        stderr = ""
        try:
            try:
                stdout, stderr = proc.communicate(timeout=cfg.subprocess_timeout)
            except subprocess.TimeoutExpired:
                # Kill the whole process group, then collect remaining output.
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except ProcessLookupError:
                    # Already dead between communicate timeout and killpg.
                    pass
                try:
                    stdout, stderr = proc.communicate(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    stdout, stderr = "", ""
                self._raise_with_artifacts(
                    f"DecontX timed out after {cfg.subprocess_timeout}s",
                    tmp_root, keep_temp, repro_cmd_str,
                    stderr=stderr or stdout or "<no stderr captured>",
                )
        finally:
            # Defensive: ensure no orphaned process group remains.
            if proc.poll() is None:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except (ProcessLookupError, OSError):
                    pass

        if proc.returncode != 0:
            self._raise_with_artifacts(
                f"DecontX R driver exited {proc.returncode}",
                tmp_root, keep_temp, repro_cmd_str,
                stderr=stderr or stdout,
            )

        # Read back corrected counts as mtx + contamination CSV + summary JSON.
        corrected_dir = out_dir / "corrected_mtx"
        contamination_csv = out_dir / "contamination.csv"
        summary_json = out_dir / "contamination_summary.json"
        missing = [p for p in (corrected_dir / "matrix.mtx",
                               corrected_dir / "barcodes.tsv",
                               corrected_dir / "features.tsv",
                               contamination_csv,
                               summary_json) if not p.exists()]
        if missing:
            self._raise_with_artifacts(
                f"DecontX produced exit-0 but missing outputs: {missing}",
                tmp_root, keep_temp, repro_cmd_str, stderr=stderr,
            )

        # Load corrected counts from mtx (rows=genes, cols=cells).
        from scipy.io import mmread as _mmread
        corrected_counts = _mmread(str(corrected_dir / "matrix.mtx")).tocsr()
        corrected_barcodes = (corrected_dir / "barcodes.tsv").read_text().splitlines()
        corrected_features = (corrected_dir / "features.tsv").read_text().splitlines()
        if len(corrected_barcodes) != adata.n_obs or len(corrected_features) != adata.n_vars:
            raise RuntimeError(
                f"DecontX returned matrix with shape "
                f"{(len(corrected_features), len(corrected_barcodes))}, "
                f"expected {(adata.n_vars, adata.n_obs)}. Refusing to replace .X."
            )
        if list(corrected_barcodes) != list(adata.obs_names):
            raise RuntimeError(
                "DecontX corrected barcodes do not match adata.obs_names; "
                "refusing to replace .X to prevent silent reordering."
            )
        summary = json.loads(summary_json.read_text())

        # Replace .X (counts) with decontXcounts; preserve raw in a layer.
        adata.layers["counts_raw_pre_decontx"] = adata.X.copy()
        # corrected_counts is genes×cells in 10X convention; transpose to cells×genes.
        adata.X = corrected_counts.T.tocsr()
        contam_df = pd.read_csv(contamination_csv)
        contam_df = contam_df.set_index("barcode").loc[list(adata.obs_names)]
        contamination = np.asarray(contam_df["decontx_contamination"].values,
                                   dtype=float)
        adata.obs["decontx_contamination"] = contamination
        if "decontx_clusters" in contam_df.columns:
            adata.obs["decontx_clusters"] = contam_df["decontx_clusters"].values

        # Re-compute key QC metrics on corrected counts for qc_post.
        # Inline (no shared qc.py helper exists yet -- if a future refactor
        # introduces one, this is the place to switch).
        qc_post = self._inline_qc_post(adata)

        # Clean up temp dir on success
        if not keep_temp:
            shutil.rmtree(tmp_root, ignore_errors=True)
        else:
            logger.info("ambient_correction: temp artifacts preserved at %s",
                        tmp_root)

        return contamination, summary, qc_post

    def _inline_qc_post(self, adata: Any) -> dict[str, float]:
        """Sparse-safe minimal QC recompute for qc_post.

        Avoids densifying the full count matrix. Computes ``total_counts_median``
        and ``n_genes_median`` directly from the (sparse or dense) ``.X``;
        skips ``top50_median`` in this fallback because the per-row top-50
        metric requires a sort that has no efficient sparse equivalent --
        provenance can rely on the pre-correction T1 reading in ``qc_pre``.
        """
        X = adata.X
        if hasattr(X, "sum") and hasattr(X, "indices"):
            # scipy sparse path -- use sparse aggregations, no densification.
            total = np.asarray(X.sum(axis=1)).ravel()
            n_genes = np.asarray((X != 0).sum(axis=1)).ravel()
        else:
            arr = np.asarray(X)
            total = arr.sum(axis=1)
            n_genes = (arr > 0).sum(axis=1)
        return {
            "total_counts_median": float(np.median(total)),
            "n_genes_median": float(np.median(n_genes)),
        }

    def _record_skip(
        self, adata: Any, cfg: AmbientCorrectionConfig,
        reason: str, triggers_fired: list[str], qc_pre: dict | None,
        ctx: PipelineContext | None = None,
    ) -> None:
        adata.uns["ambient_correction"] = {
            "engine": "decontx",
            "tool_version": "celda 1.26.0 (not invoked)",
            "decision": reason,
            "triggers_fired": triggers_fired,
            "trigger_thresholds": self._snapshot_thresholds(cfg),
            "qc_pre": qc_pre,
            "qc_post": None,
            "contamination_summary": None,
            "decontx_runtime_seconds": None,
        }
        if ctx is not None:
            ctx.metadata["ambient_correction_engine"] = "decontx"
            ctx.metadata["ambient_correction_decision"] = reason
            ctx.metadata["ambient_correction_triggers"] = ",".join(triggers_fired)
            if qc_pre:
                for key, value in qc_pre.items():
                    if (
                        isinstance(value, (int, float, np.integer, np.floating))
                        and np.isfinite(value)
                    ):
                        ctx.metadata[f"ambient_qc_pre_{key}"] = float(value)

    @staticmethod
    def _snapshot_thresholds(cfg: AmbientCorrectionConfig) -> dict[str, float]:
        return {
            "T1_top50": cfg.trigger_top50,
            "T2_mt": cfg.trigger_mt,
            "T4_count_correlation": cfg.trigger_count_correlation,
            "T5_cross_sample_cv": cfg.trigger_cross_sample_cv,
        }

    def _raise_with_artifacts(
        self, msg: str, tmp_root: Path, keep_temp: bool,
        repro_cmd: str, stderr: str | None,
    ) -> NoReturn:
        """Raise RuntimeError preserving temp dir + repro hint."""
        if keep_temp:
            preserved = f"temp dir preserved at {tmp_root}"
        else:
            shutil.rmtree(tmp_root, ignore_errors=True)
            preserved = "temp dir cleaned (set keep_temp_on_failure=True to retain)"
        stderr_tail = ""
        if stderr:
            tail = "\n".join(stderr.splitlines()[-20:])
            stderr_tail = f"\n  R stderr (last 20 lines):\n{tail}"
        raise RuntimeError(
            f"ambient_correction: {msg}.\n"
            f"  repro: {repro_cmd}\n"
            f"  {preserved}{stderr_tail}"
        )
