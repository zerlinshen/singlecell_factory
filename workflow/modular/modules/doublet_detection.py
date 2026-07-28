from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from ._gpu_utils import bind_cuda_context, gpu_available
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext


__references__ = {
    "Wolock_Scrublet_2019": {
        "title": "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data",
        "authors": "Wolock, Lopez, Klein",
        "journal": "Cell Systems",
        "year": "2019",
        "doi": "10.1016/j.cels.2018.11.005",
        "description": "Reference implementation; supports whole-dataset and per-sample (grouped) doublet rate calibration.",
    },
    "McGinnis_DoubletFinder_2019": {
        "title": "DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors",
        "authors": "McGinnis, Murrow, Gartner",
        "journal": "Cell Systems",
        "year": "2019",
        "doi": "10.1016/j.cels.2019.03.003",
        "description": (
            "Second-opinion backend invoked via the r_multiomics conda env using "
            "subprocess. DoubletFinder is sourced directly from a local clone "
            "(DOUBLETFINDER_R_PATH env, default "
            "/home/zerlinshen/downloads/external_refs/DoubletFinder/R) rather "
            "than installed as a package, so no remotes::install_github is "
            "required. Added 2026-05-21 after a head-to-head benchmark on "
            "LUSC PS01 (4,672 cells) revealed Scrublet's bimodality "
            "auto-threshold under-called by ~40x on that slice."
        ),
    },
    "Germain_scDblFinder_2021": {
        "title": "Doublet identification in single-cell sequencing data using scDblFinder",
        "authors": "Germain, Lun, Garcia Meixide, Macnair, Robinson",
        "journal": "F1000Research",
        "year": "2021",
        "description": (
            "Second-opinion Bioconductor backend invoked via the r_multiomics "
            "conda env using the same 10X mtx bridge as DoubletFinder. Added "
            "for data-shape-specific recovery when Scrublet under-calls on "
            "tumor/tissue datasets."
        ),
    },
}


logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# GPU availability probe (evaluated once at import time)
# ---------------------------------------------------------------------------
try:
    import rapids_singlecell as _rsc  # noqa: F401
    _RSC_AVAILABLE = True
except ImportError:
    _RSC_AVAILABLE = False


class DoubletDetectionModule:
    """Mandatory module: detect and remove doublets using Scrublet.

    Should run after QC and before clustering. Doublets (two cells captured
    in one droplet) can create artificial intermediate clusters and corrupt
    downstream differential expression results.

    GPU path (preferred): uses ``rsc.pp.scrublet`` when rapids-singlecell is
    available.  Falls back to CPU ``scrublet.Scrublet`` only if rsc is not
    importable.
    """

    name = "doublet_detection"
    required = True

    @staticmethod
    def _n_prin_comps(n_obs: int, n_vars: int) -> int:
        return max(2, min(30, n_obs - 1, n_vars - 1))

    @staticmethod
    def _resolved_expected_rate(cfg, n_obs: int) -> float:
        """Cell-count-scaled expected doublet rate.

        10x Chromium multiplet rate ~= 0.8% per 1000 cells (10x v3 User Guide;
        scDblFinder, Germain 2021). When cfg.scale_expected_doublet_rate is True,
        rate = min(0.08, 0.008 * n_obs/1000) — this is the point of the fix: a
        fixed 6% over-calls on ~6-7k-cell samples. When scaling is disabled the
        fixed cfg.expected_doublet_rate (default 0.06) is used as the floor/default.
        """
        if not bool(getattr(cfg, "scale_expected_doublet_rate", False)):
            return float(getattr(cfg, "expected_doublet_rate", 0.06))
        return min(0.08, 0.008 * float(n_obs) / 1000.0)

    @classmethod
    def _resolve_expected_rate_contract(cls, ctx, cfg, adata) -> float:
        """Resolve the sample-aware prior once and persist its audit contract.

        Backends that operate on the combined matrix accept one scalar prior.
        Resolve each sample's cell-count-aware candidate first, then use their
        cell-weighted mean as the single effective prior supplied to every
        backend, rank consensus, and the under-call diagnostic.
        """
        requested = float(getattr(cfg, "expected_doublet_rate", 0.06))
        if "sample" in adata.obs.columns:
            sample_sizes = adata.obs["sample"].astype(str).value_counts(sort=False)
        else:
            sample_sizes = pd.Series({"__all__": int(adata.n_obs)}, dtype="int64")

        resolved_per_sample = {
            str(sample): cls._resolved_expected_rate(cfg, int(n_obs))
            for sample, n_obs in sample_sizes.items()
        }
        total_cells = int(sample_sizes.sum())
        if total_cells > 0:
            effective = sum(
                int(sample_sizes[sample]) * resolved_per_sample[str(sample)]
                for sample in sample_sizes.index
            ) / total_cells
        else:
            effective = requested

        effective = round(float(effective), 6)
        ctx.metadata["doublet_expected_rate_requested"] = round(requested, 6)
        ctx.metadata["doublet_expected_rate_resolved_per_sample"] = {
            sample: round(float(rate), 6)
            for sample, rate in resolved_per_sample.items()
        }
        ctx.metadata["doublet_expected_rate_sample_sizes"] = {
            str(sample): int(n_obs) for sample, n_obs in sample_sizes.items()
        }
        ctx.metadata["doublet_expected_rate_effective"] = effective
        ctx.metadata["doublet_expected_rate_aggregation_rule"] = (
            "cell_weighted_mean_of_resolved_per_sample"
        )
        ctx.metadata["doublet_rate_scaling_enabled"] = bool(
            getattr(cfg, "scale_expected_doublet_rate", False)
        )
        return effective

    @staticmethod
    def _record_backend_expected_rate(ctx, backend: str, mode: str = "aggregate"):
        """Record and return the already-resolved prior used by ``backend``."""
        if "doublet_expected_rate_effective" not in ctx.metadata:
            raise RuntimeError(
                "doublet expected-rate contract must be resolved before backend execution"
            )
        requested = float(ctx.metadata["doublet_expected_rate_requested"])
        effective = float(ctx.metadata["doublet_expected_rate_effective"])
        if mode == "per_sample":
            resolved = dict(ctx.metadata["doublet_expected_rate_resolved_per_sample"])
            record = {
                "requested": round(requested, 6),
                "mode": "per_sample",
                "resolved_per_sample": resolved,
                "sample_sizes": dict(ctx.metadata["doublet_expected_rate_sample_sizes"]),
            }
            ctx.metadata.setdefault("doublet_expected_rate_by_backend", {})[backend] = record
            ctx.metadata["doublet_expected_rate_used_per_sample"] = resolved
            return resolved
        if mode != "aggregate":
            raise ValueError(f"unknown doublet expected-rate mode: {mode}")
        record = {
            "requested": round(requested, 6),
            "mode": "aggregate",
            "effective": round(effective, 6),
            "aggregation_rule": ctx.metadata["doublet_expected_rate_aggregation_rule"],
        }
        ctx.metadata.setdefault("doublet_expected_rate_by_backend", {})[backend] = record
        ctx.metadata["doublet_expected_rate_used"] = round(effective, 6)
        return effective

    @staticmethod
    def _materialize_counts_matrix(x):
        """Return a CPU CSR matrix, materializing lazy/dask/GPU arrays."""
        if sparse.issparse(x):
            return x.tocsr()
        if hasattr(x, "to_memory"):
            try:
                x = x.to_memory()
            except Exception:
                pass
        if sparse.issparse(x):
            return x.tocsr()
        if hasattr(x, "compute"):
            try:
                x = x.compute()
            except Exception:
                pass
        if sparse.issparse(x):
            return x.tocsr()
        return np.asarray(x)

    @staticmethod
    def _matrix_looks_integer_like(matrix, sample_size: int = 2048) -> bool:
        """Heuristic: does `matrix` resemble raw non-negative integer counts?

        Mirrors batch_correction._matrix_looks_integer_like so the doublet
        raw-counts resolver does not need to import the batch module.
        """
        if sparse.issparse(matrix):
            values = matrix.data
        else:
            values = np.asarray(matrix).ravel()
        if values.size == 0:
            return True
        if values.size > sample_size:
            idx = np.linspace(0, values.size - 1, sample_size, dtype=int)
            values = values[idx]
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return False
        if np.any(finite < 0):
            return False
        return bool(np.allclose(finite, np.round(finite), atol=1e-6))

    @classmethod
    def _resolve_raw_counts_matrix(cls, adata):
        """Return a raw-integer-UMI count matrix for Scrublet, or RAISE.

        Scrublet's synthetic-doublet model assumes raw integer UMI counts.
        ambient_correction overwrites adata.X with FRACTIONAL decontXcounts
        (raw preserved at layers['counts_raw_pre_decontx']), so feeding adata.X
        directly violates the model. Priority resolver (no silent fallback):
          1. layers['counts']                  (pseudobulk-DE raw layer)
          2. layers['counts_raw_pre_decontx']  (ambient-correction snapshot)
          3. adata.X IFF it looks integer-like  (ambient never fired)
          4. loud RAISE                         (Principle 9: no silent compromise)

        Returns (matrix, source_label). `matrix` shares dtype/sparsity with the
        chosen source (caller swaps it onto adata.X under try/finally).
        """
        layers = getattr(adata, "layers", {})
        if "counts" in layers:
            return layers["counts"], "layers.counts"
        if "counts_raw_pre_decontx" in layers:
            return layers["counts_raw_pre_decontx"], "layers.counts_raw_pre_decontx"
        if cls._matrix_looks_integer_like(adata.X):
            return adata.X, "adata.X_integer_like"
        raise RuntimeError(
            "doublet_detection: no raw integer-UMI count source for Scrublet. "
            "adata.X is not integer-like (likely fractional ambient-corrected "
            "counts) and neither layers['counts'] nor "
            "layers['counts_raw_pre_decontx'] is present. Scrublet's "
            "synthetic-doublet model requires raw counts; refusing to score on "
            "fractional input (Principle 9: no silent algorithmic compromise)."
        )

    @staticmethod
    def _resolve_doublet_strategy(ctx, adata) -> str:
        """Return the effective doublet strategy string.

        Priority: SC_DOUBLET_STRATEGY env var > cfg.doublet_strategy > "auto".
        "auto" → "grouped" when n_obs >= 100k and sample column with >= 2 labels exists,
                 otherwise "whole".
        """
        import os as _os
        import pandas as _pd
        flag = (
            _os.environ.get("SC_DOUBLET_STRATEGY", "").strip().lower()
            or getattr(ctx.cfg, "doublet_strategy", "auto")
        )
        if flag in {"grouped", "whole", "skip"}:
            return flag
        # "auto"
        has_sample = (
            "sample" in adata.obs.columns
            and _pd.Series(adata.obs["sample"]).nunique() >= 2
        )
        if adata.n_obs >= 100000 and has_sample:
            return "grouped"
        return "whole"

    # ------------------------------------------------------------------
    # GPU (rsc) paths
    # ------------------------------------------------------------------

    @staticmethod
    def _free_gpu_memory_pool() -> None:
        # HIGH fix (2026-05-19): release CuPy GPU memory after each RAPIDS
        # doublet call. Previously the GPU AnnData copy stayed resident until
        # Python GC eventually collected it, which accumulated across modules
        # (clustering, batch correction, DE) and pushed long pipelines toward
        # GPU OOM. CuPy's default mempool is process-wide so explicit free is
        # both safe and effective.
        try:
            import cupy as cp  # type: ignore
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass

    @staticmethod
    def _clone_for_gpu_lite(adata, extra_obs_cols=()):
        """Minimal AnnData carrying only X (+ requested obs cols) for GPU Scrublet.

        Scrublet (Wolock 2019) only consumes the counts matrix; a full
        ``adata.copy()`` duplicates obs/var/obsm/layers into RAM/VRAM needlessly.
        Analogous to (not identical to) clustering.ClusteringModule._clone_for_gpu_lite:
        that clone shares obsm/varm/uns/layers by reference, whereas this one carries
        ONLY a deep-copied X (rsc mutates it on device) + var_names + the obs columns
        Scrublet needs (e.g. the grouped batch_key), since Scrublet consumes nothing
        else. Scores are unchanged.
        """
        import anndata as ad
        import pandas as pd

        X = adata.X
        if X is None:
            X_clone = None
        elif sparse.issparse(X):
            X_clone = X.copy()
        else:
            X_clone = np.array(X, copy=True)
        obs = pd.DataFrame(index=adata.obs_names)
        for col in extra_obs_cols:
            if col is not None and col in adata.obs.columns:
                obs[col] = adata.obs[col].values
        var = pd.DataFrame(index=adata.var_names)
        return ad.AnnData(X=X_clone, obs=obs, var=var)

    def _run_rsc_scrublet_whole(self, adata, cfg, ctx) -> tuple[np.ndarray, np.ndarray, float | None]:
        """Run rsc.pp.scrublet on the full dataset."""
        import rapids_singlecell as rsc
        bind_cuda_context()  # ensure cuBLAS/cuSOLVER pre-warm before GPU ops (P1 fix; no-op if already warmed in probe)

        expected_rate = self._record_backend_expected_rate(
            ctx, "scrublet", mode="aggregate"
        )
        adata_gpu = self._clone_for_gpu_lite(adata)
        try:
            rsc.get.anndata_to_GPU(adata_gpu)
            rsc.pp.scrublet(
                adata_gpu,
                expected_doublet_rate=expected_rate,
                n_prin_comps=self._n_prin_comps(adata.n_obs, adata.n_vars),
                random_state=ctx.random_state,
                verbose=False,
            )
            scores = self._col_to_numpy(adata_gpu.obs["doublet_score"]).astype(np.float32)
            predicted = self._col_to_numpy(adata_gpu.obs["predicted_doublet"]).astype(bool)
            threshold = adata_gpu.uns.get("scrublet", {}).get("threshold", None)
            ctx.metadata["doublet_method"] = "scrublet_gpu"
            return scores, predicted, threshold
        finally:
            del adata_gpu
            self._free_gpu_memory_pool()

    def _run_rsc_scrublet_grouped(self, adata, cfg, ctx) -> tuple[np.ndarray, np.ndarray, None]:
        """Run rsc.pp.scrublet with batch_key for per-sample doublet detection."""
        import rapids_singlecell as rsc
        bind_cuda_context()  # ensure cuBLAS/cuSOLVER pre-warm before GPU ops (P1 fix; no-op if already warmed in probe)

        sample_key = "sample" if "sample" in adata.obs.columns else None
        if sample_key is None:
            raise ValueError("grouped scrublet requested without sample labels")

        # rsc.pp.scrublet(batch_key=...) accepts one scalar across batches. Use
        # the sample-aware, cell-weighted effective prior resolved once in run().
        sample_sizes = adata.obs[sample_key].value_counts()
        expected_rate = self._record_backend_expected_rate(
            ctx, "scrublet", mode="aggregate"
        )
        ctx.metadata["doublet_expected_rate_applied_per_sample"] = {
            str(sample): round(expected_rate, 6) for sample in sample_sizes.index
        }
        adata_gpu = self._clone_for_gpu_lite(adata, extra_obs_cols=(sample_key,))
        try:
            rsc.get.anndata_to_GPU(adata_gpu)
            rsc.pp.scrublet(
                adata_gpu,
                batch_key=sample_key,
                expected_doublet_rate=expected_rate,
                n_prin_comps=self._n_prin_comps(adata.n_obs, adata.n_vars),
                random_state=ctx.random_state,
                verbose=False,
            )
            scores = self._col_to_numpy(adata_gpu.obs["doublet_score"]).astype(np.float32)
            predicted = self._col_to_numpy(adata_gpu.obs["predicted_doublet"]).astype(bool)
            ctx.metadata["doublet_grouped_key"] = sample_key
            ctx.metadata["doublet_method"] = "scrublet_gpu_grouped"
            return scores, predicted, None
        finally:
            del adata_gpu
            self._free_gpu_memory_pool()

    # ------------------------------------------------------------------
    # CPU fallback paths (used only when rsc is unavailable)
    # ------------------------------------------------------------------

    def _run_cpu_scrublet_grouped(self, adata, scrublet_cls, cfg, ctx) -> tuple[np.ndarray, np.ndarray, None]:
        sample_key = "sample" if "sample" in adata.obs.columns else None
        if sample_key is None:
            raise ValueError("grouped scrublet requested without sample labels")
        random_state = ctx.random_state
        labels = pd.Series(adata.obs[sample_key].astype(str), index=adata.obs_names)
        scores = np.zeros(adata.n_obs, dtype=np.float32)
        predicted = np.zeros(adata.n_obs, dtype=bool)
        thresholds = {}
        scrublet_params = {
            "min_counts": 2,
            "min_cells": 3,
            "min_gene_variability_pctl": 85,
            "n_prin_comps": self._n_prin_comps(adata.n_obs, adata.n_vars),
        }
        expected_rates = self._record_backend_expected_rate(
            ctx, "scrublet", mode="per_sample"
        )
        per_sample_rates: dict[str, float] = {}
        for group in labels.unique():
            idx = np.where(labels.values == group)[0]
            if len(idx) < 20:
                continue
            group_rate = float(expected_rates[str(group)])
            per_sample_rates[str(group)] = round(float(group_rate), 6)
            subX = self._materialize_counts_matrix(adata.X[idx])
            scrub = scrublet_cls(subX, expected_doublet_rate=group_rate, random_state=random_state)
            try:
                s, p = scrub.scrub_doublets(**{**scrublet_params, "n_prin_comps": self._n_prin_comps(len(idx), adata.n_vars)})
                scores[idx] = s.astype(np.float32, copy=False)
                predicted[idx] = p.astype(bool, copy=False)
                thr = getattr(scrub, "threshold_", None)
                if thr is not None:
                    thresholds[str(group)] = float(thr)
            except Exception as exc:
                logger.warning("Grouped Scrublet failed for %s, falling back to singlets: %s", group, exc)
        ctx.metadata["doublet_grouped_key"] = sample_key
        ctx.metadata["doublet_grouped_thresholds"] = thresholds
        ctx.metadata["doublet_method"] = "scrublet_grouped"
        ctx.metadata["doublet_rate_scaling_enabled"] = bool(
            getattr(cfg, "scale_expected_doublet_rate", False)
        )
        ctx.metadata["doublet_expected_rate_applied_per_sample"] = per_sample_rates
        return scores, predicted, None

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _col_to_numpy(col) -> np.ndarray:
        """Convert pandas or cudf Series to a numpy array."""
        if hasattr(col, "to_numpy"):
            return col.to_numpy()
        return np.asarray(col)

    # Opt-in env var for failure-based all-singlets fallback. The
    # tiny-dataset path remains exempt because n_obs<20 / n_vars<50 is a
    # legitimate degenerate input where no doublet algorithm can produce
    # meaningful calls. The two FAILURE-based paths (rsc+cpu scrublet both
    # raise, or stand-alone scrublet raises) are silent algorithmic
    # compromises (Principle 9) and now require explicit acknowledgement.
    _NULL_FALLBACK_ENV = "SC_ALLOW_DOUBLET_NULL_FALLBACK"

    @staticmethod
    def _fallback_all_singlets(n_obs: int, reason: str) -> tuple[np.ndarray, np.ndarray]:
        scores = np.zeros(n_obs, dtype=np.float32)
        predicted = np.zeros(n_obs, dtype=bool)
        logger.info("Doublet detection fallback activated: %s", reason)
        return scores, predicted

    @classmethod
    def _require_null_fallback_opt_in(cls, primary_exc: Exception, secondary_exc: Exception | None = None) -> None:
        """Raise unless the operator has opted in to the all-singlets fallback on real failures.

        Tiny-dataset degeneracy is exempt and handled at its own call site.
        """
        if os.environ.get(cls._NULL_FALLBACK_ENV, "").strip() == "1":
            return
        detail = f"primary={primary_exc!r}"
        if secondary_exc is not None:
            detail += f"; secondary={secondary_exc!r}"
        raise RuntimeError(
            "doublet_detection: all backend attempts failed and the silent "
            "all-singlets fallback is banned by default (Principle 9). "
            f"Errors: {detail}. Set {cls._NULL_FALLBACK_ENV}=1 to acknowledge "
            "that the run will record zero doublets and proceed; the run "
            "manifest will mark this as fallback_all_singlets_opt_in."
        ) from primary_exc

    @staticmethod
    def _write_10x_mtx_triplet(adata, mtx_dir: Path) -> None:
        """Write AnnData counts as a 10X-style triplet for R subprocesses."""
        from scipy import io as _sio

        mtx_dir.mkdir(parents=True, exist_ok=True)
        X = adata.X
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        # mtx is genes x cells for Seurat::ReadMtx.
        _sio.mmwrite(str(mtx_dir / "matrix.mtx"), X.T.tocoo())
        with open(mtx_dir / "features.tsv", "w") as f:
            for g in adata.var_names:
                f.write(f"{g}\t{g}\tGene Expression\n")
        with open(mtx_dir / "barcodes.tsv", "w") as f:
            for b in adata.obs_names:
                f.write(f"{b}\n")

    # ------------------------------------------------------------------
    # DoubletFinder backend (R subprocess; sources local clone, no install)
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_doubletfinder_driver() -> Path:
        """Return the path to doubletfinder_run.R in r_multiomics_factory/scripts.

        Resolution priority mirrors ambient_correction.py:
          1. R_MULTIOMICS_SCRIPT_DIR env (explicit override)
          2. Repo-relative path from this file
        """
        env_dir = os.environ.get("R_MULTIOMICS_SCRIPT_DIR")
        if env_dir:
            driver = Path(env_dir) / "doubletfinder_run.R"
            if driver.exists():
                return driver
            raise FileNotFoundError(
                f"R_MULTIOMICS_SCRIPT_DIR={env_dir} does not contain doubletfinder_run.R"
            )
        # Repo-relative: this file lives at
        #   singlecell_factory/workflow/modular/modules/doublet_detection.py
        # R driver lives at
        #   r_multiomics_factory/scripts/doubletfinder_run.R
        here = Path(__file__).resolve()
        for parent in here.parents:
            if (parent / "singlecell_factory").exists() and (parent / "r_multiomics_factory").exists():
                driver = parent / "r_multiomics_factory" / "scripts" / "doubletfinder_run.R"
                if driver.exists():
                    return driver
        raise FileNotFoundError(
            "doubletfinder_run.R not located; set R_MULTIOMICS_SCRIPT_DIR or run "
            "from a checkout where singlecell_factory and r_multiomics_factory "
            "are siblings."
        )

    @staticmethod
    def _resolve_scdblfinder_driver() -> Path:
        """Return the path to scdblfinder_run.R in r_multiomics_factory/scripts."""
        env_dir = os.environ.get("R_MULTIOMICS_SCRIPT_DIR")
        if env_dir:
            driver = Path(env_dir) / "scdblfinder_run.R"
            if driver.exists():
                return driver
            raise FileNotFoundError(
                f"R_MULTIOMICS_SCRIPT_DIR={env_dir} does not contain scdblfinder_run.R"
            )
        here = Path(__file__).resolve()
        for parent in here.parents:
            if (parent / "singlecell_factory").exists() and (parent / "r_multiomics_factory").exists():
                driver = parent / "r_multiomics_factory" / "scripts" / "scdblfinder_run.R"
                if driver.exists():
                    return driver
        raise FileNotFoundError(
            "scdblfinder_run.R not located; set R_MULTIOMICS_SCRIPT_DIR or run "
            "from a checkout where singlecell_factory and r_multiomics_factory "
            "are siblings."
        )

    def _run_doubletfinder_via_r(
        self,
        adata,
        cfg,
        ctx,
    ) -> tuple[np.ndarray, np.ndarray, float | None]:
        """Run DoubletFinder via the r_multiomics conda env subprocess.

        Returns (pANN_scores, calls_bool, threshold). Threshold is None because
        DoubletFinder uses top-nExp ranking, not a score threshold.
        """
        expected_rate = self._record_backend_expected_rate(ctx, "doubletfinder")
        if not shutil.which("conda"):
            raise RuntimeError(
                "DoubletFinder backend requires conda on PATH to launch the "
                "r_multiomics env subprocess."
            )
        driver = self._resolve_doubletfinder_driver()

        pn = float(getattr(cfg, "doubletfinder_pn", 0.25))
        pk = float(getattr(cfg, "doubletfinder_pk", 0.09))
        pcs = int(getattr(cfg, "doubletfinder_pcs", 20))

        tmp = tempfile.mkdtemp(prefix="doubletfinder_")
        mtx_dir = Path(tmp) / "mtx"
        mtx_dir.mkdir(parents=True, exist_ok=True)
        out_csv = Path(tmp) / "out.csv"

        self._write_10x_mtx_triplet(adata, mtx_dir)

        cmd = [
            "conda", "run", "-n", str(getattr(cfg, "r_conda_env", "r_multiomics")),
            "Rscript", str(driver),
            "--mtx-dir", str(mtx_dir),
            "--out-csv", str(out_csv),
            "--expected-rate", f"{expected_rate}",
            "--pn", f"{pn}",
            "--pk", f"{pk}",
            "--pcs", f"{pcs}",
            "--seed", f"{ctx.random_state}",
        ]
        logger.info("DoubletFinder subprocess: %s", " ".join(cmd))

        timeout = float(getattr(cfg, "subprocess_timeout", 1800.0))
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
        try:
            try:
                stdout, stderr = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(os.getpgid(proc.pid), 9)
                stdout, stderr = proc.communicate()
                raise RuntimeError(
                    f"DoubletFinder subprocess timed out after {timeout}s"
                )
        except Exception:
            # Keep temp dir for debugging on failure.
            logger.error("DoubletFinder failed; preserving %s for debug", tmp)
            raise

        if stderr:
            for line in stderr.decode("utf-8", errors="replace").splitlines():
                if line.strip():
                    logger.warning("[DoubletFinder R] %s", line)

        if proc.returncode != 0:
            logger.error("DoubletFinder R driver exited with code %d; tmp dir preserved: %s",
                         proc.returncode, tmp)
            raise RuntimeError(
                f"DoubletFinder R driver failed (exit {proc.returncode}); see "
                f"stderr above and {tmp} for reproduction inputs."
            )

        if not out_csv.exists():
            raise RuntimeError(f"DoubletFinder did not produce {out_csv}")

        df_out = pd.read_csv(out_csv).set_index("barcode")
        # Align back to AnnData order.
        try:
            df_out = df_out.loc[list(adata.obs_names)]
        except KeyError as exc:
            raise RuntimeError(
                f"DoubletFinder output barcodes do not cover AnnData obs_names: {exc}"
            )

        scores = df_out["df_pANN"].to_numpy(dtype=np.float32)
        predicted = df_out["df_call"].to_numpy().astype(bool)

        # Clean temp on success.
        shutil.rmtree(tmp, ignore_errors=True)

        ctx.metadata["doublet_method"] = "doubletfinder"
        ctx.metadata["doublet_doubletfinder_params"] = {
            "expected_rate": expected_rate, "pn": pn, "pk": pk, "pcs": pcs,
        }
        return scores, predicted, None

    def _run_scdblfinder_via_r(
        self,
        adata,
        cfg,
        ctx,
    ) -> tuple[np.ndarray, np.ndarray, float | None]:
        """Run scDblFinder via the r_multiomics conda env subprocess."""
        expected_rate = self._record_backend_expected_rate(ctx, "scdblfinder")
        if not shutil.which("conda"):
            raise RuntimeError(
                "scDblFinder backend requires conda on PATH to launch the "
                "r_multiomics env subprocess."
            )
        driver = self._resolve_scdblfinder_driver()

        samples_col = getattr(cfg, "scdblfinder_samples_col", None)

        tmp = tempfile.mkdtemp(prefix="scdblfinder_")
        mtx_dir = Path(tmp) / "mtx"
        out_csv = Path(tmp) / "out.csv"
        self._write_10x_mtx_triplet(adata, mtx_dir)

        cmd = [
            "conda", "run", "-n", str(getattr(cfg, "r_conda_env", "r_multiomics")),
            "Rscript", str(driver),
            "--mtx-dir", str(mtx_dir),
            "--out-csv", str(out_csv),
            "--expected-rate", f"{expected_rate}",
            "--seed", f"{ctx.random_state}",
        ]
        samples_tsv = None
        if samples_col and samples_col in adata.obs.columns:
            samples_tsv = Path(tmp) / "samples.tsv"
            pd.DataFrame(
                {
                    "barcode": adata.obs_names.astype(str),
                    "sample": adata.obs[samples_col].astype(str).to_numpy(),
                }
            ).to_csv(samples_tsv, sep="\t", index=False)
            cmd += ["--samples-tsv", str(samples_tsv)]
        elif samples_col:
            msg = f"requested scDblFinder samples column not found: {samples_col}"
            logger.warning(msg)
            ctx.metadata["doublet_scdblfinder_samples_warning"] = msg

        logger.info("scDblFinder subprocess: %s", " ".join(cmd))

        timeout = float(getattr(cfg, "subprocess_timeout", 1800.0))
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
        try:
            try:
                stdout, stderr = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(os.getpgid(proc.pid), 9)
                stdout, stderr = proc.communicate()
                raise RuntimeError(
                    f"scDblFinder subprocess timed out after {timeout}s"
                )
        except Exception:
            logger.error("scDblFinder failed; preserving %s for debug", tmp)
            raise

        if stderr:
            for line in stderr.decode("utf-8", errors="replace").splitlines():
                if line.strip():
                    logger.warning("[scDblFinder R] %s", line)

        if proc.returncode != 0:
            logger.error("scDblFinder R driver exited with code %d; tmp dir preserved: %s",
                         proc.returncode, tmp)
            raise RuntimeError(
                f"scDblFinder R driver failed (exit {proc.returncode}); see "
                f"stderr above and {tmp} for reproduction inputs."
            )

        if not out_csv.exists():
            raise RuntimeError(f"scDblFinder did not produce {out_csv}")

        df_out = pd.read_csv(out_csv).set_index("barcode")
        try:
            df_out = df_out.loc[list(adata.obs_names)]
        except KeyError as exc:
            raise RuntimeError(
                f"scDblFinder output barcodes do not cover AnnData obs_names: {exc}"
            )

        scores = df_out["scdbl_score"].to_numpy(dtype=np.float32)
        predicted = df_out["scdbl_call"].to_numpy().astype(bool)
        shutil.rmtree(tmp, ignore_errors=True)

        ctx.metadata["doublet_method"] = "scdblfinder"
        ctx.metadata["doublet_scdblfinder_params"] = {
            "expected_rate": expected_rate,
            "samples_col": samples_col if samples_tsv is not None else None,
        }
        return scores, predicted, None

    # ------------------------------------------------------------------
    # Consensus mode: require >=2 backends to agree
    # ------------------------------------------------------------------

    def _run_consensus(self, adata, cfg, ctx) -> tuple[np.ndarray, np.ndarray, dict]:
        """Run configured doublet backends and merge calls via consensus_logic."""
        effective_rate = self._record_backend_expected_rate(ctx, "consensus")
        pair = (
            os.environ.get("SC_DOUBLET_CONSENSUS_PAIR", "").strip().lower()
            or getattr(cfg, "consensus_pair", "scrublet_doubletfinder")
        ).lower()
        aliases = {
            "scrublet_scdbl": "scrublet_scdblfinder",
            "scrublet_scdb": "scrublet_scdblfinder",
            "all3": "scrublet_doubletfinder_scdblfinder",
            "scrublet_all3": "scrublet_doubletfinder_scdblfinder",
        }
        pair = aliases.get(pair, pair)
        valid_pairs = {
            "scrublet_doubletfinder",
            "scrublet_scdblfinder",
            "scrublet_doubletfinder_scdblfinder",
        }
        if pair not in valid_pairs:
            logger.warning(
                "Unknown consensus_pair '%s'; defaulting to scrublet_doubletfinder.",
                pair,
            )
            pair = "scrublet_doubletfinder"

        scrub_scores, scrub_calls, scrub_threshold = self._run_scrublet_only(adata, cfg, ctx)
        per_backend = {
            "scrublet": {
                "scores": scrub_scores,
                "calls": scrub_calls,
                "threshold": scrub_threshold,
            }
        }

        if "doubletfinder" in pair:
            df_scores, df_calls, _ = self._run_doubletfinder_via_r(adata, cfg, ctx)
            per_backend["doubletfinder"] = {
                "scores": df_scores,
                "calls": df_calls,
            }
        if "scdblfinder" in pair:
            scdbl_scores, scdbl_calls, _ = self._run_scdblfinder_via_r(adata, cfg, ctx)
            per_backend["scdblfinder"] = {
                "scores": scdbl_scores,
                "calls": scdbl_calls,
            }

        ranks = [
            pd.Series(result["scores"]).rank(pct=True).to_numpy(dtype=np.float32)
            for result in per_backend.values()
        ]
        combined_score = np.mean(np.vstack(ranks), axis=0).astype(np.float32)

        logic = (
            os.environ.get("SC_DOUBLET_CONSENSUS_LOGIC", "").strip().lower()
            or getattr(cfg, "consensus_logic", "or")
        ).lower()
        if logic not in {"and", "or", "rank"}:
            logger.warning("Unknown consensus_logic '%s'; defaulting to 'or'.", logic)
            logic = "or"

        call_matrix = np.vstack([
            result["calls"].astype(bool) for result in per_backend.values()
        ])
        if logic == "and":
            combined_calls = np.all(call_matrix, axis=0)
        elif logic == "or":
            combined_calls = np.any(call_matrix, axis=0)
        else:
            n_expected = max(1, int(round(effective_rate * adata.n_obs)))
            order = np.argsort(-combined_score, kind="stable")
            combined_calls = np.zeros(adata.n_obs, dtype=bool)
            combined_calls[order[:n_expected]] = True

        backends = list(per_backend.keys())
        pairwise_agreement = {}
        for i, left in enumerate(backends):
            for right in backends[i + 1:]:
                pairwise_agreement[f"{left}|{right}"] = float(
                    (per_backend[left]["calls"] == per_backend[right]["calls"]).mean()
                )
        mean_agreement = (
            float(np.mean(list(pairwise_agreement.values())))
            if pairwise_agreement else 1.0
        )

        ctx.metadata["doublet_method"] = f"consensus_{logic}_{pair}"
        ctx.metadata["doublet_consensus_logic"] = logic
        ctx.metadata["doublet_consensus_pair"] = pair
        ctx.metadata["doublet_consensus_backends"] = ",".join(backends)
        ctx.metadata["doublet_consensus_agreement"] = mean_agreement
        ctx.metadata["doublet_consensus_pairwise_agreement"] = pairwise_agreement
        ctx.metadata["doublet_consensus_per_backend_rates"] = {
            **{name: float(result["calls"].mean()) for name, result in per_backend.items()},
            "final": float(combined_calls.mean()),
        }
        return combined_score, combined_calls.astype(bool), per_backend

    def _run_scrublet_only(self, adata, cfg, ctx) -> tuple[np.ndarray, np.ndarray, float | None]:
        """Helper: run scrublet via existing GPU/CPU paths and return results.

        Extracted so consensus mode can call it. Mirrors the logic inside
        `run()` for selecting whole vs grouped, GPU vs CPU.
        """
        strategy = self._resolve_doublet_strategy(ctx, adata)
        use_grouped = strategy == "grouped"
        random_state = ctx.random_state

        if _RSC_AVAILABLE and gpu_available(ctx.cfg.gpu_mode):
            try:
                if use_grouped:
                    return self._run_rsc_scrublet_grouped(adata, cfg, ctx)
                return self._run_rsc_scrublet_whole(adata, cfg, ctx)
            except Exception as exc:
                logger.warning(
                    "rsc.pp.scrublet failed inside consensus, retrying CPU: %s", exc
                )
        # CPU path.
        import scrublet as scr
        if use_grouped:
            return self._run_cpu_scrublet_grouped(adata, scr.Scrublet, cfg, ctx)
        expected_rate = self._record_backend_expected_rate(ctx, "scrublet")
        scrub = scr.Scrublet(
            self._materialize_counts_matrix(adata.X),
            expected_doublet_rate=expected_rate,
            random_state=random_state,
        )
        scores, predicted = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=self._n_prin_comps(adata.n_obs, adata.n_vars),
        )
        return scores.astype(np.float32), predicted.astype(bool), getattr(scrub, "threshold_", None)

    # ------------------------------------------------------------------
    # Under-call diagnostic
    # ------------------------------------------------------------------

    @staticmethod
    def _check_undercall(
        call_rate: float,
        expected_rate: float,
        ctx,
        backend: str,
        scope: str | None = None,
    ) -> None:
        """Emit WARNING when call rate < 0.5 * expected.

        A 0.13% call on a 5%-expected sample (40x undercall) is what surfaced
        on LUSC PS01 when Scrublet's bimodality threshold failed to find a
        clean separation. Surfacing this loudly avoids silent quality
        regressions on tumor cohorts.
        """
        if scope is None:
            ctx.metadata["doublet_undercall_expected_rate"] = round(
                float(expected_rate), 6
            )
        else:
            ctx.metadata.setdefault(
                "doublet_undercall_expected_rate_per_sample", {}
            )[scope] = round(float(expected_rate), 6)
        if expected_rate <= 0:
            return
        if call_rate < 0.5 * expected_rate:
            ratio = expected_rate / max(call_rate, 1e-6)
            msg = (
                f"{backend} call rate {call_rate*100:.2f}% is < 0.5x expected "
                f"({expected_rate*100:.2f}%); under-called by {ratio:.1f}x. "
                "Consider re-running with backend='consensus' or "
                "backend='scdblfinder'/'doubletfinder' for a second opinion."
            )
            logger.warning("DOUBLET_UNDERCALL_DIAGNOSTIC: %s", msg)
            if scope is None:
                ctx.metadata["doublet_undercall_warning"] = msg
                ctx.metadata["doublet_undercall_ratio"] = ratio
            else:
                ctx.metadata.setdefault("doublet_undercall_warning_per_sample", {})[
                    scope
                ] = msg
                ctx.metadata.setdefault("doublet_undercall_ratio_per_sample", {})[
                    scope
                ] = ratio

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Doublet detection requires loaded AnnData.")

        cfg = ctx.cfg.doublet
        random_state = ctx.random_state
        threshold = None
        self._resolve_expected_rate_contract(ctx, cfg, adata)

        # Backend selection (default = scrublet, preserves prior behavior).
        # cfg.doublet.backend overridden by SC_DOUBLET_BACKEND env var if set.
        requested_backend = (
            os.environ.get("SC_DOUBLET_BACKEND", "").strip().lower()
            or getattr(cfg, "backend", "scrublet")
        ).lower()
        backend = {"scdbl": "scdblfinder"}.get(requested_backend, requested_backend)
        fallback_reason = None
        if backend not in {"scrublet", "doubletfinder", "scdblfinder", "consensus"}:
            logger.warning(
                "Unknown doublet backend '%s'; falling back to 'scrublet'.", backend
            )
            fallback_reason = f"unknown_backend:{backend}"
            backend = "scrublet"
        ctx.metadata["doublet_backend_requested"] = requested_backend
        ctx.metadata["doublet_backend"] = backend
        if fallback_reason:
            ctx.metadata["doublet_backend_fallback_reason"] = fallback_reason

        _computed = False
        # Scrublet can fail on tiny/degenerate datasets; keep pipeline usable by
        # falling back to "all singlets" rather than aborting mandatory stage.
        # Tiny-dataset path never invokes a count-based caller, so it runs
        # outside the raw-counts .X swap below.
        if adata.n_obs < 20 or adata.n_vars < 50:
            doublet_scores, predicted_doublets = self._fallback_all_singlets(
                adata.n_obs,
                reason=f"dataset too small (n_obs={adata.n_obs}, n_vars={adata.n_vars})",
            )
            ctx.metadata["doublet_method"] = "fallback_all_singlets"
            _computed = True

        # Raw-counts .X swap (try/finally) around ALL count-based doublet calls.
        # ambient_correction may have overwritten adata.X with fractional
        # decontXcounts; Scrublet/DoubletFinder/scDblFinder all model raw
        # integer UMIs. Swap the resolved raw counts onto adata.X for the
        # compute block and ALWAYS restore the original .X afterward (before the
        # post-doublet row-subset). GPU paths use adata.copy(), so the swapped
        # .X propagates to device.
        _original_x = None
        _x_swapped = False
        if _computed:
            doublet_scores_local = doublet_scores
            predicted_local = predicted_doublets
            backend_local = backend
            threshold_local = threshold
        else:
            raw_counts, raw_source = self._resolve_raw_counts_matrix(adata)
            ctx.metadata["doublet_input_source"] = raw_source
            if raw_source != "adata.X_integer_like":
                _original_x = adata.X
                adata.X = raw_counts
                _x_swapped = True
            try:
                (
                    doublet_scores_local,
                    predicted_local,
                    threshold_local,
                    backend_local,
                ) = self._compute_doublets_with_raw_x(
                    ctx, adata, cfg, backend, random_state,
                    computed=_computed,
                    scores=None,
                    predicted=None,
                    threshold=threshold,
                )
            finally:
                if _x_swapped:
                    adata.X = _original_x

        self._finalize_doublets(
            ctx, adata, cfg, backend_local, random_state,
            doublet_scores_local, predicted_local, threshold_local,
        )
        return

    def _compute_doublets_with_raw_x(
        self, ctx, adata, cfg, backend, random_state, *,
        computed, scores, predicted, threshold,
    ) -> None:
        """Count-based doublet compute block, run with raw-integer .X in place.

        Split out of run() so the raw-counts .X swap in run() can wrap every
        Scrublet/DoubletFinder/scDblFinder/consensus path in a single
        try/finally, guaranteeing .X is restored before the row-subset.
        """
        _computed = computed
        doublet_scores = scores
        predicted_doublets = predicted
        if _computed:
            # Tiny-dataset all-singlets path already produced results in run().
            pass
        elif backend == "doubletfinder":
            try:
                doublet_scores, predicted_doublets, threshold = self._run_doubletfinder_via_r(
                    adata, cfg, ctx,
                )
                _computed = True
            except Exception as exc:
                logger.warning(
                    "DoubletFinder backend failed (%s); falling back to scrublet.", exc,
                )
                ctx.metadata["doublet_backend_fallback_reason"] = f"doubletfinder_failed:{exc}"
                backend = "scrublet"  # fall through to scrublet path
        elif backend == "scdblfinder":
            try:
                doublet_scores, predicted_doublets, threshold = self._run_scdblfinder_via_r(
                    adata, cfg, ctx,
                )
                _computed = True
            except Exception as exc:
                logger.warning(
                    "scDblFinder backend failed (%s); falling back to scrublet.", exc,
                )
                ctx.metadata["doublet_backend_fallback_reason"] = f"scdblfinder_failed:{exc}"
                backend = "scrublet"  # fall through to scrublet path
        elif backend == "consensus":
            try:
                doublet_scores, predicted_doublets, per_backend = self._run_consensus(
                    adata, cfg, ctx,
                )
                # Persist per-backend annotations on adata for downstream inspection.
                for name, result in per_backend.items():
                    adata.obs[f"{name}_score"] = result["scores"]
                    adata.obs[f"{name}_call"] = result["calls"].astype(bool)
                _computed = True
            except Exception as exc:
                logger.warning(
                    "Consensus backend failed (%s); falling back to scrublet.", exc,
                )
                ctx.metadata["doublet_backend_fallback_reason"] = f"consensus_failed:{exc}"
                backend = "scrublet"  # fall through

        if not _computed and backend == "scrublet":
            strategy = self._resolve_doublet_strategy(ctx, adata)
            use_grouped = strategy == "grouped"

            # `--gpu-mode off` is a user-facing control; gating GPU dispatch on
            # library availability alone made it unenforceable here, so a run
            # requested as CPU-only still initialised CUDA and ran GPU Scrublet.
            if _RSC_AVAILABLE and gpu_available(ctx.cfg.gpu_mode):
                # --- GPU path ---
                try:
                    if use_grouped:
                        doublet_scores, predicted_doublets, threshold = self._run_rsc_scrublet_grouped(adata, cfg, ctx)
                    else:
                        doublet_scores, predicted_doublets, threshold = self._run_rsc_scrublet_whole(adata, cfg, ctx)
                    if threshold is not None:
                        logger.info("rsc.pp.scrublet auto-threshold: %.4f", threshold)
                except Exception as exc:
                    # HIGH-1 fix (2026-05-19): retry CPU scrublet before
                    # giving up. The previous behaviour silently dropped
                    # doublet detection entirely whenever the RAPIDS path
                    # raised (e.g. rsc.pp.scrublet sparse-dtype mismatch
                    # against int counts), keeping ~2% of cells that should
                    # have been removed. CPU scrublet is the canonical
                    # reference implementation (Wolock et al. 2019, Cell
                    # Systems) so retrying it is scientifically equivalent
                    # to the GPU path's intent.
                    logger.warning("rsc.pp.scrublet failed, retrying with CPU scrublet: %s", exc)
                    try:
                        import scrublet as scr  # type: ignore

                        if use_grouped:
                            doublet_scores, predicted_doublets, threshold = self._run_cpu_scrublet_grouped(
                                adata, scr.Scrublet, cfg, ctx
                            )
                        else:
                            expected_rate = self._record_backend_expected_rate(
                                ctx, "scrublet"
                            )
                            scrub = scr.Scrublet(
                                self._materialize_counts_matrix(adata.X),
                                expected_doublet_rate=expected_rate,
                                random_state=random_state,
                            )
                            scrublet_params = {
                                "min_counts": 2,
                                "min_cells": 3,
                                "min_gene_variability_pctl": 85,
                                "n_prin_comps": self._n_prin_comps(adata.n_obs, adata.n_vars),
                            }
                            doublet_scores, predicted_doublets = scrub.scrub_doublets(**scrublet_params)
                            threshold = getattr(scrub, "threshold_", None)
                        ctx.metadata["doublet_method"] = "cpu_scrublet_fallback_from_rsc"
                        ctx.metadata["doublet_rsc_failure_reason"] = str(exc)
                    except Exception as cpu_exc:
                        logger.error(
                            "CPU scrublet retry also failed; gating on %s. rsc=%s cpu=%s",
                            self._NULL_FALLBACK_ENV, exc, cpu_exc,
                        )
                        self._require_null_fallback_opt_in(exc, cpu_exc)
                        doublet_scores, predicted_doublets = self._fallback_all_singlets(
                            adata.n_obs, reason=f"rsc={exc}; cpu_scrublet={cpu_exc}",
                        )
                        ctx.metadata["doublet_method"] = "fallback_all_singlets_opt_in"
                        ctx.metadata["doublet_method_actually_used"] = "fallback_all_singlets_opt_in"
                        ctx.metadata["doublet_null_fallback_opt_in_acknowledged"] = True
                        ctx.metadata["doublet_rsc_failure_reason"] = str(exc)
                        ctx.metadata["doublet_cpu_failure_reason"] = str(cpu_exc)
            else:
                # --- CPU fallback path (rsc not installed) ---
                import scrublet as scr

                if use_grouped:
                    doublet_scores, predicted_doublets, threshold = self._run_cpu_scrublet_grouped(adata, scr.Scrublet, cfg, ctx)
                else:
                    expected_rate = self._record_backend_expected_rate(
                        ctx, "scrublet"
                    )
                    scrub = scr.Scrublet(
                        self._materialize_counts_matrix(adata.X),
                        expected_doublet_rate=expected_rate,
                        random_state=random_state,
                    )
                    scrublet_params = {
                        "min_counts": 2,
                        "min_cells": 3,
                        "min_gene_variability_pctl": 85,
                        "n_prin_comps": self._n_prin_comps(adata.n_obs, adata.n_vars),
                    }
                    try:
                        doublet_scores, predicted_doublets = scrub.scrub_doublets(**scrublet_params)
                        threshold = getattr(scrub, "threshold_", None)
                        if threshold is not None:
                            logger.info("Scrublet auto-threshold: %.4f", threshold)
                        ctx.metadata["doublet_method"] = "scrublet"
                    except Exception as exc:
                        logger.error(
                            "Scrublet failed; gating on %s. error=%s",
                            self._NULL_FALLBACK_ENV, exc,
                        )
                        self._require_null_fallback_opt_in(exc)
                        doublet_scores, predicted_doublets = self._fallback_all_singlets(
                            adata.n_obs, reason=str(exc),
                        )
                        ctx.metadata["doublet_method"] = "fallback_all_singlets_opt_in"
                        ctx.metadata["doublet_method_actually_used"] = "fallback_all_singlets_opt_in"
                        ctx.metadata["doublet_null_fallback_opt_in_acknowledged"] = True
                        ctx.metadata["doublet_scrublet_failure_reason"] = str(exc)

        # Return raw results to run(); metadata assembly, visualization, and the
        # row-subset happen in run() AFTER the .X swap is restored.
        return doublet_scores, predicted_doublets, threshold, backend

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------

    # Colour roles for the doublet panel, taken from the theme's qualitative
    # ramp so the diagnostic and the manuscript figures speak one palette:
    # index 0 = the retained population, index 1 = the cutoff and the cells it
    # removes. The literal pair is only a fallback for a missing theme file.
    _PANEL_FALLBACK_COLOURS = ("#0072B2", "#D55E00")
    _PANEL_MUTED = "#666666"
    # A solid histogram leaves no clear space, so in-panel notes get a
    # translucent white plate rather than being drawn over the bars.
    _PANEL_NOTE_PLATE = {"facecolor": "white", "alpha": 0.85, "edgecolor": "none",
                         "boxstyle": "square,pad=0.25"}

    @staticmethod
    def _fmt_pct(pct: float) -> str:
        """Percentage with enough precision to stay truthful near zero.

        The interesting case is a pre-filtered atlas: 4 of 87,380 cells is
        0.0046%, which ``%.2f`` renders as "0.00%" — i.e. as exactly none, which
        is a different claim from "a handful".
        """
        if pct <= 0:
            return "0%"
        return f"{pct:.2g}%" if pct < 0.01 else f"{pct:.2f}%"

    @classmethod
    def _plot_doublet_scores(
        cls, ctx, adata, scores, predicted, threshold, *,
        n_doublets: int, call_rate: float, backend: str,
    ) -> None:
        """Draw the score distribution against the prior the caller applied.

        rcParams cannot fix this panel. It used to be an 8x4 inch canvas (no
        journal column fits it) holding a 50-bin linear histogram titled with
        the backend name, which on a real cohort is unreadable in three
        independent ways. The score axis spans three decades (LUSC atlas:
        median 0.010, max 0.58), so every singlet collapses into the first bin.
        The count axis is linear, so the handful of cells the threshold
        actually removes is a flat line at zero. And the panel never states
        what the caller expected, so a 529x under-call on a pre-filtered atlas
        looks exactly like a clean run with a genuinely low doublet rate.

        So the panel reports the decision rather than only the scores: the
        distribution on log axes with the applied cutoff drawn (a), and the
        resulting call rate beside the prior it must be judged against (b).
        Every number is read from ``ctx.metadata`` or from arrays already
        computed upstream; nothing here is recomputed.
        """
        from .._figure_theme import journal_figure_size, qualitative_colors

        palette = qualitative_colors("group_qualitative", 2)
        keep_colour, call_colour = (
            palette if len(palette) == 2 else list(cls._PANEL_FALLBACK_COLOURS)
        )

        scores = np.asarray(scores, dtype=float)
        called = np.asarray(predicted, dtype=bool)
        finite = np.isfinite(scores)
        method = str(ctx.metadata.get("doublet_method", backend))
        expected_rate = ctx.metadata.get("doublet_expected_rate_effective")
        observed_pct = float(call_rate) * 100.0

        values = scores[finite]
        # A constant score has no distribution to draw, so panel (a) degrades to a
        # one-line statement. Shrink the canvas with it instead of shipping a
        # journal-width figure that is mostly empty.
        degenerate = bool(values.size == 0 or np.nanmax(values) == np.nanmin(values))
        fig, axes = plt.subplots(
            1, 2,
            figsize=journal_figure_size("double", height_mm=42.0 if degenerate else 62.0),
            gridspec_kw={"width_ratios": [1.4, 1.0] if degenerate else [3.0, 1.0]},
            constrained_layout=True,
        )
        ax, ax_rate = axes[0], axes[1]
        for idx, panel_ax in enumerate((ax, ax_rate)):
            panel_ax.text(-0.02, 1.12, chr(ord("a") + idx), transform=panel_ax.transAxes,
                          ha="right", va="top", fontweight="bold", fontsize=8)

        # --- (a) score distribution with the applied cutoff ---------------
        drawable_threshold = False
        if degenerate:
            # A constant score says nothing about doublets. Declare that, rather
            # than drawing one bar on an autoscaled axis that implies a model was
            # fitted — the all-singlets fallback lands exactly here.
            constant = float(values[0]) if values.size else float("nan")
            ax.text(0.5, 0.5,
                    f"no score variation\nevery cell scored {constant:g}\n"
                    f"no doublet model was fitted ({method})",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=7, color=cls._PANEL_MUTED)
            ax.set_axis_off()
        else:
            # Bin geometrically when the scores span decades, so the singlet mode
            # is resolved instead of being packed into the first linear bin. Log
            # bins need strictly positive values, and a zero score is a real
            # backend output, so fall back to linear rather than dropping cells.
            use_log_x = bool(values.min() > 0 and values.max() / values.min() >= 100.0)
            edges = (np.geomspace(values.min(), values.max(), 51) if use_log_x
                     else np.linspace(values.min(), values.max(), 51))
            keep_mask = finite & ~called
            ax.hist(scores[keep_mask], bins=edges, color=keep_colour, linewidth=0,
                    label=f"Retained (n = {int(keep_mask.sum()):,})")
            if n_doublets:
                # Called cells sit inside the same tail as the retained ones;
                # without this overlay they cannot be told apart.
                ax.hist(scores[finite & called], bins=edges, color=call_colour,
                        linewidth=0, zorder=3,
                        label=f"Called doublet (n = {n_doublets:,})")
            if use_log_x:
                ax.set_xscale("log")
            # Bin counts span ~4 decades between the singlet mode and the tail the
            # threshold acts on. On a linear count axis that tail — i.e. every cell
            # this module removes — is indistinguishable from zero.
            ax.set_yscale("log")
            ax.set_ylim(bottom=0.5)
            ax.set_xlabel("Doublet score" + (" (log scale)" if use_log_x else ""))
            ax.set_ylabel("Cells per bin (log scale)")

            drawable_threshold = (
                threshold is not None
                and np.isfinite(float(threshold))
                and (float(threshold) > 0 or not use_log_x)
            )
            notes: list[tuple[str, str]] = []
            if drawable_threshold:
                ax.axvline(float(threshold), color=call_colour, linewidth=0.6,
                           linestyle="--", label=f"Threshold = {float(threshold):.3f}")
                ax.set_title("Score distribution and applied cutoff")
            else:
                # Rank-based backends (DoubletFinder top-nExp) and the per-sample
                # grouped lanes expose no single global cutoff. Say so in the title
                # too, so a missing line is not read as one that was forgotten.
                ax.set_title("Score distribution (no global cutoff)")
                notes.append(("cutoff is rank- or per-sample-calibrated",
                              cls._PANEL_MUTED))
            if not n_doublets:
                # The pre-filtered-atlas case. Without this the panel is just a
                # tidy distribution and the reader has to notice the absent
                # overlay to learn that the module removed nothing.
                notes.append(("no cell was called a doublet", call_colour))
            handles, _ = ax.get_legend_handles_labels()
            if len(handles) > 1:
                # Same plate as the notes: the theme's borderless legend is the
                # right publication default, but the threshold rule can land on
                # top of the labels, so give it an opaque (still borderless) bed.
                ax.legend(loc="upper right", fontsize=5.5, handlelength=1.2,
                          frameon=True, framealpha=0.85, edgecolor="none",
                          borderpad=0.3)
            for row, (text, colour) in enumerate(notes):
                ax.text(0.02, 0.97 - 0.08 * row, text, transform=ax.transAxes,
                        ha="left", va="top", fontsize=5.5, color=colour,
                        bbox=cls._PANEL_NOTE_PLATE)

        # --- (b) realised call rate against the prior that was applied ----
        ax_rate.set_title("Call rate vs prior")
        if expected_rate is None:
            ax_rate.text(0.5, 0.5, "expected rate\nnot recorded", ha="center", va="center",
                         transform=ax_rate.transAxes, fontsize=6, color=cls._PANEL_MUTED)
            ax_rate.set_axis_off()
        else:
            expected_pct = float(expected_rate) * 100.0
            ax_rate.bar([0], [observed_pct], width=0.6, color=call_colour)
            # The prior is an assumption, not a measurement: drawn open so it can
            # never be read as a second observed quantity.
            ax_rate.bar([1], [expected_pct], width=0.6, facecolor="none",
                        edgecolor=cls._PANEL_MUTED, linewidth=0.6, linestyle="--")
            ax_rate.set_xticks([0, 1])
            ax_rate.set_xticklabels(["called", "expected"])
            ax_rate.set_ylabel("% of cells")
            headroom = max(observed_pct, expected_pct)
            ax_rate.set_ylim(0, headroom * 1.35 if headroom > 0 else 1.0)
            for x, pct in ((0, observed_pct), (1, expected_pct)):
                # Both bars carry their value because the interesting case is the
                # one where the called bar is ~500x shorter than the prior and
                # would otherwise render as an unreadable zero-height stub.
                ax_rate.text(x, pct + max(headroom, 1e-9) * 0.04, cls._fmt_pct(pct),
                             ha="center", va="bottom", fontsize=5.5)
            note = cls._undercall_note(ctx)
            if note:
                ax_rate.text(0.03, 0.97, note, transform=ax_rate.transAxes,
                             ha="left", va="top", fontsize=5.5, color=call_colour)

        # Second suptitle line defines only the encodings this run actually drew,
        # so it never explains a threshold line or a prior bar that is absent.
        key = []
        if drawable_threshold:
            key.append("a: dashed line = applied threshold")
        if expected_rate is not None:
            key.append("b: open bar = expected rate (a prior, not a measurement)")
        fig.suptitle(
            f"Doublet detection — {method}; n = {adata.n_obs:,} cells, "
            f"{n_doublets:,} called ({cls._fmt_pct(observed_pct)})"
            + ("\n" + ".   ".join(key) if key else ""),
            fontsize=7,
        )
        fig.savefig(ctx.figure_dir / "doublet_scores.png")
        plt.close(fig)

    @staticmethod
    def _undercall_note(ctx) -> str | None:
        """One-line restatement of the under-call diagnostic already in metadata.

        The diagnostic is the reason this panel exists on tumour cohorts, so it
        belongs on the figure and not only in the manifest. Reads whichever scope
        ``_check_undercall`` populated; computes nothing.
        """
        ratio = ctx.metadata.get("doublet_undercall_ratio")
        if ratio:
            return f"under-called {float(ratio):.0f}-fold vs prior"
        per_sample = ctx.metadata.get("doublet_undercall_ratio_per_sample") or {}
        if per_sample:
            resolved = ctx.metadata.get("doublet_expected_rate_resolved_per_sample") or {}
            total = len(resolved) or len(per_sample)
            return f"{len(per_sample)}/{total} samples under-called"
        return None

    def _finalize_doublets(
        self, ctx, adata, cfg, backend, random_state,
        doublet_scores, predicted_doublets, threshold,
    ) -> None:
        """Record metadata, plot, and apply the row-subset (raw .X restored)."""
        ctx.metadata["doublet_backend"] = backend
        ctx.metadata.setdefault("doublet_method", backend)
        # Mirror doublet_method into a manifest-stable field so every run
        # records which engine actually produced the predicted_doublet column,
        # including the degenerate-input tiny-dataset case and the opt-in
        # null fallback. This lets reviewers gate publication claims without
        # parsing per-path metadata keys.
        ctx.metadata.setdefault("doublet_method_actually_used", ctx.metadata["doublet_method"])
        ctx.metadata["doublet_random_state"] = random_state

        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets

        n_doublets = int(predicted_doublets.sum())
        call_rate = n_doublets / adata.n_obs if adata.n_obs > 0 else 0.0
        ctx.metadata["doublets_detected"] = n_doublets
        ctx.metadata["doublet_rate_pct"] = round(call_rate * 100, 2)

        # Diagnostic: warn loudly when call_rate is < 0.5x expected. This
        # surfaced as a real failure on LUSC PS01 where Scrublet's
        # bimodality threshold collapsed at 0.13% on a slice that should
        # have yielded ~6% doublets.
        backend_prior = ctx.metadata.get("doublet_expected_rate_by_backend", {}).get(
            backend, {}
        )
        if backend_prior.get("mode") == "per_sample" and "sample" in adata.obs:
            sample_labels = adata.obs["sample"].astype(str)
            for sample, expected_rate in backend_prior["resolved_per_sample"].items():
                sample_mask = sample_labels == sample
                sample_call_rate = float(predicted_doublets[sample_mask].mean())
                self._check_undercall(
                    call_rate=sample_call_rate,
                    expected_rate=float(expected_rate),
                    ctx=ctx,
                    backend=backend,
                    scope=sample,
                )
        else:
            self._check_undercall(
                call_rate=call_rate,
                expected_rate=float(ctx.metadata["doublet_expected_rate_effective"]),
                ctx=ctx,
                backend=backend,
            )

        # Visualize doublet score distribution
        self._plot_doublet_scores(
            ctx, adata, doublet_scores, predicted_doublets, threshold,
            n_doublets=n_doublets, call_rate=call_rate, backend=backend,
        )

        if cfg.remove_doublets:
            before = adata.n_obs
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            ctx.adata = adata
            ctx.metadata["cells_after_doublet_removal"] = int(adata.n_obs)
            ctx.metadata["doublets_removed"] = int(before - adata.n_obs)
