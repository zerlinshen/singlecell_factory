"""Multimodal joint-embedding integration module (EXPERIMENTAL).

Third consumer of the bundle v2.1 extension API
(``scripts.export_singlecell_r_bundle.add_extension``); registers under the
``multimodal_obsm`` extension key.

EXPERIMENTAL: this module is OFF by default. Embeddings produced here are
visualization aids -- the cross-modality validation logic that would let a
user make new quantitative claims is not yet wired up. The exported
``multimodal_obsm`` extension carries the same plotting-only ``CLAIM_GUARD``
as the rest of the bundle.

Operating modes (selected by ``MultimodalIntegrationConfig.engine`` or the
``SC_MULTIMODAL_ENGINE`` env var):

- ``off`` (default): no-op; module records ``adata.uns["multimodal_status"]
  = {"engine": "off", "status": "skipped"}`` and exits.
- ``wnn``: writes the two latent obsm slices to a temp parquet pair, calls
  ``scripts/run_wnn.R`` via Rscript subprocess, reads the WNN UMAP back as
  ``adata.obsm["X_wnn"]``. Skips cleanly (no exception) if Rscript is
  missing or the R driver returns a non-zero exit code; mirrors the
  squidpy-skip pattern used by ``spatial_neighborhoods``.
- ``mofa``: imports ``mofapy2``/``muon`` and runs MOFA over the two latent
  modalities, writing ``adata.obsm["X_mofa"]``. Skips cleanly when neither
  package is importable.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "hao_seurat_wnn_2021": {
        "title": "Integrated analysis of multimodal single-cell data",
        "authors": "Hao et al.",
        "journal": "Cell",
        "year": "2021",
        "doi": "10.1016/j.cell.2021.04.048",
        "description": "Seurat WNN: weighted nearest neighbours for multimodal joint embeddings.",
    },
    "argelaguet_mofa_2020": {
        "title": "MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data",
        "authors": "Argelaguet et al.",
        "journal": "Genome Biology",
        "year": "2020",
        "doi": "10.1186/s13059-020-02015-1",
        "description": "MOFA factor model for multi-omics integration.",
    },
}


# Default driver-script path on disk. Mirrors the protein/spatial
# path-resolution chain: env override -> repo-relative canonical fallback.
_DEFAULT_RUN_WNN_R = (
    "/home/zerlinshen/multiomics_r_factory/scripts/run_wnn.R"
)


@dataclass(frozen=True)
class MultimodalIntegrationConfig:
    """Local config for the multimodal integration module (EXPERIMENTAL).

    Off by default. To enable, set ``engine`` to ``"wnn"`` or ``"mofa"`` (or
    set the ``SC_MULTIMODAL_ENGINE`` env var to one of those literals).
    """

    engine: str = "off"  # "off" | "wnn" | "mofa"
    rna_obsm_key: str = "X_pca"
    second_obsm_key: str = "protein_clr"
    wnn_obsm_output_key: str = "X_wnn"
    mofa_obsm_output_key: str = "X_mofa"
    n_factors: int = 10  # MOFA only
    rscript_bin: str | None = None  # falls back to PATH lookup / RSCRIPT_BIN
    run_wnn_r_path: str | None = None  # override path to scripts/run_wnn.R
    subprocess_timeout: int = 600  # seconds


class MultimodalIntegrationModule:
    """EXPERIMENTAL: multimodal joint embedding (WNN via R, MOFA via muon).

    Default engine is ``"off"`` -- the module is registered in the catalog
    but never runs unless explicitly enabled. Enabling it requires both an
    RNA latent (``adata.obsm["X_pca"]`` by default) AND a second-modality
    latent (``adata.obsm["protein_clr"]`` by default) to be present.
    """

    name = "multimodal_integration"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "obsm": ["X_wnn", "X_mofa"],
        "uns": ["multimodal_status", "mofa_factors"],
    }

    def __init__(self, config: MultimodalIntegrationConfig | None = None) -> None:
        self.config = config or MultimodalIntegrationConfig()

    # ------------------------------------------------------------------
    # Public entrypoint
    # ------------------------------------------------------------------

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        cfg = self._merge_runtime_config(ctx)
        engine = self._resolve_engine(cfg)

        if engine == "off":
            logger.info("%s: engine=off; module is a no-op (EXPERIMENTAL).", self.name)
            adata.uns["multimodal_status"] = {
                "engine": "off",
                "status": "skipped",
                "reason": "engine=off (default)",
            }
            ctx.metadata["multimodal_status"] = "off"
            ctx.status(self.name, "skipped", "engine=off (default)")
            return

        # Pre-condition: both latent obsm keys must be present.
        rna_key = cfg.rna_obsm_key
        sec_key = cfg.second_obsm_key
        if rna_key not in adata.obsm or sec_key not in adata.obsm:
            reason = (
                f"missing latents (rna='{rna_key}' present={rna_key in adata.obsm}, "
                f"second='{sec_key}' present={sec_key in adata.obsm})"
            )
            logger.info("%s: %s; skipping engine=%s.", self.name, reason, engine)
            adata.uns["multimodal_status"] = {
                "engine": engine,
                "status": "skipped",
                "reason": reason,
            }
            ctx.metadata["multimodal_status"] = "skipped_no_inputs"
            ctx.status(self.name, "skipped", reason)
            return

        if engine == "wnn":
            self._run_wnn(ctx, cfg, rna_key, sec_key)
        elif engine == "mofa":
            self._run_mofa(ctx, cfg, rna_key, sec_key)
        else:
            raise ValueError(
                f"{self.name}: unsupported engine '{engine}'; expected 'off' | 'wnn' | 'mofa'."
            )

    # ------------------------------------------------------------------
    # Engines
    # ------------------------------------------------------------------

    def _run_wnn(
        self,
        ctx: PipelineContext,
        cfg: MultimodalIntegrationConfig,
        rna_key: str,
        sec_key: str,
    ) -> None:
        adata = ctx.adata

        rscript = self._resolve_rscript(cfg)
        if rscript is None:
            reason = "Rscript not found on PATH and RSCRIPT_BIN not set"
            logger.info("%s: %s; WNN skipped.", self.name, reason)
            adata.uns["multimodal_status"] = {
                "engine": "wnn", "status": "skipped", "reason": reason,
            }
            ctx.metadata["multimodal_status"] = "skipped_no_rscript"
            ctx.status(self.name, "skipped", reason)
            return

        run_wnn_r = self._resolve_run_wnn_r_path(cfg)
        if not run_wnn_r or not Path(run_wnn_r).exists():
            reason = f"run_wnn.R driver not found at {run_wnn_r}"
            logger.info("%s: %s; WNN skipped.", self.name, reason)
            adata.uns["multimodal_status"] = {
                "engine": "wnn", "status": "skipped", "reason": reason,
            }
            ctx.metadata["multimodal_status"] = "skipped_no_driver"
            ctx.status(self.name, "skipped", reason)
            return

        cells = adata.obs_names.astype(str)
        rna_arr = np.asarray(adata.obsm[rna_key])
        sec_arr = np.asarray(adata.obsm[sec_key])

        with tempfile.TemporaryDirectory(prefix="wnn_") as tmpdir:
            tmp = Path(tmpdir)
            rna_pq = tmp / "rna_latent.parquet"
            sec_pq = tmp / "second_latent.parquet"
            out_pq = tmp / "wnn_umap.parquet"

            self._write_latent_parquet(rna_arr, cells, rna_pq, prefix="RNA_")
            self._write_latent_parquet(sec_arr, cells, sec_pq, prefix="MOD2_")

            cmd = [
                rscript, str(run_wnn_r),
                "--rna-parquet", str(rna_pq),
                "--protein-parquet", str(sec_pq),
                "--out-parquet", str(out_pq),
            ]
            logger.info("%s: invoking WNN driver: %s", self.name, " ".join(cmd))
            try:
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=cfg.subprocess_timeout,
                )
            except FileNotFoundError as exc:
                reason = f"Rscript invocation failed: {exc}"
                logger.info("%s: %s; WNN skipped.", self.name, reason)
                adata.uns["multimodal_status"] = {
                    "engine": "wnn", "status": "skipped", "reason": reason,
                }
                ctx.metadata["multimodal_status"] = "skipped_subprocess"
                ctx.status(self.name, "skipped", reason)
                return
            except subprocess.TimeoutExpired:
                reason = f"WNN driver exceeded timeout={cfg.subprocess_timeout}s"
                logger.warning("%s: %s; WNN skipped.", self.name, reason)
                adata.uns["multimodal_status"] = {
                    "engine": "wnn", "status": "skipped", "reason": reason,
                }
                ctx.metadata["multimodal_status"] = "skipped_timeout"
                ctx.status(self.name, "skipped", reason)
                return

            if proc.returncode != 0 or not out_pq.exists():
                reason = (
                    f"WNN driver exit_code={proc.returncode}; out_exists={out_pq.exists()}; "
                    f"stderr={proc.stderr.strip()[:400]}"
                )
                logger.info("%s: %s; WNN skipped.", self.name, reason)
                adata.uns["multimodal_status"] = {
                    "engine": "wnn", "status": "skipped", "reason": reason,
                }
                ctx.metadata["multimodal_status"] = "skipped_driver_failed"
                ctx.status(self.name, "skipped", reason)
                return

            df = pd.read_parquet(out_pq)

        # Align WNN UMAP rows to adata.obs_names. The driver writes columns
        # cell, WNN_1, WNN_2.
        if "cell" not in df.columns:
            raise ValueError(
                f"{self.name}: WNN driver output missing 'cell' column; got {list(df.columns)}"
            )
        coord_cols = [c for c in df.columns if c != "cell"]
        if len(coord_cols) < 2:
            raise ValueError(
                f"{self.name}: WNN driver output must have >= 2 coord columns; got {coord_cols}"
            )

        df_indexed = df.set_index("cell").reindex(cells)
        # densify-allowed: WNN UMAP is O(n_cells x 2), small dense output from Seurat
        coords = df_indexed[coord_cols[:2]].to_numpy(dtype=np.float32, copy=True)
        if np.isnan(coords).any():
            n_missing = int(np.isnan(coords).any(axis=1).sum())
            logger.warning(
                "%s: %d cells missing in WNN driver output; filling with 0.",
                self.name, n_missing,
            )
            coords = np.nan_to_num(coords, nan=0.0)

        adata.obsm[cfg.wnn_obsm_output_key] = coords
        adata.uns["multimodal_status"] = {
            "engine": "wnn",
            "status": "ok",
            "n_cells": int(coords.shape[0]),
            "rna_obsm_key": rna_key,
            "second_obsm_key": sec_key,
        }
        ctx.metadata["multimodal_status"] = "ok"
        ctx.metadata["multimodal_engine"] = "wnn"

    def _run_mofa(
        self,
        ctx: PipelineContext,
        cfg: MultimodalIntegrationConfig,
        rna_key: str,
        sec_key: str,
    ) -> None:
        adata = ctx.adata

        # Try mofapy2 first; fall back to muon.tools.mofa wrappers.
        try:
            import mofapy2  # type: ignore
            from mofapy2.run.entry_point import entry_point  # type: ignore
        except ImportError:
            try:
                import muon  # type: ignore  # noqa: F401
                # muon is installed but we still need mofapy2 underneath; if
                # that fails the import above will have raised. If muon
                # imports, mofa probably runs through it.
                from mofapy2.run.entry_point import entry_point  # type: ignore
            except ImportError:
                reason = "neither mofapy2 nor muon are importable"
                logger.info("%s: %s; MOFA skipped.", self.name, reason)
                adata.uns["multimodal_status"] = {
                    "engine": "mofa", "status": "skipped", "reason": reason,
                }
                ctx.metadata["multimodal_status"] = "skipped_no_mofa"
                ctx.status(self.name, "skipped", reason)
                return

        rna_arr = np.asarray(adata.obsm[rna_key]).astype(np.float64, copy=False)
        sec_arr = np.asarray(adata.obsm[sec_key]).astype(np.float64, copy=False)

        try:
            ent = entry_point()
            # mofapy2 expects a list-of-lists: [view][group] -> ndarray.
            ent.set_data_matrix([[rna_arr], [sec_arr]],
                                views_names=["rna_latent", "modality_2_latent"],
                                groups_names=["group0"])
            ent.set_model_options(factors=cfg.n_factors)
            ent.set_train_options(iter=200, convergence_mode="medium", verbose=False)
            ent.build()
            ent.run()

            # Extract factor matrix Z (n_cells x n_factors). API varies by version;
            # try the documented public surfaces and fall back to the model object.
            factors = None
            try:
                factors = np.asarray(ent.model.getExpectations()["Z"]["E"])
            except Exception:  # pragma: no cover
                pass
            if factors is None:
                try:
                    factors = np.asarray(ent.model.getNodes()["Z"].getExpectation())
                except Exception:  # pragma: no cover
                    pass
            if factors is None:
                raise RuntimeError(
                    "mofapy2 ran but no factor matrix could be extracted from the model"
                )

            explained_variance: list | None = None
            try:
                explained_variance = list(
                    np.asarray(ent.model.calculate_variance_explained())
                    .astype(float)
                    .ravel()
                    .tolist()
                )
            except Exception:  # pragma: no cover
                explained_variance = None
        except Exception as exc:
            reason = f"mofapy2 training failed: {exc}"
            logger.warning("%s: %s; MOFA skipped.", self.name, reason)
            adata.uns["multimodal_status"] = {
                "engine": "mofa", "status": "skipped", "reason": reason,
            }
            ctx.metadata["multimodal_status"] = "skipped_mofa_failed"
            ctx.status(self.name, "skipped", reason)
            return

        adata.obsm[cfg.mofa_obsm_output_key] = factors.astype(np.float32, copy=False)
        adata.uns["mofa_factors"] = {
            "n_factors": int(factors.shape[1]),
            "explained_variance": explained_variance,
        }
        adata.uns["multimodal_status"] = {
            "engine": "mofa",
            "status": "ok",
            "n_cells": int(factors.shape[0]),
            "n_factors": int(factors.shape[1]),
            "rna_obsm_key": rna_key,
            "second_obsm_key": sec_key,
        }
        ctx.metadata["multimodal_status"] = "ok"
        ctx.metadata["multimodal_engine"] = "mofa"

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _merge_runtime_config(self, ctx: PipelineContext) -> MultimodalIntegrationConfig:
        runtime = getattr(ctx.cfg, "multimodal_integration", None)
        if runtime is None:
            return self.config
        if isinstance(runtime, MultimodalIntegrationConfig):
            return runtime
        if isinstance(runtime, dict):
            base = self.config
            return MultimodalIntegrationConfig(
                engine=str(runtime.get("engine", base.engine)),
                rna_obsm_key=str(runtime.get("rna_obsm_key", base.rna_obsm_key)),
                second_obsm_key=str(runtime.get("second_obsm_key", base.second_obsm_key)),
                wnn_obsm_output_key=str(runtime.get("wnn_obsm_output_key", base.wnn_obsm_output_key)),
                mofa_obsm_output_key=str(runtime.get("mofa_obsm_output_key", base.mofa_obsm_output_key)),
                n_factors=int(runtime.get("n_factors", base.n_factors)),
                rscript_bin=runtime.get("rscript_bin", base.rscript_bin),
                run_wnn_r_path=runtime.get("run_wnn_r_path", base.run_wnn_r_path),
                subprocess_timeout=int(runtime.get("subprocess_timeout", base.subprocess_timeout)),
            )
        return self.config

    @staticmethod
    def _resolve_engine(cfg: MultimodalIntegrationConfig) -> str:
        env = os.environ.get("SC_MULTIMODAL_ENGINE")
        if env:
            return env.strip().lower()
        return cfg.engine.strip().lower()

    @staticmethod
    def _resolve_rscript(cfg: MultimodalIntegrationConfig) -> str | None:
        if cfg.rscript_bin and os.access(cfg.rscript_bin, os.X_OK):
            return cfg.rscript_bin
        env_override = os.environ.get("RSCRIPT_BIN")
        if env_override and os.access(env_override, os.X_OK):
            return env_override
        return shutil.which("Rscript")

    @staticmethod
    def _resolve_run_wnn_r_path(cfg: MultimodalIntegrationConfig) -> str:
        if cfg.run_wnn_r_path:
            return cfg.run_wnn_r_path
        env_override = os.environ.get("RUN_WNN_R_PATH")
        if env_override:
            return env_override
        return _DEFAULT_RUN_WNN_R

    @staticmethod
    def _write_latent_parquet(
        arr: np.ndarray, cells: pd.Index, path: Path, *, prefix: str,
    ) -> None:
        if arr.ndim != 2:
            raise ValueError(f"latent matrix must be 2D; got shape {arr.shape}")
        cols = [f"{prefix}{i + 1}" for i in range(arr.shape[1])]
        df = pd.DataFrame(arr.astype(np.float32, copy=False), index=cells, columns=cols)
        df.index.name = "cell"
        df = df.reset_index()
        df.to_parquet(path, index=False)
