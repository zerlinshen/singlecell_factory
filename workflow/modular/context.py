from __future__ import annotations

import gc
import io
import json
import logging
import os
import shutil
import threading
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

from .config import PipelineConfig

logger = logging.getLogger(__name__)

SC_CHECKPOINT_SCHEMA_VERSION = 2


class CheckpointSchemaMismatch(Exception):
    pass


@dataclass
class PipelineContext:
    """Runtime state shared by all workflow modules."""

    cfg: PipelineConfig
    run_dir: Path
    figure_dir: Path
    table_dir: Path
    adata: ad.AnnData | None = None
    metadata: dict = field(default_factory=dict)
    module_status: list[dict] = field(default_factory=list)
    _module_dirs: dict[str, Path] = field(default_factory=dict, repr=False)
    _figure_pool: ThreadPoolExecutor | None = field(default=None, repr=False)
    _figure_futures: list[Future] = field(default_factory=list, repr=False)
    _status_lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    @property
    def random_state(self) -> int:
        """Global random seed for all stochastic modules. Source of truth for AC-10 reproducibility."""
        return getattr(self.cfg, "random_state", 42)

    def status(self, module: str, ok: bool | str, message: str) -> None:
        """Record module execution status (thread-safe)."""
        if isinstance(ok, str):
            status = ok
        else:
            status = "ok" if ok else "failed"
        entry = {"module": module, "status": status, "message": message}
        with self._status_lock:
            self.module_status.append(entry)

    def set_module_dir(self, module_name: str) -> None:
        """Point figure_dir and table_dir to a per-module subfolder."""
        mod_dir = self.run_dir / module_name
        mod_dir.mkdir(parents=True, exist_ok=True)
        self.figure_dir = mod_dir
        self.table_dir = mod_dir
        self._module_dirs[module_name] = mod_dir

    def module_output_dir(self, module_name: str) -> Path | None:
        """Return the output directory of a previously-run module."""
        return self._module_dirs.get(module_name)

    # --- Checkpointing ---

    @property
    def _checkpoint_dir(self) -> Path:
        return self.run_dir / ".checkpoints"

    def _use_zarr_checkpoints(self) -> bool:
        """Check if zarr is available for faster checkpoint I/O."""
        try:
            import zarr  # noqa: F401
            return True
        except ImportError:
            return False

    def _adata_checkpoint_policy(self) -> str:
        return str(
            os.environ.get("SC_CHECKPOINT_POLICY")
            or os.environ.get("SCF_MASSIVE_CHECKPOINT_POLICY")
            or getattr(self.cfg, "checkpoint_policy", "full")
            or "full"
        ).lower()

    def _should_save_adata_checkpoint(self, module_name: str) -> bool:
        policy = self._adata_checkpoint_policy()
        return policy not in {
            "metadata",
            "metadata-only",
            "metadata_only",
            "sidecar",
            "json",
            "json-only",
            "json_only",
            "none",
            "skip-adata",
            "skip_adata",
        }

    @staticmethod
    def _is_unsupported_lazy_checkpoint_error(exc: BaseException) -> bool:
        message = str(exc)
        if isinstance(exc, NotImplementedError) and "Dataset2D" in message:
            return True
        return "Writing AnnData objects with a Dataset2D not supported" in message

    def _record_skipped_adata_checkpoint(self, module_name: str, reason: str) -> None:
        logger.warning("Skipping AnnData checkpoint after %s: %s", module_name, reason)
        self.metadata.setdefault("checkpoint_warnings", []).append(
            {"module": module_name, "reason": reason}
        )

    @staticmethod
    def _remove_partial_checkpoint(path: Path) -> None:
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        elif path.exists():
            path.unlink()

    def _write_adata_checkpoint(self, module_name: str, cp_dir: Path, adata: ad.AnnData) -> None:
        if adata is None:
            return

        h5ad_path = cp_dir / f"after_{module_name}.h5ad"

        def _write_h5ad() -> None:
            try:
                adata.write(h5ad_path)
            except (NotImplementedError, ValueError, TypeError) as exc:
                if self._is_unsupported_lazy_checkpoint_error(exc):
                    self._remove_partial_checkpoint(h5ad_path)
                    self._record_skipped_adata_checkpoint(module_name, str(exc))
                    return
                raise

        if not self._use_zarr_checkpoints():
            _write_h5ad()
            return

        zarr_path = cp_dir / f"after_{module_name}.zarr"
        try:
            adata.write_zarr(zarr_path)
        except NotImplementedError as exc:
            if self._is_unsupported_lazy_checkpoint_error(exc):
                self._remove_partial_checkpoint(zarr_path)
                self._record_skipped_adata_checkpoint(module_name, str(exc))
                return
            raise
        except (ValueError, TypeError):
            # Zarr rejects keys with forward slashes; fall back to h5ad.
            _write_h5ad()

    def _compact_adata_for_checkpoint(self) -> ad.AnnData | None:
        """Return a float32-compacted COPY of ``self.adata`` for on-disk checkpoints.

        Copy-on-write: the live ``self.adata`` is never mutated. Downcasting the
        checkpoint in place would silently change the dtypes the *next* module
        computes on, so with ``--checkpoint`` a run would diverge from a
        non-checkpoint run after the first checkpoint. We build a shallow AnnData
        copy and downcast only its float64 arrays, so subsequent modules keep
        operating on the original (float64) live object.
        """
        if self.adata is None:
            return None
        adata = self.adata

        def _compact_value(value):
            if sparse.issparse(value):
                return value.astype(np.float32) if value.dtype == np.float64 else value
            if isinstance(value, np.ndarray) and value.dtype == np.float64:
                # copy=True: never alias the live array; the live object must
                # retain float64 even though the checkpoint is downcast.
                return value.astype(np.float32, copy=True)
            return value

        def _compact_uns(value):
            if isinstance(value, dict):
                return {k: _compact_uns(v) for k, v in value.items()}
            if isinstance(value, list):
                return [_compact_uns(v) for v in value]
            return _compact_value(value)

        # Build the checkpoint object directly rather than ``adata.copy()`` +
        # overwrite: a full copy would deep-copy every array to float64 first and
        # then allocate the float32 downcasts on top, transiently doubling memory.
        # Constructing in place downcasts each large numeric array exactly once;
        # non-float64 arrays are shared by reference (read-only on the write path).
        compacted = ad.AnnData(
            X=_compact_value(adata.X),
            obs=adata.obs,
            var=adata.var,
            uns=_compact_uns(dict(adata.uns)),
            # AnnData 0.13 exposes X as ``layers[None]``. X is already passed
            # explicitly above, so only named layers belong in this mapping.
            layers={
                k: _compact_value(adata.layers[k])
                for k in adata.layers
                if k is not None
            },
            obsm={k: _compact_value(adata.obsm[k]) for k in adata.obsm},
            varm={k: _compact_value(adata.varm[k]) for k in adata.varm},
            obsp={k: _compact_value(adata.obsp[k]) for k in adata.obsp},
            varp={k: _compact_value(adata.varp[k]) for k in adata.varp},
            raw=adata.raw,
        )
        return compacted

    def save_checkpoint(self, module_name: str) -> None:
        """Save adata + metadata after a module completes successfully.

        Uses zarr format when available (3-5x faster I/O), falls back to h5ad.
        """
        if not self.cfg.checkpoint:
            return
        cp_dir = self._checkpoint_dir
        cp_dir.mkdir(parents=True, exist_ok=True)
        if self.adata is not None and self._should_save_adata_checkpoint(module_name):
            # Copy-on-write: compaction returns a float32 COPY; self.adata keeps
            # its original dtypes so later modules do not diverge onto float32.
            compacted = self._compact_adata_for_checkpoint()
            self._write_adata_checkpoint(module_name, cp_dir, compacted)
        elif self.adata is not None:
            self._record_skipped_adata_checkpoint(
                module_name,
                f"adata checkpoint policy is {self._adata_checkpoint_policy()}",
            )
        sidecar = {
            "schema_version": SC_CHECKPOINT_SCHEMA_VERSION,
            "module": module_name,
            "adata_checkpoint_policy": self._adata_checkpoint_policy(),
            "metadata": self.metadata,
            "module_status": self.module_status,
            "module_dirs": {k: str(v) for k, v in self._module_dirs.items()},
        }
        (cp_dir / f"after_{module_name}.json").write_text(
            json.dumps(sidecar, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        gc.collect()

    def load_checkpoint(self, module_name: str) -> bool:
        """Load checkpoint saved after *module_name*. Returns True if loaded."""
        cp_zarr = self._checkpoint_dir / f"after_{module_name}.zarr"
        cp_h5ad = self._checkpoint_dir / f"after_{module_name}.h5ad"
        cp_json = self._checkpoint_dir / f"after_{module_name}.json"
        if cp_zarr.exists():
            self.adata = ad.read_zarr(cp_zarr)
        elif cp_h5ad.exists():
            self.adata = ad.read_h5ad(cp_h5ad)
        else:
            return False
        if cp_json.exists():
            sidecar = json.loads(cp_json.read_text(encoding="utf-8"))
            schema_v = sidecar.get("schema_version", 1)
            if schema_v < SC_CHECKPOINT_SCHEMA_VERSION:
                raise CheckpointSchemaMismatch(
                    f"checkpoint at {cp_json} is v{schema_v} but pipeline is v{SC_CHECKPOINT_SCHEMA_VERSION}; "
                    f"regenerate prepared zarr OR rerun from cellranger. "
                    f"Use scripts/dev/regenerate_prepared_zarr.py to migrate."
                )
            self.metadata = sidecar.get("metadata", {})
            self.module_status = sidecar.get("module_status", [])
            raw_dirs = sidecar.get("module_dirs", {})
            restored_dirs: dict[str, Path] = {}
            if isinstance(raw_dirs, dict):
                for mod, path_str in raw_dirs.items():
                    p = Path(path_str)
                    if p.exists():
                        restored_dirs[str(mod)] = p
            # Backward compatibility: older checkpoints do not include module_dirs.
            if not restored_dirs:
                for entry in self.module_status:
                    mod = entry.get("module")
                    if not mod:
                        continue
                    mod_dir = self.run_dir / str(mod)
                    if mod_dir.exists():
                        restored_dirs[str(mod)] = mod_dir
            self._module_dirs = restored_dirs
        return True

    # --- Async figure saving ---

    def save_figure(self, fig, path: Path, dpi: int = 160, **kwargs) -> None:
        """Save a matplotlib figure, optionally offloading render+write to a background thread."""
        if self._figure_pool is not None:
            def _render_and_write(f, p, _dpi, _kw):
                buf = io.BytesIO()
                f.savefig(buf, format="png", dpi=_dpi, bbox_inches="tight", **_kw)
                plt.close(f)
                p.write_bytes(buf.getvalue())

            self._figure_futures.append(
                self._figure_pool.submit(_render_and_write, fig, path, dpi, kwargs)
            )
        else:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
            plt.close(fig)

    def flush_figures(self) -> None:
        """Wait for all pending async figure saves to complete."""
        errors = 0
        for f in self._figure_futures:
            try:
                f.result()
            except Exception as exc:
                errors += 1
                logger.warning("Async figure save failed: %s", exc)
        self._figure_futures.clear()
        if errors:
            logger.warning("%d figure save(s) failed", errors)
