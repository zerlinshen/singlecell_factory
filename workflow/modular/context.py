from __future__ import annotations

import gc
import io
import json
import logging
import os
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

    def _should_save_adata_checkpoint(self, module_name: str) -> bool:
        massive_policy = os.environ.get("SCF_MASSIVE_CHECKPOINT_POLICY", "").strip().lower()
        if (
            self.cfg.scale_mode == "massive"
            and massive_policy in {"metadata_only", "json_only", "sidecar_only"}
        ):
            return False
        if self.cfg.scale_mode != "massive":
            return True
        return module_name not in {"cellranger", "qc", "doublet_detection"}

    def _compact_adata_for_checkpoint(self) -> None:
        if self.adata is None:
            return
        adata = self.adata

        def _compact_value(value):
            if sparse.issparse(value):
                return value.astype(np.float32) if value.dtype == np.float64 else value
            if isinstance(value, np.ndarray) and value.dtype == np.float64:
                return value.astype(np.float32, copy=False)
            return value

        def _compact_uns(value):
            if isinstance(value, dict):
                return {k: _compact_uns(v) for k, v in value.items()}
            if isinstance(value, list):
                return [_compact_uns(v) for v in value]
            return _compact_value(value)

        adata.X = _compact_value(adata.X)
        for key in list(adata.layers.keys()):
            adata.layers[key] = _compact_value(adata.layers[key])
        for attr in ("obsm", "varm", "obsp", "varp"):
            store = getattr(adata, attr, None)
            if store is None:
                continue
            for key in list(store.keys()):
                store[key] = _compact_value(store[key])
        adata.uns = _compact_uns(dict(adata.uns))

    def save_checkpoint(self, module_name: str) -> None:
        """Save adata + metadata after a module completes successfully.

        Uses zarr format when available (3-5x faster I/O), falls back to h5ad.
        """
        if not self.cfg.checkpoint:
            return
        cp_dir = self._checkpoint_dir
        cp_dir.mkdir(parents=True, exist_ok=True)
        if self.adata is not None and self._should_save_adata_checkpoint(module_name):
            self._compact_adata_for_checkpoint()
            if self._use_zarr_checkpoints():
                try:
                    zarr_path = cp_dir / f"after_{module_name}.zarr"
                    self.adata.write_zarr(zarr_path)
                except (ValueError, TypeError):
                    # Zarr rejects keys with forward slashes; fall back to h5ad
                    self.adata.write(cp_dir / f"after_{module_name}.h5ad")
            else:
                self.adata.write(cp_dir / f"after_{module_name}.h5ad")
        elif self.adata is not None:
            logger.info("Skipping full AnnData checkpoint after %s in massive mode to reduce memory/IO pressure", module_name)
        sidecar = {
            "module": module_name,
            "adata_checkpoint_policy": os.environ.get("SCF_MASSIVE_CHECKPOINT_POLICY", ""),
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
