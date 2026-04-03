from __future__ import annotations

import io
import json
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt

from .config import PipelineConfig


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

    def status(self, module: str, ok: bool, message: str) -> None:
        """Record module execution status."""
        self.module_status.append(
            {"module": module, "status": "ok" if ok else "failed", "message": message}
        )

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

    def save_checkpoint(self, module_name: str) -> None:
        """Save adata + metadata after a module completes successfully."""
        if not self.cfg.checkpoint:
            return
        cp_dir = self._checkpoint_dir
        cp_dir.mkdir(parents=True, exist_ok=True)
        if self.adata is not None:
            self.adata.write(cp_dir / f"after_{module_name}.h5ad")
        sidecar = {
            "module": module_name,
            "metadata": self.metadata,
            "module_status": self.module_status,
        }
        (cp_dir / f"after_{module_name}.json").write_text(
            json.dumps(sidecar, indent=2, ensure_ascii=False), encoding="utf-8"
        )

    def load_checkpoint(self, module_name: str) -> bool:
        """Load checkpoint saved after *module_name*. Returns True if loaded."""
        cp_h5ad = self._checkpoint_dir / f"after_{module_name}.h5ad"
        cp_json = self._checkpoint_dir / f"after_{module_name}.json"
        if not cp_h5ad.exists():
            return False
        self.adata = ad.read_h5ad(cp_h5ad)
        if cp_json.exists():
            sidecar = json.loads(cp_json.read_text(encoding="utf-8"))
            self.metadata = sidecar.get("metadata", {})
            self.module_status = sidecar.get("module_status", [])
        return True

    # --- Async figure saving ---

    def save_figure(self, fig, path: Path, dpi: int = 160, **kwargs) -> None:
        """Save a matplotlib figure, optionally offloading disk I/O to a background thread."""
        if self._figure_pool is not None:
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", **kwargs)
            plt.close(fig)
            data = buf.getvalue()
            self._figure_futures.append(self._figure_pool.submit(path.write_bytes, data))
        else:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
            plt.close(fig)

    def flush_figures(self) -> None:
        """Wait for all pending async figure saves to complete."""
        for f in self._figure_futures:
            f.result()
        self._figure_futures.clear()
