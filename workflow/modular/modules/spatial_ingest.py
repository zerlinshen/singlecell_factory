"""Spatial transcriptomics ingest module.

Second consumer of the bundle v2.1 extension API
(``scripts.export_singlecell_r_bundle.add_extension``).

Inputs
------
The module attaches a per-cell/spot 2D coordinate matrix to ``ctx.adata``
from one of two sources, in this order:

1. ``adata.obsm["spatial"]`` already populated -- treat as a no-op for the
   matrix itself but still backfill metadata.
2. A spaceranger output dir (``cfg.spaceranger_dir``) containing
   ``spatial/tissue_positions_list.csv``. We accept the canonical Visium
   layout (one CSV per library; ``in_tissue`` flag is honoured when
   ``cfg.in_tissue_only`` is True).
3. A free-form CSV / parquet file (``cfg.coords_path``) with columns
   ``[barcode, x, y]`` plus optional ``sample`` and ``library_id`` columns.
   The barcode column joins to ``adata.obs_names``.

If none of the three are available the module records
``spatial_status="skipped_no_input"`` and returns without raising -- spatial
is a fail-soft modality in this pipeline.

Outputs
-------
- ``adata.obsm["spatial"]`` -- (n_cells x 2) float32 array of (x, y) coords.
- ``adata.uns["spatial_coord_system"]`` -- dict with ``platform``,
  ``units``, ``in_tissue_only`` keys.
- ``adata.uns["spatial"]["library_id"]`` -- per-library image-path pointers
  (STRING paths only; we never serialize image bytes here).
- ``adata.obs["spatial_in_tissue"]`` -- boolean column when an in_tissue
  flag was provided by the source CSV.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "10x_visium_2020": {
        "title": "Visium Spatial Gene Expression",
        "authors": "10x Genomics",
        "year": "2020",
        "description": "Spatial barcoded array; tissue_positions_list.csv schema.",
    },
    "merfish_chen_2015": {
        "title": "RNA imaging. Spatially resolved, highly multiplexed RNA profiling in single cells",
        "authors": "Chen et al.",
        "journal": "Science",
        "year": "2015",
        "doi": "10.1126/science.aaa6090",
    },
    "xenium_10x_2023": {
        "title": "Xenium In Situ platform",
        "authors": "10x Genomics",
        "year": "2023",
        "description": "Subcellular in-situ assay; same (x,y) per-cell schema.",
    },
}


_SUPPORTED_PLATFORMS = ("visium", "merfish", "xenium")
_DEFAULT_UNITS = {
    "visium": "pixel",
    "merfish": "micron",
    "xenium": "micron",
}


@dataclass(frozen=True)
class SpatialIngestConfig:
    """Local config for the spatial ingest module. Module-local per project rules."""

    platform: str = "visium"
    spaceranger_dir: Path | None = None
    coords_path: Path | None = None
    barcode_column: str = "barcode"
    x_column: str = "x"
    y_column: str = "y"
    sample_column: str | None = "sample"
    library_id_column: str | None = "library_id"
    in_tissue_only: bool = True
    library_image_paths: dict[str, str] | None = None
    obsm_output_key: str = "spatial"


class SpatialIngestModule:
    """Optional module: attach 2D spatial coordinates and metadata to adata.

    Registers ``adata.obsm["spatial"]`` plus ``adata.uns["spatial_coord_system"]``
    and (when configured) per-library image-path pointers in
    ``adata.uns["spatial"]["library_id"]``. The bundle exporter publishes this
    data via the v2.1 ``spatial`` extension when ``include_spatial=True``.
    """

    name = "spatial_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "obsm": ["spatial"],
        "uns": ["spatial_coord_system", "spatial"],
    }

    def __init__(self, config: SpatialIngestConfig | None = None) -> None:
        self.config = config or SpatialIngestConfig()

    # ------------------------------------------------------------------
    # Public entrypoint
    # ------------------------------------------------------------------

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        cfg = self._merge_runtime_config(ctx)
        if cfg.platform not in _SUPPORTED_PLATFORMS:
            raise ValueError(
                f"{self.name}: unsupported platform '{cfg.platform}'; "
                f"expected one of {_SUPPORTED_PLATFORMS}."
            )

        coords, in_tissue, sample, library_id = self._load_coords(adata, cfg)
        if coords is None:
            logger.info("%s: no spatial input found; skipping.", self.name)
            ctx.metadata["spatial_status"] = "skipped_no_input"
            ctx.status(self.name, "skipped", "no spatial input found")
            return

        n_cells = adata.n_obs
        if coords.shape[0] != n_cells:
            raise ValueError(
                f"{self.name}: coords row count ({coords.shape[0]}) != n_obs ({n_cells}); "
                "barcode alignment failed."
            )

        adata.obsm[cfg.obsm_output_key] = coords.astype(np.float32, copy=False)
        if in_tissue is not None:
            adata.obs["spatial_in_tissue"] = in_tissue.astype(bool)
        if sample is not None:
            adata.obs["spatial_sample"] = pd.Categorical(sample.astype(str))
        if library_id is not None:
            adata.obs["spatial_library_id"] = pd.Categorical(library_id.astype(str))

        adata.uns["spatial_coord_system"] = {
            "platform": cfg.platform,
            "units": _DEFAULT_UNITS[cfg.platform],
            "in_tissue_only": bool(cfg.in_tissue_only),
        }

        # Image paths are stored as plain strings; we never serialize bytes.
        spatial_uns = adata.uns.get("spatial")
        if not isinstance(spatial_uns, dict):
            spatial_uns = {}
        if cfg.library_image_paths:
            spatial_uns["library_id"] = {
                str(lib): str(path) for lib, path in cfg.library_image_paths.items()
            }
        adata.uns["spatial"] = spatial_uns

        ctx.metadata["spatial_n_spots"] = int(n_cells)
        ctx.metadata["spatial_platform"] = cfg.platform
        ctx.metadata["spatial_status"] = "ok"
        n_libs = 0
        if isinstance(spatial_uns.get("library_id"), dict):
            n_libs = len(spatial_uns["library_id"])
        ctx.metadata["spatial_n_libraries"] = int(n_libs)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _merge_runtime_config(self, ctx: PipelineContext) -> SpatialIngestConfig:
        runtime = getattr(ctx.cfg, "spatial_ingest", None)
        if runtime is None:
            return self.config
        if isinstance(runtime, SpatialIngestConfig):
            return runtime
        if isinstance(runtime, dict):
            base = self.config
            return SpatialIngestConfig(
                platform=str(runtime.get("platform", base.platform)),
                spaceranger_dir=runtime.get("spaceranger_dir", base.spaceranger_dir),
                coords_path=runtime.get("coords_path", base.coords_path),
                barcode_column=str(runtime.get("barcode_column", base.barcode_column)),
                x_column=str(runtime.get("x_column", base.x_column)),
                y_column=str(runtime.get("y_column", base.y_column)),
                sample_column=runtime.get("sample_column", base.sample_column),
                library_id_column=runtime.get("library_id_column", base.library_id_column),
                in_tissue_only=bool(runtime.get("in_tissue_only", base.in_tissue_only)),
                library_image_paths=runtime.get("library_image_paths", base.library_image_paths),
                obsm_output_key=str(runtime.get("obsm_output_key", base.obsm_output_key)),
            )
        return self.config

    def _load_coords(
        self, adata, cfg: SpatialIngestConfig,
    ) -> tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
        """Resolve spatial coords from one of three sources.

        Returns ``(coords, in_tissue, sample, library_id)``. Any of the
        ancillary arrays may be ``None`` if absent. Coords are always aligned
        to ``adata.obs_names`` ordering on return.
        """
        # Path 1: coords already present on adata.
        if cfg.obsm_output_key in adata.obsm:
            arr = np.asarray(adata.obsm[cfg.obsm_output_key])
            if arr.ndim != 2 or arr.shape[1] < 2:
                raise ValueError(
                    f"{self.name}: existing adata.obsm['{cfg.obsm_output_key}'] must be "
                    f"(n_cells x 2); got shape {arr.shape}"
                )
            coords = arr[:, :2].astype(np.float32, copy=False)
            in_tissue = None
            if "spatial_in_tissue" in adata.obs.columns:
                in_tissue = np.asarray(adata.obs["spatial_in_tissue"].values)
            return coords, in_tissue, None, None

        # Path 2: spaceranger directory.
        if cfg.spaceranger_dir is not None:
            sr_dir = Path(cfg.spaceranger_dir)
            csv_path = sr_dir / "spatial" / "tissue_positions_list.csv"
            if csv_path.exists():
                return self._load_visium_csv(adata, csv_path, cfg)
            logger.warning("%s: spaceranger_dir provided but %s not found.",
                           self.name, csv_path)

        # Path 3: free-form CSV or parquet.
        if cfg.coords_path is not None:
            cp = Path(cfg.coords_path)
            if cp.exists():
                return self._load_table(adata, cp, cfg)
            logger.warning("%s: coords_path %s not found.", self.name, cp)

        return None, None, None, None

    def _load_visium_csv(
        self, adata, csv_path: Path, cfg: SpatialIngestConfig,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
        """Read the canonical Visium tissue_positions_list.csv layout.

        Two known variants exist in the wild:
          * Legacy (no header): columns are
            barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
          * 2022+ (with header): same columns but with a header line.

        We sniff by checking whether the first row's second field parses as
        an int; the legacy layout has no header.
        """
        # Sniff for header by reading first line.
        with csv_path.open("r") as fh:
            first_line = fh.readline().strip()
        has_header = True
        if first_line:
            parts = first_line.split(",")
            if len(parts) >= 2:
                try:
                    int(parts[1])
                    has_header = False  # second column was already an int -> no header
                except ValueError:
                    has_header = True

        legacy_cols = [
            "barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres",
        ]
        if has_header:
            df = pd.read_csv(csv_path)
        else:
            df = pd.read_csv(csv_path, header=None, names=legacy_cols)

        # Tolerate column-name variants seen in the wild.
        col_map = {c.lower(): c for c in df.columns}
        bc_col = col_map.get("barcode") or df.columns[0]
        x_col = (col_map.get("pxl_col_in_fullres")
                 or col_map.get("x")
                 or col_map.get("array_col"))
        y_col = (col_map.get("pxl_row_in_fullres")
                 or col_map.get("y")
                 or col_map.get("array_row"))
        in_tissue_col = col_map.get("in_tissue")

        if x_col is None or y_col is None:
            raise ValueError(
                f"{self.name}: cannot find x/y columns in {csv_path}; "
                f"got columns {list(df.columns)}"
            )

        df[bc_col] = df[bc_col].astype(str)
        if cfg.in_tissue_only and in_tissue_col is not None:
            df = df[df[in_tissue_col].astype(int) == 1]

        return self._align_to_obs(adata, df, cfg, bc_col, x_col, y_col,
                                  in_tissue_col=in_tissue_col,
                                  sample_col=None, library_col=None)

    def _load_table(
        self, adata, path: Path, cfg: SpatialIngestConfig,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
        """Load a CSV or parquet with [barcode, x, y, ...] columns."""
        suffix = path.suffix.lower()
        if suffix == ".parquet":
            df = pd.read_parquet(path)
        else:
            df = pd.read_csv(path)

        for required in (cfg.barcode_column, cfg.x_column, cfg.y_column):
            if required not in df.columns:
                raise ValueError(
                    f"{self.name}: column '{required}' not found in {path}; "
                    f"got columns {list(df.columns)}"
                )

        sample_col = cfg.sample_column if cfg.sample_column in df.columns else None
        library_col = cfg.library_id_column if cfg.library_id_column in df.columns else None
        in_tissue_col = "in_tissue" if "in_tissue" in df.columns else None
        return self._align_to_obs(
            adata, df, cfg, cfg.barcode_column, cfg.x_column, cfg.y_column,
            in_tissue_col=in_tissue_col,
            sample_col=sample_col, library_col=library_col,
        )

    def _align_to_obs(
        self,
        adata,
        df: pd.DataFrame,
        cfg: SpatialIngestConfig,
        bc_col: str,
        x_col: str,
        y_col: str,
        *,
        in_tissue_col: str | None,
        sample_col: str | None,
        library_col: str | None,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
        """Reorder the spatial table so its rows align with adata.obs_names."""
        df = df.copy()
        df[bc_col] = df[bc_col].astype(str)
        df = df.drop_duplicates(subset=[bc_col], keep="first")
        df = df.set_index(bc_col)

        obs_names = pd.Index(adata.obs_names.astype(str))
        present = obs_names.intersection(df.index)
        if len(present) == 0:
            raise ValueError(
                f"{self.name}: no overlap between spatial barcodes and adata.obs_names "
                f"(adata={len(obs_names)} cells; table={len(df)} rows)"
            )
        df = df.reindex(obs_names)

        coords = df[[x_col, y_col]].to_numpy(dtype=np.float32, copy=True)
        # Cells without coords end up as NaN; clamp to 0 and warn.
        n_missing = int(np.isnan(coords).any(axis=1).sum())
        if n_missing > 0:
            logger.warning(
                "%s: %d / %d cells missing spatial coords; filling with 0.",
                self.name, n_missing, adata.n_obs,
            )
            coords = np.nan_to_num(coords, nan=0.0)

        in_tissue = None
        if in_tissue_col is not None and in_tissue_col in df.columns:
            in_tissue = df[in_tissue_col].fillna(0).astype(int).to_numpy()
        sample = (
            df[sample_col].astype(str).fillna("").to_numpy()
            if sample_col is not None and sample_col in df.columns
            else None
        )
        library_id = (
            df[library_col].astype(str).fillna("").to_numpy()
            if library_col is not None and library_col in df.columns
            else None
        )
        return coords, in_tissue, sample, library_id
