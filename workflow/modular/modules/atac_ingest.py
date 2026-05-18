"""ATAC-seq ingest (Wave 2 / P2.S11): sparse peak matrix + TF-IDF/LSI.

Memory-safe per Wave 1 lesson (no full sparse-matrix materialization).
TF-IDF + TruncatedSVD both accept sparse input natively.

Inputs (via cfg, both optional):
  --atac-peak-matrix-path PATH   sparse peak matrix in mtx.gz or h5 format
  --atac-peaks-bed-path PATH     accompanying peaks.bed (chrom, start, end)
  --atac-fragment-path PATH      fragments.tsv.gz for FRiP/TSS (Wave 2B)

Outputs:
  adata.obsm["atac_peaks"]            sparse peak count matrix (n_cells x n_peaks)
  adata.obsm["X_atac"]                compatibility LSI embedding (n_cells x n_components)
  adata.uns["atac_peaks"]             peak coords DataFrame
  adata.uns["atac_var"]               peak coords DataFrame with peak_id for linkage code
  adata.uns["atac_ingest_metadata"]   shape, dtype, n_components, density
  runs/<run-id>/atac_ingest/atac_summary.json
"""
from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

from ..context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "SnapATAC2": {
        "title": "A fast, scalable and versatile tool for analysis of single-cell omics data",
        "authors": "Zhang et al.",
        "journal": "Nature Methods",
        "year": "2024",
        "doi": "10.1038/s41592-023-02139-9",
        "description": "Python ATAC-seq analysis framework (used when --atac-fragment-path provided)",
    },
    "Signac": {
        "title": "Multimodal single-cell chromatin analysis with Signac",
        "authors": "Stuart et al.",
        "journal": "Nature Methods",
        "year": "2021",
        "doi": "10.1038/s41592-021-01282-5",
        "description": "R-side ATAC analysis and visualization (consumes bundle extensions/atac/)",
    },
    "Cusanovich2018": {
        "title": "A single-cell atlas of in vivo mammalian chromatin accessibility",
        "authors": "Cusanovich et al.",
        "journal": "Cell",
        "year": "2018",
        "doi": "10.1016/j.cell.2018.06.052",
        "description": "TF-IDF + LSI methodology for sparse ATAC peak matrices",
    },
}

_DEFAULT_N_COMPONENTS = 30


def _load_peak_matrix(path: Path) -> sp.csr_matrix:
    """Load a sparse peak matrix. Supports .mtx / .mtx.gz / .h5 (cellranger-atac)."""
    path = Path(path)
    if path.suffix == ".gz" and ".mtx" in path.name:
        with gzip.open(path, "rb") as f:
            mat = sio.mmread(f).tocsr()
    elif path.suffix == ".mtx":
        mat = sio.mmread(str(path)).tocsr()
    elif path.suffix == ".h5":
        try:
            import h5py
        except ImportError as exc:
            raise RuntimeError("h5 ATAC matrix requires h5py") from exc
        with h5py.File(path, "r") as f:
            # cellranger-atac filtered_peak_bc_matrix.h5 layout
            grp = f["matrix"] if "matrix" in f else f
            data = grp["data"][:]
            indices = grp["indices"][:]
            indptr = grp["indptr"][:]
            shape = tuple(grp["shape"][:])
            mat = sp.csr_matrix((data, indices, indptr), shape=shape)
    else:
        raise ValueError(f"unsupported peak matrix format: {path.suffix}")
    return mat


def _load_peaks_bed(path: Path) -> pd.DataFrame:
    """Read peaks.bed (chrom, start, end)."""
    peaks = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end"],
        usecols=[0, 1, 2],
        dtype={"chrom": str, "start": int, "end": int},
    )
    peaks["peak_id"] = (
        peaks["chrom"].astype(str) + ":" +
        peaks["start"].astype(str) + "-" +
        peaks["end"].astype(str)
    )
    return peaks


def _tfidf_lsi(X: sp.csr_matrix, n_components: int, random_state: int) -> np.ndarray:
    """Sparse TF-IDF + TruncatedSVD (LSI) following Cusanovich 2018.

    X: cells x peaks (CSR). NEVER densified.
    Returns dense (cells, n_components) LSI embedding.
    """
    from sklearn.feature_extraction.text import TfidfTransformer
    from sklearn.decomposition import TruncatedSVD

    # cellranger-atac peak matrices are peaks x cells; convention here is cells x peaks.
    # We accept input in cells x peaks (caller responsibility).
    n_cells = X.shape[0]
    if n_cells == 0:
        return np.zeros((0, n_components), dtype=np.float32)

    # Binarize (any non-zero peak count -> 1)
    X_bin = X.copy()
    X_bin.data = (X_bin.data > 0).astype(np.float32)

    tfidf = TfidfTransformer(norm="l1", sublinear_tf=True)
    X_tfidf = tfidf.fit_transform(X_bin)  # stays sparse

    n_components = min(n_components, X_tfidf.shape[1] - 1, X_tfidf.shape[0] - 1)
    n_components = max(2, n_components)
    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    embedding = svd.fit_transform(X_tfidf).astype(np.float32)
    return embedding


class ATACIngestModule:
    """ATAC-seq ingest: TF-IDF + LSI embedding from a sparse peak matrix.

    Fail-soft skip when no ATAC input is provided. Memory-safe (sparse-only).
    """

    name = "atac_ingest"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "obsm": ["atac_peaks", "X_atac"],
        "uns": ["atac_peaks", "atac_var"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata
        peak_matrix_path = getattr(ctx.cfg, "atac_peak_matrix_path", None)
        peaks_bed_path = getattr(ctx.cfg, "atac_peaks_bed_path", None)
        n_components = int(getattr(ctx.cfg, "atac_n_components", _DEFAULT_N_COMPONENTS))

        if not peak_matrix_path:
            ctx.status(self.name, "skipped", "atac_peak_matrix_path not set")
            ctx.metadata["atac_ingest_status"] = "skipped_no_input"
            return

        peak_matrix_path = Path(peak_matrix_path)
        if not peak_matrix_path.exists():
            ctx.status(self.name, "skipped", f"peak matrix path missing: {peak_matrix_path}")
            ctx.metadata["atac_ingest_status"] = "skipped_path_missing"
            return

        logger.info("%s: loading sparse peak matrix from %s", self.name, peak_matrix_path)
        peak_mat = _load_peak_matrix(peak_matrix_path)

        # Heuristic: align orientation to cells x peaks (n_rows == n_cells).
        # If the matrix has obs_names.size rows, it's already cells x peaks.
        # Otherwise transpose (cellranger-atac default is peaks x cells).
        if peak_mat.shape[0] == adata.n_obs:
            cells_x_peaks = peak_mat
        elif peak_mat.shape[1] == adata.n_obs:
            cells_x_peaks = peak_mat.T.tocsr()
        else:
            ctx.status(
                self.name,
                "skipped",
                f"peak matrix shape {peak_mat.shape} does not match adata.n_obs={adata.n_obs}",
            )
            ctx.metadata["atac_ingest_status"] = "skipped_shape_mismatch"
            return

        # Optional peaks.bed
        peaks_df: pd.DataFrame | None = None
        if peaks_bed_path and Path(peaks_bed_path).exists():
            peaks_df = _load_peaks_bed(Path(peaks_bed_path))
            if len(peaks_df) != cells_x_peaks.shape[1]:
                logger.warning(
                    "%s: peaks.bed has %d rows but matrix has %d peaks — discarding peak metadata",
                    self.name, len(peaks_df), cells_x_peaks.shape[1],
                )
                peaks_df = None

        random_state = int(getattr(ctx, "random_state", 42))
        logger.info(
            "%s: TF-IDF + LSI (%d cells x %d peaks, n_components=%d, random_state=%d)",
            self.name, cells_x_peaks.shape[0], cells_x_peaks.shape[1], n_components, random_state,
        )
        cells_x_peaks = cells_x_peaks.tocsr()
        embedding = _tfidf_lsi(cells_x_peaks, n_components=n_components, random_state=random_state)

        adata.obsm["atac_peaks"] = cells_x_peaks
        adata.obsm["X_atac"] = embedding
        if peaks_df is not None:
            adata.uns["atac_peaks"] = peaks_df
            adata.uns["atac_var"] = peaks_df.copy()
        adata.uns["atac_ingest_metadata"] = {
            "n_peaks": int(cells_x_peaks.shape[1]),
            "n_components": int(embedding.shape[1]),
            "peak_matrix_obsm_key": "atac_peaks",
            "lsi_compat_obsm_key": "X_atac",
            "density": float(cells_x_peaks.nnz) / (cells_x_peaks.shape[0] * cells_x_peaks.shape[1]),
            "random_state": random_state,
        }

        out_dir = ctx.run_dir / "atac_ingest"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "atac_summary.json").write_text(
            json.dumps(adata.uns["atac_ingest_metadata"], indent=2), encoding="utf-8"
        )

        ctx.metadata["atac_ingest_status"] = "ok"
        ctx.metadata["atac_n_components"] = int(embedding.shape[1])
        logger.info(
            "%s: wrote obsm['atac_peaks'] shape=%s, X_atac shape=%s, and atac_summary.json",
            self.name, cells_x_peaks.shape, embedding.shape,
        )
