"""ATAC LSI (Wave 5 / US-W5-1): TF-IDF normalisation + truncated SVD embedding.

TF-IDF+SVD (LSI) is preferred over standard PCA for ATAC peak matrices because
the data are binary/count sparse (not log-normal): TF-IDF re-weights peaks by
inverse document frequency before SVD, yielding a low-dimensional embedding that
captures open-chromatin cell-type structure while remaining memory-safe (operates
on the sparse peak matrix without densification). cisTopic (LDA) would give
probabilistic topic scores but requires orders-of-magnitude more compute and a
Bayesian inference step ill-suited to real-time pipeline runs.

Inputs (via adata):
  adata.obsm["atac_peaks"]   sparse CSR/CSC peak count matrix (n_obs × n_peaks)

Outputs:
  adata.obsm["X_lsi"]                 ndarray (n_obs, n_components - 1) — first component dropped
  adata.uns["lsi_variance_explained"]  list[float] of per-component ratios
  <module_dir>/lsi_variance_explained.png  variance-explained bar chart
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy import sparse

from ..context import PipelineContext
from .._contract_violation import ModuleContractError

logger = logging.getLogger(__name__)


__references__ = {
    "Cusanovich_2018_Cell": {
        "title": "A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility",
        "authors": "Cusanovich et al.",
        "journal": "Cell",
        "year": "2018",
        "doi": "10.1016/j.cell.2018.06.052",
        "description": (
            "Motivates TF-IDF normalisation of the binary peak matrix followed by SVD "
            "as the canonical dimensionality-reduction pipeline for single-cell ATAC "
            "(the TF-IDF step implemented here at run())."
        ),
    },
    "Stuart_2021_Nat_Methods_Signac": {
        "title": "Multimodal single-cell chromatin analysis with Signac",
        "authors": "Stuart et al.",
        "journal": "Nat Methods",
        "year": "2021",
        "doi": "10.1038/s41592-021-01282-5",
        "description": (
            "Establishes the convention of storing the LSI embedding in obsm['X_lsi'] "
            "and motivates dropping the first SVD component (depth-correlated) to remove "
            "sequencing-depth confounding — implemented via drop_first=True default."
        ),
    },
    "Granja_2021_Nat_Genet_ArchR": {
        "title": "ArchR is a scalable software package for integrative single-cell "
                 "chromatin accessibility analysis",
        "authors": "Granja et al.",
        "journal": "Nat Genet",
        "year": "2021",
        "doi": "10.1038/s41588-021-00790-6",
        "description": (
            "Confirms drop_first=True as the default for LSI (their iterativeLSI "
            "implementation discards the first LSI component on every iteration); "
            "motivates our choice of 50 components before dropping the first."
        ),
    },
}


def _tfidf_transform(X: sparse.spmatrix) -> sparse.csr_matrix:
    """Apply TF-IDF to a sparse peak count matrix.

    TF  = X / row_sum  (normalise each cell by total counts)
    IDF = log(1 + n_obs / col_sum)  (inverse document frequency across cells)
    """
    X = X.astype(np.float64)
    if not sparse.issparse(X):
        raise TypeError("_tfidf_transform requires a sparse matrix")
    X = X.tocsr()
    n_obs = X.shape[0]

    # TF: divide each row by its sum
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    tf = X.multiply(1.0 / row_sums[:, np.newaxis])

    # IDF: log(1 + N / df) per peak
    col_sums = np.asarray(X.sum(axis=0)).ravel()
    idf = np.log1p(n_obs / np.maximum(col_sums, 1.0))

    tfidf = tf.multiply(idf[np.newaxis, :])
    return tfidf.tocsr()


class ATACLSIModule:
    """ATAC LSI: TF-IDF normalisation + truncated SVD for open-chromatin embedding."""

    name = "atac_lsi"
    required = False
    mutates_structure = True
    requires_keys: dict[str, list[str]] = {"obsm": ["atac_peaks"]}
    provides_keys: dict[str, list[str]] = {
        "obsm": ["X_lsi"],
        "uns": ["lsi_variance_explained"],
    }

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        adata = ctx.adata

        # Contract: atac_peaks must exist in obsm
        if "atac_peaks" not in adata.obsm:
            raise ModuleContractError(
                f"{self.name}: required key 'atac_peaks' is absent from adata.obsm. "
                "Run atac_ingest before atac_lsi."
            )

        peak_matrix = adata.obsm["atac_peaks"]

        # Contract: atac_peaks must be sparse (densified matrices violate memory discipline)
        if not sparse.issparse(peak_matrix):
            raise ModuleContractError(
                f"{self.name}: adata.obsm['atac_peaks'] is dense (got {type(peak_matrix).__name__}). "
                "Peak matrices must remain sparse to avoid memory blowup. "
                "Re-ingest with atac_ingest, which stores CSR."
            )

        logger.info(
            "%s: peak matrix shape=%s nnz=%d", self.name, peak_matrix.shape, peak_matrix.nnz
        )

        # TF-IDF normalisation
        tfidf = _tfidf_transform(peak_matrix)

        # Truncated SVD (up to 50 components; first component dropped below).
        # Small smoke fixtures may have fewer peaks/cells than the production
        # default, so bound components to the available matrix dimensions.
        from sklearn.decomposition import TruncatedSVD

        n_components = min(50, tfidf.shape[1])
        if n_components < 2:
            raise ModuleContractError(
                f"{self.name}: at least 2 ATAC peaks are required for LSI after dropping "
                f"the first component; got {tfidf.shape[1]} peak(s)."
            )
        svd = TruncatedSVD(n_components=n_components, random_state=42)
        embedding = svd.fit_transform(tfidf)

        # Drop first (depth-correlated) component — Signac / ArchR convention
        X_lsi = embedding[:, 1:]

        adata.obsm["X_lsi"] = X_lsi
        variance_explained = svd.explained_variance_ratio_.tolist()
        adata.uns["lsi_variance_explained"] = variance_explained

        logger.info(
            "%s: X_lsi shape=%s, variance explained (all %d comps) top-5=%s",
            self.name,
            X_lsi.shape,
            n_components,
            [round(v, 4) for v in variance_explained[:5]],
        )

        # Variance-explained plot
        self._save_variance_plot(variance_explained, ctx)
        ctx.status(self.name, True, "completed")

    def _save_variance_plot(self, variance_explained: list[float], ctx: PipelineContext) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(8, 4))
            components = list(range(1, len(variance_explained) + 1))
            ax.bar(components, [v * 100 for v in variance_explained], color="steelblue")
            ax.set_xlabel("SVD Component")
            ax.set_ylabel("Variance Explained (%)")
            ax.set_title("LSI Variance Explained (component 1 dropped from X_lsi)")
            ax.axvline(x=1.5, color="red", linestyle="--", linewidth=1, label="drop_first boundary")
            ax.legend(fontsize=8)
            fig.tight_layout()

            # figure_dir is set to the module's output dir by the pipeline before run() is called
            out_dir: Path = ctx.figure_dir
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / "lsi_variance_explained.png"
            fig.savefig(str(out_path))
            plt.close(fig)
            logger.info("%s: variance plot saved to %s", self.name, out_path)
        except Exception as exc:
            logger.warning("%s: could not save variance plot: %s", self.name, exc)
