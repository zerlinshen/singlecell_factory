from __future__ import annotations

import logging
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ._scanpy_compat import import_scanpy_or_stub
from scipy import sparse

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from ._gpu_utils import gpu_available
from .._neighbors_cache import get_or_compute_neighbors


__references__ = {
    "Korsunsky_Harmony_2019": {
        "title": "Fast, sensitive and accurate integration of single-cell data with Harmony",
        "authors": "Korsunsky et al.",
        "journal": "Nature Methods",
        "year": "2019",
        "doi": "10.1038/s41592-019-0619-0",
        "description": "Default backend (rsc.pp.harmony_integrate GPU port, US-B3 parity ARI=1.000).",
    },
    "Polanski_BBKNN_2020": {
        "title": "BBKNN: fast batch alignment of single cell transcriptomes",
        "authors": "Polanski et al.",
        "journal": "Bioinformatics",
        "year": "2020",
        "doi": "10.1093/bioinformatics/btz625",
        "description": "BBKNN backend.",
    },
    "Johnson_ComBat_2007": {
        "title": "Adjusting batch effects in microarray expression data using empirical Bayes methods",
        "authors": "Johnson, Li, Rabinovic",
        "journal": "Biostatistics",
        "year": "2007",
        "doi": "10.1093/biostatistics/kxj037",
        "description": "ComBat backend (sc.pp.combat).",
    },
    "Hie_Scanorama_2019": {
        "title": "Efficient integration of heterogeneous single-cell transcriptomes using Scanorama",
        "authors": "Hie, Bryson, Berger",
        "journal": "Nature Biotechnology",
        "year": "2019",
        "doi": "10.1038/s41587-019-0113-3",
        "description": "Scanorama backend.",
    },
    "Lopez_scVI_2018": {
        "title": "Deep generative modeling for single-cell transcriptomics",
        "authors": "Lopez et al.",
        "journal": "Nature Methods",
        "year": "2018",
        "doi": "10.1038/s41592-018-0229-2",
        "description": "scVI backend via scvi-tools.",
    },
    "Haghverdi_MNN_2018": {
        "title": "Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
        "authors": "Haghverdi et al.",
        "journal": "Nature Biotechnology",
        "year": "2018",
        "doi": "10.1038/nbt.4091",
        "description": "MNN / fastMNN backends.",
    },
}


logger = logging.getLogger(__name__)


def _to_cpu_value(value):
    if sparse.issparse(value):
        data = getattr(value, "data", None)
        if hasattr(data, "get"):
            value.data = data.get()
        return value
    if hasattr(value, "get"):
        return value.get()
    return value


def _ensure_cpu_batch_inputs(adata) -> None:
    adata.X = _to_cpu_value(adata.X)
    for layer_key in list(adata.layers.keys()):
        adata.layers[layer_key] = _to_cpu_value(adata.layers[layer_key])
    for rep_key in ["X_pca", "X_pca_harmony", "X_scanorama", "X_scvi", "X_mnn", "X_fastmnn"]:
        if rep_key in adata.obsm:
            adata.obsm[rep_key] = _to_cpu_value(adata.obsm[rep_key])


class BatchCorrectionModule:
    """Optional module: batch effect correction for multi-sample integration.

    Supports multiple backends:
    1. Harmony (default) — fast, PCA-based correction
    2. BBKNN — batch-balanced k-nearest-neighbors
    3. Combat — linear model-based correction
    4. Scanorama — manifold alignment integration
    5. scVI — variational latent integration

    Should run after clustering (specifically after PCA), before downstream analysis.
    The batch key must be present in adata.obs.
    """

    name = "batch_correction"
    mutates_structure = True
    requires_keys = {"obs": ["leiden"]}
    provides_keys = {"obs": ["leiden"]}

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Batch correction requires AnnData.")

        cfg = ctx.cfg.batch
        # Sync canonical seed into clustering cfg so post-correction UMAP/Leiden are deterministic
        ctx.cfg.clustering.random_state = ctx.random_state
        batch_key = cfg.batch_key

        if batch_key not in adata.obs.columns:
            ctx.metadata["batch_correction_status"] = "skipped_missing_batch_key"
            ctx.status(
                self.name,
                "skipped",
                f"Batch key '{batch_key}' not found in adata.obs.",
            )
            return

        n_batches = adata.obs[batch_key].nunique()
        if n_batches < 2:
            ctx.metadata["batch_correction_status"] = "skipped_single_batch"
            ctx.status(self.name, "skipped", f"Only one batch found in '{batch_key}'.")
            return

        ctx.metadata["n_batches"] = n_batches
        ctx.metadata["batch_method"] = cfg.method
        ctx.metadata["gpu_mode"] = ctx.cfg.gpu_mode

        # Pre-correction UMAP (for comparison)
        if "X_umap" in adata.obsm:
            sc.pl.umap(adata, color=[batch_key], show=False)
            plt.savefig(
                ctx.figure_dir / "umap_batch_before.png", dpi=160, bbox_inches="tight",
            )
            plt.close()

        if cfg.method == "harmony":
            backend = self._run_harmony
        elif cfg.method == "bbknn":
            backend = self._run_bbknn
        elif cfg.method == "combat":
            backend = self._run_combat
        elif cfg.method == "scanorama":
            backend = self._run_scanorama
        elif cfg.method == "scvi":
            backend = self._run_scvi
        elif cfg.method == "mnn":
            backend = self._run_mnn
        elif cfg.method == "fastmnn":
            backend = self._run_fastmnn
        else:
            raise ValueError(f"Unknown batch correction method: {cfg.method}")
        _ensure_cpu_batch_inputs(adata)
        try:
            backend(adata, batch_key, ctx)
        except Exception as exc:
            if cfg.method not in {"scvi", "mnn", "fastmnn"}:
                raise
            msg = f"{cfg.method} backend unavailable or failed: {exc}"
            logger.warning("%s", msg)
            ctx.metadata["batch_correction_status"] = f"skipped_{cfg.method}_unavailable_or_failed"
            ctx.metadata["batch_correction_skip_reason"] = str(exc)
            ctx.status(self.name, "skipped", msg)
            return

        # Recompute UMAP + Leiden on corrected representation
        use_rep = "X_pca"
        if cfg.method == "harmony":
            use_rep = "X_pca_harmony"
        elif cfg.method == "scanorama":
            use_rep = "X_scanorama"
        elif cfg.method == "scvi":
            use_rep = "X_scvi"
        elif cfg.method == "mnn":
            use_rep = "X_mnn"
        elif cfg.method == "fastmnn":
            use_rep = "X_fastmnn"

        use_gpu = gpu_available(ctx.cfg.gpu_mode)
        if use_gpu:
            try:
                import rapids_singlecell as rsc
                logger.info("GPU batch post-processing: neighbors/UMAP/Leiden")
                if cfg.method != "bbknn":
                    # bbknn bypasses cache: incompatible neighbors API signature
                    get_or_compute_neighbors(
                        adata,
                        backend_id="rapids",
                        method="umap",
                        metric="euclidean",
                        n_pcs=None,
                        n_neighbors=ctx.cfg.clustering.n_neighbors,
                        use_rep=use_rep,
                        knn=True,
                        random_state=ctx.cfg.clustering.random_state,
                        compute_fn=lambda: rsc.pp.neighbors(adata, use_rep=use_rep),
                    )
                rsc.tl.umap(adata, random_state=ctx.cfg.clustering.random_state)
                rsc.tl.leiden(adata, resolution=ctx.cfg.clustering.leiden_resolution,
                              random_state=ctx.cfg.clustering.random_state)
            except Exception as exc:
                if ctx.cfg.gpu_mode == "force":
                    raise RuntimeError(f"GPU batch post-processing failed in force mode: {exc}") from exc
                logger.warning("GPU batch post-processing failed (%s), falling back to CPU", exc)
                ctx.metadata["batch_gpu_fallback_reason"] = str(exc)
                use_gpu = False

        if not use_gpu:
            try:
                if cfg.method != "bbknn":
                    # bbknn bypasses cache: incompatible neighbors API signature
                    get_or_compute_neighbors(
                        adata,
                        backend_id="scanpy",
                        method="umap",
                        metric="euclidean",
                        n_pcs=None,
                        n_neighbors=ctx.cfg.clustering.n_neighbors,
                        use_rep=use_rep,
                        knn=True,
                        random_state=ctx.cfg.clustering.random_state,
                        compute_fn=lambda: sc.pp.neighbors(adata, use_rep=use_rep),
                    )
                sc.tl.umap(adata, random_state=ctx.cfg.clustering.random_state)
                sc.tl.leiden(
                    adata,
                    resolution=ctx.cfg.clustering.leiden_resolution,
                    flavor="igraph",
                    directed=False,
                    random_state=ctx.cfg.clustering.random_state,
                )
            except Exception as exc:
                # Keep prior clustering outputs rather than failing the module.
                logger.warning(
                    "CPU batch post-processing failed (%s); keeping existing UMAP/leiden.",
                    exc,
                )
                ctx.metadata["batch_post_error"] = str(exc)
        ctx.metadata["batch_post_backend"] = "gpu" if use_gpu else "cpu"
        ctx.metadata["n_clusters_after_batch"] = int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs else 0
        ctx.metadata["batch_correction_status"] = "completed"
        ctx.metadata["batch_correction_random_state"] = ctx.random_state

        # Post-correction UMAP
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        sc.pl.umap(adata, color=[batch_key], ax=axes[0], show=False)
        axes[0].set_title(f"After {cfg.method} — batch")
        sc.pl.umap(adata, color=["leiden"], ax=axes[1], show=False, legend_loc="on data")
        axes[1].set_title(f"After {cfg.method} — clusters")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "umap_batch_after.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _run_harmony(adata, batch_key: str, ctx) -> None:
        if "X_pca" not in adata.obsm:
            raise ValueError("Harmony requires PCA (run clustering first).")
        try:
            import rapids_singlecell as rsc
        except ImportError:
            raise ImportError(
                "rapids_singlecell is required for Harmony batch correction. "
                "Activate the sc_gpu_stable conda environment."
            )
        # rsc.pp.harmony_integrate writes corrected embeddings to
        # adata.obsm["X_pca_harmony"] in place (US-B3, ARI=1.000 parity gate passed).
        ctx.metadata["harmony_device"] = "gpu"
        rsc.pp.harmony_integrate(adata, key=batch_key)
        # Ensure float32 for downstream consistency
        adata.obsm["X_pca_harmony"] = np.asarray(
            adata.obsm["X_pca_harmony"], dtype=np.float32
        )

    @staticmethod
    def _run_bbknn(adata, batch_key: str, ctx) -> None:
        try:
            import bbknn
        except ImportError:
            raise ImportError(
                "bbknn is required for BBKNN batch correction. "
                "Install with: pip install bbknn"
            )
        if "X_pca" not in adata.obsm:
            raise ValueError("BBKNN requires PCA (run clustering first).")
        bbknn.bbknn(adata, batch_key=batch_key)

    @staticmethod
    def _run_combat(adata, batch_key: str, ctx) -> None:
        sc.pp.combat(adata, key=batch_key)

    @staticmethod
    def _run_scanorama(adata, batch_key: str, ctx) -> None:
        """Scanorama integration (Hie et al., Nature Biotechnology 2019)."""
        try:
            import scanorama
        except ImportError:
            raise ImportError(
                "scanorama is required for Scanorama batch correction. "
                "Install with: pip install scanorama"
            )

        if "X_pca" not in adata.obsm:
            raise ValueError("Scanorama requires PCA (run clustering first).")
        batches = adata.obs[batch_key].unique().tolist()
        adatas = [adata[adata.obs[batch_key] == b].copy() for b in batches]
        scanorama.integrate_scanpy(adatas)
        # Reconstruct corrected embedding into original object
        adata.obsm["X_scanorama"] = np.zeros(
            (adata.n_obs, adatas[0].obsm["X_scanorama"].shape[1])
        )
        for b, sub in zip(batches, adatas):
            mask = adata.obs[batch_key] == b
            adata.obsm["X_scanorama"][mask] = sub.obsm["X_scanorama"]

    @staticmethod
    def _run_mnn(adata, batch_key: str, ctx) -> None:
        """MNN correction via scanpy.external.pp.mnn_correct (Haghverdi et al., 2018)."""
        if "X_pca" not in adata.obsm:
            raise ValueError("MNN requires PCA (run clustering first).")
        try:
            _ = sc.external.pp.mnn_correct
        except (AttributeError, ImportError):
            raise ImportError(
                "mnn_correct is unavailable. Install mnnpy and ensure scanpy.external is available. "
                "Example: pip install mnnpy"
            )

        batches = adata.obs[batch_key].astype(str)
        batch_levels = list(dict.fromkeys(batches.tolist()))
        if len(batch_levels) < 2:
            raise ValueError("MNN requires at least two batches.")
        adatas = [adata[batches == b].copy() for b in batch_levels]

        try:
            corrected = sc.external.pp.mnn_correct(*adatas, batch_key=None, index_unique=None)[0]
        except TypeError:
            # Fallback for older mnnpy/scanpy signatures.
            corrected = sc.external.pp.mnn_correct(*adatas)[0]
        corrected = corrected[adata.obs_names, :].copy()
        sc.pp.pca(corrected, n_comps=ctx.cfg.clustering.n_pcs)
        adata.obsm["X_mnn"] = corrected.obsm["X_pca"].copy()

    @staticmethod
    def _run_fastmnn(adata, batch_key: str, ctx) -> None:
        """Fast-MNN style approximation using mnn_correct with reduced SVD dimension."""
        if "X_pca" not in adata.obsm:
            raise ValueError("fastMNN requires PCA (run clustering first).")
        try:
            _ = sc.external.pp.mnn_correct
        except (AttributeError, ImportError):
            raise ImportError(
                "mnn_correct is unavailable. Install mnnpy and ensure scanpy.external is available. "
                "Example: pip install mnnpy"
            )

        batches = adata.obs[batch_key].astype(str)
        batch_levels = list(dict.fromkeys(batches.tolist()))
        if len(batch_levels) < 2:
            raise ValueError("fastMNN requires at least two batches.")
        adatas = [adata[batches == b].copy() for b in batch_levels]

        svd_dim = max(10, int(ctx.cfg.clustering.n_pcs))
        try:
            corrected = sc.external.pp.mnn_correct(
                *adatas, batch_key=None, index_unique=None, svd_dim=svd_dim
            )[0]
        except TypeError:
            corrected = sc.external.pp.mnn_correct(*adatas)[0]
        corrected = corrected[adata.obs_names, :].copy()
        sc.pp.pca(corrected, n_comps=ctx.cfg.clustering.n_pcs)
        adata.obsm["X_fastmnn"] = corrected.obsm["X_pca"].copy()
        ctx.metadata["fastmnn_backend"] = "mnn_correct_svd_approx"
        ctx.metadata["fastmnn_svd_dim"] = svd_dim

    @staticmethod
    def _matrix_looks_integer_like(matrix, sample_size: int = 2048) -> bool:
        """Heuristic: check if matrix resembles raw integer counts."""
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

    @staticmethod
    def _run_scvi(adata, batch_key: str, ctx) -> None:
        """scVI latent integration (Lopez et al., Nature Methods 2018)."""
        try:
            import scvi
        except ImportError:
            raise ImportError(
                "scvi-tools is required for scVI batch correction. "
                "Install with: pip install scvi-tools"
            )

        batch_cfg = ctx.cfg.batch
        if "counts" in adata.layers:
            layer = "counts"
            input_source = "layers.counts"
        else:
            layer = None
            if not BatchCorrectionModule._matrix_looks_integer_like(adata.X):
                raise ValueError(
                    "scVI requires raw count-like input. 'counts' layer is missing and adata.X "
                    "does not appear to be integer UMI counts."
                )
            input_source = "adata.X"

        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, layer=layer)
        model = scvi.model.SCVI(adata, n_latent=batch_cfg.scvi_n_latent)
        model.train(
            max_epochs=batch_cfg.scvi_max_epochs,
            early_stopping=batch_cfg.scvi_early_stopping,
        )
        adata.obsm["X_scvi"] = model.get_latent_representation()
        ctx.metadata["scvi_input_source"] = input_source
        ctx.metadata["scvi_train_config"] = {
            "max_epochs": batch_cfg.scvi_max_epochs,
            "n_latent": batch_cfg.scvi_n_latent,
            "early_stopping": batch_cfg.scvi_early_stopping,
        }
