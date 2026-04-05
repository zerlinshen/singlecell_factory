from __future__ import annotations

import hashlib
import gzip
import importlib
import logging
import os
import warnings
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
from time import perf_counter
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from ..context import PipelineContext

logger = logging.getLogger(__name__)
_SCVELO_NUMPY2_PATCH_DONE = False
_WORKER_STATE: dict[str, Any] = {}


# ---------------------------------------------------------------------------
# numpy 2.x + scVelo 0.3.x compatibility
# ---------------------------------------------------------------------------


def _patch_scvelo_numpy2() -> bool:
    """Patch scVelo 0.3.x for numpy 2.x compatibility.

    numpy 2.x changed matrix multiplication return shapes, causing
    'setting an array element with a sequence' in
    scvelo.tools.optimization.leastsq_generalized.

    Returns True if the patch was applied.
    """
    global _SCVELO_NUMPY2_PATCH_DONE

    if _SCVELO_NUMPY2_PATCH_DONE:
        return False

    if not _version_at_least(np.__version__, (2, 0)):
        return False

    try:
        import scvelo

        scv_version = getattr(scvelo, "__version__", "")
        # leastsq_generalized assignment bug still occurs in scVelo 0.4.x with numpy 2.x.
        # Keep patch enabled for <0.5, and skip for newer versions unless needed again.
        if _version_at_least(scv_version, (0, 5)):
            return False
    except ImportError:
        logger.debug("scVelo not installed, skipping numpy 2.x patch")
        return False

    try:
        import scvelo.tools.optimization as opt

        def _leastsq_fixed(
            x,
            y,
            x2,
            y2,
            res_std=None,
            res2_std=None,
            fit_offset=False,
            fit_offset2=False,
            perc=None,
        ):
            """numpy2-safe replacement for scvelo.tools.optimization.leastsq_generalized."""
            if perc is not None:
                if not fit_offset and isinstance(perc, (list, tuple)):
                    perc = perc[1]
                weights = opt.csr_matrix(
                    opt.get_weight(x, y, perc=perc) | opt.get_weight(x, perc=perc)
                ).astype(bool)
                x, y = weights.multiply(x).tocsr(), weights.multiply(y).tocsr()

            n_obs, n_var = x.shape
            offset = np.zeros(n_var, dtype="float32")
            offset_ss = np.zeros(n_var, dtype="float32")
            gamma = np.ones(n_var, dtype="float32")

            if (res_std is None) or (res2_std is None):
                res_std, res2_std = np.ones(n_var), np.ones(n_var)
            ones, zeros = np.ones(n_obs), np.zeros(n_obs)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                x, y = (
                    np.vstack((opt.make_dense(x) / res_std, x2 / res2_std)),
                    np.vstack((opt.make_dense(y) / res_std, y2 / res2_std)),
                )

            def _solve(A, b):
                vals = np.asarray(np.linalg.pinv(A.T.dot(A)).dot(A.T.dot(b))).reshape(-1)
                return vals

            if fit_offset and fit_offset2:
                for i in range(n_var):
                    A = np.c_[
                        np.vstack(
                            (
                                np.c_[ones / res_std[i], zeros],
                                np.c_[zeros, ones / res2_std[i]],
                            )
                        ),
                        x[:, i],
                    ]
                    vals = _solve(A, y[:, i])
                    if len(vals) > 0:
                        offset[i] = float(vals[0])
                    if len(vals) > 1:
                        offset_ss[i] = float(vals[1])
                    if len(vals) > 2:
                        gamma[i] = float(vals[2])
            elif fit_offset:
                for i in range(n_var):
                    A = np.c_[np.hstack((ones / res_std[i], zeros)), x[:, i]]
                    vals = _solve(A, y[:, i])
                    if len(vals) > 0:
                        offset[i] = float(vals[0])
                    if len(vals) > 1:
                        gamma[i] = float(vals[1])
            elif fit_offset2:
                for i in range(n_var):
                    A = np.c_[np.hstack((zeros, ones / res2_std[i])), x[:, i]]
                    vals = _solve(A, y[:, i])
                    if len(vals) > 0:
                        offset_ss[i] = float(vals[0])
                    if len(vals) > 1:
                        gamma[i] = float(vals[1])
            else:
                for i in range(n_var):
                    A = np.c_[x[:, i]]
                    vals = _solve(A, y[:, i])
                    if len(vals) > 0:
                        gamma[i] = float(vals[-1])

            offset[np.isnan(offset)] = 0
            offset_ss[np.isnan(offset_ss)] = 0
            gamma[np.isnan(gamma)] = 0
            return offset, offset_ss, gamma

        opt.leastsq_generalized = _leastsq_fixed
        try:
            vel_mod = importlib.import_module("scvelo.tools.velocity")
            vel_mod.leastsq_generalized = _leastsq_fixed
        except (AttributeError, ImportError) as exc:
            logger.warning("Unable to patch scvelo.tools.velocity.leastsq_generalized: %s", exc)
        try:
            ssm_mod = importlib.import_module("scvelo.tools._steady_state_model")
            ssm_mod.leastsq_generalized = _leastsq_fixed
        except (AttributeError, ImportError) as exc:
            logger.warning(
                "Unable to patch scvelo.tools._steady_state_model.leastsq_generalized: %s",
                exc,
            )

        logger.debug("Applied scVelo numpy 2.x compatibility patch")
        _SCVELO_NUMPY2_PATCH_DONE = True
        return True

    except AttributeError as exc:
        logger.warning("scVelo internal API changed, numpy 2.x patch skipped: %s", exc)
        return False


def _version_at_least(version_str: str, target: tuple[int, ...]) -> bool:
    """Compare version strings numerically (safe for e.g. 0.10 vs 0.4)."""
    if not version_str:
        return False
    current = _extract_numeric_version(version_str)
    if not current:
        return False
    pad = max(len(current), len(target))
    current = current + (0,) * (pad - len(current))
    goal = target + (0,) * (pad - len(target))
    return current >= goal


def _extract_numeric_version(version_str: str) -> tuple[int, ...]:
    parts: list[int] = []
    for chunk in version_str.split("."):
        digits = []
        for ch in chunk:
            if ch.isdigit():
                digits.append(ch)
            else:
                break
        if not digits:
            break
        parts.append(int("".join(digits)))
    return tuple(parts)


# ---------------------------------------------------------------------------
# Transfer helper
# ---------------------------------------------------------------------------


def _transfer_velocity_results(adata_v, adata) -> None:
    """Transfer velocity results from the scVelo working copy back to the main adata.

    Only transfers cell-indexed data (obs columns, obsm entries) and
    scalar/dict uns entries. Does NOT transfer layers or var columns
    because the gene index on adata_v differs after filter_and_normalize.
    """
    velocity_obs_cols = [
        "velocity_confidence",
        "velocity_length",
        "velocity_self_transition",
        "latent_time",
        "root_cells",
        "end_points",
    ]
    for col in velocity_obs_cols:
        if col in adata_v.obs.columns:
            adata.obs[col] = adata_v.obs[col].values

    if "velocity_umap" in adata_v.obsm:
        adata.obsm["velocity_umap"] = adata_v.obsm["velocity_umap"]

    for key in ("velocity_params",):
        if key in adata_v.uns:
            adata.uns[key] = adata_v.uns[key]


# ---------------------------------------------------------------------------
# Module
# ---------------------------------------------------------------------------


class RNAVelocityModule:
    """Optional module: RNA velocity analysis using scVelo.

    Runs all scVelo operations on an internal adata copy to avoid mutating
    the shared gene index, then transfers cell-level results back.

    Supports three data sources (tried in order):
    1. Pre-computed loom file (--velocity-loom)
    2. BAM + GTF extraction via pysam (--velocity-bam + --velocity-gtf)
    3. Falls back with actionable error if no data available.

    References:
    - Bergen et al., Nature Biotechnology, 2020 (scVelo)
    - La Manno et al., Nature, 2018 (RNA velocity concept)
    """

    name = "rna_velocity"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "X_umap" not in adata.obsm:
            raise ValueError("RNA velocity requires AnnData with UMAP embedding.")

        import scvelo as scv

        # Apply numpy 2.x patch BEFORE any scVelo operations
        _patch_scvelo_numpy2()

        cfg = ctx.cfg.velocity
        n_vars_original = adata.n_vars
        step_times: dict[str, float] = {}

        # Work on a copy to preserve the shared gene index
        adata_v = adata.copy()

        # --- Data loading ---
        t0 = perf_counter()
        self._load_layers(adata_v, ctx)
        step_times["load_layers"] = perf_counter() - t0

        # --- scVelo preprocessing (mutates gene index -- only on the copy) ---
        t0 = perf_counter()
        scv.pp.filter_and_normalize(adata_v, min_shared_counts=cfg.min_shared_counts)
        step_times["filter_and_normalize"] = perf_counter() - t0
        if adata_v.n_vars == 0:
            raise ValueError(
                f"filter_and_normalize removed all genes "
                f"(min_shared_counts={cfg.min_shared_counts}). "
                "Lower --velocity-min-shared-counts or check spliced/unspliced layer quality."
            )
        logger.info(
            "After filter_and_normalize: %d genes remain (from %d)",
            adata_v.n_vars,
            n_vars_original,
        )

        t0 = perf_counter()
        scv.pp.moments(adata_v, n_pcs=cfg.n_pcs, n_neighbors=cfg.n_neighbors)
        step_times["moments"] = perf_counter() - t0

        # --- Velocity computation ---
        t0 = perf_counter()
        self._compute_velocity(adata_v, cfg, scv)
        step_times["compute_velocity"] = perf_counter() - t0

        # --- Metrics ---
        t0 = perf_counter()
        scv.tl.velocity_confidence(adata_v)
        step_times["velocity_confidence"] = perf_counter() - t0

        # --- Transfer cell-level results to the original adata ---
        _transfer_velocity_results(adata_v, adata)

        # --- Tables ---
        t0 = perf_counter()
        self._write_tables(adata_v, ctx, scv)
        step_times["write_tables"] = perf_counter() - t0

        # --- Visualizations (use adata_v which has velocity data) ---
        t0 = perf_counter()
        self._visualize(adata_v, ctx, cfg, scv)
        step_times["visualize"] = perf_counter() - t0

        # --- Metadata ---
        if "velocity_confidence" in adata_v.obs:
            ctx.metadata["velocity_mean_confidence"] = float(
                adata_v.obs["velocity_confidence"].mean()
            )
        ctx.metadata["velocity_mode"] = cfg.mode
        ctx.metadata["velocity_step_seconds"] = {
            k: round(v, 3) for k, v in step_times.items()
        }

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_gtf_path(ctx: PipelineContext) -> Path | None:
        """Resolve GTF path for BAM extraction.

        Priority:
        1) explicit --velocity-gtf
        2) auto-discover from --transcriptome-dir (Cell Ranger reference)
        """
        vcfg = ctx.cfg.velocity
        ccfg = ctx.cfg.cellranger

        if vcfg.gtf_path:
            if vcfg.gtf_path.exists():
                return vcfg.gtf_path
            logger.warning("velocity gtf file not found: %s", vcfg.gtf_path)
            return None

        ref_dir = ccfg.transcriptome_dir
        if ref_dir is None:
            return None

        candidates = [
            ref_dir / "genes" / "genes.gtf",
            ref_dir / "genes" / "genes.gtf.gz",
            ref_dir / "genes.gtf",
            ref_dir / "genes.gtf.gz",
        ]
        for path in candidates:
            if path.exists():
                logger.info("Auto-detected velocity GTF from transcriptome dir: %s", path)
                return path
        return None

    @classmethod
    def _load_layers(cls, adata_v, ctx: PipelineContext) -> None:
        """Load spliced/unspliced layers from loom or BAM."""
        import scvelo as scv

        cfg = ctx.cfg.velocity
        if cfg.loom_path and cfg.loom_path.exists():
            ldata = scv.read(str(cfg.loom_path))
            merged = scv.utils.merge(adata_v, ldata)
            # scv.utils.merge returns a new object; propagate back
            adata_v.__dict__.update(merged.__dict__)
        elif cfg.loom_path:
            logger.warning("velocity loom file not found: %s", cfg.loom_path)

        has_layers = "spliced" in adata_v.layers and "unspliced" in adata_v.layers

        resolved_gtf = cls._resolve_gtf_path(ctx)
        if resolved_gtf is not None:
            ctx.metadata["velocity_gtf_resolved"] = str(resolved_gtf)

        if not has_layers and cfg.bam_path:
            if cfg.bam_path.exists() and resolved_gtf and resolved_gtf.exists():
                logger.info("Extracting spliced/unspliced counts from BAM + GTF...")
                extract_jobs = max(1, int(cfg.n_jobs))
                # Default cfg.n_jobs=4 is conservative; auto-bump extraction workers
                # on high-core hosts while preserving explicit user overrides.
                if extract_jobs == 4 and (os.cpu_count() or 1) >= 16:
                    extract_jobs = 8
                ctx.metadata["velocity_extract_n_jobs"] = int(extract_jobs)
                t0 = perf_counter()
                spliced, unspliced = _count_spliced_unspliced(
                    cfg.bam_path,
                    resolved_gtf,
                    adata_v.obs_names.tolist(),
                    adata_v.var_names.tolist(),
                    n_jobs=extract_jobs,
                )
                ctx.metadata["velocity_extract_seconds"] = round(perf_counter() - t0, 3)
                adata_v.layers["spliced"] = spliced
                adata_v.layers["unspliced"] = unspliced
                has_layers = True
                logger.info(
                    "Extracted layers: spliced nnz=%d, unspliced nnz=%d",
                    spliced.nnz,
                    unspliced.nnz,
                )
            else:
                logger.warning(
                    "velocity BAM/GTF path missing (bam=%s exists=%s, gtf=%s exists=%s)",
                    cfg.bam_path,
                    cfg.bam_path.exists(),
                    resolved_gtf,
                    bool(resolved_gtf and resolved_gtf.exists()),
                )

        if not has_layers:
            auto_gtf_note = (
                "Auto-detection from --transcriptome-dir failed.\n"
                "  Expected one of:\n"
                "    <transcriptome-dir>/genes/genes.gtf(.gz)\n"
                "    <transcriptome-dir>/genes.gtf(.gz)\n"
            )
            raise ValueError(
                "No spliced/unspliced layers found. RNA velocity requires these layers.\n"
                "Options:\n"
                "  1. Provide a loom file:  --velocity-loom /path/to/file.loom\n"
                "  2. Provide BAM + GTF:    --velocity-bam possorted_genome_bam.bam "
                "--velocity-gtf genes.gtf.gz\n"
                "  2b. Or set --transcriptome-dir to auto-discover genes.gtf(.gz)\n"
                "  3. Pre-process with velocyto:  velocyto run10x <sample_dir> <gtf>\n"
                "  4. Use pseudo_velocity instead (works from pseudotime, no splicing data needed)\n"
                + auto_gtf_note
            )

    # ------------------------------------------------------------------
    # Velocity computation
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_velocity(adata_v, cfg, scv) -> None:
        """Run scVelo velocity model."""
        if cfg.mode == "dynamical":
            scv.tl.recover_dynamics(adata_v, n_jobs=cfg.n_jobs)
            scv.tl.velocity(adata_v, mode="dynamical")
            if hasattr(scv.tl, "latent_time"):
                try:
                    scv.tl.latent_time(adata_v)
                except Exception as exc:
                    logger.warning("Latent time computation failed: %s", exc)
        else:
            scv.tl.velocity(adata_v, mode="stochastic")

        graph_kwargs = {
            "n_jobs": max(1, int(cfg.n_jobs)),
            "backend": "threading",
            "show_progress_bar": False,
        }

        def _run_velocity_graph():
            try:
                scv.tl.velocity_graph(adata_v, **graph_kwargs)
            except TypeError:
                # Older/mocked signatures may not accept backend/show_progress_bar or n_jobs.
                try:
                    scv.tl.velocity_graph(adata_v, n_jobs=graph_kwargs["n_jobs"])
                except TypeError:
                    scv.tl.velocity_graph(adata_v)

        try:
            _run_velocity_graph()
        except ValueError:
            # scVelo 0.3.x + numpy 2.x: force float64 on velocity layer
            if "velocity" in adata_v.layers:
                adata_v.layers["velocity"] = np.asarray(
                    adata_v.layers["velocity"], dtype=np.float64
                )
            _run_velocity_graph()

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------

    @staticmethod
    def _write_tables(adata_v, ctx: PipelineContext, scv) -> None:
        """Write velocity result tables."""
        if "velocity_confidence" in adata_v.obs:
            cols = ["velocity_confidence", "velocity_length"]
            if "latent_time" in adata_v.obs:
                cols.append("latent_time")
            adata_v.obs[cols].to_csv(ctx.table_dir / "velocity_confidence.csv")

        try:
            groupby = "leiden" if "leiden" in adata_v.obs else None
            scv.tl.rank_velocity_genes(adata_v, groupby=groupby, n_genes=10)
            velocity_genes = pd.DataFrame(adata_v.uns["rank_velocity_genes"]["names"]).head(10)
            velocity_genes.to_csv(ctx.table_dir / "velocity_top_genes.csv", index=False)
        except Exception as exc:
            logger.warning("rank_velocity_genes failed: %s", exc)

    # ------------------------------------------------------------------
    # Visualizations
    # ------------------------------------------------------------------

    def _visualize(self, adata_v, ctx: PipelineContext, cfg, scv) -> None:
        """Generate all velocity plots."""
        color_col = "cell_type" if "cell_type" in adata_v.obs else (
            "leiden" if "leiden" in adata_v.obs else None
        )

        self._plot_stream(adata_v, ctx, scv, color_col)
        self._plot_grid(adata_v, ctx, scv, color_col)
        self._plot_velocity_length_distribution(adata_v, ctx)

        if cfg.mode == "dynamical":
            self._plot_latent_time(adata_v, ctx)
            self._plot_phase_portraits(adata_v, ctx, scv)

    @staticmethod
    def _plot_stream(adata_v, ctx: PipelineContext, scv, color_col) -> None:
        """Stream plot (main velocity visualization)."""
        try:
            scv.pl.velocity_embedding_stream(
                adata_v, basis="umap", color=color_col, show=False,
            )
            fig = plt.gcf()
            ctx.save_figure(fig, ctx.figure_dir / "velocity_stream_umap.png", dpi=160)
        except Exception as exc:
            logger.warning("velocity_embedding_stream failed: %s", exc)
            plt.close("all")
            try:
                scv.pl.velocity_embedding(
                    adata_v, basis="umap", color=color_col,
                    arrow_length=3, arrow_size=2, show=False,
                )
                fig = plt.gcf()
                ctx.save_figure(fig, ctx.figure_dir / "velocity_stream_umap.png", dpi=160)
            except Exception as exc2:
                logger.warning("velocity_embedding fallback also failed: %s", exc2)
                plt.close("all")

    @staticmethod
    def _plot_grid(adata_v, ctx: PipelineContext, scv, color_col) -> None:
        """Grid-based velocity field."""
        try:
            scv.pl.velocity_embedding_grid(
                adata_v, basis="umap", color=color_col, show=False,
            )
            fig = plt.gcf()
            ctx.save_figure(fig, ctx.figure_dir / "velocity_grid_umap.png", dpi=160)
        except Exception as exc:
            logger.warning("velocity_embedding_grid failed: %s", exc)
            plt.close("all")

    @staticmethod
    def _plot_latent_time(adata_v, ctx: PipelineContext) -> None:
        """UMAP colored by latent time (dynamical mode only)."""
        if "latent_time" not in adata_v.obs:
            return
        fig, ax = plt.subplots(figsize=(8, 7))
        sc = ax.scatter(
            adata_v.obsm["X_umap"][:, 0],
            adata_v.obsm["X_umap"][:, 1],
            c=adata_v.obs["latent_time"],
            cmap="viridis",
            s=3,
            alpha=0.7,
        )
        plt.colorbar(sc, ax=ax, label="Latent Time")
        ax.set_title("scVelo Latent Time (Dynamical Mode)")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        plt.tight_layout()
        ctx.save_figure(fig, ctx.figure_dir / "velocity_latent_time_umap.png", dpi=160)

    @staticmethod
    def _plot_phase_portraits(adata_v, ctx: PipelineContext, scv) -> None:
        """Phase portraits for top velocity genes (dynamical mode)."""
        try:
            top_genes = None
            if "fit_likelihood" in adata_v.var:
                top_genes = adata_v.var_names[
                    np.argsort(adata_v.var["fit_likelihood"].values)[-6:]
                ].tolist()
            elif "rank_velocity_genes" in adata_v.uns:
                names = pd.DataFrame(adata_v.uns["rank_velocity_genes"]["names"])
                top_genes = names.iloc[0, :min(6, names.shape[1])].tolist()

            if not top_genes:
                return

            scv.pl.scatter(adata_v, basis=top_genes, ncols=3, show=False)
            fig = plt.gcf()
            ctx.save_figure(fig, ctx.figure_dir / "velocity_phase_portraits.png", dpi=160)
        except Exception as exc:
            logger.warning("Phase portrait generation failed: %s", exc)
            plt.close("all")

    @staticmethod
    def _plot_velocity_length_distribution(adata_v, ctx: PipelineContext) -> None:
        """Histogram of velocity length as QC diagnostic."""
        if "velocity_length" not in adata_v.obs:
            return
        vlen = adata_v.obs["velocity_length"].dropna()
        if len(vlen) < 10:
            return
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(vlen, bins=50, color="steelblue", alpha=0.7, edgecolor="white")
        ax.axvline(
            vlen.median(), color="red", linestyle="--",
            label=f"Median: {vlen.median():.3f}",
        )
        ax.set_xlabel("Velocity Length")
        ax.set_ylabel("Number of Cells")
        ax.set_title("RNA Velocity Length Distribution")
        ax.legend()
        plt.tight_layout()
        ctx.save_figure(fig, ctx.figure_dir / "velocity_length_distribution.png", dpi=160)


# ---------------------------------------------------------------------------
# BAM-based spliced/unspliced counting
# ---------------------------------------------------------------------------


def _parse_gtf_exons(gtf_path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """Parse GTF to extract exon intervals per gene symbol.

    Returns: {gene_symbol: [(chrom, start, end), ...]}
    """
    exons: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GTF is 1-based, convert to 0-based
            end = int(parts[4])
            attrs = parts[8]
            gene_name = _extract_attr(attrs, "gene_name")
            if gene_name:
                exons[gene_name].append((chrom, start, end))

    return dict(exons)


def _extract_attr(attr_str: str, key: str) -> str | None:
    """Extract attribute value from GTF attribute string."""
    for part in attr_str.split(";"):
        part = part.strip()
        if part.startswith(key):
            val = part.split('"')
            if len(val) >= 2:
                return val[1]
    return None


def _build_exon_lookup(
    exons: dict[str, list[tuple[str, int, int]]],
) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """Build {chrom: {gene: [(start, end), ...]}} for fast lookup."""
    lookup: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for gene, intervals in exons.items():
        for chrom, start, end in intervals:
            lookup[chrom][gene].append((start, end))
    for chrom in lookup:
        for gene in lookup[chrom]:
            lookup[chrom][gene].sort()
    return {chrom: dict(genes) for chrom, genes in lookup.items()}


def _is_exonic(pos: int, end: int, gene_exons: list[tuple[int, int]]) -> bool:
    """Check if a read region [pos, end) is fully within exons."""
    for estart, eend in gene_exons:
        if pos >= estart and end <= eend:
            return True
    return False


def _init_count_worker(
    bam_path: Path,
    exon_lookup: dict[str, dict[str, list[tuple[int, int]]]],
    barcode_idx: dict[str, int],
    gene_idx: dict[str, int],
    min_mapq: int,
) -> None:
    """Initialize worker state once to avoid per-task pickling overhead."""
    _WORKER_STATE["bam_path"] = str(bam_path)
    _WORKER_STATE["exon_lookup"] = exon_lookup
    _WORKER_STATE["barcode_idx"] = barcode_idx
    _WORKER_STATE["gene_idx"] = gene_idx
    _WORKER_STATE["min_mapq"] = int(min_mapq)


def _count_chromosome(chrom: str) -> tuple[list[int], list[int], list[int], list[int], int, int]:
    """Count spliced/unspliced reads for a single chromosome.

    Returns (spliced_rows, spliced_cols, unspliced_rows, unspliced_cols,
             n_reads, n_counted).
    """
    import pysam

    if not _WORKER_STATE:
        raise RuntimeError("RNA velocity BAM worker state not initialized")

    bam_path = _WORKER_STATE["bam_path"]
    exon_lookup = _WORKER_STATE["exon_lookup"]
    barcode_idx = _WORKER_STATE["barcode_idx"]
    gene_idx = _WORKER_STATE["gene_idx"]
    min_mapq = _WORKER_STATE["min_mapq"]

    chrom_exons = exon_lookup.get(chrom, {})
    spliced_rows: list[int] = []
    spliced_cols: list[int] = []
    unspliced_rows: list[int] = []
    unspliced_cols: list[int] = []
    n_reads = 0
    n_counted = 0

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for read in bam.fetch(contig=chrom):
        n_reads += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue

        try:
            cb = read.get_tag("CB")
            gene = read.get_tag("GN")
        except KeyError:
            continue

        ci = barcode_idx.get(cb)
        if ci is None:
            continue
        gi = gene_idx.get(gene)
        if gi is None:
            continue

        cigartuples = read.cigartuples
        has_splice = False
        if cigartuples:
            for op, _ in cigartuples:
                if op == 3:  # CIGAR N
                    has_splice = True
                    break

        read_end = read.reference_end
        if read_end is None:
            continue

        if has_splice:
            spliced_rows.append(ci)
            spliced_cols.append(gi)
        else:
            gene_exon_list = chrom_exons.get(gene)
            if gene_exon_list and _is_exonic(
                read.reference_start, read_end, gene_exon_list
            ):
                spliced_rows.append(ci)
                spliced_cols.append(gi)
            else:
                unspliced_rows.append(ci)
                unspliced_cols.append(gi)
        n_counted += 1

    bam.close()
    return spliced_rows, spliced_cols, unspliced_rows, unspliced_cols, n_reads, n_counted


def _velocity_cache_dir() -> Path:
    """Return cache directory for extracted spliced/unspliced layers."""
    root = os.environ.get("SCF_VELOCITY_CACHE_DIR", "/tmp/singlecell_factory_velocity_cache")
    path = Path(root)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _velocity_cache_key(
    bam_path: Path,
    gtf_path: Path,
    barcodes: list[str],
    gene_names: list[str],
    min_mapq: int,
) -> str:
    """Build a stable content key for velocity extraction cache."""
    h = hashlib.blake2b(digest_size=20)
    for p in (bam_path, gtf_path):
        st = p.stat()
        h.update(str(p.resolve()).encode("utf-8"))
        h.update(str(st.st_size).encode("ascii"))
        h.update(str(st.st_mtime_ns).encode("ascii"))
    h.update(str(min_mapq).encode("ascii"))
    h.update(b"classifier=strict_exon_v1")
    h.update(str(len(barcodes)).encode("ascii"))
    h.update(str(len(gene_names)).encode("ascii"))
    for bc in barcodes:
        h.update(bc.encode("utf-8"))
        h.update(b"\0")
    for gene in gene_names:
        h.update(gene.encode("utf-8"))
        h.update(b"\0")
    return h.hexdigest()


def _velocity_cache_paths(cache_key: str) -> tuple[Path, Path]:
    cache_dir = _velocity_cache_dir()
    return (
        cache_dir / f"{cache_key}.spliced.npz",
        cache_dir / f"{cache_key}.unspliced.npz",
    )


def _count_spliced_unspliced(
    bam_path: Path,
    gtf_path: Path,
    barcodes: list[str],
    gene_names: list[str],
    min_mapq: int = 255,
    n_jobs: int = 4,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Count spliced and unspliced reads from BAM using GTF exon annotations.

    Strategy:
    - A read is 'spliced' if its CIGAR contains N (splice junction) OR
      it maps entirely within exonic regions.
    - A read is 'unspliced' if it overlaps intronic regions (no splice
      junction and not fully exonic).
    - Uses CB (cell barcode) and GN (gene name) tags from Cell Ranger BAM.
    - Parallelized by chromosome using multiprocessing.
    """
    import pysam

    logger.info("RNA velocity BAM classification: strict exon/intron mode")

    cache_key = _velocity_cache_key(
        bam_path=bam_path,
        gtf_path=gtf_path,
        barcodes=barcodes,
        gene_names=gene_names,
        min_mapq=min_mapq,
    )
    sp_cache, un_cache = _velocity_cache_paths(cache_key)
    if sp_cache.exists() and un_cache.exists():
        try:
            spliced_cached = sparse.load_npz(sp_cache).tocsr()
            unspliced_cached = sparse.load_npz(un_cache).tocsr()
            if spliced_cached.shape == (len(barcodes), len(gene_names)) and unspliced_cached.shape == (len(barcodes), len(gene_names)):
                logger.info("Loaded cached velocity layers from %s", sp_cache.parent)
                return spliced_cached, unspliced_cached
        except Exception as exc:
            logger.warning("Velocity cache load failed, recomputing: %s", exc)

    logger.info("Parsing GTF exon annotations...")
    exons = _parse_gtf_exons(gtf_path)
    exon_lookup = _build_exon_lookup(exons)

    barcode_idx = {bc: i for i, bc in enumerate(barcodes)}
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    n_cells = len(barcodes)
    n_genes = len(gene_names)

    # Get chromosomes from BAM header
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    chromosomes = [ref for ref in bam.references if ref in exon_lookup]
    bam.close()

    if not chromosomes:
        logger.warning("No overlapping chromosomes between BAM and GTF")
        return (
            sparse.csr_matrix((n_cells, n_genes), dtype=np.float32),
            sparse.csr_matrix((n_cells, n_genes), dtype=np.float32),
        )

    logger.info(
        "Scanning BAM across %d chromosomes with %d workers...",
        len(chromosomes),
        min(n_jobs, len(chromosomes)),
    )

    all_sr: list[int] = []
    all_sc: list[int] = []
    all_ur: list[int] = []
    all_uc: list[int] = []
    total_reads = 0
    total_counted = 0

    n_workers = max(1, min(max(1, n_jobs), len(chromosomes)))
    init_args = (
        bam_path,
        exon_lookup,
        barcode_idx,
        gene_idx,
        min_mapq,
    )
    if n_workers <= 1:
        # Single-worker: avoid multiprocessing overhead
        _init_count_worker(*init_args)
        try:
            for chrom in chromosomes:
                sr, sc_, ur, uc, nr, nc = _count_chromosome(chrom)
                all_sr.extend(sr)
                all_sc.extend(sc_)
                all_ur.extend(ur)
                all_uc.extend(uc)
                total_reads += nr
                total_counted += nc
        finally:
            _WORKER_STATE.clear()
    else:
        with Pool(processes=n_workers, initializer=_init_count_worker, initargs=init_args) as pool:
            for sr, sc_, ur, uc, nr, nc in pool.imap_unordered(_count_chromosome, chromosomes):
                all_sr.extend(sr)
                all_sc.extend(sc_)
                all_ur.extend(ur)
                all_uc.extend(uc)
                total_reads += nr
                total_counted += nc

    logger.info("BAM scan complete: %dM reads, %d counted.", total_reads // 1_000_000, total_counted)

    # Build sparse matrices -- COO sums duplicate entries on .tocsr()
    if all_sr:
        spliced = sparse.coo_matrix(
            (np.ones(len(all_sr), dtype=np.float32),
             (np.array(all_sr, dtype=np.int32), np.array(all_sc, dtype=np.int32))),
            shape=(n_cells, n_genes),
        ).tocsr()
    else:
        spliced = sparse.csr_matrix((n_cells, n_genes), dtype=np.float32)

    if all_ur:
        unspliced = sparse.coo_matrix(
            (np.ones(len(all_ur), dtype=np.float32),
             (np.array(all_ur, dtype=np.int32), np.array(all_uc, dtype=np.int32))),
            shape=(n_cells, n_genes),
        ).tocsr()
    else:
        unspliced = sparse.csr_matrix((n_cells, n_genes), dtype=np.float32)

    try:
        sparse.save_npz(sp_cache, spliced)
        sparse.save_npz(un_cache, unspliced)
    except Exception as exc:
        logger.warning("Velocity cache save failed: %s", exc)

    return spliced, unspliced
