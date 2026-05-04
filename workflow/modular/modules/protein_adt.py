"""CITE-seq / ADT proteomics module.

First real consumer of the bundle v2.1 extension API
(``scripts.export_singlecell_r_bundle.add_extension``).

Inputs
------
The module reads ADT counts from one of two locations on the active
``ctx.adata`` (checked in this order):

1. ``adata.obsm["protein_counts"]`` — a (cells x proteins) array of integer
   counts (sparse or dense). Protein names are taken from
   ``adata.uns["protein_names"]`` if present, otherwise auto-generated as
   ``ADT_1..ADT_K``.
2. ``adata.uns["protein_modality"]`` — a dict produced by 10x-style
   feature_bc_matrix loaders, with keys ``counts`` ((cells x proteins) array)
   and ``names`` (list of protein IDs, one per column).

If neither input is present the module records ``protein_status="skipped_no_input"``
in ctx.metadata and returns without raising — proteomics is a fail-soft modality
in this pipeline.

Outputs
-------
- ``adata.obsm["protein_clr"]`` — dense (cells x proteins) float32 CLR matrix.
- ``adata.obs["protein_total_counts"]`` — per-cell total ADT counts.
- ``adata.obs["protein_n_detected"]`` — per-cell number of proteins with > 0 counts.
- ``adata.obs["X_wnn_prep"]`` — boolean placeholder for a future multimodal
  module to detect cells with usable ADT data.
- ``adata.uns["protein_qc"]`` — dict with per-protein detection rate, isotype
  flags, panel size, normalization name.
- A QC PNG ``protein_detection_rate_violin.png`` written to the module's
  output dir.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

from ..context import PipelineContext

logger = logging.getLogger(__name__)


__references__ = {
    "stoeckius_citeseq_2017": {
        "title": "Simultaneous epitope and transcriptome measurement in single cells",
        "authors": "Stoeckius et al.",
        "journal": "Nature Methods",
        "year": "2017",
        "doi": "10.1038/nmeth.4380",
        "description": "CITE-seq protocol; CLR normalization for ADT counts.",
    },
    "mule_dsb_2022": {
        "title": "Normalizing and denoising protein expression data from droplet-based single cell profiling",
        "authors": "Mulè et al.",
        "journal": "Nature Communications",
        "year": "2022",
        "doi": "10.1038/s41467-022-29356-8",
        "description": "DSB normalization (stretch-goal stub).",
    },
}


@dataclass(frozen=True)
class ProteinADTConfig:
    """Local config for the ADT module. Kept module-local per project rules."""

    isotype_controls: tuple[str, ...] = ()
    normalization: str = "CLR"  # "CLR" | "DSB"
    obsm_input_key: str = "protein_counts"
    obsm_output_key: str = "protein_clr"
    pseudocount: float = 1.0
    detection_min_count: int = 1


class ProteinADTModule:
    """Optional module: CITE-seq / ADT QC + CLR normalization.

    Registers an ``adata.obsm["protein_clr"]`` matrix and per-cell QC columns,
    then (when the bundle exporter runs with ``include_protein=True``) the
    matrix is published as ``protein.parquet`` under the v2.1 extension API.
    """

    name = "protein_adt"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {
        "obsm": ["protein_clr"],
        "obs": ["protein_total_counts", "protein_n_detected", "X_wnn_prep"],
        "uns": ["protein_qc"],
    }

    def __init__(self, config: ProteinADTConfig | None = None) -> None:
        self.config = config or ProteinADTConfig()

    # ------------------------------------------------------------------
    # Public entrypoint
    # ------------------------------------------------------------------

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        cfg = self._merge_runtime_config(ctx)

        counts, protein_names = self._extract_protein_counts(adata, cfg)
        if counts is None:
            logger.info("%s: no protein modality found; skipping.", self.name)
            ctx.metadata["protein_status"] = "skipped_no_input"
            ctx.status(self.name, "skipped", "no protein modality on adata")
            return

        if cfg.normalization.upper() == "DSB":
            raise NotImplementedError(
                "DSB normalization requires raw + empty droplets; route via R "
                "subprocess in a later phase"
            )
        if cfg.normalization.upper() != "CLR":
            raise ValueError(
                f"{self.name}: unsupported normalization '{cfg.normalization}'; "
                "expected 'CLR' or 'DSB'."
            )

        n_cells, n_proteins = counts.shape
        logger.info(
            "%s: ADT panel size n_cells=%d n_proteins=%d isotypes=%d",
            self.name, n_cells, n_proteins, len(cfg.isotype_controls),
        )

        total_counts = self._row_totals(counts)
        n_detected = self._row_detected(counts, cfg.detection_min_count)
        detection_rate = self._col_detection_rate(counts, cfg.detection_min_count)

        clr = self._clr_normalize(counts, pseudocount=cfg.pseudocount)

        adata.obsm[cfg.obsm_output_key] = clr.astype(np.float32, copy=False)
        adata.obs["protein_total_counts"] = total_counts
        adata.obs["protein_n_detected"] = n_detected
        # Multimodal flag: cells with at least one detected protein are WNN-ready.
        adata.obs["X_wnn_prep"] = (n_detected > 0).astype(bool)

        isotype_set = set(cfg.isotype_controls)
        isotype_flags = [name in isotype_set for name in protein_names]
        adata.uns["protein_qc"] = {
            "panel_size": int(n_proteins),
            "n_cells": int(n_cells),
            "normalization": "CLR",
            "isotype_controls": list(cfg.isotype_controls),
            "isotype_flags": isotype_flags,
            "protein_names": list(protein_names),
            "detection_rate": [float(x) for x in detection_rate],
            "pseudocount": float(cfg.pseudocount),
        }

        ctx.metadata["protein_n_proteins"] = int(n_proteins)
        ctx.metadata["protein_n_isotype_controls"] = int(sum(isotype_flags))
        ctx.metadata["protein_normalization"] = "CLR"
        ctx.metadata["protein_status"] = "ok"

        # Per-protein detection-rate violin (binned).
        try:
            self._plot_detection_rate(detection_rate, isotype_flags, ctx)
        except Exception as exc:  # plotting is non-fatal
            logger.warning("%s: detection-rate plot failed: %s", self.name, exc)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _merge_runtime_config(self, ctx: PipelineContext) -> ProteinADTConfig:
        """Allow runtime overrides via ctx.cfg.protein_adt without requiring it."""
        runtime = getattr(ctx.cfg, "protein_adt", None)
        if runtime is None:
            return self.config
        # Accept either a ProteinADTConfig or a plain mapping (forward-compat).
        if isinstance(runtime, ProteinADTConfig):
            return runtime
        if isinstance(runtime, dict):
            base = self.config
            return ProteinADTConfig(
                isotype_controls=tuple(runtime.get("isotype_controls", base.isotype_controls)),
                normalization=str(runtime.get("normalization", base.normalization)),
                obsm_input_key=str(runtime.get("obsm_input_key", base.obsm_input_key)),
                obsm_output_key=str(runtime.get("obsm_output_key", base.obsm_output_key)),
                pseudocount=float(runtime.get("pseudocount", base.pseudocount)),
                detection_min_count=int(runtime.get("detection_min_count", base.detection_min_count)),
            )
        return self.config

    def _extract_protein_counts(
        self, adata, cfg: ProteinADTConfig,
    ) -> tuple[np.ndarray | sparse.spmatrix | None, list[str]]:
        """Locate the ADT count matrix and protein names on the adata."""
        # Path 1: obsm key.
        key = cfg.obsm_input_key
        if key in adata.obsm:
            counts = adata.obsm[key]
            n_cells, n_proteins = counts.shape
            names = adata.uns.get("protein_names")
            if names is None or len(names) != n_proteins:
                names = [f"ADT_{i + 1}" for i in range(n_proteins)]
            return counts, list(map(str, names))

        # Path 2: 10x-style feature_bc_matrix dict in uns.
        modality = adata.uns.get("protein_modality")
        if isinstance(modality, dict) and "counts" in modality:
            counts = modality["counts"]
            names = modality.get("names")
            n_proteins = counts.shape[1]
            if names is None or len(names) != n_proteins:
                names = [f"ADT_{i + 1}" for i in range(n_proteins)]
            return counts, list(map(str, names))

        return None, []

    @staticmethod
    def _row_totals(counts) -> np.ndarray:
        if sparse.issparse(counts):
            return np.asarray(counts.sum(axis=1)).ravel().astype(np.float32)
        return np.asarray(counts).sum(axis=1).astype(np.float32)

    @staticmethod
    def _row_detected(counts, min_count: int) -> np.ndarray:
        if sparse.issparse(counts):
            return np.asarray((counts >= min_count).sum(axis=1)).ravel().astype(np.int32)
        return (np.asarray(counts) >= min_count).sum(axis=1).astype(np.int32)

    @staticmethod
    def _col_detection_rate(counts, min_count: int) -> np.ndarray:
        n_cells = counts.shape[0]
        if n_cells == 0:
            return np.zeros(counts.shape[1], dtype=np.float32)
        if sparse.issparse(counts):
            return (
                np.asarray((counts >= min_count).sum(axis=0)).ravel().astype(np.float32)
                / float(n_cells)
            )
        return (
            (np.asarray(counts) >= min_count).sum(axis=0).astype(np.float32) / float(n_cells)
        )

    @staticmethod
    def _clr_normalize(counts, *, pseudocount: float = 1.0) -> np.ndarray:
        """Centered log-ratio across proteins for each cell.

        clr(x) = log(x) - mean(log(x)) computed per cell. We add a pseudocount
        before taking the log so zero-count proteins do not produce -inf.
        """
        if sparse.issparse(counts):
            # densify-allowed: protein panels are O(100) features per cell, dense is appropriate
            dense = counts.toarray()
        else:
            dense = np.asarray(counts)
        dense = dense.astype(np.float64, copy=True)
        np.add(dense, pseudocount, out=dense)
        np.log(dense, out=dense)
        # Per-cell mean across proteins (axis=1).
        row_mean = dense.mean(axis=1, keepdims=True)
        dense -= row_mean
        return dense

    @staticmethod
    def _plot_detection_rate(
        detection_rate: Sequence[float],
        isotype_flags: Sequence[bool],
        ctx: PipelineContext,
    ) -> None:
        if len(detection_rate) == 0:
            return
        rates = np.asarray(detection_rate, dtype=np.float32)
        # Bucket per-protein detection rate into quartiles for the violin.
        # When the panel is tiny we fall back to two buckets (low/high).
        n_buckets = 4 if len(rates) >= 4 else 2
        # quantile cut points.
        qs = np.linspace(0, 1, n_buckets + 1)
        cuts = np.quantile(rates, qs)
        # Make cuts strictly increasing so digitize buckets correctly.
        for i in range(1, len(cuts)):
            if cuts[i] <= cuts[i - 1]:
                cuts[i] = cuts[i - 1] + 1e-6
        bucket_idx = np.clip(
            np.digitize(rates, cuts[1:-1], right=True), 0, n_buckets - 1
        )
        buckets: list[np.ndarray] = []
        labels: list[str] = []
        for b in range(n_buckets):
            sel = rates[bucket_idx == b]
            if len(sel) == 0:
                # matplotlib.violinplot requires non-empty datasets; fill with the
                # bucket midpoint so empty quantile slices still draw a flat line.
                lo = cuts[b]
                hi = cuts[b + 1]
                sel = np.array([(lo + hi) / 2.0], dtype=np.float32)
            buckets.append(sel)
            labels.append(f"Q{b + 1}\n[{cuts[b]:.2f},{cuts[b + 1]:.2f}]")

        fig, ax = plt.subplots(figsize=(max(5, n_buckets * 1.6), 4))
        ax.violinplot(buckets, showmeans=True, showmedians=True)
        ax.set_xticks(range(1, n_buckets + 1))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("per-protein detection rate")
        ax.set_title(
            f"ADT detection-rate buckets (n_proteins={len(rates)}, "
            f"isotypes={int(sum(isotype_flags))})"
        )
        ax.set_ylim(-0.05, 1.05)
        fig.tight_layout()
        fig.savefig(
            ctx.figure_dir / "protein_detection_rate_violin.png",
            dpi=160,
            bbox_inches="tight",
        )
        plt.close(fig)
