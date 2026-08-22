from __future__ import annotations

import json
import hashlib
import logging
import matplotlib
from pathlib import Path

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from ._scanpy_compat import import_scanpy_or_stub

sc = import_scanpy_or_stub()

from ..context import PipelineContext
from .._reference_mapping import (
    ASSIGNMENT_STATUS_ACCEPTED,
    align_reference_genes,
    calibrate_reference_ood_threshold,
    map_knn_reference,
    resolve_reference_device,
)


__references__ = {
    "Tirosh_marker_scoring_2016": {
        "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
        "authors": "Tirosh et al.",
        "journal": "Science",
        "year": "2016",
        "doi": "10.1126/science.aad0501",
        "description": "Cluster-vote marker-mean scoring approach for cell type assignment.",
    },
    "Stuart_label_transfer_2019": {
        "title": "Comprehensive Integration of Single-Cell Data",
        "authors": "Stuart et al.",
        "journal": "Cell",
        "year": "2019",
        "doi": "10.1016/j.cell.2019.05.031",
        "description": "Reference-based label transfer principles; this module uses explicit CPU/GPU KNN transfer with reference-only OOD calibration.",
    },
}



# Curated marker catalogue used ONLY when neither ctx.cfg.markers nor the
# marker_db_loader provided a marker map. The authoritative path for any
# scientific run is to populate `adata.uns["marker_db_index"]` via
# `marker_db_loader` (which pins CellMarker 2.0 / PanglaoDB / CellTypist /
# scTypeDB version strings into the run manifest). This dict is a tumor /
# immune / stromal starter set assembled by hand from common single-cell
# atlas conventions; it is NOT a versioned DB.
#
# When tuning markers for a new tissue or condition, prefer (a) configuring
# `marker_db_loader` for that (tissue, condition) pair or (b) passing
# `--markers-json` to the CLI rather than editing this dict.
DEFAULT_MARKERS = {
    "Tumor epithelial": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1"],
    "T cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "IL7R"],
    "NK cell": ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "B cell": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Myeloid/Macro": ["LYZ", "CD68", "CST3", "FCGR3A", "AIF1"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL3A1", "FAP"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "Plasma cell": ["IGHG1", "IGKC", "MZB1", "XBP1"],
    "Mast cell": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "Dendritic cell": ["CD1C", "CLEC9A", "FCER1A", "IRF8"],
}

# DEFAULT_MARKERS is a hand-curated tumor/immune/stromal *starter* pack (NSCLC-
# flavoured). It is never a tissue-validated atlas. Scientific audit 2026-08-05
# (C1 Trevino brain smoke) showed majority "NK cell" labels on neural clusters
# when this pack ran without a tissue-specific marker map.
DEFAULT_MARKER_SOURCE = "default_tumor_immune_starter"
# Tissues for which the starter pack is *plausible* (still non-claimable without
# an explicit operator-supplied marker map or marker_db). Used only for messaging.
_STARTER_PACK_PLAUSIBLE_TISSUES = frozenset(
    {
        "lung",
        "nsclc",
        "lusc",
        "luad",
        "lung_cancer",
        "lung-cancer",
        "tumor_lung",
    }
)

EPITHELIAL_QC_MARKERS = ("EPCAM", "KRT8", "KRT18")

# Label written by the confidence gates for cells whose marker evidence did not
# clear ``annotation_confidence_threshold``.
UNASSIGNED_LABEL = "Unknown"

# Theme token ``defaults.missing_value_color`` (nature_high_impact.yaml).
# "Unknown" is the annotation declining to assign an identity, not a cell type:
# giving it a categorical hue alongside the real types both consumes a colour
# and implies a biological identity that was never called, so it is held out of
# the qualitative allocation and drawn in the theme's missing-value grey.
UNASSIGNED_COLOR = "#CFCFCF"

# Theme token ``palettes.continuous_viridis`` (low #440154 / mid #21918C /
# high #FDE725) is viridis. It is the theme's only continuous ramp usable on a
# scatter panel: ``sequential_expression`` and ``qc_sequential`` start at
# near-white (#F7FBFF / #F7FCF5), which would make the lowest-value cells
# invisible against the white canvas instead of merely low.
SEQUENTIAL_CMAP = "viridis"

logger = logging.getLogger(__name__)


class AnnotationModule:
    """Optional module: marker-scoring based cell-type annotation with confidence scores."""

    name = "annotation"
    requires_keys = {"obs": ["leiden"]}
    provides_keys = {"obs": ["cell_type"]}

    @staticmethod
    def _resolve_leiden_key(adata) -> str:
        """Pick the cluster column to bind cell_type against.

        Prefer ``leiden_corrected`` (the post-batch-correction snapshot written
        by batch_correction) so cell_type is a deterministic function of the
        FINAL post-correction clustering regardless of execution mode. Fall
        back to ``leiden`` when batch_correction did not run.
        """
        if "leiden_corrected" in adata.obs:
            return "leiden_corrected"
        return "leiden"

    def run(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None or "leiden" not in adata.obs:
            raise ValueError("Annotation requires clustered AnnData.")
        leiden_key = self._resolve_leiden_key(adata)
        ctx.metadata["annotation_leiden_source"] = leiden_key

        using_default_markers = not bool(ctx.cfg.markers)
        marker_map = ctx.cfg.markers or DEFAULT_MARKERS
        available = {
            cell_type: [gene for gene in genes if gene in adata.var_names]
            for cell_type, genes in marker_map.items()
        }
        available = {cell_type: genes for cell_type, genes in available.items() if genes}
        if not available:
            raise ValueError("No valid marker genes found in dataset.")

        strategy = getattr(ctx.cfg, "annotation_strategy", "cluster_voting")
        ctx.metadata["annotation_strategy"] = strategy

        # --- Marker provenance + claim honesty (audit 2026-08-05 W1.2) ---
        tissue = str(getattr(ctx.cfg, "tissue", "") or "").strip().lower()
        marker_db_present = bool(
            isinstance(getattr(adata, "uns", None), dict)
            and adata.uns.get("marker_db_index")
        )
        if using_default_markers:
            marker_source = DEFAULT_MARKER_SOURCE
            # Default starter pack is never confirmatory cell-type evidence.
            # Operators must supply --markers-json or run marker_db_loader for
            # claimable annotation (tissue-validated panels).
            annotation_claimable = False
            claim_reason = "default_tumor_immune_starter_nonclaimable"
            if tissue and tissue not in _STARTER_PACK_PLAUSIBLE_TISSUES:
                claim_reason = (
                    "default_starter_pack_tissue_mismatch:"
                    f"tissue={tissue or 'unspecified'}"
                )
                logger.warning(
                    "ANNOTATION_TISSUE_MISMATCH: using DEFAULT_MARKERS (tumor/"
                    "immune starter pack) while cfg.tissue=%r. Labels are "
                    "exploratory and NOT claimable. Supply --markers-json or "
                    "run marker_db_loader for tissue-validated markers "
                    "(scientific audit 2026-08-05 C1: brain cells labeled NK).",
                    tissue,
                )
            else:
                logger.warning(
                    "ANNOTATION_DEFAULT_STARTER_PACK: using DEFAULT_MARKERS "
                    "(tumor/immune starter; NSCLC-flavoured). "
                    "annotation_claimable=False. For claimable labels pass "
                    "--markers-json or populate marker_db_loader "
                    "(cfg.tissue=%r).",
                    tissue or "unspecified",
                )
        elif marker_db_present:
            marker_source = "marker_db_or_operator"
            annotation_claimable = True
            claim_reason = "operator_or_marker_db_panel"
        else:
            marker_source = "operator_markers_json"
            annotation_claimable = True
            claim_reason = "operator_supplied_markers"

        ctx.metadata["annotation_marker_source"] = marker_source
        ctx.metadata["annotation_tissue"] = tissue or "unspecified"
        ctx.metadata["annotation_claimable"] = annotation_claimable
        ctx.metadata["annotation_claim_reason"] = claim_reason
        if "annotation" not in adata.uns or not isinstance(adata.uns.get("annotation"), dict):
            adata.uns["annotation"] = {}
        adata.uns["annotation"].update(
            {
                "marker_source": marker_source,
                "tissue": tissue or "unspecified",
                "claimable": annotation_claimable,
                "claim_reason": claim_reason,
            }
        )

        # Score genes on full adata (per cell) — populates score_* obs columns and
        # provides per-cell confidence used by reference mapping conservative gate.
        # random_state: control-bin sampling in score_genes is RNG-dependent
        # (Tirosh 2016); propagate the canonical seed for reproducibility.
        score_seed = getattr(ctx, "random_state", 0)
        for cell_type, genes in available.items():
            sc.tl.score_genes(
                adata, genes, score_name=f"score_{cell_type}",
                use_raw=False, random_state=score_seed,
            )

        score_cols = [f"score_{cell_type}" for cell_type in available]
        score_mat = adata.obs[score_cols].copy()
        score_mat.columns = list(available.keys())

        if strategy == "cluster_voting":
            # Cluster-level assignment uses raw mean expression per marker gene set,
            # not aggregated score_genes values. score_genes background subtraction is
            # unreliable when absolute marker expression is low (background genes in the
            # same expression bin dominate the signal). Raw mean expression correctly
            # captures relative marker enrichment per cluster.
            cluster_score_matrix = self._cluster_mean_expression(adata, available, leiden_key)
            cluster_label_map = cluster_score_matrix.idxmax(axis=1).to_dict()
            adata.obs["cell_type"] = adata.obs[leiden_key].astype(str).map(cluster_label_map).astype(str)
            # Confidence: per-cell score_genes margin (max - second), for reference mapping gate
            max_scores = score_mat.max(axis=1)
            vals = score_mat.values
            if vals.shape[1] >= 2:
                second_scores = pd.Series(
                    np.partition(vals, -2, axis=1)[:, -2], index=score_mat.index,
                )
            else:
                second_scores = pd.Series(np.zeros(vals.shape[0], dtype=float), index=score_mat.index)
            adata.obs["annotation_confidence"] = (max_scores - second_scores).values
            # Unknown gate: cluster-level raw mean max score
            cluster_max_score = cluster_score_matrix.max(axis=1)
            low_conf_clusters = set(
                cluster_max_score[cluster_max_score < ctx.cfg.annotation_confidence_threshold].index.astype(str)
            )
            low_conf = adata.obs[leiden_key].astype(str).isin(low_conf_clusters)
            adata.obs.loc[low_conf, "cell_type"] = "Unknown"
            # Write audit trail
            cluster_score_matrix.to_csv(ctx.table_dir / "cluster_score_matrix.csv")
        else:
            # cell_argmax: original per-cell path
            max_scores = score_mat.max(axis=1)
            vals = score_mat.values
            if vals.shape[1] >= 2:
                second_scores = pd.Series(
                    np.partition(vals, -2, axis=1)[:, -2], index=score_mat.index,
                )
            else:
                second_scores = pd.Series(np.zeros(vals.shape[0], dtype=float), index=score_mat.index)
            adata.obs["cell_type"] = score_mat.idxmax(axis=1).values
            adata.obs["annotation_confidence"] = (max_scores - second_scores).values
            low_conf = max_scores < ctx.cfg.annotation_confidence_threshold
            adata.obs.loc[low_conf, "cell_type"] = "Unknown"
            # Write cluster_score_matrix.csv for consistency
            cluster_score_matrix = score_mat.groupby(
                adata.obs[leiden_key].astype(str), observed=True
            ).mean()
            cluster_score_matrix.index.name = leiden_key
            cluster_score_matrix.to_csv(ctx.table_dir / "cluster_score_matrix.csv")

        adata.obs["cell_type_marker"] = adata.obs["cell_type"].astype(str)

        self._try_reference_mapping(adata, ctx)

        low_conf_mask = adata.obs["cell_type"] == "Unknown"
        ctx.metadata["annotation_unknown_pct"] = round(
            float(low_conf_mask.sum()) / max(len(low_conf_mask), 1) * 100, 2
        )

        score_cols = [f"score_{ct}" for ct in available]
        base_cols = ["leiden", "cell_type", "annotation_confidence", "cell_type_marker"] + score_cols
        extra_cols = [
            c for c in [
                "reference_predicted_label",
                "reference_confidence",
                "reference_distance",
                "reference_assignment_status",
                "reference_cell_type",
                "reference_ood",
            ] if c in adata.obs.columns
        ]
        adata.obs[base_cols + extra_cols].to_csv(ctx.table_dir / "cell_type_annotation.csv")

        cluster_summary = (
            adata.obs.groupby("leiden", observed=True)["cell_type"]
            .agg(lambda x: x.value_counts().idxmax())
        ).rename("majority_cell_type")
        cluster_summary.to_csv(ctx.table_dir / "cluster_majority_cell_type.csv")
        self._write_epithelial_marker_qc(adata, ctx)

        # --- Visualizations ---
        cell_type_what = "Marker-based cell type"
        if not annotation_claimable:
            cell_type_what = "Marker-based cell type (EXPLORATORY / non-claimable)"
        self._plot_label_umap(
            adata, ctx,
            key="cell_type",
            filename="umap_cell_type.png",
            what=cell_type_what,
            palette_token="cell_type_qualitative",
            label_noun="cell types",
        )
        self._plot_score_umap(
            adata, ctx,
            key="annotation_confidence",
            filename="umap_annotation_confidence.png",
            what="Marker annotation confidence",
            # The number is a margin between two marker scores, so it means
            # nothing without saying which two; an unlabelled ramp is the
            # single most common reviewer complaint about a colour-coded panel.
            cbar_label="Marker-score margin (best − runner-up, a.u.)",
        )

        if "reference_cell_type" in adata.obs:
            self._plot_label_umap(
                adata, ctx,
                key="reference_cell_type",
                filename="umap_reference_cell_type.png",
                what="Reference label transfer",
                palette_token="cell_type_qualitative",
                label_noun="reference labels",
            )
        if "reference_confidence" in adata.obs:
            neighbours = ctx.metadata.get("reference_mapping_k", getattr(ctx.cfg, "reference_k", None))
            k_text = f"k = {int(neighbours)}" if neighbours is not None else "k nearest"
            self._plot_score_umap(
                adata, ctx,
                key="reference_confidence",
                filename="umap_reference_confidence.png",
                what="Reference label transfer confidence",
                cbar_label=f"Neighbour vote fraction ({k_text})",
                # A vote fraction is bounded by construction. Autoscaling it to
                # the observed range would paint a 0.95-1.00 panel across the
                # whole ramp and imply a spread that is not there.
                vlimits=(0.0, 1.0),
            )

        # Stacked bar: cell type composition per cluster
        self._plot_composition(adata, ctx)

    @staticmethod
    def _cluster_mean_expression(
        adata, available: dict[str, list[str]], leiden_key: str = "leiden"
    ) -> pd.DataFrame:
        """Compute mean expression of each marker gene set per leiden cluster.

        Returns DataFrame with rows=cluster, cols=cell_type, values=mean log-expr of gene set.
        Uses raw matrix slicing to avoid score_genes background subtraction artefacts.
        """
        X = adata.X
        leiden_vals = adata.obs[leiden_key].astype(str).values
        clusters = sorted(set(leiden_vals), key=lambda c: (not c.isdigit(), int(c) if c.isdigit() else c))
        gene_index = {g: i for i, g in enumerate(adata.var_names)}

        rows: dict[str, dict[str, float]] = {}
        for cluster in clusters:
            mask = leiden_vals == cluster
            x_clust = X[mask]
            if sparse.issparse(x_clust):
                mean_expr = np.asarray(x_clust.mean(axis=0)).flatten()
            else:
                mean_expr = np.asarray(x_clust, dtype=np.float64).mean(axis=0)
            row: dict[str, float] = {}
            for cell_type, genes in available.items():
                idxs = [gene_index[g] for g in genes if g in gene_index]
                row[cell_type] = float(mean_expr[idxs].mean()) if idxs else 0.0
            rows[cluster] = row

        df = pd.DataFrame(rows).T
        df.index.name = leiden_key
        return df

    @staticmethod
    def _try_reference_mapping(adata, ctx: PipelineContext) -> None:
        """Reference KNN label transfer with formal OOD calibration and fail-closed contracts."""
        ref_path = ctx.cfg.reference_adata
        if ref_path is None:
            ctx.metadata["reference_mapping_status"] = "skipped_no_reference"
            return

        mode = str(ctx.cfg.reference_override_mode).strip().lower()
        if mode not in {"all", "conservative"}:
            raise ValueError(
                f"Unsupported reference_override_mode {ctx.cfg.reference_override_mode!r}; "
                "expected 'all' or 'conservative'."
            )

        ref_path = Path(ref_path)
        if not ref_path.exists():
            raise FileNotFoundError(f"Requested reference H5AD file does not exist: {ref_path}")

        ref_stat_before = ref_path.stat()

        try:
            ref = ad.read_h5ad(ref_path)
        except Exception as exc:
            raise RuntimeError(f"Failed to read requested reference file {ref_path}: {exc}") from exc
        ref_stat_after = ref_path.stat()
        if (
            ref_stat_before.st_size != ref_stat_after.st_size
            or ref_stat_before.st_mtime_ns != ref_stat_after.st_mtime_ns
            or ref_stat_before.st_ino != ref_stat_after.st_ino
        ):
            raise RuntimeError(f"Requested reference file changed while it was being read: {ref_path}")

        ref_hasher = hashlib.sha256()
        with ref_path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                ref_hasher.update(chunk)
        reference_sha256 = ref_hasher.hexdigest()

        label_key = ctx.cfg.reference_label_key
        if label_key not in ref.obs.columns:
            raise ValueError(f"Requested reference label key {label_key!r} not found in {ref_path} .obs columns.")

        requested_device = getattr(ctx.cfg, "reference_device", "auto")
        device, device_policy = resolve_reference_device(
            requested_device,
            validation_domain=getattr(ctx.cfg, "reference_validation_domain", None),
        )
        ood_mode = getattr(ctx.cfg, "reference_ood_mode", "reference_quantile")
        dist_quantile = getattr(ctx.cfg, "reference_distance_quantile", 0.95)
        fixed_threshold = getattr(ctx.cfg, "reference_fixed_distance_threshold", None)
        group_key = getattr(ctx.cfg, "reference_calibration_group_key", None)
        cal_fraction = getattr(ctx.cfg, "reference_calibration_fraction", 0.2)
        min_shared = getattr(ctx.cfg, "reference_min_shared_genes", 50)
        k = int(ctx.cfg.reference_k)
        min_conf = float(ctx.cfg.reference_min_confidence)
        seed = int(getattr(ctx, "random_state", 42))

        # Resolve the exact one-to-one feature space first. Calibration and
        # query mapping must use the same reference columns so their cosine
        # distances are comparable.
        _, aligned_ref_genes, alignment_route, n_shared_genes, _, _ = align_reference_genes(
            adata,
            ref,
            min_shared=min_shared,
        )

        # Calibrate distance threshold on reference-only data.
        dist_thresh, cal_receipt = calibrate_reference_ood_threshold(
            ref_adata=ref,
            ref_genes=aligned_ref_genes,
            mode=ood_mode,
            quantile=dist_quantile,
            fixed_threshold=fixed_threshold,
            group_key=group_key,
            cal_fraction=cal_fraction,
            seed=seed,
            k=k,
        )
        cal_receipt["alignment_route"] = alignment_route
        cal_receipt["n_shared_genes"] = int(n_shared_genes)

        # Execute deterministic KNN mapping
        map_df, map_meta = map_knn_reference(
            query_adata=adata,
            ref_adata=ref,
            label_key=label_key,
            k=k,
            min_confidence=min_conf,
            distance_threshold=dist_thresh,
            device=device,
            seed=seed,
            min_shared_genes=min_shared,
        )

        # Populate adata.obs
        adata.obs["reference_predicted_label"] = map_df["reference_predicted_label"]
        adata.obs["reference_confidence"] = map_df["reference_confidence"]
        adata.obs["reference_distance"] = map_df["reference_distance"]
        adata.obs["reference_assignment_status"] = map_df["reference_assignment_status"]
        adata.obs["reference_cell_type"] = map_df["reference_cell_type"]
        adata.obs["reference_ood"] = map_df["reference_ood"]

        # Label override logic: ONLY accepted cells may override cell_type
        accepted_mask = adata.obs["reference_assignment_status"] == ASSIGNMENT_STATUS_ACCEPTED
        if mode == "all":
            override_mask = accepted_mask
        else:
            marker_unknown = adata.obs["cell_type_marker"].astype(str) == "Unknown"
            marker_low = adata.obs["annotation_confidence"] < float(ctx.cfg.annotation_confidence_threshold)
            override_mask = accepted_mask & (marker_unknown | marker_low)

        adata.obs.loc[override_mask, "cell_type"] = adata.obs.loc[override_mask, "reference_cell_type"].values

        applied = int(override_mask.sum())
        total = int(adata.n_obs)
        ctx.metadata["reference_mapping_status"] = "completed"
        ctx.metadata["reference_mapping_source"] = str(ref_path)
        ctx.metadata["reference_mapping_source_sha256"] = reference_sha256
        ctx.metadata["reference_mapping_source_size_bytes"] = int(ref_stat_after.st_size)
        ctx.metadata["reference_mapping_label_key"] = label_key
        ctx.metadata["reference_mapping_backend"] = "knn"
        ctx.metadata["reference_mapping_scanvi_status"] = "not_run_missing_real_compatible_model_artifact"
        ctx.metadata["reference_mapping_shared_genes"] = int(map_meta["n_shared_genes"])
        ctx.metadata["reference_mapping_k"] = int(map_meta["effective_k"])
        ctx.metadata["reference_mapping_requested_device"] = requested_device
        ctx.metadata["reference_mapping_device"] = device
        ctx.metadata["reference_mapping_device_policy"] = device_policy
        ctx.metadata["reference_mapping_ood_mode"] = ood_mode
        ctx.metadata["reference_mapping_distance_threshold"] = float(dist_thresh)
        ctx.metadata["reference_mapping_override_mode"] = mode
        ctx.metadata["reference_mapping_overridden_cells"] = applied
        ctx.metadata["reference_mapping_overridden_pct"] = round(100.0 * applied / max(total, 1), 2)
        ctx.metadata["reference_mapping_accepted_cells"] = int(map_meta["accepted_cells"])
        ctx.metadata["reference_mapping_rejected_cells"] = int(map_meta["rejected_cells"])
        ctx.metadata["reference_mapping_rejected_low_confidence"] = int(map_meta["rejected_low_confidence"])
        ctx.metadata["reference_mapping_rejected_ood_distance"] = int(map_meta["rejected_ood_distance"])
        ctx.metadata["reference_mapping_cal_receipt"] = cal_receipt
        ctx.metadata["reference_mapping_timings"] = map_meta.get("timings", {})

        logger.info(
            "Reference mapping applied: %d/%d cells overridden (mode=%s, k=%d, min_conf=%.3f, dist_thresh=%.4f, accepted=%d/%d).",
            applied,
            total,
            ctx.metadata["reference_mapping_override_mode"],
            k,
            min_conf,
            dist_thresh,
            ctx.metadata["reference_mapping_accepted_cells"],
            total,
        )

    @staticmethod
    def _label_palette(adata, key: str, palette_token: str) -> tuple[list[str], list[str]]:
        """Return (categories, one distinct colour per category) for a label column.

        The colours are keyed by category so the SAME cell type keeps the SAME
        colour in every panel of the module (UMAP and composition bars); a type
        that changes colour between two panels of one run reads as two types.

        ``qualitative_colors`` is used rather than indexing a fixed palette
        modulo its length, which would hand two cell types the identical colour
        while the legend still lists both.
        """
        from .._figure_theme import qualitative_colors

        series = adata.obs[key]
        counts = series.astype(str).value_counts()
        # scanpy plots categories in ``cat.categories`` order, so the palette has
        # to be built in that order for colour i to land on category i.
        categories = (
            [str(c) for c in series.cat.categories]
            if hasattr(series, "cat")
            else sorted(str(c) for c in counts.index)
        )
        # Only labels that some cell carries take a hue. A categorical column
        # keeps levels no cell carries (a marker panel entry that never won a
        # cluster), and spending palette entries on them pushes real types
        # towards the ambiguous tail of the ramp for nothing.
        assignable = [c for c in categories if c != UNASSIGNED_LABEL and int(counts.get(c, 0)) > 0]
        hues = qualitative_colors(palette_token, len(assignable))
        if len(hues) != len(assignable):
            # Theme unavailable (sibling repo absent) — let scanpy pick, rather
            # than pass a short palette it would silently wrap.
            return categories, []
        hue_of = dict(zip(assignable, hues))
        colours = [hue_of.get(c, UNASSIGNED_COLOR) for c in categories]
        return categories, colours

    @classmethod
    def _plot_label_umap(
        cls, adata, ctx: PipelineContext, *,
        key: str, filename: str, what: str, palette_token: str,
        label_noun: str = "labels",
    ) -> None:
        """Draw a categorical annotation UMAP at journal geometry.

        rcParams cannot reach this panel: ``sc.pl.umap`` supplies its own
        palette, figure size and title, so a themed run still produced scanpy's
        default hues on an arbitrary canvas titled with the obs key. Two further
        defects are specific to a label panel — scanpy's right-margin legend
        lists every level of the categorical, including levels no cell carries
        (a legend row pointing at nothing), and it reports no cell counts. The
        legend is therefore rebuilt here from the observed labels, carrying n
        per label.
        """
        from .._figure_theme import journal_figure_size
        from matplotlib.lines import Line2D

        counts = adata.obs[key].astype(str).value_counts()
        categories, palette = cls._label_palette(adata, key, palette_token)

        # Same geometry as the cluster UMAP (clustering._plot_umap_clusters) so
        # the label panel and the cluster panel register point-for-point when a
        # reader compares them.
        fig, ax = plt.subplots(
            figsize=journal_figure_size("double", height_mm=110.0), constrained_layout=True
        )
        sc.pl.umap(
            adata,
            color=[key],
            palette=palette or None,
            show=False,
            legend_loc=None,        # rebuilt below with counts and live labels only
            frameon=False,          # UMAP axes carry no units; a box implies scale
            size=3,
            ax=ax,
        )
        present = [c for c in categories if int(counts.get(c, 0)) > 0]
        ax.set_title(f"{what} (n = {adata.n_obs:,} cells, {len(present)} {label_noun})")

        # scanpy records the colours it actually used, which is the only safe
        # source for the swatches when the theme palette was not applied.
        used = [str(c) for c in adata.uns.get(f"{key}_colors", [])] or palette
        if len(used) >= len(categories):
            colour_of = dict(zip(categories, used))
            handles = [
                Line2D(
                    [], [], marker="o", linestyle="none", markersize=3,
                    markerfacecolor=colour_of[label], markeredgecolor="none",
                    label=f"{label} ({int(counts[label]):,})",
                )
                for label in present
            ]
            # "outside" reserves the legend's space inside the canvas, so the
            # saved figure stays at the journal column width instead of being
            # widened by the tight bounding box.
            fig.legend(
                handles=handles, loc="outside right upper", frameon=False,
                handletextpad=0.3, labelspacing=0.35, borderaxespad=0.0,
                title="cells",
            )
        fig.savefig(ctx.figure_dir / filename)
        plt.close(fig)

    @staticmethod
    def _plot_score_umap(
        adata, ctx: PipelineContext, *,
        key: str, filename: str, what: str, cbar_label: str,
        vlimits: tuple[float, float] | None = None,
    ) -> None:
        """Draw a continuous annotation UMAP with a colourbar that states its units."""
        from .._figure_theme import journal_figure_size

        fig, ax = plt.subplots(
            figsize=journal_figure_size("double", height_mm=110.0), constrained_layout=True
        )
        limits = {} if vlimits is None else {"vmin": vlimits[0], "vmax": vlimits[1]}
        pre_existing = set(fig.axes)
        sc.pl.umap(
            adata,
            color=[key],
            cmap=SEQUENTIAL_CMAP,
            show=False,
            frameon=False,
            size=3,
            ax=ax,
            **limits,
        )
        ax.set_title(f"{what} (n = {adata.n_obs:,} cells)")
        # scanpy creates the colourbar as an extra axes; label whichever axes the
        # call added, since an unlabelled ramp does not say what the number is.
        for cax in (axis for axis in fig.axes if axis not in pre_existing):
            cax.set_ylabel(cbar_label)
        fig.savefig(ctx.figure_dir / filename)
        plt.close(fig)

    @classmethod
    def _plot_composition(cls, adata, ctx: PipelineContext) -> None:
        """Stacked cell-type composition per cluster.

        Previously a 12x5 inch pandas bar plot in ``tab20`` (which wraps after 20
        categories and is not the theme palette), titled "Fraction" against bare
        cluster ids: a reader could not tell whether a bar stood for 47 cells or
        13,000, and cell types kept different colours than the same types on the
        UMAP. Composition is a per-group proportion, so the group size is the
        first thing needed to judge it.
        """
        from .._figure_theme import journal_figure_size

        ct_counts = adata.obs.groupby(["leiden", "cell_type"], observed=True).size().unstack(fill_value=0)
        ct_frac = ct_counts.div(ct_counts.sum(axis=1), axis=0)

        # A cell type no cell carries must not become a legend entry with nothing
        # to point at. Whether the unstack keeps such a level depends on the
        # pandas version's handling of unused categoricals, so drop it here
        # explicitly. This cannot change ct_frac: the denominator is the row sum,
        # to which an all-zero column contributes nothing.
        live = [column for column in ct_counts.columns if int(ct_counts[column].sum()) > 0]
        empty = [str(column) for column in ct_counts.columns if column not in live]
        ct_counts, ct_frac = ct_counts[live], ct_frac[live]
        if empty:
            logger.info("Composition plot: omitting cell types with no cells: %s", ", ".join(empty))

        categories, palette = cls._label_palette(adata, "cell_type", "cell_type_qualitative")
        colour_of = dict(zip(categories, palette)) if palette else {}

        cluster_sizes = ct_counts.sum(axis=1)
        clusters = [str(c) for c in ct_frac.index]
        positions = np.arange(len(clusters), dtype=float)

        fig, ax = plt.subplots(
            figsize=journal_figure_size("double", height_mm=74.0), constrained_layout=True
        )
        bottom = np.zeros(len(clusters))
        for column in ct_frac.columns:
            values = ct_frac[column].to_numpy(dtype=float)
            ax.bar(
                positions, values, bottom=bottom, width=0.82, linewidth=0,
                color=colour_of.get(str(column)), label=str(column),
            )
            bottom += values

        ax.set_xticks(positions)
        # n per cluster sits under its own bar: the fraction axis alone hides a
        # 47-cell cluster and a 13,000-cell cluster behind identical bars.
        ax.set_xticklabels(
            [f"{cluster}\n{int(cluster_sizes.iloc[i]):,}" for i, cluster in enumerate(clusters)],
            fontsize=5, linespacing=1.1,
        )
        ax.set_xlim(-0.6, len(clusters) - 0.4)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel("Leiden cluster (cells per cluster below the label)")
        ax.set_ylabel("Fraction of cluster")
        title = (
            f"Cell-type composition per cluster "
            f"(n = {adata.n_obs:,} cells in {len(clusters)} clusters)"
        )
        # Under cluster-level voting every cell in a cluster carries the cluster's
        # label, so every bar is one full-height block. Say that, instead of
        # leaving a reader to wonder whether a uniform panel is a finding.
        if len(clusters) and float(ct_frac.max(axis=1).min()) >= 1.0:
            title += "\nevery cluster resolves to a single type — composition is the cluster-level assignment"
        ax.set_title(title)
        # "outside" keeps the legend inside the canvas, so the saved figure stays
        # at the journal column width instead of being widened by the tight
        # bounding box around a legend hanging off the axes.
        fig.legend(
            loc="outside center right", frameon=False,
            handlelength=0.9, handletextpad=0.4, labelspacing=0.35, borderaxespad=0.0,
        )
        fig.savefig(ctx.figure_dir / "cell_type_composition.png")
        plt.close(fig)

    @staticmethod
    def _marker_expression_frame(adata, markers: list[str]) -> pd.DataFrame:
        expr = adata[:, markers].X
        if sparse.issparse(expr):
            # densify-allowed: marker panel is O(n_cells x <=3 genes)
            expr = expr.toarray()
        else:
            expr = np.asarray(expr)
        return pd.DataFrame(expr, index=adata.obs_names, columns=markers)

    @staticmethod
    def _summarize_marker_distribution(expr: pd.DataFrame, labels, group_name: str) -> pd.DataFrame:
        work = expr.copy()
        work[group_name] = pd.Series(labels, index=expr.index).astype(str)
        rows = []
        for group, sub in work.groupby(group_name, observed=True):
            marker_values = sub.drop(columns=[group_name])
            for marker in marker_values.columns:
                values = marker_values[marker].to_numpy(dtype=float)
                rows.append(
                    {
                        group_name: str(group),
                        "marker": marker,
                        "n_cells": int(values.shape[0]),
                        "mean": float(np.mean(values)),
                        "median": float(np.median(values)),
                        "pct_positive": round(float(np.mean(values > 0.0)) * 100.0, 3),
                    }
                )
        return pd.DataFrame(rows)

    @classmethod
    def _write_epithelial_marker_qc(cls, adata, ctx: PipelineContext) -> None:
        present = [marker for marker in EPITHELIAL_QC_MARKERS if marker in adata.var_names]
        ctx.metadata["epithelial_marker_qc_present_markers"] = present
        if not present:
            ctx.metadata["epithelial_marker_qc_status"] = "skipped_no_epithelial_markers"
            return

        expr = cls._marker_expression_frame(adata, present)
        if "leiden" in adata.obs:
            cls._summarize_marker_distribution(
                expr,
                adata.obs["leiden"],
                "leiden",
            ).to_csv(ctx.table_dir / "epithelial_marker_summary_by_leiden.csv", index=False)
        if "cell_type" in adata.obs:
            cls._summarize_marker_distribution(
                expr,
                adata.obs["cell_type"],
                "cell_type",
            ).to_csv(ctx.table_dir / "epithelial_marker_summary_by_cell_type.csv", index=False)

        epithelial_mask = adata.obs["cell_type"].astype(str).str.contains(
            "epithelial",
            case=False,
            regex=False,
        ) if "cell_type" in adata.obs else pd.Series(False, index=adata.obs_names)
        marker_score = expr.mean(axis=1)
        qc_payload = {
            "present_markers": present,
            "epithelial_cells": int(epithelial_mask.sum()),
            "non_epithelial_cells": int((~epithelial_mask).sum()),
            "epithelial_marker_score_median": (
                float(marker_score[epithelial_mask].median()) if epithelial_mask.any() else None
            ),
            "non_epithelial_marker_score_median": (
                float(marker_score[~epithelial_mask].median()) if (~epithelial_mask).any() else None
            ),
            "status": "completed",
        }
        ctx.metadata["epithelial_marker_qc_status"] = "completed"
        ctx.metadata["epithelial_marker_score_median"] = qc_payload[
            "epithelial_marker_score_median"
        ]
        (ctx.table_dir / "epithelial_marker_qc.json").write_text(
            json.dumps(qc_payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        if "X_umap" in adata.obsm:
            cls._plot_marker_umaps(adata, ctx, present)

    @staticmethod
    def _plot_marker_umaps(adata, ctx: PipelineContext, markers: list[str]) -> None:
        """One UMAP panel per epithelial QC marker, on one journal-width canvas.

        ``sc.pl.umap(color=[...])`` builds its own multi-panel grid at its own
        size and titles each panel with the bare gene symbol, leaving the reader
        no way to know what the colour encodes — counts, scaled values, which
        matrix. The panel exists to check that the epithelial call is supported
        by epithelial markers, so the matrix behind the colour has to be named.
        """
        from .._figure_theme import journal_figure_size

        # A lone panel belongs in a single column; two or three need the double.
        column = "single" if len(markers) == 1 else "double"
        fig, axes = plt.subplots(
            1, len(markers),
            figsize=journal_figure_size(column, height_mm=66.0),
            constrained_layout=True, squeeze=False,
        )
        for ax, gene in zip(axes[0], markers):
            sc.pl.umap(
                adata,
                color=[gene],
                cmap=SEQUENTIAL_CMAP,
                show=False,
                frameon=False,      # UMAP axes carry no units; a box implies scale
                size=3,
                ax=ax,
            )
            ax.set_title(gene, style="italic")   # gene symbols are set in italic

        # sc.pl.umap reads adata.raw when it exists (scanpy's use_raw default), so
        # the caption must name that matrix rather than assume adata.X. The
        # clustering module records what adata.raw holds; fall back to the
        # unqualified word when it did not.
        source = "adata.raw" if adata.raw is not None else "adata.X"
        semantics = str(ctx.metadata.get("raw_semantics", "")) if adata.raw is not None else ""
        unit = "log-normalised expression" if semantics == "normalized_log1p_all_genes" else "expression"
        # Each panel keeps its own colour range (scanpy's per-gene default), so
        # say so: identical yellow in two panels is not identical expression.
        scaling = ", scaled per panel" if len(markers) > 1 else ""
        fig.suptitle(
            f"Epithelial QC markers — colour = {unit} ({source}){scaling}; "
            f"n = {adata.n_obs:,} cells",
            fontsize=7,
        )
        fig.savefig(ctx.figure_dir / "umap_epithelial_markers.png")
        plt.close(fig)
