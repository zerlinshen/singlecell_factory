#!/usr/bin/env python3
"""Render public-matrix RNA velocity reproduction outputs for Trevino Figure 3A.

Inputs are project-root artifacts:
- a layer-aligned public RNA velocity AnnData with spliced/unspliced layers
- the public RNA pipeline viability final AnnData with UMAP + annotation

The script transfers pipeline annotations/UMAP onto the velocity object by cell
barcode, samples a bounded stratified subset, computes scVelo stochastic
velocity, and writes auditable figure outputs plus a summary JSON.  This is a
public processed-matrix reproduction lane, not FASTQ-level reproduction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scvelo as scv  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sparse_colsum(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).ravel().astype(float, copy=False)
    return np.asarray(matrix).sum(axis=0).astype(float, copy=False)


def choose_cells(obs: pd.DataFrame, max_cells: int, seed: int) -> list[str]:
    if len(obs) <= max_cells:
        return list(map(str, obs.index))
    rng = np.random.default_rng(seed)
    if "cell_type" in obs.columns:
        chosen: list[str] = []
        groups = obs.groupby("cell_type", observed=True).indices
        # Preserve rare cell types while capping total cells.
        floor = min(250, max(30, max_cells // max(1, len(groups) * 6)))
        remaining = max_cells
        for _, idx in sorted(groups.items(), key=lambda kv: len(kv[1])):
            idx = np.asarray(idx)
            take = min(len(idx), floor, remaining)
            if take > 0:
                chosen.extend(obs.index[rng.choice(idx, size=take, replace=False)].astype(str).tolist())
                remaining -= take
        if remaining > 0:
            pool = np.setdiff1d(np.arange(len(obs)), obs.index.get_indexer(chosen), assume_unique=False)
            if len(pool) > 0:
                take = min(remaining, len(pool))
                chosen.extend(obs.index[rng.choice(pool, size=take, replace=False)].astype(str).tolist())
        return chosen[:max_cells]
    return obs.index[rng.choice(len(obs), size=max_cells, replace=False)].astype(str).tolist()


def save_figure_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": out_dir / f"{stem}.png",
        "svg": out_dir / f"{stem}.svg",
        "pdf": out_dir / f"{stem}.pdf",
    }
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    pdf_mode = "native"
    try:
        fig.savefig(paths["pdf"], dpi=300, bbox_inches="tight")
    except ValueError as exc:
        # Matplotlib's PDF backend can reject non-finite path vertices emitted by
        # scVelo stream plots even when the PNG/SVG render succeeds.  Preserve a
        # PDF deliverable by embedding the already-rendered PNG and record the
        # fallback mode in the summary.
        if "finite numbers" not in str(exc):
            raise
        pdf_mode = "png_embedded_fallback_after_pdf_nonfinite"
        img = plt.imread(paths["png"])
        fallback_fig, ax = plt.subplots(figsize=(max(4, img.shape[1] / 300), max(4, img.shape[0] / 300)))
        ax.imshow(img)
        ax.axis("off")
        fallback_fig.savefig(paths["pdf"], dpi=300, bbox_inches="tight")
        plt.close(fallback_fig)
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    } | {"pdf_mode": pdf_mode}


def expression_vector(adata: ad.AnnData, gene: str) -> np.ndarray | None:
    if gene not in adata.var_names:
        return None
    values = adata[:, gene].X
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel().astype(float, copy=False)


def marker_pseudotime_metrics(adata: ad.AnnData, markers: list[str]) -> tuple[dict[str, dict], plt.Figure | None]:
    if "velocity_pseudotime" not in adata.obs.columns:
        return {}, None
    pt = np.asarray(adata.obs["velocity_pseudotime"], dtype=float)
    finite = np.isfinite(pt)
    if finite.sum() < 10:
        return {}, None
    bins = np.linspace(float(np.nanmin(pt[finite])), float(np.nanmax(pt[finite])), 21)
    bin_mid = (bins[:-1] + bins[1:]) / 2
    metrics: dict[str, dict] = {}
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for marker in markers:
        expr = expression_vector(adata, marker)
        if expr is None:
            metrics[marker] = {"present": False}
            continue
        valid = finite & np.isfinite(expr)
        if valid.sum() < 10:
            metrics[marker] = {"present": True, "valid_cells": int(valid.sum())}
            continue
        rho, pval = spearmanr(pt[valid], expr[valid])
        binned = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            mask = valid & (pt >= lo) & (pt < hi)
            binned.append(float(np.nanmedian(expr[mask])) if mask.any() else np.nan)
        ax.plot(bin_mid, binned, marker="o", linewidth=1.2, markersize=3, label=marker)
        metrics[marker] = {
            "present": True,
            "valid_cells": int(valid.sum()),
            "spearman_rho": None if np.isnan(rho) else float(rho),
            "spearman_pvalue": None if np.isnan(pval) else float(pval),
            "binned_median_expression": binned,
        }
    ax.set_xlabel("velocity pseudotime")
    ax.set_ylabel("normalized expression (median per bin)")
    ax.set_title("Figure 3A quantitative check — marker trends along velocity pseudotime")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    return metrics, fig


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.family": "DejaVu Sans",
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "legend.fontsize": 7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    scv.settings.verbosity = 1


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--velocity-h5ad", required=True, type=Path)
    p.add_argument("--pipeline-h5ad", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--max-cells", type=int, default=8000)
    p.add_argument("--top-genes", type=int, default=2500)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()

    configure_plotting()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    velocity = ad.read_h5ad(args.velocity_h5ad)
    pipeline = ad.read_h5ad(args.pipeline_h5ad)
    common = pipeline.obs_names.intersection(velocity.obs_names)
    if len(common) == 0:
        raise ValueError("No overlapping cells between velocity and pipeline AnnData")
    velocity = velocity[common, :].copy()
    pipeline = pipeline[common, :].copy()

    for col in [
        "leiden",
        "cell_type",
        "cell_type_marker",
        "annotation_confidence",
        "sample",
        "development_stage",
        "trevino_seurat_clusters",
    ]:
        if col in pipeline.obs.columns:
            velocity.obs[col] = pipeline.obs[col].astype(str) if col != "annotation_confidence" else pipeline.obs[col]
    if "X_umap" not in pipeline.obsm:
        raise KeyError("pipeline AnnData missing obsm['X_umap']")
    velocity.obsm["X_umap"] = np.asarray(pipeline.obsm["X_umap"])
    del pipeline

    selected_cells = choose_cells(velocity.obs, args.max_cells, args.seed)
    velocity = velocity[selected_cells, :].copy()

    trend_markers = ["SOX2", "PAX6", "EOMES", "NEUROD1", "NEUROD6", "SLC17A7", "SATB2", "GAD2"]

    score = sparse_colsum(velocity.X)
    for layer in ["spliced", "unspliced"]:
        if layer not in velocity.layers:
            raise KeyError(f"Velocity input missing layer {layer!r}")
        score += sparse_colsum(velocity.layers[layer])
    top_n = min(args.top_genes, velocity.n_vars)
    top = np.argsort(score)[-top_n:]
    marker_indices = np.asarray([velocity.var_names.get_loc(g) for g in trend_markers if g in velocity.var_names])
    top = np.sort(np.unique(np.concatenate([top, marker_indices])))
    velocity = velocity[:, top].copy()

    retain_genes = [g for g in trend_markers if g in velocity.var_names]
    scv.pp.filter_and_normalize(velocity, min_shared_counts=20, retain_genes=retain_genes)
    # scvelo 0.3.4 computes PCA/neighbors here; keep explicit params recorded in summary.
    scv.pp.moments(velocity, n_pcs=20, n_neighbors=20)
    scv.tl.velocity(velocity, mode="stochastic")
    scv.tl.velocity_graph(velocity)
    vel_layer = velocity.layers["velocity"]
    if sparse.issparse(vel_layer):
        velocity.obs["velocity_length"] = np.sqrt(np.asarray(vel_layer.multiply(vel_layer).sum(axis=1)).ravel())
    else:
        velocity.obs["velocity_length"] = np.linalg.norm(np.asarray(vel_layer), axis=1)
    try:
        scv.tl.velocity_confidence(velocity)
    except Exception as exc:  # confidence is useful but non-essential for Figure 3A rendering
        velocity.uns["velocity_confidence_warning"] = str(exc)
    try:
        scv.tl.velocity_pseudotime(velocity)
    except Exception as exc:
        velocity.uns["velocity_pseudotime_warning"] = str(exc)

    computed_h5ad = args.out_dir / "figure3a_velocity_computed_subset.h5ad"
    velocity.write_h5ad(computed_h5ad, compression="gzip")

    figure_outputs: dict[str, dict[str, str]] = {}
    # Stream plot on pipeline UMAP, colored by marker annotation.
    fig, ax = plt.subplots(figsize=(8, 6))
    scv.pl.velocity_embedding_stream(
        velocity,
        basis="umap",
        color="cell_type" if "cell_type" in velocity.obs.columns else "leiden",
        title="Trevino Fig. 3A reproduction — public RNA velocity (matrix-level)",
        legend_loc="right margin",
        density=1.4,
        linewidth=1.0,
        show=False,
        ax=ax,
    )
    figure_outputs["velocity_stream_cell_type"] = save_figure_triple(
        fig, args.out_dir, "figure3a_velocity_stream_cell_type"
    )

    fig = sc.pl.embedding(
        velocity,
        basis="umap",
        color="cell_type" if "cell_type" in velocity.obs.columns else "leiden",
        title="Public RNA pipeline annotations on Figure 3A velocity subset",
        show=False,
        return_fig=True,
        size=12,
    )
    figure_outputs["annotation_umap"] = save_figure_triple(fig, args.out_dir, "figure3a_annotation_umap")

    fig = sc.pl.embedding(
        velocity,
        basis="umap",
        color="velocity_length",
        title="RNA velocity length on pipeline UMAP",
        show=False,
        return_fig=True,
        size=12,
        cmap="magma",
    )
    figure_outputs["velocity_length_umap"] = save_figure_triple(
        fig, args.out_dir, "figure3a_velocity_length_umap"
    )

    if "velocity_pseudotime" in velocity.obs.columns:
        fig = sc.pl.embedding(
            velocity,
            basis="umap",
            color="velocity_pseudotime",
            title="RNA velocity pseudotime on pipeline UMAP",
            show=False,
            return_fig=True,
            size=12,
            cmap="viridis",
        )
        figure_outputs["velocity_pseudotime_umap"] = save_figure_triple(
            fig, args.out_dir, "figure3a_velocity_pseudotime_umap"
        )

    marker_metrics, marker_fig = marker_pseudotime_metrics(velocity, trend_markers)
    if marker_fig is not None:
        figure_outputs["marker_pseudotime_trends"] = save_figure_triple(
            marker_fig, args.out_dir, "figure3a_marker_pseudotime_trends"
        )

    cell_type_counts = (
        {str(k): int(v) for k, v in velocity.obs["cell_type"].value_counts().items()}
        if "cell_type" in velocity.obs.columns
        else {}
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "truth_boundary": (
            "Public FASTQ/SRA raw reads not found; Figure 3A output is public "
            "processed-matrix RNA velocity reproduction from GSE162170 matrices."
        ),
        "velocity_h5ad": str(args.velocity_h5ad),
        "pipeline_h5ad": str(args.pipeline_h5ad),
        "out_dir": str(args.out_dir),
        "computed_subset_h5ad": str(computed_h5ad),
        "parameters": {
            "max_cells": args.max_cells,
            "top_genes": args.top_genes,
            "retained_marker_genes": retain_genes,
            "seed": args.seed,
            "min_shared_counts": 20,
            "moments_n_pcs": 20,
            "moments_n_neighbors": 20,
            "velocity_mode": "stochastic",
            "basis": "pipeline X_umap transferred by cell barcode",
        },
        "input_overlap_cells": int(len(common)),
        "computed_shape": [int(velocity.n_obs), int(velocity.n_vars)],
        "layers_after_compute": list(map(str, velocity.layers.keys())),
        "has_velocity": bool("velocity" in velocity.layers),
        "has_velocity_graph": bool("velocity_graph" in velocity.uns),
        "has_velocity_pseudotime": bool("velocity_pseudotime" in velocity.obs.columns),
        "velocity_graph_shape": list(velocity.uns["velocity_graph"].shape)
        if "velocity_graph" in velocity.uns
        else None,
        "velocity_pseudotime_summary": (
            {
                "min": float(np.nanmin(velocity.obs["velocity_pseudotime"])),
                "max": float(np.nanmax(velocity.obs["velocity_pseudotime"])),
                "mean": float(np.nanmean(velocity.obs["velocity_pseudotime"])),
                "finite_cells": int(np.isfinite(velocity.obs["velocity_pseudotime"]).sum()),
            }
            if "velocity_pseudotime" in velocity.obs.columns
            else None
        ),
        "marker_pseudotime_metrics": marker_metrics,
        "cell_type_counts": cell_type_counts,
        "figure_outputs": figure_outputs,
        "software": {
            "scvelo": scv.__version__,
            "scanpy": sc.__version__,
        },
    }
    (args.out_dir / "figure3a_velocity_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
