#!/usr/bin/env python3
"""Render Trevino Figure 3F public-matrix TF expression / motif activity heatmaps.

This reconstructs a Figure-3F-like panel from public processed matrices:

- TF/motif pair candidates from author TF_MotifExpressionCorrelation.RDS
  exported to TSV
- RNA expression from the full public RNA matrix binned by Figure 3A velocity
  pseudotime
- chromVAR motif activity from the public ATAC matrix binned by Figure 3B
  transferred pseudotime

Boundary: not exact author 363-pseudobulk or motif-cluster object parity.
"""
from __future__ import annotations

import argparse
import gzip
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
from scipy import sparse  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    paths = {ext: out_dir / f"{stem}.{ext}" for ext in ("png", "svg", "pdf")}
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    }


def make_bins(values: np.ndarray, n_bins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(values)
    bins_sorted = np.floor(np.arange(len(values)) * n_bins / len(values)).astype(int)
    bins_sorted = np.minimum(bins_sorted, n_bins - 1)
    bins = np.empty_like(bins_sorted)
    bins[order] = bins_sorted
    counts = np.bincount(bins, minlength=n_bins).astype(float)
    centers = np.array([np.nanmean(values[bins == i]) if np.any(bins == i) else np.nan for i in range(n_bins)])
    return bins, counts, centers


def binned_means(matrix, bins: np.ndarray, n_bins: int) -> np.ndarray:
    out = []
    for b in range(n_bins):
        mask = bins == b
        sub = matrix[mask, :]
        mean = sub.mean(axis=0)
        if sparse.issparse(mean):
            mean = mean.A1
        else:
            mean = np.asarray(mean).ravel()
        out.append(mean.astype(np.float32, copy=False))
    return np.vstack(out).T


def row_zscore(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    mean = matrix.mean(axis=1, keepdims=True)
    std = matrix.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    out = (matrix - mean) / std
    out[~np.isfinite(out)] = 0.0
    return out


def read_chromvar_rows(path: Path, wanted_motif_ids: list[str], selected_col_idx: np.ndarray) -> dict[str, np.ndarray]:
    wanted = set(wanted_motif_ids)
    rows: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        max_idx = int(selected_col_idx.max())
        for line in fh:
            if not line.strip():
                continue
            motif_id, rest = line.split(maxsplit=1)
            if motif_id not in wanted:
                continue
            values = np.fromstring(rest, sep=" ", dtype=np.float32)
            if values.size <= max_idx:
                raise ValueError(f"chromVAR row {motif_id} has {values.size} columns, need {max_idx}")
            rows[motif_id] = values[selected_col_idx]
            if len(rows) == len(wanted):
                break
    missing = wanted.difference(rows)
    if missing:
        raise ValueError(f"Missing chromVAR motif rows: {sorted(missing)[:10]}")
    return rows


def read_chromvar_index_and_header(path: Path) -> tuple[list[str], set[str]]:
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        motifs = {line.split(maxsplit=1)[0] for line in fh if line.strip()}
    return header, motifs


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tf-motif-correlation-tsv", required=True, type=Path)
    p.add_argument("--link-summary", required=True, type=Path)
    p.add_argument("--chromvar", required=True, type=Path)
    p.add_argument("--atac-pseudotime-h5ad", required=True, type=Path)
    p.add_argument("--rna-full-h5ad", required=True, type=Path)
    p.add_argument("--rna-pseudotime-h5ad", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--min-rho", type=float, default=0.4)
    p.add_argument("--n-pairs", type=int, default=31)
    p.add_argument("--n-bins", type=int, default=60)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    corr = pd.read_csv(args.tf_motif_correlation_tsv, sep="\t")
    required = {"motif.name", "gene.symbol", "motif.id", "gene.id", "glia.pb.rho"}
    missing = sorted(required.difference(corr.columns))
    if missing:
        raise ValueError(f"correlation TSV missing columns: {missing}")
    links = pd.read_csv(args.link_summary, sep="\t", usecols=["Gene symbol"])
    linked_genes = set(links["Gene symbol"].astype(str))
    chrom_header, chrom_motifs = read_chromvar_index_and_header(args.chromvar)

    rna_full = ad.read_h5ad(args.rna_full_h5ad, backed="r")
    rna_genes = set(map(str, rna_full.var_names))
    candidates = corr[
        (corr["glia.pb.rho"] >= args.min_rho)
        & corr["gene.symbol"].astype(str).isin(linked_genes)
        & corr["gene.symbol"].astype(str).isin(rna_genes)
        & corr["motif.id"].astype(str).isin(chrom_motifs)
    ].copy()
    candidates = candidates.sort_values("glia.pb.rho", ascending=False)
    selected = candidates.drop_duplicates("gene.symbol", keep="first").head(args.n_pairs).copy()
    if len(selected) < min(args.n_pairs, 10):
        raise ValueError(f"Too few selected TF/motif pairs: {len(selected)}")
    selected_genes = selected["gene.symbol"].astype(str).tolist()
    selected_motifs = selected["motif.id"].astype(str).tolist()

    rna_pt = ad.read_h5ad(args.rna_pseudotime_h5ad, backed="r")
    rna_cells = rna_pt.obs_names.to_numpy()
    rna_pseudotime = rna_pt.obs["velocity_pseudotime"].to_numpy(dtype=float)
    finite = np.isfinite(rna_pseudotime)
    rna_cells = rna_cells[finite]
    rna_pseudotime = rna_pseudotime[finite]
    rna_pt.file.close()
    obs_idx = rna_full.obs_names.get_indexer(rna_cells)
    keep = obs_idx >= 0
    obs_idx = obs_idx[keep]
    rna_pseudotime = rna_pseudotime[keep]
    var_idx = rna_full.var_names.get_indexer(selected_genes)
    rna_subset = rna_full[obs_idx, var_idx].to_memory()
    rna_full.file.close()
    rna_bins, rna_bin_counts, rna_bin_centers = make_bins(rna_pseudotime, args.n_bins)
    expr = rna_subset.X
    if sparse.issparse(expr):
        expr = expr.copy()
        expr.data = np.log1p(expr.data)
    else:
        expr = np.log1p(np.maximum(np.asarray(expr, dtype=np.float32), 0.0))
    expr_profile = binned_means(expr, rna_bins, args.n_bins)

    atac = ad.read_h5ad(args.atac_pseudotime_h5ad, backed="r")
    atac_cells = atac.obs_names.to_numpy()
    atac_pt = atac.obs["transferred_pseudotime"].to_numpy(dtype=float)
    finite_atac = np.isfinite(atac_pt)
    atac_cells = atac_cells[finite_atac]
    atac_pt = atac_pt[finite_atac]
    atac.file.close()
    cell_to_col = {cell: i for i, cell in enumerate(chrom_header)}
    selected_col_idx = np.array([cell_to_col.get(str(c), -1) for c in atac_cells], dtype=int)
    keep = selected_col_idx >= 0
    selected_col_idx = selected_col_idx[keep]
    atac_pt = atac_pt[keep]
    atac_bins, atac_bin_counts, atac_bin_centers = make_bins(atac_pt, args.n_bins)
    chrom_rows = read_chromvar_rows(args.chromvar, selected_motifs, selected_col_idx)
    motif_matrix = np.vstack([chrom_rows[m] for m in selected_motifs]).astype(np.float32)
    motif_profile = binned_means(motif_matrix.T, atac_bins, args.n_bins)

    expr_z = row_zscore(expr_profile)
    motif_z = row_zscore(motif_profile)
    order = np.argsort(np.argmax(expr_z, axis=1))
    expr_z = expr_z[order, :]
    motif_z = motif_z[order, :]
    selected = selected.iloc[order].reset_index(drop=True)
    labels = [f"{r['gene.symbol']} / {r['motif.name']}" for _, r in selected.iterrows()]

    matrix_npz = args.out_dir / "figure3f_tf_motif_heatmap_matrices.npz"
    np.savez_compressed(
        matrix_npz,
        expression_z=expr_z,
        motif_activity_z=motif_z,
        rna_bin_centers=rna_bin_centers,
        atac_bin_centers=atac_bin_centers,
    )
    selected_tsv = args.out_dir / "figure3f_selected_tf_motif_pairs.tsv"
    selected.to_csv(selected_tsv, sep="\t", index=False)

    fig_h = max(6.0, 0.22 * len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(10, fig_h), sharey=True, constrained_layout=True)
    vmax = 2.5
    im0 = axes[0].imshow(expr_z, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    im1 = axes[1].imshow(motif_z, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[0].set_title("TF expression (RNA)")
    axes[1].set_title("motif activity (chromVAR)")
    for ax in axes:
        ax.set_xlabel("pseudotime bin")
        ax.set_xticks([0, args.n_bins - 1])
        ax.set_xticklabels(["early", "late"])
    axes[0].set_yticks(range(len(labels)))
    axes[0].set_yticklabels(labels, fontsize=7)
    axes[0].set_ylabel("selected TF / motif")
    fig.colorbar(im1, ax=axes, location="right", shrink=0.75, label="row z-score")
    fig.suptitle("Trevino Fig. 3F public-matrix TF/motif dynamics")
    figure_outputs = {"tf_motif_heatmap": save_triple(fig, args.out_dir, "figure3f_tf_motif_heatmap")}

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "GENERATED_PUBLIC_MATRIX_TF_MOTIF_HEATMAP_DERIVED_PAIRS",
        "truth_boundary": (
            "TF/motif heatmap from public RNA/chromVAR matrices binned by Figure 3A/3B pseudotime and "
            "TF_MotifExpressionCorrelation-derived pairs. Not exact author 363-pseudobulk, motif-cluster, or original 31/24 object parity."
        ),
        "tf_motif_correlation_tsv": str(args.tf_motif_correlation_tsv),
        "link_summary": str(args.link_summary),
        "chromvar": str(args.chromvar),
        "atac_pseudotime_h5ad": str(args.atac_pseudotime_h5ad),
        "rna_full_h5ad": str(args.rna_full_h5ad),
        "rna_pseudotime_h5ad": str(args.rna_pseudotime_h5ad),
        "out_dir": str(args.out_dir),
        "parameters": {"min_rho": args.min_rho, "n_pairs": args.n_pairs, "n_bins": args.n_bins},
        "candidate_pairs_after_filters": int(len(candidates)),
        "selected_pairs": int(len(selected)),
        "selected_unique_genes": int(selected["gene.symbol"].nunique()),
        "selected_unique_motifs": int(selected["motif.id"].nunique()),
        "rna_cells_binned": int(len(rna_pseudotime)),
        "atac_cells_binned": int(len(atac_pt)),
        "rna_bin_count_min": int(rna_bin_counts.min()),
        "rna_bin_count_max": int(rna_bin_counts.max()),
        "atac_bin_count_min": int(atac_bin_counts.min()),
        "atac_bin_count_max": int(atac_bin_counts.max()),
        "selected_pairs_preview": selected.head(10).to_dict("records"),
        "outputs": {
            "selected_pairs_tsv": str(selected_tsv),
            "matrix_npz": str(matrix_npz),
            "figure_outputs": figure_outputs,
        },
        "checks": {
            "selected_pairs_eq_requested": int(len(selected)) == min(args.n_pairs, int(selected["gene.symbol"].nunique())),
            "selected_pairs_ge_10": int(len(selected)) >= 10,
            "finite_expression_matrix": bool(np.isfinite(expr_z).all()),
            "finite_motif_activity_matrix": bool(np.isfinite(motif_z).all()),
            "rna_bins_nonempty": bool(rna_bin_counts.min() > 0),
            "atac_bins_nonempty": bool(atac_bin_counts.min() > 0),
        },
    }
    summary_path = args.out_dir / "figure3f_tf_motif_heatmap_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
