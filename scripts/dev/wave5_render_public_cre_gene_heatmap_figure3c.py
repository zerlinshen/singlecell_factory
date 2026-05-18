#!/usr/bin/env python3
"""Render Trevino Figure 3C public-matrix CRE-gene heatmap evidence.

Inputs are the earliest public computable artifacts currently available locally:

- Table S3B GluN peak-gene links (13,989 rows)
- public ATAC counts matrix + consensus peak BED
- Figure 3B ATAC cells with transferred pseudotime
- full public RNA velocity-prepared matrix + Figure 3A RNA cells with velocity
  pseudotime

Boundary: this is not the author's exact 363-pseudobulk / CCA / fragment-level
pipeline. It reconstructs a public processed-matrix heatmap by pseudotime bins.
The local Table S3B extraction contains 10 link clusters although the paper
caption/methods state k=5, so the script derives a k=5 clustering from the
current RNA expression profiles and records the Table S3B discrepancy.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402
from sklearn.cluster import KMeans  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": out_dir / f"{stem}.png",
        "svg": out_dir / f"{stem}.svg",
        "pdf": out_dir / f"{stem}.pdf",
    }
    fig.savefig(paths["png"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {f"{k}_path": str(v) for k, v in paths.items()} | {
        f"{k}_sha256": sha256_file(v) for k, v in paths.items()
    }


def coord_key(chrom: object, start: object, end: object) -> str:
    return f"{chrom}:{int(start)}:{int(end)}"


def read_gene_activity_header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().strip().split()


def build_peak_row_index(consensus_bed: Path) -> dict[str, int]:
    mapping: dict[str, int] = {}
    with gzip.open(consensus_bed, "rt") as fh:
        for idx, line in enumerate(fh):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            mapping[coord_key(fields[0], fields[1], fields[2])] = idx
    return mapping


def make_bins(values: np.ndarray, n_bins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(values)
    bins_sorted = np.floor(np.arange(len(values)) * n_bins / len(values)).astype(int)
    bins_sorted = np.minimum(bins_sorted, n_bins - 1)
    bins = np.empty_like(bins_sorted)
    bins[order] = bins_sorted
    bin_counts = np.bincount(bins, minlength=n_bins).astype(float)
    bin_centers = np.array(
        [np.nanmean(values[bins == i]) if np.any(bins == i) else np.nan for i in range(n_bins)],
        dtype=float,
    )
    return bins, bin_counts, bin_centers


def binned_means(matrix, bins: np.ndarray, n_bins: int) -> np.ndarray:
    """Return genes/features x bins means."""
    out = []
    for b in range(n_bins):
        mask = bins == b
        if not np.any(mask):
            out.append(np.zeros(matrix.shape[1], dtype=np.float32))
            continue
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
    mean = np.nanmean(matrix, axis=1, keepdims=True)
    std = np.nanstd(matrix, axis=1, keepdims=True)
    std[std == 0] = 1.0
    out = (matrix - mean) / std
    out[~np.isfinite(out)] = 0.0
    return out.astype(np.float32, copy=False)


def parse_atac_count_profiles(
    *,
    counts_tsv_gz: Path,
    wanted_peak_rows: Iterable[int],
    selected_col_idx: np.ndarray,
    atac_bins: np.ndarray,
    atac_bin_counts: np.ndarray,
    n_bins: int,
    progress_every: int,
) -> dict[int, np.ndarray]:
    wanted = set(map(int, wanted_peak_rows))
    profiles: dict[int, np.ndarray] = {}
    if not wanted:
        return profiles
    max_wanted = max(wanted)
    with gzip.open(counts_tsv_gz, "rt") as fh:
        for row_idx, line in enumerate(fh):
            if row_idx > max_wanted and len(profiles) == len(wanted):
                break
            if row_idx not in wanted:
                continue
            values = np.fromstring(line, sep=" ", dtype=np.float32)
            if values.size <= int(selected_col_idx.max()):
                raise ValueError(
                    f"ATAC counts row {row_idx} has {values.size} columns; "
                    f"need index {int(selected_col_idx.max())}"
                )
            selected = values[selected_col_idx]
            sums = np.bincount(atac_bins, weights=selected, minlength=n_bins).astype(np.float32)
            means = sums / np.maximum(atac_bin_counts, 1.0)
            profiles[row_idx] = np.log1p(means).astype(np.float32, copy=False)
            if progress_every and len(profiles) % progress_every == 0:
                print(
                    f"parsed_peak_profiles={len(profiles)}/{len(wanted)} "
                    f"at_counts_row={row_idx}",
                    flush=True,
                )
    missing = wanted.difference(profiles)
    if missing:
        raise ValueError(f"Missing {len(missing)} requested peak rows from ATAC counts matrix")
    return profiles


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table-s3b", required=True, type=Path)
    p.add_argument("--consensus-bed", required=True, type=Path)
    p.add_argument("--atac-counts", required=True, type=Path)
    p.add_argument("--atac-gene-activity", required=True, type=Path, help="Used only for ATAC cell-order header.")
    p.add_argument("--atac-pseudotime-h5ad", required=True, type=Path)
    p.add_argument("--rna-full-h5ad", required=True, type=Path)
    p.add_argument("--rna-pseudotime-h5ad", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--n-bins", type=int, default=60)
    p.add_argument("--k", type=int, default=5)
    p.add_argument("--max-links", type=int, default=0, help="Optional deterministic cap for smoke/debug runs.")
    p.add_argument("--seed", type=int, default=13)
    p.add_argument("--progress-every", type=int, default=1000)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.table_s3b, sep="\t")
    expected_cols = {
        "Peak chromosome",
        "Peak start",
        "Peak end",
        "Gene symbol",
        "correlation",
        "Link cluster",
    }
    missing_cols = sorted(expected_cols.difference(table.columns))
    if missing_cols:
        raise ValueError(f"Table S3B missing columns: {missing_cols}")
    table["coord_key"] = [
        coord_key(c, s, e)
        for c, s, e in zip(table["Peak chromosome"], table["Peak start"], table["Peak end"], strict=True)
    ]
    table["input_row"] = np.arange(len(table))
    if args.max_links and args.max_links > 0:
        table = table.sort_values(["Link cluster", "correlation"], ascending=[True, False]).head(args.max_links).copy()
        table = table.sort_values("input_row").reset_index(drop=True)

    peak_row_index = build_peak_row_index(args.consensus_bed)
    table["peak_row_idx"] = table["coord_key"].map(peak_row_index)
    peak_mapped = table["peak_row_idx"].notna()

    rna_pt = ad.read_h5ad(args.rna_pseudotime_h5ad, backed="r")
    if "velocity_pseudotime" not in rna_pt.obs:
        raise KeyError("RNA pseudotime H5AD missing obs['velocity_pseudotime']")
    rna_cells = rna_pt.obs_names.to_numpy()
    rna_pseudotime = rna_pt.obs["velocity_pseudotime"].to_numpy(dtype=float)
    finite_rna = np.isfinite(rna_pseudotime)
    rna_cells = rna_cells[finite_rna]
    rna_pseudotime = rna_pseudotime[finite_rna]
    rna_pt.file.close()

    rna_genes_requested = list(dict.fromkeys(table.loc[peak_mapped, "Gene symbol"].astype(str)))
    rna_full = ad.read_h5ad(args.rna_full_h5ad, backed="r")
    rna_obs_idx = rna_full.obs_names.get_indexer(rna_cells)
    keep_rna_cells = rna_obs_idx >= 0
    rna_cells = rna_cells[keep_rna_cells]
    rna_pseudotime = rna_pseudotime[keep_rna_cells]
    rna_obs_idx = rna_obs_idx[keep_rna_cells]
    rna_var_idx = rna_full.var_names.get_indexer(rna_genes_requested)
    present_gene_mask = rna_var_idx >= 0
    present_genes = [g for g, ok in zip(rna_genes_requested, present_gene_mask, strict=True) if ok]
    rna_var_idx = rna_var_idx[present_gene_mask]
    if len(present_genes) < 10:
        raise ValueError(f"Too few linked genes present in RNA matrix: {len(present_genes)}")
    rna_subset = rna_full[rna_obs_idx, rna_var_idx].to_memory()
    rna_full.file.close()
    rna_bins, rna_bin_counts, rna_bin_centers = make_bins(rna_pseudotime, args.n_bins)
    rna_expr = rna_subset.X
    if sparse.issparse(rna_expr):
        rna_expr = rna_expr.copy()
        rna_expr.data = np.log1p(rna_expr.data)
    else:
        rna_expr = np.log1p(np.maximum(np.asarray(rna_expr, dtype=np.float32), 0.0))
    gene_profiles_matrix = binned_means(rna_expr, rna_bins, args.n_bins)
    gene_profiles = {
        gene: gene_profiles_matrix[i, :]
        for i, gene in enumerate(present_genes)
    }

    atac = ad.read_h5ad(args.atac_pseudotime_h5ad, backed="r")
    if "transferred_pseudotime" not in atac.obs:
        raise KeyError("ATAC pseudotime H5AD missing obs['transferred_pseudotime']")
    atac_cells = atac.obs_names.to_numpy()
    atac_pseudotime = atac.obs["transferred_pseudotime"].to_numpy(dtype=float)
    finite_atac = np.isfinite(atac_pseudotime)
    atac_cells = atac_cells[finite_atac]
    atac_pseudotime = atac_pseudotime[finite_atac]
    atac.file.close()
    atac_bins, atac_bin_counts, atac_bin_centers = make_bins(atac_pseudotime, args.n_bins)
    atac_header = read_gene_activity_header(args.atac_gene_activity)
    cell_to_col = {cell: i for i, cell in enumerate(atac_header)}
    selected_col_idx = np.array([cell_to_col.get(str(cell), -1) for cell in atac_cells], dtype=int)
    keep_atac = selected_col_idx >= 0
    atac_cells = atac_cells[keep_atac]
    atac_pseudotime = atac_pseudotime[keep_atac]
    atac_bins = atac_bins[keep_atac]
    selected_col_idx = selected_col_idx[keep_atac]
    atac_bin_counts = np.bincount(atac_bins, minlength=args.n_bins).astype(float)
    if len(atac_cells) < args.n_bins:
        raise ValueError(f"Too few ATAC cells with count-matrix columns: {len(atac_cells)}")

    link_filter = peak_mapped & table["Gene symbol"].astype(str).isin(present_genes)
    usable = table.loc[link_filter].copy()
    usable["peak_row_idx"] = usable["peak_row_idx"].astype(int)
    wanted_peak_rows = sorted(usable["peak_row_idx"].unique())
    peak_profiles = parse_atac_count_profiles(
        counts_tsv_gz=args.atac_counts,
        wanted_peak_rows=wanted_peak_rows,
        selected_col_idx=selected_col_idx,
        atac_bins=atac_bins,
        atac_bin_counts=atac_bin_counts,
        n_bins=args.n_bins,
        progress_every=args.progress_every,
    )

    accessibility = np.vstack([peak_profiles[int(idx)] for idx in usable["peak_row_idx"]]).astype(np.float32)
    expression = np.vstack([gene_profiles[str(g)] for g in usable["Gene symbol"]]).astype(np.float32)
    acc_z = row_zscore(accessibility)
    expr_z = row_zscore(expression)

    kmeans = KMeans(n_clusters=args.k, n_init=50, random_state=args.seed)
    raw_clusters = kmeans.fit_predict(expr_z)
    bin_axis = np.linspace(0.0, 1.0, args.n_bins)
    cluster_score = {
        c: float(np.average(bin_axis, weights=np.maximum(kmeans.cluster_centers_[c] - kmeans.cluster_centers_[c].min(), 1e-6)))
        for c in range(args.k)
    }
    ordered_clusters = sorted(cluster_score, key=cluster_score.get)
    cluster_rename = {old: f"k{i + 1}" for i, old in enumerate(ordered_clusters)}
    usable["derived_k5_cluster"] = [cluster_rename[c] for c in raw_clusters]
    usable["expr_peak_bin"] = np.argmax(expr_z, axis=1)
    usable["acc_peak_bin"] = np.argmax(acc_z, axis=1)
    cluster_order_rank = {name: i for i, name in enumerate([cluster_rename[c] for c in ordered_clusters])}
    usable["_cluster_rank"] = usable["derived_k5_cluster"].map(cluster_order_rank)
    sort_order = np.lexsort((usable["expr_peak_bin"].to_numpy(), usable["_cluster_rank"].to_numpy()))
    usable_sorted = usable.iloc[sort_order].reset_index(drop=True)
    acc_sorted = acc_z[sort_order, :]
    expr_sorted = expr_z[sort_order, :]
    boundaries = np.flatnonzero(
        usable_sorted["derived_k5_cluster"].to_numpy()[1:] != usable_sorted["derived_k5_cluster"].to_numpy()[:-1]
    ) + 1

    matrix_npz = args.out_dir / "figure3c_cre_gene_heatmap_matrices.npz"
    np.savez_compressed(
        matrix_npz,
        accessibility_z=acc_sorted,
        expression_z=expr_sorted,
        atac_bin_centers=atac_bin_centers,
        rna_bin_centers=rna_bin_centers,
    )
    link_summary = args.out_dir / "figure3c_link_row_summary.tsv"
    usable_sorted.drop(columns=["_cluster_rank"]).to_csv(link_summary, sep="\t", index=False)

    height = max(6.0, min(18.0, len(usable_sorted) / 900.0))
    fig, axes = plt.subplots(1, 2, figsize=(10, height), sharey=True, constrained_layout=True)
    vmax = 2.5
    im0 = axes[0].imshow(acc_sorted, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    im1 = axes[1].imshow(expr_sorted, aspect="auto", interpolation="nearest", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[0].set_title("CRE accessibility (ATAC counts)")
    axes[1].set_title("Linked gene expression (RNA)")
    for ax in axes:
        ax.set_xlabel("pseudotime bin")
        ax.set_xticks([0, args.n_bins - 1])
        ax.set_xticklabels(["early", "late"])
        for b in boundaries:
            ax.axhline(b - 0.5, color="black", linewidth=0.4)
    axes[0].set_ylabel(f"linked CRE-gene rows (n={len(usable_sorted)})")
    # Label derived clusters at their row midpoints.
    label_positions = []
    labels = []
    start = 0
    cluster_vals = usable_sorted["derived_k5_cluster"].to_numpy()
    for b in list(boundaries) + [len(cluster_vals)]:
        label_positions.append((start + b - 1) / 2)
        labels.append(cluster_vals[start])
        start = b
    axes[0].set_yticks(label_positions)
    axes[0].set_yticklabels(labels, fontsize=8)
    fig.colorbar(im1, ax=axes, location="right", shrink=0.7, label="row z-score")
    fig.suptitle("Trevino Fig. 3C public-matrix reconstruction — derived k=5")
    figure_outputs = {
        "cre_gene_heatmap": save_triple(fig, args.out_dir, "figure3c_cre_gene_heatmap_public_matrix")
    }

    table_cluster_counts = table["Link cluster"].astype(str).value_counts().sort_index().to_dict()
    derived_cluster_counts = usable_sorted["derived_k5_cluster"].value_counts().sort_index().to_dict()
    status = (
        "GENERATED_PUBLIC_MATRIX_SUBSET_DERIVED_K5_WITH_TABLE_CLUSTER_DISCREPANCY"
        if len(usable_sorted) >= 0.9 * len(table)
        else "PARTIAL_PUBLIC_MATRIX_SUBSET_DERIVED_K5"
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "truth_boundary": (
            "Public processed-matrix reconstruction using Table S3B links, public ATAC counts, "
            "Figure 3B transferred ATAC pseudotime, and Figure 3A RNA velocity pseudotime. "
            "Not the author's exact 363 pseudobulks, CCA transfer, fragment-level processing, or original heatmap object."
        ),
        "table_s3b": str(args.table_s3b),
        "consensus_bed": str(args.consensus_bed),
        "atac_counts": str(args.atac_counts),
        "atac_gene_activity_header_source": str(args.atac_gene_activity),
        "atac_pseudotime_h5ad": str(args.atac_pseudotime_h5ad),
        "rna_full_h5ad": str(args.rna_full_h5ad),
        "rna_pseudotime_h5ad": str(args.rna_pseudotime_h5ad),
        "out_dir": str(args.out_dir),
        "parameters": {"n_bins": args.n_bins, "k": args.k, "seed": args.seed, "max_links": args.max_links},
        "input_links": int(len(table)),
        "usable_links": int(len(usable_sorted)),
        "usable_link_fraction": float(len(usable_sorted) / max(len(table), 1)),
        "unique_link_genes_requested": int(len(rna_genes_requested)),
        "unique_link_genes_present_in_rna": int(len(present_genes)),
        "peak_coordinate_mapped_links": int(peak_mapped.sum()),
        "unique_peak_rows_loaded": int(len(wanted_peak_rows)),
        "rna_cells_with_finite_pseudotime": int(len(rna_cells)),
        "atac_cells_with_finite_transferred_pseudotime_and_counts": int(len(atac_cells)),
        "rna_bin_count_min": int(rna_bin_counts.min()),
        "rna_bin_count_max": int(rna_bin_counts.max()),
        "atac_bin_count_min": int(atac_bin_counts.min()),
        "atac_bin_count_max": int(atac_bin_counts.max()),
        "table_s3b_cluster_count": int(table["Link cluster"].astype(str).nunique()),
        "table_s3b_cluster_counts": {str(k): int(v) for k, v in table_cluster_counts.items()},
        "paper_declared_k": 5,
        "derived_cluster_counts": {str(k): int(v) for k, v in derived_cluster_counts.items()},
        "known_discrepancy": (
            "Paper caption/methods state k=5, but local Table S3B extraction has "
            f"{table['Link cluster'].astype(str).nunique()} Link cluster labels."
        ),
        "outputs": {
            "link_row_summary_tsv": str(link_summary),
            "matrix_npz": str(matrix_npz),
            "figure_outputs": figure_outputs,
        },
        "checks": {
            "peak_mapping_complete_for_input_table": bool(peak_mapped.all()),
            "usable_link_fraction_ge_0_9": bool(len(usable_sorted) >= 0.9 * len(table)),
            "derived_k_matches_paper_declared_k": bool(args.k == 5),
            "finite_accessibility_matrix": bool(np.isfinite(acc_sorted).all()),
            "finite_expression_matrix": bool(np.isfinite(expr_sorted).all()),
        },
    }
    summary_path = args.out_dir / "figure3c_cre_gene_heatmap_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
