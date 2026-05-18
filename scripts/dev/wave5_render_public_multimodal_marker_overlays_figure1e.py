#!/usr/bin/env python3
"""Render Trevino Figure 1E public-matrix multimodal marker overlays.

For SOX9/EOMES/NEUROD2/DLX2, plot RNA expression on the public RNA UMAP,
ATAC gene activity on the public ATAC subset UMAP, and chromVAR motif activity
on the same ATAC UMAP.

Boundary: this uses processed public matrices and a 12k ATAC subset embedding.
DLX2 has no exact motif row in the local chromVAR matrix; DLX6 is recorded as a
DLX-family motif proxy rather than claimed as exact DLX2 motif parity.
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


MARKERS = ["SOX9", "EOMES", "NEUROD2", "DLX2"]
MOTIF_ROWS = {
    "SOX9": {"motif_id": "MA0077.1_SOX9", "mapping": "exact"},
    "EOMES": {"motif_id": "MA0800.1_EOMES", "mapping": "exact"},
    "NEUROD2": {"motif_id": "MA0668.1_NEUROD2", "mapping": "exact"},
    "DLX2": {"motif_id": "MA0882.1_DLX6", "mapping": "DLX-family proxy; exact DLX2 row absent"},
}


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


def to_1d(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def read_rows(path: Path, wanted_rows: dict[str, str], selected_cells: list[str]) -> dict[str, np.ndarray]:
    wanted = set(wanted_rows.values())
    out: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        col_index = pd.Index(header)
        idx = col_index.get_indexer(selected_cells)
        if np.any(idx < 0):
            missing = [selected_cells[i] for i in np.where(idx < 0)[0][:5]]
            raise ValueError(f"selected cells missing from {path}: {missing}")
        for line in fh:
            if not line.strip():
                continue
            row_id, rest = line.split(maxsplit=1)
            if row_id not in wanted:
                continue
            values = np.fromstring(rest, sep=" ", dtype=np.float32)
            out[row_id] = values[idx]
            if len(out) == len(wanted):
                break
    missing = wanted.difference(out)
    if missing:
        raise ValueError(f"missing rows from {path}: {sorted(missing)}")
    return {key: out[row_id] for key, row_id in wanted_rows.items()}


def scatter_signal(ax: plt.Axes, coords: np.ndarray, values: np.ndarray, title: str, cmap: str) -> None:
    finite = np.isfinite(values)
    vals = values.copy()
    if finite.any():
        lo, hi = np.nanpercentile(vals[finite], [1, 99])
        if hi > lo:
            vals = np.clip(vals, lo, hi)
    ax.scatter(coords[:, 0], coords[:, 1], s=1.2, c="#d9d9d9", linewidths=0)
    sc = ax.scatter(coords[finite, 0], coords[finite, 1], s=2.0, c=vals[finite], cmap=cmap, linewidths=0)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna-h5ad", required=True, type=Path)
    p.add_argument("--atac-umap-h5ad", required=True, type=Path)
    p.add_argument("--atac-gene-activities", required=True, type=Path)
    p.add_argument("--chromvar", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rna = ad.read_h5ad(args.rna_h5ad, backed="r")
    rna_coords = np.asarray(rna.obsm["X_umap"])
    rna_present = [m for m in MARKERS if m in set(map(str, rna.var_names))]
    rna_values = {}
    for marker in rna_present:
        idx = rna.var_names.get_loc(marker)
        rna_values[marker] = np.log1p(np.maximum(to_1d(rna[:, idx].X).astype(np.float32), 0))
    rna.file.close()

    atac = ad.read_h5ad(args.atac_umap_h5ad, backed="r")
    atac_coords = np.asarray(atac.obsm["X_umap"])
    atac_cells = list(map(str, atac.obs_names))
    atac.file.close()

    atac_gene_rows = read_rows(args.atac_gene_activities, {m: m for m in MARKERS}, atac_cells)
    atac_gene_values = {m: np.log1p(np.maximum(v.astype(np.float32), 0)) for m, v in atac_gene_rows.items()}
    motif_rows = {m: MOTIF_ROWS[m]["motif_id"] for m in MARKERS}
    motif_values = read_rows(args.chromvar, motif_rows, atac_cells)

    fig, axes = plt.subplots(len(MARKERS), 3, figsize=(12, 13), constrained_layout=True)
    signal_rows = []
    for row, marker in enumerate(MARKERS):
        motif_info = MOTIF_ROWS[marker]
        scatter_signal(axes[row, 0], rna_coords, rna_values[marker], f"{marker} RNA expression", "viridis")
        scatter_signal(axes[row, 1], atac_coords, atac_gene_values[marker], f"{marker} ATAC gene activity", "magma")
        scatter_signal(
            axes[row, 2],
            atac_coords,
            motif_values[marker],
            f"{marker} motif activity\\n{motif_info['motif_id']} ({motif_info['mapping']})",
            "coolwarm",
        )
        for modality, values in [
            ("rna_expression", rna_values[marker]),
            ("atac_gene_activity", atac_gene_values[marker]),
            ("chromvar_motif_activity", motif_values[marker]),
        ]:
            finite = np.isfinite(values)
            signal_rows.append(
                {
                    "marker": marker,
                    "modality": modality,
                    "motif_id": motif_info["motif_id"] if modality == "chromvar_motif_activity" else "",
                    "motif_mapping": motif_info["mapping"] if modality == "chromvar_motif_activity" else "",
                    "n_values": int(values.size),
                    "finite_values": int(finite.sum()),
                    "mean": float(np.nanmean(values)),
                    "std": float(np.nanstd(values)),
                    "nonzero_fraction": float(np.nanmean(values > 0)),
                }
            )
    figure_outputs = {"figure1e_multimodal_marker_overlays": save_triple(fig, args.out_dir, "figure1e_multimodal_marker_overlays")}
    signal_tsv = args.out_dir / "figure1e_marker_signal_summary.tsv"
    pd.DataFrame(signal_rows).to_csv(signal_tsv, sep="\t", index=False)

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "GENERATED_PUBLIC_MATRIX_FIGURE1E_WITH_DLX_FAMILY_MOTIF_PROXY",
        "truth_boundary": (
            "RNA expression, ATAC gene activity, and chromVAR motif activity were rendered from public processed matrices. "
            "ATAC panels use the 12k public ATAC subset embedding. DLX2 uses a DLX6 motif row as a DLX-family proxy because no exact DLX2 chromVAR row was found."
        ),
        "inputs": {
            "rna_h5ad": str(args.rna_h5ad),
            "atac_umap_h5ad": str(args.atac_umap_h5ad),
            "atac_gene_activities": str(args.atac_gene_activities),
            "chromvar": str(args.chromvar),
        },
        "out_dir": str(args.out_dir),
        "markers": MARKERS,
        "motif_rows": MOTIF_ROWS,
        "metrics": {
            "rna_cells": int(rna_coords.shape[0]),
            "atac_cells": int(atac_coords.shape[0]),
            "rna_markers_present": rna_present,
            "atac_gene_activity_markers_present": sorted(atac_gene_values),
            "chromvar_motif_markers_present": sorted(motif_values),
            "dlx2_exact_motif_present": False,
        },
        "outputs": {"signal_summary_tsv": str(signal_tsv), "figure_outputs": figure_outputs},
        "checks": {
            "rna_markers_all_present": sorted(rna_present) == sorted(MARKERS),
            "atac_gene_activity_markers_all_present": sorted(atac_gene_values) == sorted(MARKERS),
            "chromvar_motif_rows_all_mapped": sorted(motif_values) == sorted(MARKERS),
            "rna_cells_ge_50000": int(rna_coords.shape[0]) >= 50000,
            "atac_cells_ge_10000": int(atac_coords.shape[0]) >= 10000,
            "rna_umap_finite": bool(np.isfinite(rna_coords).all()),
            "atac_umap_finite": bool(np.isfinite(atac_coords).all()),
            "dlx2_proxy_recorded": MOTIF_ROWS["DLX2"]["mapping"].startswith("DLX-family proxy"),
            "all_signal_vectors_finite": all(np.isfinite(v).all() for v in list(rna_values.values()) + list(atac_gene_values.values()) + list(motif_values.values())),
        },
    }
    summary_path = args.out_dir / "figure1e_multimodal_marker_overlay_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
