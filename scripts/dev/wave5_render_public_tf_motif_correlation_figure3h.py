#!/usr/bin/env python3
"""Render Trevino Figure 3H public-resource TF motif/gene correlation heatmap.

This produces an honest Figure-3H-like panel from the public author resource
``TF_MotifExpressionCorrelation.RDS`` after export to TSV and the Figure 3F
selected TF/motif pairs.

Boundary: the paper panel is described at the motif-cluster level. The local
public resource contains motif-level correlations and no motif-cluster mapping
object was found, so this is not exact 24 motif-cluster parity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tf-motif-correlation-tsv", required=True, type=Path)
    p.add_argument("--selected-pairs-tsv", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--min-display-rho", type=float, default=0.4)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    corr = pd.read_csv(args.tf_motif_correlation_tsv, sep="\t")
    selected = pd.read_csv(args.selected_pairs_tsv, sep="\t")
    required = {"motif.name", "gene.symbol", "motif.id", "gene.id", "glia.pb.rho"}
    missing_corr = sorted(required.difference(corr.columns))
    missing_sel = sorted(required.difference(selected.columns))
    if missing_corr or missing_sel:
        raise ValueError(f"missing columns: corr={missing_corr}, selected={missing_sel}")

    selected_genes = list(dict.fromkeys(selected["gene.symbol"].astype(str).tolist()))
    motif_label_by_id = (
        selected.assign(motif_label=selected["motif.name"].astype(str) + " (" + selected["motif.id"].astype(str) + ")")
        .drop_duplicates("motif.id")
        .set_index("motif.id")["motif_label"]
        .to_dict()
    )
    selected_motif_ids = list(motif_label_by_id)
    sub = corr[
        corr["gene.symbol"].astype(str).isin(selected_genes)
        & corr["motif.id"].astype(str).isin(selected_motif_ids)
    ].copy()
    # If duplicated motif/gene combinations appear, keep the strongest absolute
    # correlation so the displayed matrix remains one value per cell.
    sub["abs_rho"] = sub["glia.pb.rho"].abs()
    sub = sub.sort_values("abs_rho", ascending=False).drop_duplicates(["gene.symbol", "motif.id"])
    matrix = sub.pivot(index="gene.symbol", columns="motif.id", values="glia.pb.rho")
    matrix = matrix.reindex(index=selected_genes, columns=selected_motif_ids)
    matrix.columns = [motif_label_by_id[c] for c in matrix.columns]

    matrix_tsv = args.out_dir / "figure3h_tf_motif_gene_correlation_matrix.tsv"
    matrix.to_csv(matrix_tsv, sep="\t", na_rep="NA")
    long_tsv = args.out_dir / "figure3h_tf_motif_gene_correlation_long.tsv"
    sub.drop(columns=["abs_rho"]).to_csv(long_tsv, sep="\t", index=False)

    fig_w = max(7.0, 0.32 * matrix.shape[1])
    fig_h = max(7.0, 0.22 * matrix.shape[0])
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    arr = matrix.to_numpy(dtype=float)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e5e5")
    im = ax.imshow(np.ma.masked_invalid(arr), cmap=cmap, vmin=-1.0, vmax=1.0, aspect="auto")
    ax.set_title("Trevino Fig. 3H public-resource TF motif/gene correlations")
    ax.set_xlabel("selected motif from Figure 3F")
    ax.set_ylabel("selected TF gene from Figure 3F")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=6)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.75, label="glia.pb.rho")
    figure_outputs = {"tf_motif_gene_correlation_heatmap": save_triple(fig, args.out_dir, "figure3h_tf_motif_gene_correlation_heatmap")}

    finite = np.isfinite(arr)
    strong = finite & (np.abs(arr) >= args.min_display_rho)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "GENERATED_AUTHOR_RESOURCE_MOTIF_GENE_CORRELATION_HEATMAP_NOT_CLUSTER_PARITY",
        "truth_boundary": (
            "Heatmap from TF_MotifExpressionCorrelation motif-level rows and Figure 3F selected pairs. "
            "The paper panel is motif-cluster-level; no local/public 24 motif-cluster mapping object was found, "
            "so exact Figure 3H cluster parity is not claimed."
        ),
        "tf_motif_correlation_tsv": str(args.tf_motif_correlation_tsv),
        "selected_pairs_tsv": str(args.selected_pairs_tsv),
        "out_dir": str(args.out_dir),
        "selected_genes": len(selected_genes),
        "selected_motif_ids": len(selected_motif_ids),
        "observed_gene_motif_pairs": int(finite.sum()),
        "possible_gene_motif_pairs": int(matrix.shape[0] * matrix.shape[1]),
        "coverage_fraction": float(finite.sum() / max(1, matrix.shape[0] * matrix.shape[1])),
        "strong_abs_rho_pairs": int(strong.sum()),
        "matrix_shape": list(matrix.shape),
        "outputs": {
            "matrix_tsv": str(matrix_tsv),
            "long_tsv": str(long_tsv),
            "figure_outputs": figure_outputs,
        },
        "checks": {
            "selected_genes_eq_31": len(selected_genes) == 31,
            "selected_motif_ids_ge_20": len(selected_motif_ids) >= 20,
            "observed_pairs_gt_0": bool(finite.sum() > 0),
            "strong_pairs_gt_0": bool(strong.sum() > 0),
            "finite_values_within_correlation_bounds": bool(np.nanmin(arr) >= -1.0 and np.nanmax(arr) <= 1.0),
        },
    }
    summary_path = args.out_dir / "figure3h_tf_motif_gene_correlation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
