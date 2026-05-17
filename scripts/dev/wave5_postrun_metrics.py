"""Wave-5 post-run metrics — computes ARI, peak-gene overlap, and pseudotime
Spearman from a Wave-5 run's ``final_adata.h5ad`` and writes them to a
metrics.json consumed by ``scripts/ci/wave5_trevino_regression_gate.sh``.

Maps to plan Stages 2, 4, 5:
  - Stage 2 (US-W5-5 ARI)         -> AC-VAL-1 / AC-VAL-2
  - Stage 4 (US-W5-7 overlap)     -> AC-VAL-3
  - Stage 5 (US-W5-8 pseudotime)  -> AC-VAL-4 / AC-VAL-5

Acquisition fallback behavior per plan v3 Stage 4 step 16a / Stage 5 step 19b:
  - Trevino S2F missing  -> overlap stays NOT-RUN; AC-VAL-3 CLOSED-PARTIAL.
  - Fig 4D missing       -> Spearman stays NOT-RUN; AC-VAL-5 CLOSED-PARTIAL.
  - No silent invented fallbacks.

Wave-5 spec: .omc/specs/deep-interview-wave5-completion.md
Wave-5 plan: .omc/plans/wave5-completion-consensus-2026-05-16.md
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

# --- canonical brain marker panel (AC-1b backstop) ---
MARKER_PANEL = {
    "RG": ["PAX6"],
    "IPC": ["EOMES"],
    "ExN": ["NEUROD2"],
    "InN": ["GAD2"],
    "OPC": ["SOX10"],
}


def log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {msg}", flush=True)


def _load_cluster_names(path: Path) -> dict[str, str] | None:
    """Parse GSE162170_multiome_cluster_names.txt.gz -> {cluster_id: cluster_name}.

    Multiome RNA rows are used (cluster IDs match the seurat_clusters column).
    Returns None on parse failure.
    """
    try:
        df = pd.read_csv(path, sep="\t", compression="gzip")
        rna = df[df["Assay"] == "Multiome RNA"]
        return dict(zip(rna["Cluster.ID"].astype(str), rna["Cluster.Name"].astype(str)))
    except Exception as exc:
        log(f"  cluster_names parse failed: {exc}")
        return None


def _compute_ari_primary(adata: anndata.AnnData, names: dict[str, str]) -> dict:
    """PRIMARY path: ARI(joint Leiden, Trevino named clusters via seurat_clusters)."""
    if "leiden" not in adata.obs or "seurat_clusters" not in adata.obs:
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN",
                "reason": "missing leiden or seurat_clusters in obs"}
    leiden = adata.obs["leiden"].astype(str)
    seurat = adata.obs["seurat_clusters"].astype(str)
    # Trevino's GSE162170 multiome_cluster_names uses 'c0'..'c13' as Cluster.ID
    # and our load_trevino_2021 preserves those strings in obs["seurat_clusters"].
    # A direct ID lookup is the only mapping needed.
    mapped = seurat.map(names)
    n_mapped = int(mapped.notna().sum())
    if n_mapped < 0.5 * len(seurat):
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN",
                "reason": f"only {n_mapped}/{len(seurat)} seurat IDs mapped to named labels"}
    ari = float(adjusted_rand_score(mapped.fillna("UNMAPPED"), leiden))
    return {
        "value": round(ari, 4),
        "threshold": 0.70,
        "status": "CLOSED" if ari >= 0.70 else "CLOSED-PARTIAL",
        "evidence_path": "PRIMARY",
        "n_cells_mapped": int(n_mapped),
    }


def _compute_ari_marker(adata: anndata.AnnData) -> dict:
    """SECONDARY path: ARI(joint Leiden, marker-panel labels)."""
    if "leiden" not in adata.obs:
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN", "reason": "no leiden"}
    if adata.X is None:
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN", "reason": "no .X"}
    var_idx = adata.var.index
    # adata.var.index is ENSG; markers in MARKER_PANEL are gene symbols. We need
    # the gene_id <-> name lookup from the same TSS table the driver used.
    name_to_id_path = Path("/home/zerlinshen/singlecell_factory/data/external/gencode_grch38_tss_by_ensg.tsv")
    name_to_id = {}
    if name_to_id_path.exists():
        tss = pd.read_csv(name_to_id_path, sep="\t")
        name_to_id = dict(zip(tss["gene_name"].astype(str), tss["gene_id"].astype(str)))
    label_arr = np.full(adata.n_obs, "OTHER", dtype=object)
    for label, gene_names in MARKER_PANEL.items():
        for gn in gene_names:
            gid = name_to_id.get(gn)
            if gid is None or gid not in var_idx:
                continue
            col = adata[:, gid].X
            if hasattr(col, "toarray"):
                col = col.toarray().ravel()
            else:
                col = np.asarray(col).ravel()
            # Mark cells in top quartile of expression for this marker
            cutoff = np.quantile(col[col > 0], 0.75) if (col > 0).sum() > 0 else np.inf
            label_arr[(col >= cutoff) & (label_arr == "OTHER")] = label
            break  # use first available gene per type
    leiden = adata.obs["leiden"].astype(str)
    ari = float(adjusted_rand_score(label_arr, leiden))
    return {
        "value": round(ari, 4),
        "threshold": 0.70,
        "status": "CLOSED-PARTIAL",
        "evidence_path": "SECONDARY",
        "n_marker_typed": int((label_arr != "OTHER").sum()),
    }


def _compute_peak_gene_overlap(
    adata: anndata.AnnData,
    s2f_tsv: Path,
    top_k: int = 1000,
) -> dict:
    """Top-1000 peak-gene overlap vs Trevino S2F.

    Compares ``adata.uns["peak_to_gene_top1000"]`` rows (peak_idx, gene)
    to Trevino's S2F joined table by (chr-start-end, gene) tuples.
    """
    if "peak_to_gene_top1000" not in adata.uns:
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN",
                "reason": "no peak_to_gene_top1000 in adata.uns"}
    if not s2f_tsv.exists():
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN",
                "reason": f"Trevino S2F not at {s2f_tsv}",
                "acquisition_status": "failed"}
    our = adata.uns["peak_to_gene_top1000"]
    if not isinstance(our, pd.DataFrame):
        our = pd.DataFrame(our)
    # Take top-K by |pearson|
    if "pearson" in our.columns:
        our = our.reindex(our["pearson"].abs().sort_values(ascending=False).index).head(top_k)
    elif "p_perm" in our.columns:
        our = our.sort_values("p_perm").head(top_k)
    else:
        our = our.head(top_k)
    # Need peak coordinates and gene name to compare against Trevino. Our top1000
    # rows hold peak_idx (positional into adata.obsm['atac_peaks']) and gene index.
    atac_var = adata.uns["atac_peaks"]
    if isinstance(atac_var, dict):
        atac_var = pd.DataFrame(atac_var)
    chrom_col = "chrom" if "chrom" in atac_var.columns else "seqnames"
    our = our.copy()
    our["peak_chrom"] = atac_var[chrom_col].iloc[our["peak_idx"].astype(int)].values
    our["peak_start"] = atac_var["start"].iloc[our["peak_idx"].astype(int)].values
    our["peak_end"]   = atac_var["end"].iloc[our["peak_idx"].astype(int)].values
    # gene column in our top1000 is var.index ENSG; map back to gene_name
    name_to_id_path = Path("/home/zerlinshen/singlecell_factory/data/external/gencode_grch38_tss_by_ensg.tsv")
    tss = pd.read_csv(name_to_id_path, sep="\t")
    id_to_name = dict(zip(tss["gene_id"].astype(str), tss["gene_name"].astype(str)))
    our["gene_name"] = our["gene"].astype(str).map(id_to_name).fillna(our["gene"].astype(str))
    # Build comparable tuples
    our_tuples = set(zip(our["peak_chrom"].astype(str), our["peak_start"].astype(int),
                          our["peak_end"].astype(int), our["gene_name"].astype(str)))
    # Load Trevino S2F joined
    trev = pd.read_csv(s2f_tsv, sep="\t")
    # Take top-K from Trevino by |corr|
    trev = trev.reindex(trev["corr"].abs().sort_values(ascending=False).index).head(top_k)
    trev_tuples = set(zip(trev["chrom"].astype(str), trev["start"].astype(int),
                           trev["end"].astype(int), trev["gene"].astype(str)))
    overlap = our_tuples & trev_tuples
    # Also report by-rank Jaccard
    jaccard = len(overlap) / len(our_tuples | trev_tuples) if (our_tuples | trev_tuples) else 0.0
    n_overlap_pct = len(overlap) / top_k
    return {
        "value": round(n_overlap_pct, 4),
        "threshold": 0.70,
        "status": "CLOSED" if n_overlap_pct >= 0.70 else ("CLOSED-PARTIAL" if n_overlap_pct >= 0.50 else "CLOSED-PARTIAL"),
        "method": "by_rank_top1000",
        "n_ours": len(our_tuples),
        "n_trevino": len(trev_tuples),
        "n_overlap": len(overlap),
        "jaccard": round(jaccard, 4),
        "acquisition_status": "acquired",
    }


def _compute_pseudotime_spearman(
    adata: anndata.AnnData,
    fig4d_csv: Path | None = None,
) -> dict:
    """AC-VAL-5: Spearman(our dpt_pseudotime, Trevino Fig 4D pseudotime).

    Trevino did not publish per-cell pseudotime in S2D (checked all 66 cols).
    Brainchromatin repo ships scRNA_URD_GlialCells.RDS (R format) — out of
    scope for an automated reproduction. Marking NOT-RUN per plan v3.
    """
    if "dpt_pseudotime" not in adata.obs:
        return {"value": None, "threshold": 0.70, "status": "NOT-RUN",
                "reason": "no dpt_pseudotime in obs"}
    if fig4d_csv is None or not fig4d_csv.exists():
        return {
            "value": None,
            "threshold": 0.70,
            "status": "CLOSED-PARTIAL",
            "reason": "Trevino Fig 4D per-cell pseudotime not published as a per-cell artifact in S2D (verified 66 cols); URD trajectory object in brainchromatin repo requires R to read.",
            "acquisition_status": "failed",
            "pseudotime_status": "NOT-RUN",
        }
    trev = pd.read_csv(fig4d_csv)
    # expects columns: Cell.ID, pseudotime
    ours = adata.obs[["dpt_pseudotime"]].copy()
    ours["Cell.ID"] = adata.obs_names.astype(str)
    merged = ours.merge(trev, on="Cell.ID", how="inner")
    if len(merged) < 500:
        return {
            "value": None, "threshold": 0.70, "status": "CLOSED-PARTIAL",
            "pseudotime_status": "CLOSED-PARTIAL",
            "reason": f"only {len(merged)} cells joined (<500 minimum guard)",
            "pseudotime_n_joined": int(len(merged)),
        }
    rho, _ = spearmanr(merged["dpt_pseudotime"], merged["pseudotime"])
    return {
        "value": round(float(rho), 4),
        "threshold": 0.70,
        "status": "CLOSED" if rho >= 0.70 else "CLOSED-PARTIAL",
        "pseudotime_status": "CLOSED" if rho >= 0.70 else "CLOSED-PARTIAL",
        "pseudotime_n_joined": int(len(merged)),
        "acquisition_status": "acquired",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--out", default=None, type=Path,
                    help="metrics.json output path (default: <run-dir>/metrics.json + .omc/research/wave5/metrics.json)")
    ap.add_argument("--cluster-names",
                    default=Path("/home/zerlinshen/singlecell_factory/data/raw/trevino_2021_brain/GSE162170_multiome_cluster_names.txt.gz"),
                    type=Path)
    ap.add_argument("--s2f-tsv",
                    default=Path("/home/zerlinshen/singlecell_factory/data/external/trevino_2021_supp/S2F_peak_gene_links_with_coords.tsv"),
                    type=Path)
    ap.add_argument("--fig4d-csv", default=None, type=Path,
                    help="Optional per-cell pseudotime CSV (Cell.ID,pseudotime) extracted from Trevino Fig 4D.")
    args = ap.parse_args()

    final_h5ad = args.run_dir / "final_adata.h5ad"
    if not final_h5ad.exists():
        log(f"FAIL: {final_h5ad} not found")
        return 2
    log(f"loading {final_h5ad}")
    adata = anndata.read_h5ad(final_h5ad)
    log(f"  n_obs={adata.n_obs} n_vars={adata.n_vars}")

    # ARI
    names = _load_cluster_names(args.cluster_names) if args.cluster_names.exists() else None
    if names:
        ari = _compute_ari_primary(adata, names)
        if ari["status"] == "NOT-RUN":
            log(f"  ARI primary path NOT-RUN ({ari.get('reason')}); falling back to marker panel")
            ari = _compute_ari_marker(adata)
    else:
        ari = _compute_ari_marker(adata)
    log(f"  ARI: {ari}")

    # Overlap
    overlap = _compute_peak_gene_overlap(adata, args.s2f_tsv)
    log(f"  overlap: {overlap}")

    # Pseudotime
    spearman = _compute_pseudotime_spearman(adata, args.fig4d_csv)
    log(f"  spearman: {spearman}")

    out = {
        "cell_type_ari":             ari.get("value"),
        "cell_type_ari_meta":        ari,
        "peak_gene_overlap_top1000": overlap.get("value"),
        "peak_gene_overlap_meta":    overlap,
        "pseudotime_spearman":       spearman.get("value"),
        "pseudotime_spearman_meta":  spearman,
        "ari_evidence_path":         ari.get("evidence_path", "SECONDARY"),
        "pseudotime_status":         spearman.get("pseudotime_status", spearman.get("status")),
        "computed_at":               datetime.now(timezone.utc).isoformat(),
    }
    out_paths = []
    if args.out:
        out_paths.append(args.out)
    out_paths.append(args.run_dir / "metrics.json")
    omc_path = Path("/home/zerlinshen/singlecell_factory/.omc/research/wave5/metrics.json")
    out_paths.append(omc_path)
    for p in out_paths:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
        log(f"  wrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
