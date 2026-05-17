"""Wave-5.5 v5.0.4 Stage C4 — quantitative similarity metrics + reproduction-rate.

# STATUS: one-off, promote-on-reuse

Per AC-V5-COMP-QUANT-1 + AC-V5-REPRO-RATE-1:
- Per-figure metrics computed using LITERAL formulas pinned in plan §4.1.
- Per-metric normalization to [0,1]:
    Procrustes RMS:  score = 1 - min(rms/0.2, 1.0)
    ARI:             score = max(0.0, ari)
    NMI:             score = nmi
    diag_frac:       score = np.trace(C_normalized) / n_classes
    Spearman:        score = max(0.0, rho)
    overlap:         score = overlap
- Per-figure aggregation: arithmetic mean of applicable metrics (per FIGURE_SPEC.md §reproduction_rate).
- Overall: equal-weighted mean across Stream A main + Stream A2 (excludes Supplementary F-5
  and NOT-RUN placeholders from numerator AND denominator).

Outputs:
- <canonical-run>/python/figures/v5/quantitative_metrics/<fig-id>.json
- <canonical-run>/python/figures/v5/reproduction_rate.json
"""

from __future__ import annotations

import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_RUN = Path("/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0")
ADATA_PATH = CANONICAL_RUN / "python" / "final_adata.h5ad"
V42_LEDGER = REPO_ROOT / "ops" / "run_ledger" / "wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json"
THREE_WAY = CANONICAL_RUN / "python" / "figures" / "v5" / "rna_only" / "three_way_comparison.json"
OUT_DIR = CANONICAL_RUN / "python" / "figures" / "v5" / "quantitative_metrics"
REPRO_RATE_OUT = CANONICAL_RUN / "python" / "figures" / "v5" / "reproduction_rate.json"

# Applicable metrics per figure per FIGURE_SPEC.md §reproduction_rate.per_figure_metrics
FIGURE_METRICS = {
    "F-1B": ["procrustes", "ari", "nmi", "diag_frac"],
    "F-1C": ["spearman"],
    "F-1D": ["spearman"],
    "F-2A": ["procrustes"],
    "F-2BD": ["spearman"],
    "F-4A": ["spearman"],
    "F-4D": ["spearman"],
    "F-4EG": ["spearman"],
    "F-5": ["overlap", "spearman"],  # supplementary
    "F-RNA-1B": ["ari", "nmi", "diag_frac"],
    "F-RNA-1C": ["spearman"],
    "F-RNA-4A": ["spearman"],
}
NOT_RUN_FIGURES = {"F-3", "F-6", "F-7"}
SUPPLEMENTARY_FIGURES = {"F-5"}


# ---------- per-metric normalization ----------

def norm_procrustes(rms: float) -> float:
    return 1.0 - min(rms / 0.2, 1.0)


def norm_ari(ari: float) -> float:
    return max(0.0, ari)


def norm_nmi(nmi: float) -> float:
    return nmi


def norm_diag_frac(diag_frac: float) -> float:
    return diag_frac


def norm_spearman(rho: float) -> float:
    return max(0.0, rho)


def norm_overlap(overlap: float) -> float:
    return overlap


# ---------- metric computers ----------

def compute_confusion_diag_frac(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    """Per AC-V5-REPRO-RATE-1: np.trace(C_normalized) / n_classes (row-normalized confusion)."""
    from sklearn.metrics import confusion_matrix
    classes = sorted(set(labels_a) | set(labels_b))
    cm = confusion_matrix(labels_a, labels_b, labels=classes).astype(float)
    row_sums = cm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    cm_norm = cm / row_sums
    n_classes = len(classes)
    if n_classes == 0:
        return 0.0
    return float(np.trace(cm_norm) / n_classes)


def compute_metrics_for_figure(fig_id: str, ctx: dict) -> dict:
    """Compute applicable metrics for a single figure. ctx contains adata, v42_metrics, three_way."""
    metrics_needed = FIGURE_METRICS.get(fig_id, [])
    metrics_raw: dict = {}
    metrics_scored: dict = {}
    rationale = {}

    adata = ctx["adata"]
    v42_metrics = ctx["v42_metrics"]
    three_way = ctx["three_way"]

    for m in metrics_needed:
        if m == "procrustes":
            # Trevino's published UMAP/LSI coordinates not available locally.
            # NOT-COMPUTED with rationale per FIGURE_SPEC.md NOT-RUN precedent.
            metrics_raw[m] = None
            metrics_scored[m] = None
            rationale[m] = "Trevino's published UMAP/LSI coordinates not available locally (paper supplementary tables would need download); Procrustes alignment NOT-COMPUTED. Per per-figure aggregation: excluded from mean."
        elif m == "ari":
            if fig_id.startswith("F-RNA-"):
                ari_val = three_way["RNA-only-vs-Trevino"]["ari"]
                metrics_raw[m] = ari_val
                metrics_scored[m] = norm_ari(ari_val)
                rationale[m] = "ARI(RNA-only Leiden vs Trevino seurat_clusters) from three_way_comparison.json"
            else:
                # Multiome figures: ARI(WNN-leiden vs Trevino seurat_clusters)
                # v4.2 ledger stores cell_type_ari as {'value': ..., 'threshold': ..., 'status': ...}
                v42_ari = v42_metrics.get("cell_type_ari")
                if isinstance(v42_ari, dict):
                    ari_val = float(v42_ari.get("value", three_way["WNN-vs-Trevino"]["ari"]))
                elif v42_ari is not None:
                    ari_val = float(v42_ari)
                else:
                    ari_val = float(three_way["WNN-vs-Trevino"]["ari"])
                metrics_raw[m] = ari_val
                metrics_scored[m] = norm_ari(ari_val)
                rationale[m] = "ARI(WNN leiden vs Trevino seurat_clusters) from three_way_comparison.json / v4.2 ledger"
        elif m == "nmi":
            if fig_id.startswith("F-RNA-"):
                nmi_val = three_way["RNA-only-vs-Trevino"]["nmi"]
            else:
                nmi_val = three_way["WNN-vs-Trevino"]["nmi"]
            metrics_raw[m] = nmi_val
            metrics_scored[m] = norm_nmi(nmi_val)
            rationale[m] = "NMI from three_way_comparison.json"
        elif m == "diag_frac":
            if fig_id.startswith("F-RNA-"):
                # Need to load RNA-only Leiden labels — store via three_way computation.
                # Approximation: compute on-the-fly between WNN-leiden and Trevino as proxy.
                trev = adata.obs["seurat_clusters"].astype(str).values
                wnn = adata.obs["leiden"].astype(str).values
                diag = compute_confusion_diag_frac(trev, wnn)
                metrics_raw[m] = diag
                metrics_scored[m] = norm_diag_frac(diag)
                rationale[m] = "diag_frac(WNN leiden, Trevino seurat_clusters) — RNA-only Leiden labels not persisted; using WNN as proxy for confusion-matrix structure"
            else:
                trev = adata.obs["seurat_clusters"].astype(str).values
                wnn = adata.obs["leiden"].astype(str).values
                diag = compute_confusion_diag_frac(trev, wnn)
                metrics_raw[m] = diag
                metrics_scored[m] = norm_diag_frac(diag)
                rationale[m] = "diag_frac(WNN leiden, Trevino seurat_clusters) — row-normalized confusion matrix trace per AC-V5-REPRO-RATE-1"
        elif m == "spearman":
            # Trevino's per-cell measurements (gene expression / pseudotime / peak accessibility)
            # not available externally with per-cell alignment to ours.
            metrics_raw[m] = None
            metrics_scored[m] = None
            rationale[m] = "Trevino per-cell measurements (gene expression / pseudotime / peak accessibility) not available in cell-aligned form; Spearman NOT-COMPUTED. Per per-figure aggregation: excluded from mean."
        elif m == "overlap":
            # F-5: from v4.2 ledger peak_gene_overlap_top1000_strict
            ov = v42_metrics.get("peak_gene_overlap_top1000_strict")
            if ov is None:
                # try alternate location
                ov = ctx["v42_full"].get("peak_gene_overlap_top1000_strict")
            if ov is None:
                metrics_raw[m] = None
                metrics_scored[m] = None
                rationale[m] = "peak_gene_overlap_top1000_strict not present in v4.2 ledger; NOT-COMPUTED"
            else:
                metrics_raw[m] = float(ov)
                metrics_scored[m] = norm_overlap(float(ov))
                rationale[m] = "from v4.2 ledger peak_gene_overlap_top1000_strict"
    # per-figure aggregation: arithmetic mean of available scored metrics
    available = [v for v in metrics_scored.values() if v is not None]
    if available:
        per_fig_score = float(np.mean(available))
    else:
        per_fig_score = None
    return {
        "fig_id": fig_id,
        "metrics_raw": metrics_raw,
        "metrics_scored": metrics_scored,
        "rationale": rationale,
        "per_figure_score": per_fig_score,
        "n_scored_metrics": len(available),
        "n_total_metrics_declared": len(metrics_needed),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(ADATA_PATH)

    v42_full = json.loads(V42_LEDGER.read_text())
    v42_metrics = v42_full.get("metrics", {})
    three_way = json.loads(THREE_WAY.read_text())

    ctx = {
        "adata": adata,
        "v42_metrics": v42_metrics,
        "v42_full": v42_full,
        "three_way": three_way,
    }

    per_figure_results = {}
    all_figs = list(FIGURE_METRICS.keys()) + list(NOT_RUN_FIGURES)
    for fig_id in all_figs:
        if fig_id in NOT_RUN_FIGURES:
            res = {
                "fig_id": fig_id,
                "status": "NOT-RUN",
                "rationale": f"{fig_id} is a NOT-RUN placeholder per plan §1.3 (chromVAR/TF-network/GWAS out of v4.2 scope)",
                "per_figure_score": None,
                "n_scored_metrics": 0,
                "n_total_metrics_declared": 0,
            }
        else:
            res = compute_metrics_for_figure(fig_id, ctx)
        per_figure_results[fig_id] = res
        (OUT_DIR / f"{fig_id}.json").write_text(json.dumps(res, indent=2) + "\n")

    # ---- reproduction-rate aggregation ----
    # Overall: equal-weighted mean over Stream A main + Stream A2, EXCLUDE supplementary + NOT-RUN
    excluded_from_overall = SUPPLEMENTARY_FIGURES | NOT_RUN_FIGURES
    reproducible = {
        fid: r for fid, r in per_figure_results.items()
        if fid not in excluded_from_overall
    }
    reproducible_with_scores = {fid: r for fid, r in reproducible.items() if r.get("per_figure_score") is not None}
    if reproducible_with_scores:
        overall = float(np.mean([r["per_figure_score"] for r in reproducible_with_scores.values()]))
    else:
        overall = None

    supplementary_score = per_figure_results.get("F-5", {}).get("per_figure_score")

    repro = {
        "schema_version": "v5.0.4",
        "computed_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall": overall,
        "not_run_count": len(NOT_RUN_FIGURES),
        "supplementary_score": supplementary_score,
        "per_figure": {fid: r.get("per_figure_score") for fid, r in per_figure_results.items()},
        "per_figure_n_scored_over_n_declared": {
            fid: f"{r.get('n_scored_metrics', 0)}/{r.get('n_total_metrics_declared', 0)}"
            for fid, r in per_figure_results.items()
        },
        "aggregation_formula_source": "plan v5.0.4 §4.1 AC-V5-REPRO-RATE-1 (literal) + FIGURE_SPEC.md §reproduction_rate (duplicate)",
        "exclusion_policy": "Supplementary (F-5) and NOT-RUN (F-3/F-6/F-7) excluded from overall numerator AND denominator per AC-V5-REPRO-RATE-1",
        "limitations": [
            "Procrustes (UMAP/LSI coord alignment): NOT-COMPUTED — Trevino published coords not locally available; affects F-1B, F-2A.",
            "Per-cell Spearman (gene expr / pseudotime / peak accessibility): NOT-COMPUTED — Trevino per-cell measurements not in cell-aligned form; affects F-1C, F-1D, F-2BD, F-4A, F-4D, F-4EG, F-5, F-RNA-1C, F-RNA-4A.",
            "F-RNA diag_frac uses WNN-leiden as proxy since RNA-only Leiden labels were not persisted post-Stage A2.",
        ],
        "interpretation_note": "Score is a DIAGNOSTIC AID, not a binary pass/fail. Per-metric breakdown rendered alongside aggregate on every comparison.pdf page. NOT-RUN figures called out separately on cover page. Many figures have only partial metric coverage due to Trevino data accessibility limitations; per-figure n_scored/n_declared indicates this.",
    }
    REPRO_RATE_OUT.write_text(json.dumps(repro, indent=2) + "\n")
    print(f"Overall reproduction-rate: {overall}")
    print(f"Per-figure scores: {repro['per_figure']}")
    print(f"Per-figure coverage (scored/declared): {repro['per_figure_n_scored_over_n_declared']}")
    print(f"Output: {REPRO_RATE_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
