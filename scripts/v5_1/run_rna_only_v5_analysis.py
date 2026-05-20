"""Wave-5.5 v5.0.4 Stream A2 RNA-only analysis (Stage A2).

# STATUS: one-off, promote-on-reuse

Per plan v5.0.4 §5 Stage A2:
- Load <canonical-run>/python/final_adata.h5ad, restrict to RNA modality only.
- BEFORE recompute (Architect iter-3 Fix-4 anti-contamination): delete
  obsm.X_wnn, X_umap, X_pca, X_diffmap, obsp.connectivities, distances,
  obs.dpt_pseudotime, obsm.atac_peaks. Assert post-delete set is empty.
- Recompute on RNA-only matrix with literal params pinned in plan §5 Stage A2 step 2 +
  FIGURE_SPEC.md §rna_only (sourced from Trevino 2021 Methods + scanpy defaults):
    normalize_total(1e4) -> log1p
    HVG(3000, flavor=seurat_v3)
    PCA(50)
    neighbors(n_neighbors=30, use_rep=X_pca)
    Leiden(resolution=0.3)
    UMAP(n_components=2, min_dist=0.5, random_state=42)
- Write processed h5ad with obsm["X_pca"], obsm["X_umap"], obs["leiden_rna_only"].
- Emit rna_only_metrics.json with ARI/NMI numbers:
    {"rna_only_vs_trevino_ari": ..., "wnn_vs_trevino_ari": ..., ...NMI...}
    plus parameter dict and v4.2 published baselines.

Trevino reference clusters: `obs.seurat_clusters` (14 categories, published).
WNN reference clusters:     `obs.leiden`         (53 categories, v4.2 WNN-based).

Outputs (written to OUT_ROOT):
  processed_rna_only.h5ad       — h5ad with X_pca, X_umap, leiden_rna_only
  rna_only_metrics.json         — ARI/NMI + params
  _anti_contamination_report.json
"""

from __future__ import annotations

import hashlib
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import scanpy as sc

CANONICAL_RUN = Path("/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0")
ADATA_PATH = CANONICAL_RUN / "python" / "final_adata.h5ad"

REPO_ROOT = Path(__file__).resolve().parents[3]
FIGURE_SPEC = REPO_ROOT / ".omc" / "research" / "wave5" / "FIGURE_SPEC.md"
FIGURE_SPEC_SIGNOFF_SHA = "6d0b1881448bf8ca832e9dd7d0373c4f03d897ce3f027751eb48a20743f501af"

OUT_ROOT = CANONICAL_RUN / "python" / "rna_only_analysis"

# Plan §5 Stage A2 step 2: LITERAL params (no Stage-0 author discretion)
RNA_ONLY_PARAMS = {
    "normalize_target_sum": 1e4,
    "hvg_n_top": 3000,
    "hvg_flavor": "seurat_v3",
    "n_pcs": 50,
    "n_neighbors": 30,
    "leiden_resolution": 0.3,
    "umap_n_components": 2,
    "umap_min_dist": 0.5,
    "random_state": 42,
}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def assert_spec_hash() -> str:
    cur = sha256_file(FIGURE_SPEC)
    if cur != FIGURE_SPEC_SIGNOFF_SHA:
        raise RuntimeError(
            f"FIGURE_SPEC.md SHA mismatch! signoff={FIGURE_SPEC_SIGNOFF_SHA} cur={cur}"
        )
    return cur


def strip_wnn_contamination(adata) -> dict:
    """Architect iter-3 Fix-4: explicit anti-contamination deletes BEFORE any RNA-only recompute.

    Enumerate actual keys first; assert post-delete set is empty.
    """
    obsm_to_delete = [k for k in adata.obsm.keys()
                      if k in {"X_wnn", "X_umap", "X_pca", "X_diffmap", "atac_peaks"}]
    obsp_to_delete = [k for k in adata.obsp.keys()
                      if k in {"connectivities", "distances"}]
    obs_to_delete = [k for k in adata.obs.columns
                     if k in {"dpt_pseudotime"}]
    uns_to_delete = [k for k in adata.uns.keys()
                     if k in {"neighbors", "umap", "pca", "diffmap_evals", "iroot", "paga",
                              "peak_to_gene", "peak_to_gene_linkages", "peak_to_gene_top1000",
                              "multimodal_status"}]

    before = {
        "obsm": list(adata.obsm.keys()),
        "obsp": list(adata.obsp.keys()),
        "obs": [c for c in adata.obs.columns
                if c in {"dpt_pseudotime", "Sample.Age", "seurat_clusters", "leiden"}],
        "uns": list(adata.uns.keys()),
    }

    for k in obsm_to_delete:
        del adata.obsm[k]
    for k in obsp_to_delete:
        del adata.obsp[k]
    for k in obs_to_delete:
        del adata.obs[k]
    for k in uns_to_delete:
        del adata.uns[k]

    # Post-delete assertion per plan §5 Stage A2 step 1: WNN-derived/ATAC keys gone
    forbidden_after = {"X_wnn", "X_umap", "X_pca", "X_diffmap", "atac_peaks"}
    remaining_forbidden = forbidden_after & set(adata.obsm.keys())
    assert not remaining_forbidden, f"anti-contamination FAILED: {remaining_forbidden} still in obsm"
    forbidden_obsp = {"connectivities", "distances"} & set(adata.obsp.keys())
    assert not forbidden_obsp, f"anti-contamination FAILED: {forbidden_obsp} still in obsp"

    return {"before": before, "deleted": {
        "obsm": obsm_to_delete, "obsp": obsp_to_delete,
        "obs": obs_to_delete, "uns": uns_to_delete,
    }}


def recompute_rna_only(adata):
    """Recompute embeddings from RNA-only counts per FIGURE_SPEC.md §rna_only."""
    # Restore counts layer if available; otherwise X is already normalized.
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=RNA_ONLY_PARAMS["normalize_target_sum"])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=RNA_ONLY_PARAMS["hvg_n_top"],
        flavor=RNA_ONLY_PARAMS["hvg_flavor"],
        layer="counts" if "counts" in adata.layers else None,
    )
    sc.tl.pca(adata, n_comps=RNA_ONLY_PARAMS["n_pcs"], random_state=RNA_ONLY_PARAMS["random_state"])
    sc.pp.neighbors(adata, n_neighbors=RNA_ONLY_PARAMS["n_neighbors"],
                    use_rep="X_pca", random_state=RNA_ONLY_PARAMS["random_state"])
    sc.tl.leiden(adata, resolution=RNA_ONLY_PARAMS["leiden_resolution"],
                 random_state=RNA_ONLY_PARAMS["random_state"],
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, n_components=RNA_ONLY_PARAMS["umap_n_components"],
               min_dist=RNA_ONLY_PARAMS["umap_min_dist"],
               random_state=RNA_ONLY_PARAMS["random_state"])
    # Rename leiden -> leiden_rna_only to avoid collision with WNN leiden in original
    adata.obs["leiden_rna_only"] = adata.obs["leiden"].copy()
    return adata


def compute_three_way_comparison(adata_rna, adata_orig) -> dict:
    """ARI/NMI for {RNA-only-vs-Trevino, WNN-vs-Trevino, RNA-only-vs-WNN}."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    rna_only_clusters = adata_rna.obs["leiden_rna_only"].astype(str).values
    trevino_clusters = adata_orig.obs["seurat_clusters"].astype(str).values
    wnn_clusters = adata_orig.obs["leiden"].astype(str).values

    metrics = {
        "rna_only_vs_trevino_ari": float(adjusted_rand_score(trevino_clusters, rna_only_clusters)),
        "rna_only_vs_trevino_nmi": float(
            normalized_mutual_info_score(trevino_clusters, rna_only_clusters)
        ),
        "wnn_vs_trevino_ari": float(adjusted_rand_score(trevino_clusters, wnn_clusters)),
        "wnn_vs_trevino_nmi": float(
            normalized_mutual_info_score(trevino_clusters, wnn_clusters)
        ),
        "rna_only_vs_wnn_ari": float(adjusted_rand_score(wnn_clusters, rna_only_clusters)),
        "rna_only_vs_wnn_nmi": float(
            normalized_mutual_info_score(wnn_clusters, rna_only_clusters)
        ),
        "n_clusters_rna_only": int(len(set(rna_only_clusters))),
        "n_clusters_trevino": int(len(set(trevino_clusters))),
        "n_clusters_wnn": int(len(set(wnn_clusters))),
        # v4.2 published baselines for cross-check
        "v4_2_published_baselines": {
            "RNA-only-vs-Trevino_ARI_baseline": 0.6736,
            "WNN-vs-Trevino_ARI_baseline": 0.168,
            "source": "ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json",
        },
        "params": RNA_ONLY_PARAMS,
    }
    return metrics


def main() -> int:
    spec_sha = assert_spec_hash()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata_orig = ad.read_h5ad(ADATA_PATH)
        adata = adata_orig.copy()

    # Architect iter-3 Fix-4: anti-contamination deletes BEFORE recompute
    contamination_report = strip_wnn_contamination(adata)
    print(f"anti-contamination deletes: {contamination_report['deleted']}")
    (OUT_ROOT / "_anti_contamination_report.json").write_text(
        json.dumps(contamination_report, indent=2) + "\n", encoding="utf-8"
    )

    print("recompute (normalize -> HVG -> PCA -> neighbors -> Leiden -> UMAP)...")
    adata = recompute_rna_only(adata)
    print(f"RNA-only Leiden produced {adata.obs['leiden_rna_only'].nunique()} clusters")

    print("compute three-way comparison metrics (ARI/NMI)...")
    metrics = compute_three_way_comparison(adata, adata_orig)
    (OUT_ROOT / "rna_only_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
    )
    print("metrics:")
    for k, v in metrics.items():
        if k not in ("v4_2_published_baselines", "params"):
            print(f"  {k}: {v}")

    print("writing processed h5ad...")
    out_h5ad = OUT_ROOT / "processed_rna_only.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")

    # Manifest
    manifest = {
        "figure_spec_sha": spec_sha,
        "produced_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "analysis_script_sha": sha256_file(Path(__file__)),
        "outputs": {
            "processed_h5ad": str(out_h5ad),
            "metrics_json": str(OUT_ROOT / "rna_only_metrics.json"),
            "anti_contamination_report": str(OUT_ROOT / "_anti_contamination_report.json"),
        },
        "obsm_populated": ["X_pca", "X_umap"],
        "obs_populated": ["leiden_rna_only"],
    }
    (OUT_ROOT / "_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\nRNA-only analysis complete. Outputs at: {OUT_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
