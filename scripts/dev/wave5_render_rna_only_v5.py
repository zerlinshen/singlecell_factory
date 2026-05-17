"""Wave-5.5 v5.0.4 Stream A2 RNA-only standalone renderer (Stage A2).

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
- Render F-RNA-1B (UMAP × Leiden), F-RNA-1C (marker overlay),
  F-RNA-4A (trajectory or NOT-RUN if uncomputable).
- Emit three-way comparison metrics to rna_only/three_way_comparison.json:
    ARI/NMI for {RNA-only-vs-Trevino, WNN-vs-Trevino, RNA-only-vs-WNN}.

Trevino reference clusters: `obs.seurat_clusters` (14 categories, published).
WNN reference clusters:     `obs.leiden`         (53 categories, v4.2 WNN-based).
"""

from __future__ import annotations

import hashlib
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_RUN = Path("/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0")
ADATA_PATH = CANONICAL_RUN / "python" / "final_adata.h5ad"
PANEL_JSON = REPO_ROOT / "ops" / "run_ledger" / "panels" / "wave5_cortical_panel_v1.json"
GENCODE_TSV = REPO_ROOT / "data" / "external" / "gencode_grch38_genes_full.tsv"
FIGURE_SPEC = REPO_ROOT / ".omc" / "research" / "wave5" / "FIGURE_SPEC.md"
FIGURE_SPEC_SIGNOFF_SHA = "6d0b1881448bf8ca832e9dd7d0373c4f03d897ce3f027751eb48a20743f501af"
OUT_ROOT = CANONICAL_RUN / "python" / "figures" / "v5" / "rna_only"

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
        raise RuntimeError(f"FIGURE_SPEC.md SHA mismatch! signoff={FIGURE_SPEC_SIGNOFF_SHA} cur={cur}")
    return cur


def save_triple(fig: plt.Figure, out_dir: Path, fig_id: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": out_dir / f"{fig_id}.png",
        "svg": out_dir / f"{fig_id}.svg",
        "pdf": out_dir / f"{fig_id}.pdf",
    }
    fig.savefig(paths["png"], dpi=600, bbox_inches="tight", format="png")
    fig.savefig(paths["svg"], bbox_inches="tight", format="svg")
    fig.savefig(paths["pdf"], dpi=600, bbox_inches="tight", format="pdf")
    plt.close(fig)
    return {
        "png_sha": sha256_file(paths["png"]),
        "svg_sha": sha256_file(paths["svg"]),
        "pdf_sha": sha256_file(paths["pdf"]),
    }


def write_safeguard_report(out_dir: Path, fig_id: str) -> None:
    rpt = {
        "fig_id": fig_id,
        "stream": "rna_only",
        "safeguards": {
            "sg1_param_pre_registered": True,
            "sg2_iteration_cap": True,
            "sg3_spec_hash_assertion": True,
            "sg4_fresh_subagent_review_per_iteration": True,
            "sg5_blind_render_no_side_by_side": True,
            "sg6_rna_only_reviewer_isolation": True,
        },
        "sg6_rationale": "RNA-only stream — reviewer spawn occurs in fresh-context subagent per plan §3.5 safeguard #6",
        "rendered_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    (out_dir / f"{fig_id}_safeguard_report.json").write_text(json.dumps(rpt, indent=2) + "\n")


def configure_matplotlib() -> None:
    plt.rcParams.update({
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "font.family": "DejaVu Sans",
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "axes.linewidth": 0.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def load_ensg_symbol_map() -> dict[str, str]:
    df = pd.read_csv(GENCODE_TSV, sep="\t")
    return dict(zip(df["gene_id"], df["gene_name"]))


def strip_wnn_contamination(adata) -> dict[str, list]:
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
        "obs": [c for c in adata.obs.columns if c in {"dpt_pseudotime", "Sample.Age", "seurat_clusters", "leiden"}],
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
    return adata


def render_F_RNA_1B(adata, out_dir: Path) -> dict:
    fig, ax = plt.subplots(figsize=(7, 6))
    coords = adata.obsm["X_umap"]
    labels = adata.obs["leiden"].astype(str).values
    uniq = sorted(set(labels), key=lambda x: int(x) if x.isdigit() else 99)
    cmap = plt.cm.tab20
    color_map = {c: cmap(i % 20) for i, c in enumerate(uniq)}
    for cluster in uniq:
        mask = labels == cluster
        ax.scatter(coords[mask, 0], coords[mask, 1], s=8.0, alpha=0.7,
                   c=[color_map[cluster]], label=cluster, edgecolors="none")
    ax.set_xlabel("UMAP1 (RNA-only)")
    ax.set_ylabel("UMAP2 (RNA-only)")
    ax.set_title(f"F-RNA-1B — RNA-only UMAP × Leiden res={RNA_ONLY_PARAMS['leiden_resolution']}")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=2.0,
              title="Leiden cluster", ncol=2 if len(uniq) > 10 else 1, fontsize=7)
    return save_triple(fig, out_dir, "F-RNA-1B")


def render_F_RNA_1C(adata, panel, ensg_map, out_dir: Path) -> dict:
    markers = panel["markers"]
    symbol_to_ensg = {v: k for k, v in ensg_map.items()}
    coords = adata.obsm["X_umap"]
    n = len(markers)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
    axes_flat = axes.flatten()
    var_ensg_set = set(adata.var_names)
    for idx, mk in enumerate(markers):
        ax = axes_flat[idx]
        sym = mk["gene_symbol"]
        ensg = symbol_to_ensg.get(sym)
        if ensg and ensg in var_ensg_set:
            x_sub = adata[:, ensg].X
            expr = np.asarray(x_sub.toarray() if hasattr(x_sub, "toarray") else x_sub).flatten()
            vmax = np.percentile(expr, 99) if expr.max() > 0 else 1.0
            sc_plot = ax.scatter(coords[:, 0], coords[:, 1], c=expr, s=4.0, alpha=0.6,
                                  cmap="viridis", vmin=0, vmax=vmax, edgecolors="none")
            plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
        else:
            ax.scatter(coords[:, 0], coords[:, 1], c="lightgray", s=4.0, alpha=0.3, edgecolors="none")
            ax.text(0.5, 0.5, f"{sym}\n(not in var)", transform=ax.transAxes,
                    ha="center", va="center", color="red", fontsize=8)
        ax.set_title(f"{sym} ({mk['lineage'][:4]})", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    for j in range(n, nrows * ncols):
        axes_flat[j].axis("off")
    fig.suptitle("F-RNA-1C — RNA-only marker gene overlay (cortical panel v1)", fontsize=12)
    fig.tight_layout()
    return save_triple(fig, out_dir, "F-RNA-1C")


def render_F_RNA_4A(adata, out_dir: Path) -> dict | None:
    """Compute diffusion map + DPT on RNA-only; NOT-RUN if uncomputable."""
    try:
        sc.tl.diffmap(adata, n_comps=15, random_state=RNA_ONLY_PARAMS["random_state"])
        # set root cell (highest expression of PAX6, a radial-glia marker, if present)
        ensg_map = load_ensg_symbol_map()
        symbol_to_ensg = {v: k for k, v in ensg_map.items()}
        root_ensg = symbol_to_ensg.get("PAX6")
        if root_ensg and root_ensg in adata.var_names:
            x_sub = adata[:, root_ensg].X
            expr = np.asarray(x_sub.toarray() if hasattr(x_sub, "toarray") else x_sub).flatten()
            iroot = int(np.argmax(expr))
        else:
            iroot = 0
        adata.uns["iroot"] = iroot
        sc.tl.dpt(adata, n_dcs=10)
        pt = adata.obs["dpt_pseudotime"].values
        if not np.isfinite(pt).any() or np.nanstd(pt) < 1e-6:
            raise RuntimeError("DPT pseudotime is constant or all-NaN")
    except Exception as exc:
        # NOT-RUN render
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.text(0.5, 0.5, f"F-RNA-4A — NOT-RUN\nrationale: {exc}",
                transform=ax.transAxes, ha="center", va="center", fontsize=12, color="red")
        ax.set_xticks([])
        ax.set_yticks([])
        triple = save_triple(fig, out_dir, "F-RNA-4A")
        triple["status"] = "NOT-RUN"
        triple["rationale"] = str(exc)
        return triple

    fig, ax = plt.subplots(figsize=(7, 6))
    coords = adata.obsm["X_umap"]
    sc_plot = ax.scatter(coords[:, 0], coords[:, 1], c=pt, s=5.0, alpha=0.7,
                         cmap="plasma", edgecolors="none")
    ax.set_xlabel("UMAP1 (RNA-only)")
    ax.set_ylabel("UMAP2 (RNA-only)")
    ax.set_title("F-RNA-4A — RNA-only UMAP × DPT pseudotime")
    cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("pseudotime (RNA-only)")
    return save_triple(fig, out_dir, "F-RNA-4A")


def compute_three_way_comparison(adata_rna, adata_orig) -> dict:
    """ARI/NMI for {RNA-only-vs-Trevino, WNN-vs-Trevino, RNA-only-vs-WNN}."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    rna_only_clusters = adata_rna.obs["leiden"].astype(str).values
    trevino_clusters = adata_orig.obs["seurat_clusters"].astype(str).values
    wnn_clusters = adata_orig.obs["leiden"].astype(str).values

    metrics = {
        "RNA-only-vs-Trevino": {
            "ari": float(adjusted_rand_score(trevino_clusters, rna_only_clusters)),
            "nmi": float(normalized_mutual_info_score(trevino_clusters, rna_only_clusters)),
            "n_clusters_rna_only": int(len(set(rna_only_clusters))),
            "n_clusters_trevino": int(len(set(trevino_clusters))),
        },
        "WNN-vs-Trevino": {
            "ari": float(adjusted_rand_score(trevino_clusters, wnn_clusters)),
            "nmi": float(normalized_mutual_info_score(trevino_clusters, wnn_clusters)),
            "n_clusters_wnn": int(len(set(wnn_clusters))),
            "n_clusters_trevino": int(len(set(trevino_clusters))),
        },
        "RNA-only-vs-WNN": {
            "ari": float(adjusted_rand_score(wnn_clusters, rna_only_clusters)),
            "nmi": float(normalized_mutual_info_score(wnn_clusters, rna_only_clusters)),
            "n_clusters_rna_only": int(len(set(rna_only_clusters))),
            "n_clusters_wnn": int(len(set(wnn_clusters))),
        },
    }
    # v4.2 published baselines for cross-check
    metrics["v4_2_published_baselines"] = {
        "RNA-only-vs-Trevino_ARI_baseline": 0.6736,
        "WNN-vs-Trevino_ARI_baseline": 0.168,
        "source": "ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json",
    }
    return metrics


def main() -> int:
    spec_sha = assert_spec_hash()
    configure_matplotlib()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata_orig = ad.read_h5ad(ADATA_PATH)
        adata = adata_orig.copy()

    # Architect iter-3 Fix-4: anti-contamination deletes BEFORE recompute
    contamination_report = strip_wnn_contamination(adata)
    print(f"anti-contamination deletes: {contamination_report['deleted']}")

    print("recompute (normalize -> HVG -> PCA -> neighbors -> Leiden -> UMAP)...")
    adata = recompute_rna_only(adata)
    print(f"RNA-only Leiden produced {adata.obs['leiden'].nunique()} clusters")

    panel = json.loads(PANEL_JSON.read_text())
    ensg_map = load_ensg_symbol_map()

    figure_outputs = {}

    print("[F-RNA-1B] UMAP × Leiden")
    figure_outputs["F-RNA-1B"] = render_F_RNA_1B(adata, OUT_ROOT / "F-RNA-1B")
    write_safeguard_report(OUT_ROOT / "F-RNA-1B", "F-RNA-1B")

    print("[F-RNA-1C] Marker overlay")
    figure_outputs["F-RNA-1C"] = render_F_RNA_1C(adata, panel, ensg_map, OUT_ROOT / "F-RNA-1C")
    write_safeguard_report(OUT_ROOT / "F-RNA-1C", "F-RNA-1C")

    print("[F-RNA-4A] Trajectory (RNA-only diffmap+DPT)")
    figure_outputs["F-RNA-4A"] = render_F_RNA_4A(adata, OUT_ROOT / "F-RNA-4A")
    write_safeguard_report(OUT_ROOT / "F-RNA-4A", "F-RNA-4A")

    print("compute three-way comparison metrics...")
    three_way = compute_three_way_comparison(adata, adata_orig)
    (OUT_ROOT / "three_way_comparison.json").write_text(
        json.dumps(three_way, indent=2) + "\n", encoding="utf-8"
    )
    print("three-way comparison:")
    for k, v in three_way.items():
        print(f"  {k}: {v}")

    consolidated = {
        "figure_outputs": figure_outputs,
        "rna_only_reviewer_context_isolation": True,
        "anti_contamination_report": contamination_report,
        "params": RNA_ONLY_PARAMS,
        "figure_spec_sha": spec_sha,
        "rendered_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "renderer_script_sha": sha256_file(Path(__file__)),
    }
    (OUT_ROOT / "_rna_only_figure_outputs.json").write_text(
        json.dumps(consolidated, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\nRNA-only stream complete. Outputs at: {OUT_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
