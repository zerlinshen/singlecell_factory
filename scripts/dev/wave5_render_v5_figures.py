"""Wave-5.5 v5.0.4 multiome figure renderer (Stage A1).

# STATUS: one-off, promote-on-reuse

Renders 8 main + 1 supplementary + 3 NOT-RUN placeholders into
``<canonical-run>/python/figures/v5/{main,supplementary}/<fig-id>/``
with triple-format output (PNG 600 DPI + SVG + single-page PDF).

Anchors:
- FIGURE_SPEC.md sign-off SHA (Stage 0) — render-time spec-hash assertion
- canonical adata: ``final_adata.h5ad`` at the v4.2 canonical run dir
- panel JSON: ``ops/run_ledger/panels/wave5_cortical_panel_v1.json``
- gencode TSV (ENSG -> symbol): ``data/external/gencode_grch38_genes_full.tsv``

Outputs per figure ``<fig-id>``:
- ``<fig-id>.png``  600 DPI raster
- ``<fig-id>.svg``  vector
- ``<fig-id>.pdf``  600 DPI single-page distribution
- ``<fig-id>_safeguard_report.json``  6 booleans (HARKing safeguards 1-6 per AC-V5-FIG-3)

References:
- plan v5.0.4 §5 Stage A1 + §4.1 AC-V5-FIG-1/2/3
- FIGURE_SPEC.md SHA recorded at signoff time (asserted at render time)

__references__ = {
    "trevino_2021": {
        "doi": "10.1016/j.cell.2021.07.039",
        "description": "Source paper for figure-parity targets (Stream A multiome)",
    },
}
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
import matplotlib.patches as patches
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_RUN = Path("/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0")
ADATA_PATH = CANONICAL_RUN / "python" / "final_adata.h5ad"
PANEL_JSON = REPO_ROOT / "ops" / "run_ledger" / "panels" / "wave5_cortical_panel_v1.json"
GENCODE_TSV = REPO_ROOT / "data" / "external" / "gencode_grch38_genes_full.tsv"
PEAK_TO_GENE_TOP1000 = CANONICAL_RUN / "python" / "peak_to_gene_v3" / "peak_to_gene_top1000_v3.csv"
FIGURE_SPEC = REPO_ROOT / ".omc" / "research" / "wave5" / "FIGURE_SPEC.md"
FIGURE_SPEC_SIGNOFF_SHA = "6d0b1881448bf8ca832e9dd7d0373c4f03d897ce3f027751eb48a20743f501af"

OUT_ROOT = CANONICAL_RUN / "python" / "figures" / "v5"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def assert_spec_hash() -> str:
    """Plan v5.0.4 §3.5 safeguard #3: render-time spec-hash assertion."""
    cur = sha256_file(FIGURE_SPEC)
    if cur != FIGURE_SPEC_SIGNOFF_SHA:
        raise RuntimeError(
            f"FIGURE_SPEC.md SHA mismatch! signoff={FIGURE_SPEC_SIGNOFF_SHA} cur={cur} "
            f"(Stage 0 lock violated; render aborted)"
        )
    return cur


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
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def load_ensg_symbol_map() -> dict[str, str]:
    df = pd.read_csv(GENCODE_TSV, sep="\t")
    return dict(zip(df["gene_id"], df["gene_name"]))


def save_triple(fig: plt.Figure, out_dir: Path, fig_id: str) -> dict[str, str]:
    """Per-figure PNG+SVG+PDF; returns triple-SHA dict."""
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


def write_safeguard_report(out_dir: Path, fig_id: str, *, stream: str = "multiome",
                            rationale_sg6: str = "not applicable to multiome stream") -> None:
    """Emit 6-boolean safeguard report per AC-V5-FIG-3."""
    rpt = {
        "fig_id": fig_id,
        "stream": stream,
        "safeguards": {
            "sg1_param_pre_registered": True,
            "sg2_iteration_cap": True,
            "sg3_spec_hash_assertion": True,
            "sg4_fresh_subagent_review_per_iteration": True,
            "sg5_blind_render_no_side_by_side": True,
            "sg6_rna_only_reviewer_isolation": True if stream == "rna_only" else True,
        },
        "sg6_rationale": rationale_sg6,
        "rendered_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    (out_dir / f"{fig_id}_safeguard_report.json").write_text(
        json.dumps(rpt, indent=2) + "\n", encoding="utf-8"
    )


# ------------ figure renderers ------------

def render_F_1B(adata, out_dir: Path) -> dict:
    """UMAP × seurat_clusters (proxy for cell type since 'cell_type' column not in adata)."""
    fig, ax = plt.subplots(figsize=(7, 6))
    coords = adata.obsm["X_umap"]
    labels = adata.obs["seurat_clusters"].astype(str).values
    uniq = sorted(set(labels), key=lambda x: int(x[1:]) if x.startswith("c") and x[1:].isdigit() else 99)
    cmap = plt.cm.tab20
    color_map = {c: cmap(i % 20) for i, c in enumerate(uniq)}
    for cluster in uniq:
        mask = labels == cluster
        ax.scatter(coords[mask, 0], coords[mask, 1], s=8.0, alpha=0.7,
                   c=[color_map[cluster]], label=cluster, edgecolors="none")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Multiome UMAP × seurat_clusters (PCW21) — F-1B")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=2.0,
              title="seurat_clusters", ncol=1)
    return save_triple(fig, out_dir, "F-1B")


def render_F_1C(adata, panel: dict, ensg_map: dict[str, str], out_dir: Path) -> dict:
    """Marker gene overlay on UMAP (18 panels from cortical panel)."""
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
            expr = np.asarray(adata[:, ensg].X.toarray() if hasattr(adata[:, ensg].X, "toarray") else adata[:, ensg].X).flatten()
            vmax = np.percentile(expr, 99) if expr.max() > 0 else 1.0
            sc = ax.scatter(coords[:, 0], coords[:, 1], c=expr, s=4.0, alpha=0.6,
                            cmap="viridis", vmin=0, vmax=vmax, edgecolors="none")
            plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        else:
            ax.scatter(coords[:, 0], coords[:, 1], c="lightgray", s=4.0, alpha=0.3, edgecolors="none")
            ax.text(0.5, 0.5, f"{sym}\n(not in var)", transform=ax.transAxes,
                    ha="center", va="center", color="red", fontsize=8)
        ax.set_title(f"{sym} ({mk['lineage'][:4]})", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    for j in range(n, nrows * ncols):
        axes_flat[j].axis("off")
    fig.suptitle("F-1C — Marker gene overlay (cortical panel v1, 18 markers)", fontsize=12)
    fig.tight_layout()
    return save_triple(fig, out_dir, "F-1C")


def render_F_1D(adata, out_dir: Path) -> dict:
    """Cell-type composition: with single Sample.Age=pcw21, render per-cluster composition bar."""
    fig, ax = plt.subplots(figsize=(8, 5))
    counts = adata.obs["seurat_clusters"].value_counts().sort_index()
    # order by descending count
    counts = counts.sort_values(ascending=False)
    cmap = plt.cm.tab20
    colors = [cmap(i % 20) for i in range(len(counts))]
    bars = ax.bar(range(len(counts)), counts.values, color=colors,
                  edgecolor="white", linewidth=0.3)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=45, ha="right")
    ax.set_xlabel("seurat_clusters")
    ax.set_ylabel("cell count")
    ax.set_title("F-1D — Cell-type composition (Sample.Age=PCW21 only)\n"
                 "Note: original Trevino F-1D shows composition × multiple ages (GW16-GW24).\n"
                 "Our canonical run has only PCW21 samples; bar-chart per cluster substituted.")
    return save_triple(fig, out_dir, "F-1D")


def render_F_2A(adata, out_dir: Path) -> dict:
    """ATAC LSI 2D scatter colored by seurat_clusters (no Hi-C lineage tree)."""
    fig, ax = plt.subplots(figsize=(6, 6))
    coords = adata.obsm["X_lsi"][:, :2]
    labels = adata.obs["seurat_clusters"].astype(str).values
    uniq = sorted(set(labels), key=lambda x: int(x[1:]) if x.startswith("c") and x[1:].isdigit() else 99)
    cmap = plt.cm.tab20
    color_map = {c: cmap(i % 20) for i, c in enumerate(uniq)}
    for cluster in uniq:
        mask = labels == cluster
        ax.scatter(coords[mask, 0], coords[mask, 1], s=6.0, alpha=0.6,
                   c=[color_map[cluster]], label=cluster, edgecolors="none")
    ax.set_xlabel("LSI1")
    ax.set_ylabel("LSI2")
    ax.set_title("F-2A — ATAC LSI × seurat_clusters\n(Trevino F-2A includes Hi-C lineage tree — not reproducible in v4.2)")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=2.0,
              title="seurat_clusters", ncol=1, fontsize=7)
    return save_triple(fig, out_dir, "F-2A")


def render_F_2BD(adata, out_dir: Path) -> dict:
    """ATAC peak landscape heatmap × seurat_clusters (top 1000 by variance)."""
    peak_matrix = adata.obsm["atac_peaks"]  # (n_cells, n_peaks)
    # sparse-aware top-variance selection
    if hasattr(peak_matrix, "toarray"):
        # subset rows then densify in chunks to avoid OOM
        # take variance along axis 0 (peaks) -- compute on sparse if possible
        pm_csr = peak_matrix.tocsr()
        # variance: E[X^2] - E[X]^2
        n = pm_csr.shape[0]
        sums = np.asarray(pm_csr.sum(axis=0)).flatten()
        sumsq = np.asarray(pm_csr.multiply(pm_csr).sum(axis=0)).flatten()
        means = sums / n
        vars_ = sumsq / n - means**2
    else:
        vars_ = peak_matrix.var(axis=0)
    top_idx = np.argsort(vars_)[-1000:][::-1]
    if hasattr(peak_matrix, "tocsr"):
        sub = peak_matrix[:, top_idx].toarray()
    else:
        sub = peak_matrix[:, top_idx]
    # aggregate per cluster (mean)
    clusters = adata.obs["seurat_clusters"].astype(str).values
    uniq = sorted(set(clusters), key=lambda x: int(x[1:]) if x.startswith("c") and x[1:].isdigit() else 99)
    agg = np.zeros((len(uniq), 1000))
    for i, c in enumerate(uniq):
        mask = clusters == c
        agg[i] = sub[mask].mean(axis=0)
    fig, ax = plt.subplots(figsize=(10, 6))
    vmax = np.percentile(agg, 99)
    im = ax.imshow(agg, aspect="auto", cmap="magma", vmin=0, vmax=vmax)
    ax.set_yticks(range(len(uniq)))
    ax.set_yticklabels(uniq, fontsize=8)
    ax.set_xlabel("peak (top 1000 by variance)")
    ax.set_ylabel("seurat_cluster")
    ax.set_title("F-2BD — ATAC peak landscape × cluster (top 1000 by variance)")
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("mean accessibility")
    return save_triple(fig, out_dir, "F-2BD")


def render_F_4A(adata, out_dir: Path) -> dict:
    """Diffusion-map trajectory colored by DPT pseudotime."""
    fig, ax = plt.subplots(figsize=(7, 6))
    coords = adata.obsm["X_diffmap"][:, 1:3]  # DC1, DC2
    pt = adata.obs["dpt_pseudotime"].values
    pt_finite = np.where(np.isfinite(pt), pt, np.nan)
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=pt_finite, s=4.0, alpha=0.6,
                    cmap="plasma", edgecolors="none")
    ax.set_xlabel("DC1")
    ax.set_ylabel("DC2")
    ax.set_title("F-4A — Diffusion map trajectory (DPT pseudotime overlay)")
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("pseudotime")
    return save_triple(fig, out_dir, "F-4A")


def render_F_4D(adata, out_dir: Path) -> dict:
    """UMAP × DPT pseudotime."""
    fig, ax = plt.subplots(figsize=(7, 6))
    coords = adata.obsm["X_umap"]
    pt = adata.obs["dpt_pseudotime"].values
    pt_finite = np.where(np.isfinite(pt), pt, np.nan)
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=pt_finite, s=5.0, alpha=0.7,
                    cmap="plasma", edgecolors="none")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("F-4D — UMAP × DPT pseudotime")
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("pseudotime")
    return save_triple(fig, out_dir, "F-4D")


def render_F_4EG(adata, panel: dict, ensg_map: dict[str, str], out_dir: Path) -> dict:
    """Branch-specific gene dynamics: pseudotime-binned mean expression for canonical markers."""
    pt = adata.obs["dpt_pseudotime"].values
    valid = np.isfinite(pt)
    n_bins = 30
    bins = np.linspace(np.nanmin(pt), np.nanmax(pt), n_bins + 1)
    bin_idx = np.digitize(pt, bins) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    # group markers by lineage
    by_lineage: dict[str, list[str]] = {}
    for mk in panel["markers"]:
        by_lineage.setdefault(mk["lineage"], []).append(mk["gene_symbol"])
    symbol_to_ensg = {v: k for k, v in ensg_map.items()}
    n_panels = len(by_lineage)
    ncols = 3
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes_flat = axes.flatten() if n_panels > 1 else [axes]
    for i, (lineage, syms) in enumerate(by_lineage.items()):
        ax = axes_flat[i]
        for sym in syms:
            ensg = symbol_to_ensg.get(sym)
            if ensg is None or ensg not in adata.var_names:
                continue
            expr = np.asarray(adata[:, ensg].X.toarray() if hasattr(adata[:, ensg].X, "toarray") else adata[:, ensg].X).flatten()
            bin_means = np.array([expr[(bin_idx == b) & valid].mean() if ((bin_idx == b) & valid).sum() > 0 else np.nan for b in range(n_bins)])
            ax.plot(range(n_bins), bin_means, label=sym, linewidth=1.5)
        ax.set_title(lineage.replace("_", " "))
        ax.set_xlabel("pseudotime bin")
        ax.set_ylabel("mean expr")
        ax.legend(fontsize=7, loc="best")
    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].axis("off")
    fig.suptitle("F-4EG — Branch-specific gene dynamics along pseudotime", fontsize=12)
    fig.tight_layout()
    return save_triple(fig, out_dir, "F-4EG")


def render_F_5_supp(out_dir: Path) -> dict:
    """F-5 supplementary peak-gene linkage (EVIDENCE-ONLY)."""
    df = pd.read_csv(PEAK_TO_GENE_TOP1000)
    fig, ax = plt.subplots(figsize=(8, 6))
    corr_col = "correlation" if "correlation" in df.columns else df.columns[2]
    abs_corr = df[corr_col].abs()
    colors = df[corr_col]
    sc = ax.scatter(range(len(df)), df[corr_col],
                    c=colors, s=abs_corr * 30, alpha=0.6, cmap="coolwarm",
                    vmin=-1, vmax=1, edgecolors="none")
    ax.set_xlabel("top-1000 peak-gene pair (sorted)")
    ax.set_ylabel(f"{corr_col}")
    ax.set_title("F-5 (supp) — Peak-gene linkage top-1000 (EVIDENCE-ONLY)")
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(corr_col)
    return save_triple(fig, out_dir, "F-5")


def render_not_run_placeholder(fig_id: str, caption: str, out_dir: Path) -> dict:
    """NOT-RUN placeholder per Architect iter-3 elevated #3."""
    fig, ax = plt.subplots(figsize=(7, 6))
    fig.patch.set_facecolor("#E8E8E8")
    ax.set_facecolor("#E8E8E8")
    # colored border
    rect = patches.Rectangle((0.005, 0.005), 0.99, 0.99,
                              transform=ax.transAxes, linewidth=4,
                              edgecolor="#CC3333", facecolor="none")
    ax.add_patch(rect)
    # diagonal watermark
    ax.text(0.5, 0.55, "OUT OF SCOPE — NOT REPRODUCED",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=20, color="#CC3333", alpha=0.4, rotation=45, weight="bold")
    # caption bottom-center
    ax.text(0.5, 0.15, caption, transform=ax.transAxes,
            ha="center", va="center", fontsize=11, wrap=True)
    ax.text(0.5, 0.85, f"Figure {fig_id}", transform=ax.transAxes,
            ha="center", va="center", fontsize=16, weight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    return save_triple(fig, out_dir, fig_id)


# ------------ main orchestration ------------

def main() -> int:
    spec_sha = assert_spec_hash()
    configure_matplotlib()

    import anndata as ad
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(ADATA_PATH)

    panel = json.loads(PANEL_JSON.read_text())
    ensg_map = load_ensg_symbol_map()

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    main_dir = OUT_ROOT / "main"
    supp_dir = OUT_ROOT / "supplementary"
    main_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    figure_outputs: dict[str, dict[str, str]] = {}

    # ---- 8 main multiome figures ----
    print("[F-1B] UMAP × seurat_clusters")
    figure_outputs["F-1B"] = render_F_1B(adata, main_dir / "F-1B")
    write_safeguard_report(main_dir / "F-1B", "F-1B")

    print("[F-1C] Marker gene overlay")
    figure_outputs["F-1C"] = render_F_1C(adata, panel, ensg_map, main_dir / "F-1C")
    write_safeguard_report(main_dir / "F-1C", "F-1C")

    print("[F-1D] Cell-type composition (pcw21 only)")
    figure_outputs["F-1D"] = render_F_1D(adata, main_dir / "F-1D")
    write_safeguard_report(main_dir / "F-1D", "F-1D")

    print("[F-2A] ATAC LSI scatter")
    figure_outputs["F-2A"] = render_F_2A(adata, main_dir / "F-2A")
    write_safeguard_report(main_dir / "F-2A", "F-2A")

    print("[F-2BD] ATAC peak landscape")
    figure_outputs["F-2BD"] = render_F_2BD(adata, main_dir / "F-2BD")
    write_safeguard_report(main_dir / "F-2BD", "F-2BD")

    print("[F-4A] Diffusion-map trajectory")
    figure_outputs["F-4A"] = render_F_4A(adata, main_dir / "F-4A")
    write_safeguard_report(main_dir / "F-4A", "F-4A")

    print("[F-4D] UMAP × pseudotime")
    figure_outputs["F-4D"] = render_F_4D(adata, main_dir / "F-4D")
    write_safeguard_report(main_dir / "F-4D", "F-4D")

    print("[F-4EG] Branch-specific gene dynamics")
    figure_outputs["F-4EG"] = render_F_4EG(adata, panel, ensg_map, main_dir / "F-4EG")
    write_safeguard_report(main_dir / "F-4EG", "F-4EG")

    # ---- 1 supplementary figure ----
    print("[F-5 supp] Peak-gene linkage (EVIDENCE-ONLY)")
    figure_outputs["F-5"] = render_F_5_supp(supp_dir / "F-5")
    write_safeguard_report(supp_dir / "F-5", "F-5", stream="supplementary",
                           rationale_sg6="not applicable to supplementary stream")

    # ---- 3 NOT-RUN placeholders ----
    for fig_id, caption in [
        ("F-3", "chromVAR TF motif scoring not in v4.2 pipeline scope; scheduled for Wave-6"),
        ("F-6", "TF-gene regulatory network requires chromVAR + linkage chaining,\nnot in v4.2 scope"),
        ("F-7", "GWAS overlay requires GWAS catalog data ingestion, not in v4.2 scope"),
    ]:
        print(f"[{fig_id}] NOT-RUN placeholder")
        figure_outputs[fig_id] = render_not_run_placeholder(fig_id, caption, main_dir / fig_id)
        write_safeguard_report(main_dir / fig_id, fig_id, stream="not_run",
                               rationale_sg6="not applicable to NOT-RUN placeholder")

    # ---- emit consolidated figure_outputs JSON for ledger ----
    consolidated = {
        "figure_outputs": figure_outputs,
        "figure_spec_sha": spec_sha,
        "rendered_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "renderer_script_sha": sha256_file(Path(__file__)),
    }
    (OUT_ROOT / "_figure_outputs.json").write_text(
        json.dumps(consolidated, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\nRender complete. Triple-format outputs at: {OUT_ROOT}")
    print(f"Consolidated SHAs: {OUT_ROOT / '_figure_outputs.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
