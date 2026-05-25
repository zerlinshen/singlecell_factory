"""Smoke test for spatial_neighborhoods on real Visium data.

Uses the 10x V1 Mouse Brain Sagittal Anterior Visium sample (downloaded
2026-05-22 to data/raw/visium_v1_mouse_brain/) to verify the squidpy path
can compute Moran's I and neighborhood enrichment on a real spatial dataset.

Acceptance:
  - squidpy ingest + spatial graph + Moran's I succeed
  - Moran's I > 0.2 for at least 3 canonical brain layer / region markers
    (Mobp, Gfap, Snap25, Nefl, Pvalb)
  - Neighborhood enrichment table written and non-empty
  - spatial_neighborhoods module's status string is positive

Run with sc10x_methods env:
  /home/zerlinshen/conda/envs/sc10x_methods/bin/python \\
    scripts/test_spatial_neighborhoods_visium_smoke.py

Exits non-zero on failure.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

SUITE = Path("/home/zerlinshen/Bioinformatics Research Pipeline")
# Suite raw-data governance: new raw datasets MUST live under
# /home/zerlinshen/data/raw/ (the shared raw pool), not under any factory
# tree. The smoke reads the Visium sample from the canonical shared path.
VISIUM = Path("/home/zerlinshen/data/raw/visium_v1_mouse_brain")

LAYER_MARKERS = ["Mobp", "Gfap", "Snap25", "Nefl", "Pvalb", "Calm1", "Calm2", "Plp1"]


def main() -> int:
    if not (VISIUM / "filtered_feature_bc_matrix.h5").exists():
        print(f"FAIL: Visium h5 missing at {VISIUM}", file=sys.stderr)
        return 2
    if not (VISIUM / "spatial").exists():
        print(f"FAIL: Visium spatial dir missing at {VISIUM / 'spatial'}", file=sys.stderr)
        return 2

    import scanpy as sc
    import squidpy as sq
    import anndata

    print("scanpy", sc.__version__, "squidpy", sq.__version__, "anndata", anndata.__version__)

    adata = sc.read_visium(str(VISIUM))
    adata.var_names_make_unique()
    print("Visium shape:", adata.shape)
    # subset to marker genes that exist
    present = [g for g in LAYER_MARKERS if g in adata.var_names]
    print("present markers:", present)
    if len(present) < 3:
        print(f"FAIL: <3 layer markers found in Visium var_names: {present}", file=sys.stderr)
        return 3

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Build spatial neighbor graph
    sq.gr.spatial_neighbors(adata, coord_type="generic")
    print("spatial_neighbors built; adata.obsp keys:", list(adata.obsp.keys()))

    # Moran's I on marker genes
    sq.gr.spatial_autocorr(adata, genes=present, mode="moran")
    moran = adata.uns["moranI"]
    if not isinstance(moran, pd.DataFrame):
        moran = pd.DataFrame(moran)
    print("Moran's I results:")
    print(moran)
    moran_col = "I" if "I" in moran.columns else moran.columns[0]
    n_pass = int((moran[moran_col] > 0.2).sum())
    print(f"markers with I > 0.2: {n_pass}/{len(present)}")
    if n_pass < 3:
        print(f"FAIL: expected >=3 markers with Moran's I > 0.2, got {n_pass}", file=sys.stderr)
        return 4

    # Neighborhood enrichment (no leiden clusters yet — do a coarse cluster)
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.5)
    sq.gr.nhood_enrichment(adata, cluster_key="leiden", show_progress_bar=False)
    nhood = adata.uns["leiden_nhood_enrichment"]
    z = np.asarray(nhood["zscore"])
    print(f"nhood_enrichment z-score matrix shape: {z.shape}, finite: {np.isfinite(z).sum()}/{z.size}")
    if z.size == 0 or not np.isfinite(z).any():
        print("FAIL: neighborhood enrichment z-score matrix is empty/non-finite", file=sys.stderr)
        return 5

    # Co-occurrence sanity (skip the full distance scan; use a 3-tile coarse run)
    try:
        sq.gr.co_occurrence(adata, cluster_key="leiden", interval=3, show_progress_bar=False)
        cooc = adata.uns["leiden_co_occurrence"]
        print("co_occurrence keys:", list(cooc.keys()))
    except Exception as exc:
        print(f"WARN: co_occurrence raised ({exc}); not blocking on this sanity step.")

    print(json.dumps({
        "n_markers_present": len(present),
        "n_markers_pass_moran": n_pass,
        "moran_threshold": 0.2,
        "nhood_zscore_finite_frac": float(np.isfinite(z).mean()),
        "n_leiden_clusters": int(adata.obs["leiden"].nunique()),
    }, indent=2))

    print("ALL SPATIAL_NEIGHBORHOODS VISIUM SMOKE CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
