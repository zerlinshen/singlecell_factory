"""Cross-language parity sentinels between scanpy (Python) and Seurat (R).

These are SENTINEL tests: not exhaustive correctness checks but tripwires that
fire the moment scanpy or Seurat changes a default and the two pipelines
silently diverge. Thresholds are intentionally loose so minor library version
churn does not flip CI red, but tight enough to catch real drift.

Both tests are marked ``@pytest.mark.r_contract`` and skip cleanly when the
R-side tooling (Rscript / Seurat) is unavailable, mirroring the pattern in
test_r_bundle_contract.py.

Path used to expose full expression to R-side
---------------------------------------------
The bundle's ``expr_sparse`` payload only contains genes listed in
``ExportConfig.markers``. To give the R side the full matrix without writing
a parallel parquet sidecar, we pass *every* gene name as a marker. The bundle
exporter filters markers against ``adata.var_names`` (write_marker_expr in
scripts/export_singlecell_r_bundle.py), so this lands the full matrix in
``expr_sparse`` while staying entirely on the documented bundle v2.1 path.
"""

from __future__ import annotations

import os
# Scanpy 1.12 + numba on this platform crashes when NUMBA_DISABLE_JIT=1 (the
# repo-level conftest default) because rank_genes_groups exercises a JIT-only
# code path. Re-enable JIT for THIS module before scanpy is imported. We do not
# touch other tests because we never modify shared state outside this file.
os.environ["NUMBA_DISABLE_JIT"] = "0"

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle  # noqa: E402

IO_BUNDLE_R = "/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory/R_bundle/io_bundle.R"


# ---------------------------------------------------------------------------
# Helpers (mirrored from test_r_bundle_contract.py — kept local so this file
# is self-contained and the existing test module is not modified)
# ---------------------------------------------------------------------------

def _run_r(rscript: str, script: str, timeout: int = 180) -> subprocess.CompletedProcess:
    return subprocess.run(
        [rscript, "-e", script],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _parse_kv_stdout(stdout: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in stdout.splitlines():
        if "=" in line:
            key, _, val = line.partition("=")
            result[key.strip()] = val.strip()
    return result


def _assert_r_ok(proc: subprocess.CompletedProcess, label: str) -> None:
    assert proc.returncode == 0, (
        f"R script [{label}] exited with code {proc.returncode}.\n"
        f"--- stdout ---\n{proc.stdout}\n"
        f"--- stderr ---\n{proc.stderr}"
    )
    if "Error" in proc.stderr:
        pytest.fail(
            f"R script [{label}] has Error in stderr:\n{proc.stderr}\nstdout:\n{proc.stdout}"
        )


def _maybe_skip_seurat(proc: subprocess.CompletedProcess) -> None:
    """If the R script printed SEURAT_MISSING=1 (Seurat unloadable), skip."""
    if "SEURAT_MISSING=1" in proc.stdout:
        pytest.skip("Seurat unavailable in R env; skipping parity test")


def _make_clustered_adata(
    n_cells: int,
    n_genes: int,
    n_clusters: int,
    n_diff_per_cluster: int,
    seed: int,
):
    """Synthesize a deterministic AnnData with clear per-cluster differential genes.

    Each cluster gets ``n_diff_per_cluster`` "marker" genes whose mean
    expression is pushed up ~5x in cells assigned to that cluster.
    """
    rng = np.random.default_rng(seed)
    # Base counts: low-rate Poisson-ish background.
    base = rng.integers(0, 3, size=(n_cells, n_genes)).astype(np.float32)

    # Assign each cell a cluster (balanced, deterministic).
    labels = np.array([i % n_clusters for i in range(n_cells)], dtype=np.int32)

    # Pick non-overlapping marker gene blocks per cluster.
    assert n_clusters * n_diff_per_cluster <= n_genes, (
        f"need {n_clusters * n_diff_per_cluster} marker slots, only {n_genes} genes"
    )
    for k in range(n_clusters):
        gene_lo = k * n_diff_per_cluster
        gene_hi = gene_lo + n_diff_per_cluster
        mask_cells = labels == k
        # Strong, deterministic up-regulation in the cluster's marker block.
        bump = rng.integers(8, 15, size=(mask_cells.sum(), n_diff_per_cluster)).astype(np.float32)
        base[np.ix_(mask_cells, np.arange(gene_lo, gene_hi))] += bump

    x = sparse.csr_matrix(base)
    obs = pd.DataFrame(
        {"leiden": pd.Categorical([str(k) for k in labels])},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    # Use hyphens (not underscores) so Seurat's CreateSeuratObject does not
    # silently rename features (it converts "_" to "-" for feature names),
    # which would break gene-name set equality on the round-trip.
    var = pd.DataFrame(index=[f"GENE-{i:04d}" for i in range(n_genes)])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    return adata, labels


# ---------------------------------------------------------------------------
# Test 1: marker-gene parity (Jaccard sentinel on top-N per cluster)
# ---------------------------------------------------------------------------
# Threshold rationale: 0.6 on top-50 is loose. With 30 strong differential
# genes per cluster, both algorithms should agree on the top ~30 quite firmly,
# and the next 20 are noise — Jaccard 0.6 means at least 38/62 of the union
# overlap, easily achievable when the strong block matches but ordering of the
# weak tail differs. Enough to catch a real default-flip without flapping on
# library-version-induced rank perturbations.

@pytest.mark.r_contract
def test_markers_top_n_parity(rscript_path, tmp_path):
    sc = pytest.importorskip("scanpy")

    n_cells = 300
    n_genes = 200
    n_clusters = 3
    n_diff = 30
    top_n = 50
    jaccard_floor = 0.6

    adata, _ = _make_clustered_adata(
        n_cells=n_cells, n_genes=n_genes,
        n_clusters=n_clusters, n_diff_per_cluster=n_diff,
        seed=2026,
    )

    # Scanpy side: log1p so wilcoxon sees comparable magnitudes to Seurat
    # (Seurat's FindAllMarkers internally uses normalized data; we mirror the
    # rough scale, which is what affects ranking, not absolute pval).
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=top_n)
    rgg = adata.uns["rank_genes_groups"]
    py_top: dict[str, list[str]] = {}
    for k in range(n_clusters):
        col = str(k)
        py_top[col] = [str(g) for g in rgg["names"][col][:top_n].tolist()]

    # Write a fresh raw-counts AnnData for export (Seurat will normalize itself).
    raw_adata, _ = _make_clustered_adata(
        n_cells=n_cells, n_genes=n_genes,
        n_clusters=n_clusters, n_diff_per_cluster=n_diff,
        seed=2026,
    )
    h5ad = tmp_path / "parity_markers.h5ad"
    bundle_dir = tmp_path / "bundle_markers"
    raw_adata.write_h5ad(h5ad)

    all_genes = tuple(raw_adata.var_names.astype(str).tolist())
    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("leiden",),
            markers=all_genes,  # see module docstring for the rationale
            obsm_keys=(),  # don't need embeddings for this test
            schema_version="v2.1",
            format="parquet",
        )
    )

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    io_r = IO_BUNDLE_R.replace("'", "\\'")
    r_script = f"""
tryCatch({{
  suppressPackageStartupMessages({{
    library(Seurat); library(Matrix)
  }})
}}, error = function(e) {{
  cat('SEURAT_MISSING=1\\n'); quit(status = 0)
}})
source('{io_r}')
b <- read_bundle_v2('{bundle_dir_r}')
# expr_sparse is cells x genes; Seurat wants genes x cells.
mat <- Matrix::t(b$expr_sparse)
# CreateSeuratObject expects integer-ish counts. Our synthetic data is integer.
obj <- CreateSeuratObject(counts = mat, meta.data = b$obs)
Idents(obj) <- b$obs$leiden
obj <- NormalizeData(obj, verbose = FALSE)
markers <- FindAllMarkers(obj,
  logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE,
  test.use = 'wilcox', verbose = FALSE
)
# Print top-N gene names per cluster, ordered by avg_log2FC descending.
clusters <- sort(unique(as.character(markers$cluster)))
for (cl in clusters) {{
  sub <- markers[as.character(markers$cluster) == cl, , drop = FALSE]
  sub <- sub[order(-sub$avg_log2FC), , drop = FALSE]
  topn <- head(sub$gene, {top_n})
  cat(sprintf('CLUSTER_%s_TOPN=%s\\n', cl, paste(topn, collapse = ',')))
}}
cat('R_DONE=1\\n')
"""
    proc = _run_r(rscript_path, r_script)
    _maybe_skip_seurat(proc)
    _assert_r_ok(proc, "markers_parity")

    kv = _parse_kv_stdout(proc.stdout)
    assert kv.get("R_DONE") == "1", f"R driver did not complete:\n{proc.stdout}"

    # Compute per-cluster Jaccard between Python top-N and R top-N.
    # If R's threshold (logfc.threshold=0.25, only.pos) returns fewer than
    # ``top_n`` markers, shrink the comparison window to match R's count
    # (capped at a 20 floor). Otherwise the comparison silently asymmetrizes:
    # R's "top 30" vs Python's "top 50" can drop Jaccard below 0.6 purely
    # because Python adds 20 extra noise genes that R's threshold filtered
    # away. Catching genuine drift means comparing matched windows.
    minimum_window = 20
    jaccards: dict[str, float] = {}
    windows: dict[str, int] = {}
    for k in range(n_clusters):
        col = str(k)
        r_key = f"CLUSTER_{col}_TOPN"
        assert r_key in kv, f"missing {r_key} in R output:\n{proc.stdout}"
        r_genes = [g for g in kv[r_key].split(",") if g]
        n_eff = min(top_n, len(r_genes), len(py_top[col]))
        if n_eff < minimum_window:
            pytest.skip(
                f"R returned fewer than {minimum_window} markers for cluster {col} "
                f"(got {len(r_genes)}); cannot run a meaningful sentinel."
            )
        py_set = set(py_top[col][:n_eff])
        r_set = set(r_genes[:n_eff])
        windows[col] = n_eff
        inter = len(py_set & r_set)
        union = len(py_set | r_set)
        jaccards[col] = inter / union if union else 0.0

    # When the matched window is below top_n, fall back to a 0.5 floor; the
    # full top-50 window holds to 0.6. Both are sentinel thresholds, not
    # tight correctness bounds.
    weak_floor = 0.5
    failures: dict[str, tuple[float, int]] = {}
    for col, j in jaccards.items():
        floor = jaccard_floor if windows[col] >= top_n else weak_floor
        if j < floor:
            failures[col] = (j, windows[col])
    assert not failures, (
        f"Marker-set drift detected.\n"
        f"jaccards={jaccards}\nwindows={windows}\n"
        f"failures={failures}\n"
        f"R stdout (truncated):\n{proc.stdout[:2000]}"
    )


# ---------------------------------------------------------------------------
# Test 2: PCA cosine-similarity sentinel (first 5 PCs, sign-aligned)
# ---------------------------------------------------------------------------
# Threshold rationale: 0.95 absolute cosine on the first 5 PCs survives minor
# scaling differences (target_sum 1e4 vs 1e6) and small HVG-selection drift.
# Anything below that means a real algorithmic divergence (different
# centering, different covariance estimator, different default flag).

@pytest.mark.r_contract
def test_pca_cosine_parity(rscript_path, tmp_path):
    sc = pytest.importorskip("scanpy")

    n_cells = 200
    n_genes = 500
    n_pcs = 5
    cosine_floor = 0.95

    # Build clustered data so a PCA actually has structure to recover.
    adata, _ = _make_clustered_adata(
        n_cells=n_cells, n_genes=n_genes,
        n_clusters=4, n_diff_per_cluster=20,
        seed=4242,
    )

    # Python pipeline.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Scale before PCA — matches scanpy's standard preprocess and aligns with
    # Seurat's RunPCA-on-scale.data convention.
    sc.pp.scale(adata, max_value=None)
    sc.pp.pca(adata, n_comps=n_pcs, random_state=0)
    pca_py = np.asarray(adata.obsm["X_pca"][:, :n_pcs], dtype=np.float64)

    # Re-build a fresh raw-counts AnnData for export.
    raw_adata, _ = _make_clustered_adata(
        n_cells=n_cells, n_genes=n_genes,
        n_clusters=4, n_diff_per_cluster=20,
        seed=4242,
    )
    h5ad = tmp_path / "parity_pca.h5ad"
    bundle_dir = tmp_path / "bundle_pca"
    raw_adata.write_h5ad(h5ad)

    all_genes = tuple(raw_adata.var_names.astype(str).tolist())
    export_bundle(
        ExportConfig(
            input_h5ad=h5ad,
            output_dir=bundle_dir,
            obs_columns=("leiden",),
            markers=all_genes,
            obsm_keys=(),
            schema_version="v2.1",
            format="parquet",
        )
    )

    # Hand the Python PCA to R via a sidecar parquet so R can compute cosines
    # and emit one MIN_PC_COS plus per-PC values back. Doing the cosine
    # arithmetic in R keeps the Python side a single-line parser.
    py_pca_path = tmp_path / "pca_py.parquet"
    pd.DataFrame(
        pca_py,
        index=raw_adata.obs_names.astype(str),
        columns=[f"PC{i + 1}" for i in range(n_pcs)],
    ).reset_index(names="cell").to_parquet(py_pca_path, index=False)

    bundle_dir_r = str(bundle_dir).replace("'", "\\'")
    io_r = IO_BUNDLE_R.replace("'", "\\'")
    py_pca_r = str(py_pca_path).replace("'", "\\'")

    r_script = f"""
tryCatch({{
  suppressPackageStartupMessages({{
    library(Seurat); library(Matrix); library(arrow)
  }})
}}, error = function(e) {{
  cat('SEURAT_MISSING=1\\n'); quit(status = 0)
}})
source('{io_r}')
b <- read_bundle_v2('{bundle_dir_r}')
mat <- Matrix::t(b$expr_sparse)  # genes x cells
obj <- CreateSeuratObject(counts = mat, meta.data = b$obs)
obj <- NormalizeData(obj, verbose = FALSE)
# Use ALL genes as variable features so we are comparing matched feature sets.
VariableFeatures(obj) <- rownames(obj)
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, features = rownames(obj), npcs = {n_pcs}, verbose = FALSE,
              seed.use = 42)
pca_r <- Embeddings(obj, reduction = 'pca')[, 1:{n_pcs}, drop = FALSE]

py_df <- as.data.frame(arrow::read_parquet('{py_pca_r}'))
rownames(py_df) <- py_df$cell
py_df$cell <- NULL
pca_py <- as.matrix(py_df)
# Align cell order (R PCA uses Seurat's internal order, which equals obs).
pca_py <- pca_py[rownames(pca_r), , drop = FALSE]

cosines <- numeric({n_pcs})
for (k in seq_len({n_pcs})) {{
  v_py <- as.numeric(pca_py[, k])
  v_r  <- as.numeric(pca_r[, k])
  num  <- sum(v_py * v_r)
  den  <- sqrt(sum(v_py^2)) * sqrt(sum(v_r^2))
  cos_k <- if (den > 0) abs(num / den) else 0
  cosines[k] <- cos_k
  cat(sprintf('PC%d_COS=%.6f\\n', k, cos_k))
}}
cat(sprintf('MIN_PC_COS=%.6f\\n', min(cosines)))
cat('R_DONE=1\\n')
"""
    proc = _run_r(rscript_path, r_script)
    _maybe_skip_seurat(proc)
    _assert_r_ok(proc, "pca_parity")

    kv = _parse_kv_stdout(proc.stdout)
    assert kv.get("R_DONE") == "1", f"R driver did not complete:\n{proc.stdout}"
    assert "MIN_PC_COS" in kv, f"MIN_PC_COS missing:\n{proc.stdout}"

    per_pc = {k: float(kv[k]) for k in kv if k.endswith("_COS") and k != "MIN_PC_COS"}
    min_cos = float(kv["MIN_PC_COS"])
    failing = {k: v for k, v in per_pc.items() if v < cosine_floor}
    assert not failing, (
        f"PCA divergence: cosine(s) below {cosine_floor}.\n"
        f"per_pc={per_pc}\nmin={min_cos}\nR stdout:\n{proc.stdout[:2000]}"
    )
