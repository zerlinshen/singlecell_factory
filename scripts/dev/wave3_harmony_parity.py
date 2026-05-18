"""US-B3 Harmony Parity Harness
================================
Compares harmonypy.run_harmony (CPU, v0.0.10) against
rsc.pp.harmony_integrate (GPU) on a synthetic 50k / 4-batch fixture.

Decision gate: ARI >= 0.95 -> migration approved.

Usage:
    /home/zerlinshen/conda/envs/sc_gpu_stable/bin/python \
        scripts/dev/wave3_harmony_parity.py
"""

from __future__ import annotations

import sys
import time
import warnings

import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
RNG = np.random.default_rng(SEED)

N_CELLS = 50_000
N_GENES = 2_000
N_BATCHES = 4
N_CLUSTERS = 8
N_PCS = 30
ARI_THRESHOLD = 0.95


# ---------------------------------------------------------------------------
# 1. Build synthetic AnnData
# ---------------------------------------------------------------------------
def build_fixture() -> ad.AnnData:
    """50k cells, 4 batches, 8 true clusters defined by Gaussian mixture in PC space."""
    print("Building 50k synthetic fixture ...")
    cells_per_cluster = N_CELLS // N_CLUSTERS

    pcs_list = []
    batch_labels = []
    cluster_labels = []

    # Cluster centres in a 30-dim PC space
    centres = RNG.standard_normal((N_CLUSTERS, N_PCS)) * 5.0

    for ci in range(N_CLUSTERS):
        n = cells_per_cluster if ci < N_CLUSTERS - 1 else N_CELLS - ci * cells_per_cluster
        noise = RNG.standard_normal((n, N_PCS)) * 0.8
        pcs_list.append(centres[ci] + noise)
        cluster_labels.extend([f"c{ci}"] * n)

    X_pca_true = np.vstack(pcs_list).astype(np.float32)

    # Assign batches cyclically (each batch spans all clusters => non-trivial correction)
    batch_arr = np.array([f"batch{i % N_BATCHES}" for i in range(N_CELLS)])

    # Build a minimal count matrix (log1p-normalised synthetic)
    # Use a sparse-friendly random int matrix so sc.pp.pca can operate
    counts = RNG.poisson(lam=1.5, size=(N_CELLS, N_GENES)).astype(np.float32)

    obs = pd.DataFrame(
        {"batch": batch_arr, "true_cluster": cluster_labels},
        index=[f"cell{i}" for i in range(N_CELLS)],
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_GENES)])
    adata = ad.AnnData(X=counts, obs=obs, var=var)

    # Store the ground-truth PCA (with mild batch offsets baked in for realism)
    batch_offsets = RNG.standard_normal((N_BATCHES, N_PCS)) * 2.0
    batch_idx = np.array([int(b[-1]) for b in batch_arr])
    X_pca_batched = X_pca_true + batch_offsets[batch_idx]
    adata.obsm["X_pca"] = X_pca_batched.astype(np.float32)

    print(f"  Fixture ready: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{N_BATCHES} batches, {N_CLUSTERS} clusters, {N_PCS} PCs")
    return adata


# ---------------------------------------------------------------------------
# 2. ARM A — harmonypy CPU
# ---------------------------------------------------------------------------
def run_arm_cpu(adata: ad.AnnData) -> tuple[np.ndarray, float]:
    import harmonypy

    x = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    print("Arm A (harmonypy CPU): running ...")
    t0 = time.perf_counter()
    ho = harmonypy.run_harmony(x, adata.obs, "batch")
    wall_cpu = time.perf_counter() - t0

    z = np.asarray(ho.Z_corr)
    corrected = z.T if z.T.shape == x.shape else z
    corrected = corrected.astype(np.float32)
    print(f"  Done in {wall_cpu:.2f}s — output shape {corrected.shape}")
    return corrected, wall_cpu


# ---------------------------------------------------------------------------
# 3. ARM B — rsc.pp.harmony_integrate (GPU)
# ---------------------------------------------------------------------------
def run_arm_gpu(adata: ad.AnnData) -> tuple[np.ndarray, float]:
    import rapids_singlecell as rsc

    # rsc.pp.harmony_integrate writes to adata.obsm["X_pca_harmony"] in place
    print("Arm B (rsc.pp.harmony_integrate GPU): running ...")
    t0 = time.perf_counter()
    rsc.pp.harmony_integrate(adata, key="batch")
    wall_gpu = time.perf_counter() - t0

    corrected = np.asarray(adata.obsm["X_pca_harmony"], dtype=np.float32)
    print(f"  Done in {wall_gpu:.2f}s — output shape {corrected.shape}")
    return corrected, wall_gpu


# ---------------------------------------------------------------------------
# 4. Cluster via neighbors + leiden
# ---------------------------------------------------------------------------
def cluster_embedding(adata_base: ad.AnnData, embedding: np.ndarray, label: str) -> np.ndarray:
    """Returns leiden labels (str array) for the given corrected embedding."""
    import anndata as ad

    tmp = ad.AnnData(
        X=adata_base.X,
        obs=adata_base.obs.copy(),
        var=adata_base.var.copy(),
    )
    tmp.obsm["X_corrected"] = embedding

    sc.pp.neighbors(tmp, use_rep="X_corrected", n_neighbors=15, random_state=SEED)
    sc.tl.leiden(
        tmp,
        resolution=0.5,
        flavor="igraph",
        directed=False,
        random_state=SEED,
        n_iterations=2,
    )
    labels = tmp.obs["leiden"].values
    n_clusters = len(np.unique(labels))
    print(f"  [{label}] leiden -> {n_clusters} clusters")
    return labels


# ---------------------------------------------------------------------------
# 5. Metrics
# ---------------------------------------------------------------------------
def frobenius_ratio(A: np.ndarray, B: np.ndarray) -> float:
    diff = np.linalg.norm(A.astype(np.float64) - B.astype(np.float64), "fro")
    nrm_a = np.linalg.norm(A.astype(np.float64), "fro")
    return float(diff / nrm_a) if nrm_a > 0 else float("inf")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 60)
    print("US-B3 Harmony Parity Harness")
    print(f"ARI threshold: {ARI_THRESHOLD}")
    print("=" * 60)

    adata = build_fixture()

    # ARM A (CPU)
    corrected_cpu, wall_cpu = run_arm_cpu(adata)

    # ARM B (GPU) — rsc writes into adata.obsm["X_pca_harmony"]
    # Pass a copy so the CPU result in obsm["X_pca"] is not modified
    adata_gpu = adata.copy()
    corrected_gpu, wall_gpu = run_arm_gpu(adata_gpu)

    # Cluster both arms
    print("Clustering arm A ...")
    labels_cpu = cluster_embedding(adata, corrected_cpu, "cpu")
    print("Clustering arm B ...")
    labels_gpu = cluster_embedding(adata, corrected_gpu, "gpu")

    # Metrics
    ari = adjusted_rand_score(labels_cpu, labels_gpu)
    frob = frobenius_ratio(corrected_cpu, corrected_gpu)

    speedup = wall_cpu / wall_gpu if wall_gpu > 0 else float("inf")

    print()
    print("=" * 60)
    print("PARITY RESULTS")
    print("=" * 60)
    print(f"  ARI (cpu vs gpu labels)  : {ari:.6f}")
    print(f"  Frobenius ratio ||A-B||/||A||: {frob:.6f}")
    print(f"  Wall time CPU            : {wall_cpu:.2f}s")
    print(f"  Wall time GPU            : {wall_gpu:.2f}s")
    print(f"  GPU speedup              : {speedup:.1f}x")
    print()

    passed = ari >= ARI_THRESHOLD
    if passed:
        print(f"  DECISION: PASS (ARI {ari:.4f} >= {ARI_THRESHOLD}) — migration APPROVED")
    else:
        print(f"  DECISION: FAIL (ARI {ari:.4f} < {ARI_THRESHOLD}) — migration DEFERRED")
    print("=" * 60)

    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
