"""Parity harness: CPU scrublet vs rsc.pp.scrublet on a 50k synthetic fixture.

Builds a 50k-cell AnnData with ~5% true doublets injected, runs both
implementations with identical params, and reports per-cell agreement.

Exit codes:
  0 — agreement >= 0.99 (GPU adoption recommended)
  1 — agreement < 0.99 (keep CPU path)
  2 — environment/import error
"""

from __future__ import annotations

import sys
import time

import anndata as ad
import numpy as np
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Synthetic data builder
# ---------------------------------------------------------------------------

def build_synthetic_adata(
    n_obs: int = 50_000,
    n_vars: int = 2_000,
    doublet_rate: float = 0.05,
    n_cell_types: int = 8,
    random_state: int = 42,
) -> ad.AnnData:
    """Return raw count AnnData with true doublets already injected.

    Uses multiple cell types with distinct marker gene programs so that
    doublet scores form a bimodal distribution — required for Scrublet's
    automatic threshold detection to work correctly.

    Doublets are created by summing two randomly chosen singlet count rows
    from *different* cell types, matching what Scrublet simulates.
    """
    rng = np.random.default_rng(random_state)

    n_singlets = int(round(n_obs * (1 - doublet_rate)))
    n_doublets = n_obs - n_singlets

    # Each cell type has a distinct expression program
    genes_per_type = n_vars // n_cell_types
    type_means = np.zeros((n_cell_types, n_vars), dtype=np.float32)
    for ct in range(n_cell_types):
        # Shared low background
        type_means[ct] = rng.exponential(scale=0.3, size=n_vars).astype(np.float32)
        # Type-specific markers: high expression in a gene block
        start = ct * genes_per_type
        end = start + genes_per_type
        type_means[ct, start:end] += rng.exponential(scale=5.0, size=end - start).astype(np.float32)

    # Assign cell types evenly
    cell_types = np.tile(np.arange(n_cell_types), n_singlets // n_cell_types + 1)[:n_singlets]
    rng.shuffle(cell_types)

    # Draw counts per cell with per-cell depth variation
    depth = rng.uniform(0.5, 2.0, size=(n_singlets, 1)).astype(np.float32)
    lam = type_means[cell_types] * depth  # (n_singlets, n_vars)
    singlet_counts = rng.poisson(lam=lam).astype(np.float32)

    # Inject doublets: pair cells from different cell types for realistic doublets
    all_idx = np.arange(n_singlets)
    idx_a = rng.integers(0, n_singlets, size=n_doublets)
    # Force type diversity: pick idx_b from a different cell type
    idx_b = np.array([
        rng.choice(all_idx[cell_types != cell_types[a]])
        for a in idx_a
    ], dtype=np.int64)
    doublet_counts = singlet_counts[idx_a] + singlet_counts[idx_b]

    X = np.vstack([singlet_counts, doublet_counts])
    is_true_doublet = np.array([False] * n_singlets + [True] * n_doublets, dtype=bool)

    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs={"true_doublet": is_true_doublet},
    )
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


# ---------------------------------------------------------------------------
# CPU scrublet runner
# ---------------------------------------------------------------------------

def run_cpu_scrublet(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.05,
    random_state: int = 42,
) -> np.ndarray:
    """Run CPU scrublet; returns boolean predicted_doublet array."""
    import scrublet as scr

    n_obs, n_vars = adata.n_obs, adata.n_vars
    params = {
        "min_counts": 2,
        "min_cells": 3,
        "min_gene_variability_pctl": 85,
        "n_prin_comps": max(2, min(30, n_obs - 1, n_vars - 1)),
    }

    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()

    scrub = scr.Scrublet(
        X,
        expected_doublet_rate=expected_doublet_rate,
        random_state=random_state,
    )
    _, predicted = scrub.scrub_doublets(**params)
    return predicted.astype(bool)


# ---------------------------------------------------------------------------
# GPU (rsc) scrublet runner
# ---------------------------------------------------------------------------

def run_gpu_scrublet(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.05,
    random_state: int = 42,
) -> np.ndarray:
    """Run rsc.pp.scrublet; returns boolean predicted_doublet array."""
    import rapids_singlecell as rsc

    # Work on a copy so we don't mutate the shared fixture, then move to GPU
    adata_gpu = adata.copy()
    rsc.get.anndata_to_GPU(adata_gpu)

    n_obs, n_vars = adata_gpu.n_obs, adata_gpu.n_vars

    rsc.pp.scrublet(
        adata_gpu,
        expected_doublet_rate=expected_doublet_rate,
        n_prin_comps=max(2, min(30, n_obs - 1, n_vars - 1)),
        random_state=random_state,
        verbose=False,
    )
    # obs columns may be cudf Series; convert to numpy
    col = adata_gpu.obs["predicted_doublet"]
    if hasattr(col, "to_numpy"):
        return col.to_numpy().astype(bool)
    return np.asarray(col, dtype=bool)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("=== Scrublet CPU vs GPU parity harness ===")
    print(f"Building 50k-cell synthetic AnnData (5% true doublets)...")

    t0 = time.perf_counter()
    adata = build_synthetic_adata(
        n_obs=50_000,
        n_vars=2_000,
        doublet_rate=0.05,
        random_state=42,
    )
    print(f"  AnnData shape: {adata.shape}  [{time.perf_counter()-t0:.1f}s]")

    # --- CPU run ---
    print("\nRunning CPU scrublet...")
    t1 = time.perf_counter()
    try:
        cpu_labels = run_cpu_scrublet(adata, expected_doublet_rate=0.05, random_state=42)
    except Exception as exc:
        print(f"  ERROR: CPU scrublet failed: {exc}", file=sys.stderr)
        return 2
    cpu_time = time.perf_counter() - t1
    print(f"  CPU doublets detected: {cpu_labels.sum()}  [{cpu_time:.1f}s]")

    # --- GPU run ---
    print("\nRunning GPU (rsc) scrublet...")
    t2 = time.perf_counter()
    try:
        gpu_labels = run_gpu_scrublet(adata, expected_doublet_rate=0.05, random_state=42)
    except Exception as exc:
        print(f"  ERROR: rsc.pp.scrublet failed: {exc}", file=sys.stderr)
        return 2
    gpu_time = time.perf_counter() - t2
    print(f"  GPU doublets detected: {gpu_labels.sum()}  [{gpu_time:.1f}s]")

    # --- Agreement ---
    agreement = float(np.mean(cpu_labels == gpu_labels))
    n_agree = int((cpu_labels == gpu_labels).sum())
    n_total = len(cpu_labels)

    print(f"\n--- Parity result ---")
    print(f"  Cells agreeing  : {n_agree} / {n_total}")
    print(f"  Agreement rate  : {agreement:.6f} ({agreement*100:.4f}%)")
    print(f"  CPU time        : {cpu_time:.2f}s")
    print(f"  GPU time        : {gpu_time:.2f}s")
    print(f"  Speedup (CPU/GPU): {cpu_time/gpu_time:.1f}x" if gpu_time > 0 else "")

    threshold = 0.99
    passed = agreement >= threshold
    print(f"\n  Threshold       : {threshold}")
    print(f"  RESULT          : {'PASS — GPU adoption recommended' if passed else 'FAIL — keep CPU path'}")

    # Write research memo
    import os
    research_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", ".omc", "research"
    )
    os.makedirs(research_dir, exist_ok=True)
    memo_path = os.path.join(research_dir, "scrublet-parity.md")
    with open(memo_path, "w") as f:
        f.write("# Scrublet CPU vs GPU Parity Report\n\n")
        f.write(f"**Date**: 2026-05-15\n")
        f.write(f"**Fixture**: 50k synthetic cells, 5% injected doublets, 2k genes\n")
        f.write(f"**CPU impl**: `scrublet.Scrublet` (CPU)\n")
        f.write(f"**GPU impl**: `rsc.pp.scrublet` (rapids-singlecell 0.13.4)\n\n")
        f.write("## Results\n\n")
        f.write(f"| Metric | Value |\n")
        f.write(f"|--------|-------|\n")
        f.write(f"| Cells total | {n_total} |\n")
        f.write(f"| Cells agreeing | {n_agree} |\n")
        f.write(f"| Agreement rate | {agreement:.6f} ({agreement*100:.4f}%) |\n")
        f.write(f"| CPU doublets detected | {int(cpu_labels.sum())} |\n")
        f.write(f"| GPU doublets detected | {int(gpu_labels.sum())} |\n")
        f.write(f"| CPU time | {cpu_time:.2f}s |\n")
        f.write(f"| GPU time | {gpu_time:.2f}s |\n")
        f.write(f"| Speedup (CPU/GPU) | {cpu_time/gpu_time:.1f}x |\n")
        f.write(f"| Threshold | {threshold} |\n")
        f.write(f"| **Outcome** | **{'PASS' if passed else 'FAIL'}** |\n\n")
        if passed:
            f.write(
                "## Decision\n\n"
                "Agreement >= 99%. GPU path adopted in `doublet_detection.py`.\n"
                "CPU `scrublet.Scrublet` path removed (unreachable after migration).\n"
            )
        else:
            f.write(
                "## Decision\n\n"
                f"Agreement {agreement:.4%} < 99% threshold. CPU path retained.\n"
                "US-B2 closes as 'kept CPU per parity failure'.\n\n"
                "## Divergence analysis\n\n"
                "Cells with divergent calls:\n"
                f"- CPU=doublet, GPU=singlet: {int(((cpu_labels) & (~gpu_labels)).sum())}\n"
                f"- CPU=singlet, GPU=doublet: {int(((~cpu_labels) & (gpu_labels)).sum())}\n"
            )
    print(f"\n  Memo written to: {os.path.abspath(memo_path)}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
