#!/usr/bin/env bash
# C benchmark Lane-2 GPU determinism probe (G-C0 A6 prerequisite).
#
# Per G-C0 contract: before launching the full Lane-2 run on 877k cells, run
# Lane-2 (rapids-singlecell drop-in) twice on a small NC tumor slice with
# random_state=42, compute ARI between the two leiden labelings.
#
# Pass: ARI >= 0.999 (effectively bit-identical)
# Acceptable: ARI >= 0.95 (PRNG stream differs slightly — admissible)
# Fail:   ARI < 0.95 — rsc.tl.leiden does not honor seeds; escalate
#
# Wall clock: ~10 min on RTX 5090 + 5k cells.

set -euo pipefail

REPO_ROOT="/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory"
PROJECT_ROOT="/home/zerlinshen/projects/nc-reproduction"
RESULTS_BASE="${PROJECT_ROOT}/runs"
PROBE_DIR="${RESULTS_BASE}/2026-05-20-C-lane2-determinism-probe"
PYBIN="/home/zerlinshen/conda/envs/sc_gpu_stable/bin/python"

mkdir -p "${PROBE_DIR}"

# Inputs and pins
ZARR="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr"
N_PROBE=5000

cd "${REPO_ROOT}"

cat > "${PROBE_DIR}/probe_runner.py" <<'PY'
"""Lane-2 GPU determinism probe: load 5k cells from tumor cohort, run rsc clustering
twice with random_state=42, compute ARI between leiden labelings."""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
import rapids_singlecell as rsc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ZARR_PATH = Path(os.environ["ZARR_PATH"])
PROBE_DIR = Path(os.environ["PROBE_DIR"])
N_PROBE = int(os.environ.get("N_PROBE", "5000"))
SEED = 42

t0 = time.time()
print(f"[probe] loading zarr from {ZARR_PATH}")
adata_full = ad.read_zarr(str(ZARR_PATH))
print(f"[probe] loaded {adata_full.n_obs} cells x {adata_full.n_vars} genes")

# Subset to tumor cohort by 'condition' field (G-C0 A1)
if "condition" not in adata_full.obs.columns:
    raise SystemExit("ERROR: adata.obs['condition'] missing — cannot subset to tumor")
tumor_mask = adata_full.obs["condition"].astype(str) == "tumor"
print(f"[probe] tumor cells: {int(tumor_mask.sum())}")

# Random 5k subset (seed=42) for reproducibility
rng = np.random.default_rng(SEED)
tumor_idx = np.where(tumor_mask.values)[0]
sel = rng.choice(tumor_idx, size=min(N_PROBE, len(tumor_idx)), replace=False)
adata_sub = adata_full[sel].copy()
print(f"[probe] subset to {adata_sub.n_obs} cells; took {time.time()-t0:.1f}s")

def run_lane2(adata, label):
    """Run rapids-singlecell clustering with random_state=42; return leiden labels."""
    print(f"\n[{label}] starting GPU clustering pass")
    t = time.time()
    # Materialize to in-memory float32 sparse
    from scipy import sparse
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X.astype(np.float32))
    else:
        adata.X = adata.X.tocsr().astype(np.float32)
    # Normalize + log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # HVG via scanpy CPU (deterministic across runs)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)
    adata._inplace_subset_var(adata.var["highly_variable"].values)
    # PCA via rsc (GPU)
    rsc.pp.pca(adata, n_comps=40, random_state=SEED)
    # Neighbors + UMAP + Leiden via rsc
    rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep="X_pca", random_state=SEED)
    rsc.tl.umap(adata, random_state=SEED)
    rsc.tl.leiden(adata, resolution=0.8, random_state=SEED)
    labels = adata.obs["leiden"].astype(str).values.copy()
    print(f"[{label}] done in {time.time()-t:.1f}s; n_clusters={len(set(labels))}")
    return labels

# Two passes with the same random_state — must produce identical-ish labelings
adata_r1 = adata_sub.copy()
labels_r1 = run_lane2(adata_r1, "run1")

adata_r2 = adata_sub.copy()
labels_r2 = run_lane2(adata_r2, "run2")

ari = float(adjusted_rand_score(labels_r1, labels_r2))
nmi = float(normalized_mutual_info_score(labels_r1, labels_r2))

if ari >= 0.999:
    verdict = "PASS_STRICT"
elif ari >= 0.95:
    verdict = "PASS_ACCEPTABLE"
else:
    verdict = "FAIL"

# rapids-singlecell version capture
try:
    rsc_version = rsc.__version__
except Exception:
    rsc_version = "unknown"

report = {
    "probe_type": "lane2_gpu_determinism",
    "n_cells_probed": int(adata_sub.n_obs),
    "n_genes": int(adata_sub.n_vars),
    "random_state": SEED,
    "ari_run1_vs_run2": ari,
    "nmi_run1_vs_run2": nmi,
    "n_clusters_run1": len(set(labels_r1)),
    "n_clusters_run2": len(set(labels_r2)),
    "verdict": verdict,
    "wall_clock_seconds": time.time() - t0,
    "rapids_singlecell_version": rsc_version,
    "g_c0_thresholds": {"strict": 0.999, "acceptable": 0.95},
}
out_path = PROBE_DIR / "lane2_determinism_probe.json"
out_path.write_text(json.dumps(report, indent=2))
print(f"\n[probe] verdict={verdict} ari={ari:.6f} nmi={nmi:.6f}")
print(f"[probe] report written: {out_path}")
PY

env \
  ZARR_PATH="${ZARR}" \
  PROBE_DIR="${PROBE_DIR}" \
  N_PROBE="${N_PROBE}" \
  PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" \
  "${PYBIN}" "${PROBE_DIR}/probe_runner.py" 2>&1 | tee "${PROBE_DIR}/probe_runner.log"

cat "${PROBE_DIR}/lane2_determinism_probe.json"
echo "PROBE_DIR=${PROBE_DIR}"
