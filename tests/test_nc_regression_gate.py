"""US-009 / US-010 — verify nc_regression_gate sub-gates on synthetic fixtures.

Builds two synthetic run dirs with the canonical artifact layout
(annotation/cell_type_annotation.csv, differential_expression/marker_genes.csv,
final_adata.h5ad) and checks that nc_regression_compare.py returns the
right exit code per spec for both pass and fail paths.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


FACTORY_ROOT = Path(__file__).resolve().parent.parent
COMPARE_PY = FACTORY_ROOT / "scripts" / "ci" / "nc_regression_compare.py"
SPEC_YAML = FACTORY_ROOT / "scripts" / "ci" / "nc_regression_tolerance_spec.yaml"
GATE_SH = FACTORY_ROOT / "scripts" / "ci" / "nc_regression_gate.sh"


def _build_run_dir(tmp: Path, name: str, leiden_assignments: dict[str, str], top_de: dict[str, list[str]]) -> Path:
    run = tmp / name
    (run / "annotation").mkdir(parents=True, exist_ok=True)
    (run / "differential_expression").mkdir(parents=True, exist_ok=True)
    # annotation/cell_type_annotation.csv: index=barcode, columns include leiden
    rows = [{"": bc, "leiden": cluster} for bc, cluster in leiden_assignments.items()]
    pd.DataFrame(rows).to_csv(run / "annotation" / "cell_type_annotation.csv", index=False)
    # differential_expression/marker_genes.csv: cluster, gene, fdr
    de_rows = []
    for cid, genes in top_de.items():
        for rank, gene in enumerate(genes):
            de_rows.append({"cluster": cid, "gene": gene, "fdr": 1e-6 * (rank + 1)})
    pd.DataFrame(de_rows).to_csv(run / "differential_expression" / "marker_genes.csv", index=False)
    # final_adata.h5ad: tiny synthetic anndata
    n = len(leiden_assignments)
    X = sp.random(n, 20, density=0.3, dtype=np.float32, random_state=42).tocsr()
    obs = pd.DataFrame({"leiden": [leiden_assignments[bc] for bc in leiden_assignments]},
                       index=list(leiden_assignments.keys()))
    adata = ad.AnnData(X=X, obs=obs)
    adata.write_h5ad(run / "final_adata.h5ad")
    # bundle_manifest.json (optional gate-4)
    (run / "bundle_manifest.json").write_text(json.dumps({
        "schema_version": "singlecell_r_bundle_v2.2",
    }))
    return run


@pytest.fixture
def perfect_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Baseline and candidate with identical leiden + identical top-50 DE."""
    barcodes = {f"cell_{i:04d}": str(i % 4) for i in range(120)}
    de = {str(c): [f"GENE_{c}_{r}" for r in range(50)] for c in range(4)}
    baseline = _build_run_dir(tmp_path, "baseline", barcodes, de)
    candidate = _build_run_dir(tmp_path, "candidate", barcodes, de)
    return baseline, candidate


@pytest.fixture
def divergent_clusters_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Baseline and candidate with completely shuffled leiden (low ARI)."""
    rng = np.random.default_rng(42)
    barcodes_b = {f"cell_{i:04d}": str(i % 4) for i in range(120)}
    # Re-shuffle cluster labels for the candidate so ARI is near zero
    shuffled = rng.permutation([str(i % 4) for i in range(120)])
    barcodes_c = {f"cell_{i:04d}": shuffled[i] for i in range(120)}
    de = {str(c): [f"GENE_{c}_{r}" for r in range(50)] for c in range(4)}
    baseline = _build_run_dir(tmp_path, "baseline", barcodes_b, de)
    candidate = _build_run_dir(tmp_path, "candidate", barcodes_c, de)
    return baseline, candidate


@pytest.fixture
def diverged_de_pair(tmp_path: Path) -> tuple[Path, Path]:
    """Baseline and candidate with identical leiden but non-overlapping top-50 DE."""
    barcodes = {f"cell_{i:04d}": str(i % 4) for i in range(120)}
    de_b = {str(c): [f"GENE_B_{c}_{r}" for r in range(50)] for c in range(4)}
    de_c = {str(c): [f"GENE_C_{c}_{r}" for r in range(50)] for c in range(4)}
    baseline = _build_run_dir(tmp_path, "baseline", barcodes, de_b)
    candidate = _build_run_dir(tmp_path, "candidate", barcodes, de_c)
    return baseline, candidate


def _run_gate(baseline: Path, candidate: Path) -> tuple[int, str, str]:
    """Run nc_regression_compare.py directly (no conda wrapper)."""
    proc = subprocess.run(
        [sys.executable, str(COMPARE_PY),
         "--baseline-run-dir", str(baseline),
         "--candidate-run-dir", str(candidate)],
        capture_output=True, text=True,
    )
    return proc.returncode, proc.stdout, proc.stderr


def test_spec_yaml_loads():
    """Tolerance spec is well-formed YAML with expected top-level keys."""
    import yaml
    doc = yaml.safe_load(SPEC_YAML.read_text())
    assert "cluster_label_ari" in doc
    assert "top_de_markers_jaccard" in doc
    assert "structural_byte_identical" in doc
    assert "bundle_schema" in doc
    assert doc["cluster_label_ari"]["threshold_min"] == 0.99
    assert doc["top_de_markers_jaccard"]["top_k"] == 50


def test_gate_script_exists_and_executable():
    assert GATE_SH.exists()
    assert GATE_SH.stat().st_mode & 0o111, "nc_regression_gate.sh must be executable"


def test_pass_path_all_four_gates(perfect_pair):
    baseline, candidate = perfect_pair
    rc, stdout, stderr = _run_gate(baseline, candidate)
    assert rc == 0, f"expected exit 0; got {rc}\nstdout={stdout}\nstderr={stderr}"
    assert "ALL 4 SUB-GATES PASSED" in stdout


def test_fail_path_ari(divergent_clusters_pair):
    baseline, candidate = divergent_clusters_pair
    rc, stdout, stderr = _run_gate(baseline, candidate)
    assert rc == 10, f"expected exit 10 (ARI fail); got {rc}\nstdout={stdout}"
    assert "[gate-1 cluster_label_ari] FAIL" in stdout


def test_fail_path_de_jaccard(diverged_de_pair):
    baseline, candidate = diverged_de_pair
    rc, stdout, stderr = _run_gate(baseline, candidate)
    assert rc == 11, f"expected exit 11 (Jaccard fail); got {rc}\nstdout={stdout}"
    assert "[gate-2 top_de_markers_jaccard] FAIL" in stdout
