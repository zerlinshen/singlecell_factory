from pathlib import Path
import importlib.util
import sys

from anndata import AnnData
import numpy as np
import pandas as pd
from scipy import sparse


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "nc2024_css_fidelity_audit.py"
SPEC = importlib.util.spec_from_file_location("nc2024_css_fidelity_audit", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


def test_stratified_indices_keeps_extra_indices():
    obs = pd.DataFrame(
        {
            "cell_type": ["A"] * 8 + ["B"] * 2,
            "immune_subtype": ["x"] * 5 + ["y"] * 5,
            "disease": ["LUAD"] * 10,
        }
    )
    idx = audit.stratified_indices(obs, subset_size=5, seed=1, min_per_stratum=1, extra_indices=[8, 9])
    assert {8, 9}.issubset(set(idx))
    assert len(idx) == 5


def test_cluster_majority_metrics_reports_match_fraction():
    clusters = np.array(["0", "0", "1", "1"])
    truth = pd.Series(["T", "T", "B", "T"])
    metrics = audit.cluster_majority_metrics(clusters, truth)
    assert metrics["n_clusters"] == 2
    assert metrics["majority_match_fraction"] == 0.75


def test_rare_group_metrics_includes_low_frequency_group(tmp_path):
    cfg = audit.AuditConfig(
        input_h5ad=tmp_path / "x.h5ad",
        output_dir=tmp_path,
        subset_size=10,
        seed=1,
        n_top_genes=10,
        n_pcs=2,
        n_neighbors=2,
        leiden_resolution=0.4,
        batch_key="sample",
        rare_min_cells=2,
        rare_max_fraction=0.3,
        min_per_stratum=1,
        elf3_top_n=0,
    )
    obs = pd.DataFrame({"cell_type": ["common"] * 7 + ["rare"] * 3})
    clean = np.array(["0"] * 7 + ["1", "1", "2"])
    css = np.array(["0"] * 7 + ["1", "2", "2"])
    rare = audit.rare_group_metrics(clean, css, obs, cfg)
    assert "rare" in set(rare["label"])


def test_expression_group_summary_reports_full_distribution():
    obs = pd.DataFrame({"cell_type": ["Tumor", "Tumor", "Macro"]})
    expr = np.array([0.0, 2.0, 1.0], dtype=np.float32)
    table, metrics = audit.expression_group_summary(expr, obs, "ELF3", keys=("cell_type",))
    assert metrics["elf3_positive_n_full"] == 2
    tumor = table[table["label"] == "Tumor"].iloc[0]
    assert tumor["ELF3_positive_n"] == 1


def test_expression_group_summary_handles_absent_gene():
    obs = pd.DataFrame({"cell_type": ["Tumor", "Macro"]})
    table, metrics = audit.expression_group_summary(None, obs, "ELF3", keys=("cell_type",))
    assert table.empty
    assert metrics["elf3_present_full"] is False


def test_run_audit_end_to_end_tiny_h5ad(tmp_path):
    rng = np.random.default_rng(7)
    x = sparse.csr_matrix(rng.poisson(0.3, size=(48, 16)).astype(np.float32))
    var_names = ["ELF3"] + [f"G{i}" for i in range(15)]
    adata = AnnData(x)
    adata.var_names = var_names
    adata.obs_names = [f"cell{i}" for i in range(48)]
    adata.obs["sample"] = ["s1"] * 16 + ["s2"] * 16 + ["s3"] * 16
    adata.obs["cell_type"] = ["Tumor"] * 12 + ["Macro"] * 18 + ["T"] * 18
    adata.obs["immune_subtype"] = ["Non-immune"] * 12 + ["Macro_M2"] * 18 + ["CD4_memory"] * 18
    adata.obs["disease"] = ["LUAD"] * 24 + ["LUSC"] * 24
    adata.obs["leiden"] = [str(i % 4) for i in range(48)]
    adata.var["highly_variable"] = [True] * 12 + [False] * 4
    adata.var["dispersions_norm"] = np.linspace(0.0, 1.0, 16)
    h5ad = tmp_path / "tiny.h5ad"
    adata.write_h5ad(h5ad)

    cfg = audit.AuditConfig(
        input_h5ad=h5ad,
        output_dir=tmp_path / "audit",
        subset_size=36,
        seed=11,
        n_top_genes=8,
        n_pcs=2,
        n_neighbors=3,
        leiden_resolution=0.4,
        batch_key="sample",
        rare_min_cells=2,
        rare_max_fraction=0.5,
        min_per_stratum=1,
        elf3_top_n=5,
    )
    summary = audit.run_audit(cfg)
    assert summary.exists()
    data = __import__("json").loads(summary.read_text(encoding="utf-8"))
    assert data["subset_size_actual"] == 36
    assert "clean_css_ari" in data
    assert data["elf3_present"] is True
    assert (cfg.output_dir / "elf3_full_cohort_group_summary.csv").exists()
