from pathlib import Path

import numpy as np
from anndata import AnnData

from workflow.standard import StandardWorkflowConfig, run_standard_workflow


class DummyScanpy:
    class PP:
        @staticmethod
        def calculate_qc_metrics(adata, qc_vars=None, inplace=True):
            adata.obs["n_genes_by_counts"] = np.sum(adata.X > 0, axis=1)
            adata.obs["pct_counts_mt"] = 0.0

        @staticmethod
        def normalize_total(adata, target_sum=1e4):
            pass

        @staticmethod
        def log1p(adata):
            pass

        @staticmethod
        def highly_variable_genes(adata, flavor="seurat", n_top_genes=2000):
            adata.var["highly_variable"] = True

        @staticmethod
        def scale(adata, max_value=10):
            pass

        @staticmethod
        def neighbors(adata, n_neighbors=30, n_pcs=40):
            adata.uns["neighbors"] = {"ok": True}

        @staticmethod
        def filter_genes(adata, min_cells=3):
            pass

    class TL:
        @staticmethod
        def pca(adata, n_comps=40, svd_solver="arpack"):
            adata.obsm["X_pca"] = np.zeros((adata.n_obs, min(n_comps, 2)))

        @staticmethod
        def umap(adata, random_state=0):
            adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))

        @staticmethod
        def leiden(adata, resolution=0.5, flavor="igraph", directed=False, random_state=0):
            adata.obs["leiden"] = "0"

    pp = PP()
    tl = TL()

    @staticmethod
    def read_10x_mtx(path, var_names="gene_symbols", cache=False):
        x = np.array([[1.0, 0.0], [2.0, 3.0]], dtype=float)
        adata = AnnData(x)
        adata.var_names = ["G1", "G2"]
        adata.obs_names = ["C1", "C2"]
        return adata


def test_run_standard_workflow(monkeypatch, tmp_path):
    import workflow.standard as mod

    monkeypatch.setattr(mod, "sc", DummyScanpy)
    tenx_dir = tmp_path / "tenx"
    tenx_dir.mkdir()
    cfg = StandardWorkflowConfig(
        project="p1",
        tenx_dir=tenx_dir,
        output_dir=tmp_path / "out",
        mode="full",
    )
    out = run_standard_workflow(cfg)
    assert out.exists()
    assert out.name == "processed.h5ad"


def test_run_standard_workflow_missing_input(tmp_path):
    cfg = StandardWorkflowConfig(
        project="p1",
        tenx_dir=tmp_path / "missing",
        output_dir=tmp_path / "out",
        mode="full",
    )
    try:
        run_standard_workflow(cfg)
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("Expected FileNotFoundError")


def test_parse_args_and_main(monkeypatch, tmp_path):
    import workflow.standard as mod

    monkeypatch.setattr(mod, "sc", DummyScanpy)
    tenx_dir = tmp_path / "tenx"
    tenx_dir.mkdir()
    monkeypatch.setattr(
        mod.argparse,
        "ArgumentParser",
        lambda: __import__("argparse").ArgumentParser(prog="x"),
    )
    monkeypatch.setattr(
        mod,
        "parse_args",
        lambda: type(
            "Args",
            (),
            {
                "project": "p2",
                "tenx_dir": str(tenx_dir),
                "output_dir": str(tmp_path / "out2"),
                "mode": "full",
                "n_jobs": 2,
                "no_cache": False,
            },
        )(),
    )
    mod.main()


def test_parse_args_function(monkeypatch):
    import workflow.standard as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--project",
            "p3",
            "--tenx-dir",
            "/tmp/x",
            "--output-dir",
            "/tmp/o",
            "--mode",
            "fast",
            "--n-jobs",
            "4",
            "--no-cache",
        ],
    )
    args = mod.parse_args()
    assert args.project == "p3"
    assert args.mode == "fast"
    assert args.n_jobs == 4
    assert args.no_cache is True


def test_cache_hit(monkeypatch, tmp_path):
    import workflow.standard as mod

    out_dir = tmp_path / "o"
    p = out_dir / "p4"
    p.mkdir(parents=True, exist_ok=True)
    cached = p / "processed.h5ad"
    cached.write_text("x", encoding="utf-8")
    cfg = StandardWorkflowConfig(
        project="p4",
        tenx_dir=tmp_path,
        output_dir=out_dir,
        mode="full",
        use_cache=True,
    )
    assert mod.run_standard_workflow(cfg) == cached
