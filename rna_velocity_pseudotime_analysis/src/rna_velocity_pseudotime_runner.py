import argparse
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--run_id", required=True)
    return parser.parse_args()


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def safe_import(name: str):
    try:
        module = __import__(name)
        return module, None
    except Exception as exc:
        return None, str(exc)


def ensure_isolation() -> None:
    required = [
        "RVPT_NAMESPACE",
        "RVPT_RUNTIME_ROOT",
        "RVPT_CACHE_DIR",
        "RVPT_RESULTS_DIR",
        "RVPT_REPORTS_DIR",
    ]
    missing = [k for k in required if not os.environ.get(k)]
    if missing:
        raise RuntimeError(f"Missing required isolation variables: {','.join(missing)}")
    if os.environ["RVPT_NAMESPACE"] != "rna_velocity_pseudotime_analysis":
        raise RuntimeError("Invalid namespace")


def discover_default_input() -> tuple[str | None, str]:
    root = Path(__file__).resolve().parents[2]
    output_dir = root / "output"
    if not output_dir.exists():
        return None, "output目录不存在"
    candidates = sorted(output_dir.glob("*/data/after_qc.h5ad"), key=lambda p: p.stat().st_mtime, reverse=True)
    if candidates:
        return str(candidates[0]), "auto:after_qc.h5ad"
    candidates = sorted(output_dir.glob("*/data/processed.h5ad"), key=lambda p: p.stat().st_mtime, reverse=True)
    if candidates:
        return str(candidates[0]), "auto:processed.h5ad"
    return None, "未找到可用h5ad输入"


def ensure_obs_embeddings(adata, sc_module):
    if "X_pca" not in adata.obsm:
        sc_module.pp.normalize_total(adata, target_sum=1e4)
        sc_module.pp.log1p(adata)
        sc_module.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars))
        if "highly_variable" in adata.var.columns and int(adata.var["highly_variable"].sum()) > 50:
            adata = adata[:, adata.var["highly_variable"]].copy()
        sc_module.pp.scale(adata, max_value=10)
        sc_module.tl.pca(adata, n_comps=min(40, max(2, adata.n_vars - 1)))
    if "neighbors" not in adata.uns:
        sc_module.pp.neighbors(adata, n_neighbors=min(30, max(10, int(np.sqrt(max(adata.n_obs, 100)) // 2))))
    if "X_umap" not in adata.obsm:
        sc_module.tl.umap(adata)
    if "leiden" not in adata.obs.columns:
        try:
            sc_module.tl.leiden(adata, resolution=0.8)
        except Exception:
            adata.obs["leiden"] = pd.Categorical(["0"] * adata.n_obs)
    return adata


def run_scanpy_pseudotime(adata, sc_module, figures_dir: Path, tables_dir: Path) -> dict:
    status = {"module": "scanpy_pseudotime", "status": "ok", "message": "dpt完成"}
    try:
        if "neighbors" not in adata.uns:
            sc_module.pp.neighbors(adata, n_neighbors=min(30, max(10, int(np.sqrt(max(adata.n_obs, 100)) // 2))))
        if "X_diffmap" not in adata.obsm:
            sc_module.tl.diffmap(adata)
        root_idx = 0
        adata.uns["iroot"] = root_idx
        sc_module.tl.dpt(adata)
        fig_path = figures_dir / "pseudotime_dpt_umap.png"
        sc_module.pl.umap(adata, color=["dpt_pseudotime"], show=False)
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        table = adata.obs[["leiden", "dpt_pseudotime"]].copy()
        table.to_csv(tables_dir / "pseudotime_dpt_by_cell.csv")
        cluster_summary = table.groupby("leiden", observed=False)["dpt_pseudotime"].agg(["count", "mean", "median", "std"]).reset_index()
        cluster_summary.to_csv(tables_dir / "pseudotime_dpt_by_cluster.csv", index=False)
    except Exception as exc:
        status = {"module": "scanpy_pseudotime", "status": "failed", "message": str(exc)}
    return status


def run_velocity_visualization(adata_plot, adata_velocity, cfg: dict, figures_dir: Path, tables_dir: Path) -> dict:
    scv_module, err = safe_import("scvelo")
    if scv_module is None:
        return {"module": "scvelo_velocity", "status": "skipped", "message": f"scvelo不可用: {err}"}
    if "spliced" not in adata_velocity.layers or "unspliced" not in adata_velocity.layers:
        return {"module": "scvelo_velocity", "status": "skipped", "message": "缺少spliced/unspliced层"}
    try:
        adata_scv = adata_velocity.copy()
        if "X_umap" in adata_plot.obsm:
            adata_scv.obsm["X_umap"] = adata_plot.obsm["X_umap"].copy()
        if "leiden" in adata_plot.obs.columns:
            adata_scv.obs["leiden"] = adata_plot.obs["leiden"].astype(str).values
        else:
            adata_scv.obs["leiden"] = "0"
        mode = cfg.get("modules", {}).get("velocity", {}).get("mode", "dynamical")
        try:
            scv_module.pp.filter_and_normalize(adata_scv, min_shared_counts=20, n_top_genes=min(2000, adata_scv.n_vars))
        except TypeError:
            scv_module.pp.filter_and_normalize(adata_scv, min_shared_counts=20)
        scv_module.pp.moments(adata_scv, n_pcs=30, n_neighbors=30)
        velocity_mode = mode
        try:
            scv_module.tl.velocity(adata_scv, mode=velocity_mode)
        except Exception:
            velocity_mode = "stochastic"
            scv_module.tl.velocity(adata_scv, mode=velocity_mode)
        scv_module.tl.velocity_graph(adata_scv)
        scv_module.pl.velocity_embedding_stream(adata_scv, basis="umap", color="leiden", show=False)
        plt.savefig(figures_dir / "velocity_stream_umap.png", dpi=150, bbox_inches="tight")
        plt.close()
        scv_module.pl.velocity_embedding_grid(adata_scv, basis="umap", color="leiden", show=False)
        plt.savefig(figures_dir / "velocity_grid_umap.png", dpi=150, bbox_inches="tight")
        plt.close()
        speed = np.sqrt(np.asarray(adata_scv.layers["velocity"]).power(2).sum(axis=1)).A1 if hasattr(adata_scv.layers["velocity"], "power") else np.sqrt((np.asarray(adata_scv.layers["velocity"]) ** 2).sum(axis=1))
        adata_plot.obs["velocity_speed"] = speed
        idx = pd.Index(adata_plot.var_names).get_indexer(adata_scv.var_names)
        def project_layer(layer_name: str):
            layer = adata_scv.layers[layer_name]
            if sp.issparse(layer):
                coo = layer.tocoo()
                valid = idx[coo.col] >= 0
                return sp.csr_matrix(
                    (coo.data[valid], (coo.row[valid], idx[coo.col[valid]])),
                    shape=(adata_plot.n_obs, adata_plot.n_vars),
                )
            out = np.zeros((adata_plot.n_obs, adata_plot.n_vars), dtype=np.float32)
            valid_cols = idx >= 0
            out[:, idx[valid_cols]] = np.asarray(layer)[:, valid_cols]
            return out
        adata_plot.layers["velocity"] = project_layer("velocity")
        if "Ms" in adata_scv.layers:
            adata_plot.layers["Ms"] = project_layer("Ms")
        if "Mu" in adata_scv.layers:
            adata_plot.layers["Mu"] = project_layer("Mu")
        if "velocity_graph" in adata_scv.uns:
            adata_plot.uns["velocity_graph"] = adata_scv.uns["velocity_graph"]
        adata_plot.obs[["leiden", "velocity_speed"]].groupby("leiden", observed=False).agg(["count", "mean", "median", "std"]).to_csv(
            tables_dir / "velocity_speed_by_cluster.csv"
        )
        return {"module": "scvelo_velocity", "status": "ok", "message": f"mode={velocity_mode}"}
    except Exception as exc:
        return {"module": "scvelo_velocity", "status": "failed", "message": str(exc)}


def run_cellrank_projection(adata, figures_dir: Path) -> dict:
    cr_module, err = safe_import("cellrank")
    if cr_module is None:
        return {"module": "cellrank_projection", "status": "skipped", "message": f"cellrank不可用: {err}"}
    if "velocity" not in adata.layers:
        return {"module": "cellrank_projection", "status": "skipped", "message": "缺少velocity层"}
    try:
        from cellrank.kernels import VelocityKernel

        vk = VelocityKernel(adata)
        vk.compute_transition_matrix()
        ax = vk.plot_projection(basis="umap", color="leiden", show=False)
        fig = ax.figure if ax is not None else plt.gcf()
        fig.savefig(figures_dir / "cellrank_velocity_projection_umap.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        return {"module": "cellrank_projection", "status": "ok", "message": "transition matrix完成"}
    except Exception as exc:
        return {"module": "cellrank_projection", "status": "failed", "message": str(exc)}


def write_manifest(result_dir: Path, statuses: list[dict], input_info: dict, run_id: str):
    manifest = pd.DataFrame(statuses)
    manifest.to_csv(result_dir / "tables" / "module_status.csv", index=False)
    with open(result_dir / "rna_velocity_pseudotime_summary.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "run_id": run_id,
                "generated_at": datetime.utcnow().isoformat() + "Z",
                "input": input_info,
                "module_status": statuses,
            },
            f,
            indent=2,
            ensure_ascii=False,
        )


def write_report(report_dir: Path, run_id: str, dataset: str, statuses: list[dict], input_info: dict):
    with open(report_dir / "rna_velocity_pseudotime_report.md", "w", encoding="utf-8") as f:
        f.write("# RNA Velocity and Pseudotime Analysis Report\n\n")
        f.write("可选分析步骤，非生产必需。\n\n")
        f.write(f"- Run ID: {run_id}\n")
        f.write(f"- Dataset: {dataset}\n")
        f.write(f"- Input Source: {input_info.get('source')}\n")
        f.write(f"- Input Path: {input_info.get('path')}\n")
        f.write(f"- Generated At: {datetime.utcnow().isoformat()}Z\n\n")
        f.write("## Module Status\n\n")
        for s in statuses:
            f.write(f"- {s['module']}: {s['status']} ({s['message']})\n")


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    ensure_isolation()

    report_dir = Path(os.environ["RVPT_REPORTS_DIR"])
    result_dir = Path(os.environ["RVPT_RESULTS_DIR"])
    figures_dir = result_dir / "figures"
    tables_dir = result_dir / "tables"
    report_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    sc_module, sc_err = safe_import("scanpy")
    if sc_module is None:
        statuses = [{"module": "scanpy_core", "status": "failed", "message": f"scanpy不可用: {sc_err}"}]
        input_info = {"source": "none", "path": ""}
        write_manifest(result_dir, statuses, input_info, args.run_id)
        write_report(report_dir, args.run_id, args.dataset, statuses, input_info)
        return

    configured_h5ad = cfg.get("input", {}).get("h5ad", "")
    configured_tenx = cfg.get("input", {}).get("tenx_dir", "")
    input_path = ""
    input_source = ""
    if configured_h5ad and Path(configured_h5ad).is_file():
        input_path = configured_h5ad
        input_source = "config:h5ad"
    elif configured_tenx and Path(configured_tenx).is_dir():
        input_path = configured_tenx
        input_source = "config:tenx_dir"
    else:
        discovered, source = discover_default_input()
        if discovered:
            input_path = discovered
            input_source = source

    statuses = []
    if not input_path:
        statuses.append({"module": "input_loader", "status": "failed", "message": "未找到可用输入数据"})
        input_info = {"source": "none", "path": ""}
        write_manifest(result_dir, statuses, input_info, args.run_id)
        write_report(report_dir, args.run_id, args.dataset, statuses, input_info)
        return

    input_info = {"source": input_source, "path": input_path}
    try:
        if Path(input_path).is_dir():
            adata = sc_module.read_10x_mtx(input_path, var_names="gene_symbols", cache=False)
        else:
            adata = sc_module.read_h5ad(input_path)
        statuses.append({"module": "input_loader", "status": "ok", "message": input_source})
    except Exception as exc:
        statuses.append({"module": "input_loader", "status": "failed", "message": str(exc)})
        write_manifest(result_dir, statuses, input_info, args.run_id)
        write_report(report_dir, args.run_id, args.dataset, statuses, input_info)
        return

    loom_path = cfg.get("input", {}).get("loom", "")
    if loom_path and Path(loom_path).is_file():
        scv_module, _ = safe_import("scvelo")
        if scv_module is not None:
            try:
                ldata = sc_module.read_loom(loom_path)
                adata = scv_module.utils.merge(adata, ldata)
                statuses.append({"module": "loom_merge", "status": "ok", "message": loom_path})
            except Exception as exc:
                statuses.append({"module": "loom_merge", "status": "failed", "message": str(exc)})
        else:
            statuses.append({"module": "loom_merge", "status": "skipped", "message": "scvelo不可用"})
    else:
        statuses.append({"module": "loom_merge", "status": "skipped", "message": "未提供loom"})

    adata_velocity = adata.copy()

    try:
        adata = ensure_obs_embeddings(adata, sc_module)
        sc_module.pl.umap(adata, color=["leiden"], show=False)
        plt.savefig(figures_dir / "umap_leiden.png", dpi=150, bbox_inches="tight")
        plt.close()
        statuses.append({"module": "scanpy_embedding", "status": "ok", "message": "umap/leiden已就绪"})
    except Exception as exc:
        statuses.append({"module": "scanpy_embedding", "status": "failed", "message": str(exc)})

    statuses.append(run_scanpy_pseudotime(adata, sc_module, figures_dir, tables_dir))
    statuses.append(run_velocity_visualization(adata, adata_velocity, cfg, figures_dir, tables_dir))
    statuses.append(run_cellrank_projection(adata, figures_dir))

    try:
        adata.write(result_dir / "integrated.h5ad")
        statuses.append({"module": "result_writer", "status": "ok", "message": "integrated.h5ad已保存"})
    except Exception as exc:
        statuses.append({"module": "result_writer", "status": "failed", "message": str(exc)})

    write_manifest(result_dir, statuses, input_info, args.run_id)
    write_report(report_dir, args.run_id, args.dataset, statuses, input_info)


if __name__ == "__main__":
    main()
