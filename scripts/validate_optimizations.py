import argparse
import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.metrics import adjusted_rand_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class RunResult:
    mode: str
    wall_seconds: float
    n_cells: int
    n_genes: int
    labels: pd.Series
    top_markers: dict[str, set[str]]
    qc_audit: dict[str, float]


def _load_adata(tenx_path: str):
    adata = sc.read_10x_mtx(tenx_path, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    return adata


def _extract_top_markers(adata, key: str, n: int = 50) -> dict[str, set[str]]:
    df = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    result: dict[str, set[str]] = {}
    for group in df["group"].astype(str).unique():
        top = df[df["group"].astype(str) == group].head(n)["names"].astype(str)
        result[group] = set(top.tolist())
    return result


def _run_baseline(adata, seed: int) -> RunResult:
    t0 = time.perf_counter()
    qc_before_cells = int(adata.n_obs)
    qc_before_genes = int(adata.n_vars)

    gene_upper = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = gene_upper.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)
    sc.tl.pca(adata, n_comps=40, svd_solver="arpack", mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=0.8,
        random_state=seed,
        flavor="igraph",
        directed=False,
        key_added="leiden_baseline",
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden_baseline",
        method="wilcoxon",
        use_raw=True,
        n_genes=300,
        key_added="de_baseline",
    )

    wall = time.perf_counter() - t0
    markers = _extract_top_markers(adata, key="de_baseline", n=50)
    return RunResult(
        mode="baseline",
        wall_seconds=float(wall),
        n_cells=int(adata.n_obs),
        n_genes=int(adata.n_vars),
        labels=adata.obs["leiden_baseline"].astype(str).copy(),
        top_markers=markers,
        qc_audit={
            "cells_before_qc": qc_before_cells,
            "genes_before_qc": qc_before_genes,
            "cells_after_qc": int(adata.n_obs),
            "genes_after_qc": int(adata.n_vars),
            "mean_pct_mt": float(adata.obs["pct_counts_mt"].mean()),
        },
    )


def _run_optimized(adata, seed: int, de_method: str = "wilcoxon") -> RunResult:
    t0 = time.perf_counter()
    qc_before_cells = int(adata.n_obs)
    qc_before_genes = int(adata.n_vars)

    gene_upper = adata.var_names.astype(str).str.upper()
    adata.var["mt"] = gene_upper.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True,
        log1p=False,
        percent_top=None,
    )
    obs_drop = [
        c for c in adata.obs.columns
        if c.startswith("log1p") or c.startswith("total_counts_") or c.startswith("pct_counts_in_top_")
    ]
    adata.obs.drop(columns=obs_drop, inplace=True, errors="ignore")

    mask = adata.obs["n_genes_by_counts"] >= 200
    adata = adata[mask, :].copy()
    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)
    solver = "arpack" if sparse.issparse(adata.X) else "randomized"
    sc.tl.pca(adata, n_comps=40, svd_solver=solver, mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=0.8,
        random_state=seed,
        flavor="igraph",
        directed=False,
        key_added="leiden_optimized",
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden_optimized",
        method=de_method,
        use_raw=False,
        n_genes=300,
        key_added="de_optimized",
    )

    wall = time.perf_counter() - t0
    markers = _extract_top_markers(adata, key="de_optimized", n=50)
    return RunResult(
        mode="optimized",
        wall_seconds=float(wall),
        n_cells=int(adata.n_obs),
        n_genes=int(adata.n_vars),
        labels=adata.obs["leiden_optimized"].astype(str).copy(),
        top_markers=markers,
        qc_audit={
            "cells_before_qc": qc_before_cells,
            "genes_before_qc": qc_before_genes,
            "cells_after_qc": int(adata.n_obs),
            "genes_after_qc": int(adata.n_vars),
            "mean_pct_mt": float(adata.obs["pct_counts_mt"].mean()),
        },
    )


def _pairwise_precision_recall_f1(labels_a: pd.Series, labels_b: pd.Series) -> tuple[float, float, float]:
    a = labels_a.astype(str).to_numpy()
    b = labels_b.astype(str).to_numpy()
    n = a.shape[0]
    if n < 3:
        return 1.0, 1.0, 1.0

    same_a = a[:, None] == a[None, :]
    same_b = b[:, None] == b[None, :]
    iu = np.triu_indices(n, k=1)
    ya = same_a[iu]
    yb = same_b[iu]

    tp = np.sum(ya & yb)
    fp = np.sum(~ya & yb)
    fn = np.sum(ya & ~yb)

    precision = float(tp / max(tp + fp, 1))
    recall = float(tp / max(tp + fn, 1))
    f1 = 0.0 if (precision + recall) == 0 else float(2 * precision * recall / (precision + recall))
    return precision, recall, f1


def _marker_f1_overlap(base_markers: dict[str, set[str]], opt_markers: dict[str, set[str]]) -> tuple[float, list[float]]:
    per_cluster_f1 = []
    for _, opt_set in opt_markers.items():
        best_f1 = 0.0
        for _, base_set in base_markers.items():
            inter = len(base_set & opt_set)
            p = inter / max(len(opt_set), 1)
            r = inter / max(len(base_set), 1)
            f1 = 0.0 if (p + r) == 0 else 2 * p * r / (p + r)
            if f1 > best_f1:
                best_f1 = f1
        per_cluster_f1.append(float(best_f1))

    if not per_cluster_f1:
        return 0.0, []
    return float(np.mean(per_cluster_f1)), per_cluster_f1


def _ci95(arr: np.ndarray) -> tuple[float, float]:
    if len(arr) <= 1:
        v = float(arr[0]) if len(arr) else 0.0
        return v, v
    mean = float(np.mean(arr))
    sem = float(stats.sem(arr))
    if sem == 0.0 or np.isnan(sem):
        return mean, mean
    lo, hi = stats.t.interval(0.95, df=len(arr) - 1, loc=mean, scale=sem)
    return float(lo), float(hi)


def _scenario_drop_extreme(adata, rng: np.random.Generator):
    if sparse.issparse(adata.X):
        coo = adata.X.tocoo(copy=True)
        keep = rng.random(coo.data.shape[0]) > 0.90
        coo.data = coo.data[keep]
        coo.row = coo.row[keep]
        coo.col = coo.col[keep]
        adata.X = coo.tocsr()
    else:
        mask = rng.random(adata.X.shape) > 0.90
        adata.X = adata.X * mask
    return adata


def _scenario_low_cells(adata, rng: np.random.Generator):
    n = min(80, adata.n_obs)
    idx = rng.choice(np.arange(adata.n_obs), size=n, replace=False)
    return adata[idx].copy()


def _scenario_nan_injection(adata):
    dense = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X).copy()
    dense[:5, :5] = np.nan
    adata.X = sparse.csr_matrix(np.nan_to_num(dense, nan=0.0))
    return adata


def _run_robustness_tests(tenx_path: str, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    scenarios = {
        "extreme_dropout_90pct": lambda x: _scenario_drop_extreme(x, rng),
        "low_cell_count_80": lambda x: _scenario_low_cells(x, rng),
        "nan_injection_sanitized": _scenario_nan_injection,
    }

    rows = []
    for name, fn in scenarios.items():
        adata = _load_adata(tenx_path)
        try:
            adata = fn(adata)
            _ = _run_optimized(adata, seed=seed, de_method="wilcoxon")
            rows.append({"scenario": name, "passed": True, "error": ""})
        except Exception as exc:
            rows.append({"scenario": name, "passed": False, "error": str(exc)})
    return rows


def _save_visuals(out_dir: Path, perf_df: pd.DataFrame, metric_df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.boxplot(
        [perf_df["baseline_time_sec"].to_numpy(), perf_df["optimized_time_sec"].to_numpy()],
        tick_labels=["baseline", "optimized"],
    )
    ax.set_title("Runtime Distribution")
    ax.set_ylabel("Seconds")
    fig.tight_layout()
    fig.savefig(out_dir / "runtime_boxplot.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(metric_df["cluster_ari"], bins=10, alpha=0.7, label="ARI")
    ax.hist(metric_df["pairwise_f1"], bins=10, alpha=0.7, label="Pairwise F1")
    ax.set_title("Metric Distribution Across Repeats")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "metric_distribution.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(metric_df["marker_f1"], bins=10, alpha=0.8)
    ax.set_title("Marker F1 Distribution")
    ax.set_xlabel("Marker F1")
    fig.tight_layout()
    fig.savefig(out_dir / "marker_f1_distribution.png", dpi=150)
    plt.close(fig)


def _write_markdown_report(out_dir: Path, report: dict):
    md = []
    md.append("# LUSC Optimization Validation Report")
    md.append("")
    md.append(f"- Dataset: `{report['dataset']}`")
    md.append(f"- Repeats: `{report['repeats']}`")
    md.append(f"- Acceptance loss threshold: `{report['loss_threshold_pct']}%`")
    md.append("")
    md.append("## Performance")
    md.append(
        f"- Baseline mean time: `{report['performance']['baseline_mean_sec']:.3f}s` "
        f"(95% CI `{report['performance']['baseline_ci95_sec'][0]:.3f}` - "
        f"`{report['performance']['baseline_ci95_sec'][1]:.3f}`)"
    )
    md.append(
        f"- Optimized mean time: `{report['performance']['optimized_mean_sec']:.3f}s` "
        f"(95% CI `{report['performance']['optimized_ci95_sec'][0]:.3f}` - "
        f"`{report['performance']['optimized_ci95_sec'][1]:.3f}`)"
    )
    md.append(f"- Mean speedup: `{report['performance']['speedup_pct_mean']:.2f}%`")
    md.append(f"- Runtime significance p-value: `{report['performance']['runtime_p_value']:.6f}`")
    md.append("")
    md.append("## Accuracy & Reliability")
    md.append(f"- Cluster ARI mean: `{report['accuracy']['cluster_ari_mean']:.4f}`")
    md.append(
        f"- Pairwise Precision/Recall/F1 mean: "
        f"`{report['accuracy']['pairwise_precision_mean']:.4f}` / "
        f"`{report['accuracy']['pairwise_recall_mean']:.4f}` / "
        f"`{report['accuracy']['pairwise_f1_mean']:.4f}`"
    )
    md.append(f"- Marker F1 mean: `{report['accuracy']['marker_f1_mean']:.4f}`")
    md.append(
        f"- Max relative loss across key metrics: `{report['accuracy']['max_metric_loss_pct']:.4f}%` "
        f"(threshold `{report['loss_threshold_pct']}%`)"
    )
    md.append(f"- Accuracy significance p-value: `{report['accuracy']['accuracy_p_value']:.6f}`")
    md.append("")
    md.append("## Algorithm Shift (Wilcoxon -> t-test_overestim_var)")
    md.append(
        f"- Marker F1 (same clustering, different DE test): "
        f"`{report['algorithm_shift']['marker_f1_mean']:.4f}`"
    )
    md.append(f"- Interpretation: {report['algorithm_shift']['interpretation']}")
    md.append("")
    md.append("## Robustness")
    for row in report["robustness"]["scenarios"]:
        state = "PASS" if row["passed"] else "FAIL"
        msg = f"- `{row['scenario']}`: **{state}**"
        if row["error"]:
            msg += f" (`{row['error']}`)"
        md.append(msg)
    md.append("")
    md.append("## Rollback Decision")
    md.append(f"- Final decision: **{report['decision']['status']}**")
    md.append(f"- Reason: {report['decision']['reason']}")
    (out_dir / "optimization_validation_report.md").write_text("\n".join(md), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Validate acceleration with LUSC real run.")
    parser.add_argument(
        "--lusc-path",
        default="/home/zerlinshen/singlecell_factory/data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix",
    )
    parser.add_argument("--out-dir", default="/home/zerlinshen/singlecell_factory/reports/validation")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--loss-threshold-pct", type=float, default=0.5)
    parser.add_argument("--strict-mode", action="store_true")
    parser.add_argument("--min-ari-lb", type=float, default=0.995)
    parser.add_argument("--min-pairwise-f1-lb", type=float, default=0.995)
    parser.add_argument("--min-marker-f1-lb", type=float, default=0.995)
    parser.add_argument("--min-shift-marker-f1", type=float, default=0.85)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Warm-up to reduce one-time JIT/cache bias.
    for i in range(args.warmup_runs):
        _ = _run_baseline(_load_adata(args.lusc_path), seed=args.seed + i)
        _ = _run_optimized(_load_adata(args.lusc_path), seed=args.seed + i, de_method="wilcoxon")

    perf_rows = []
    metric_rows = []
    audit_rows = []

    for i in range(args.repeats):
        seed = args.seed + i
        logger.info("Repeat %d/%d (seed=%d)", i + 1, args.repeats, seed)

        base = _run_baseline(_load_adata(args.lusc_path), seed=seed)
        opt = _run_optimized(_load_adata(args.lusc_path), seed=seed, de_method="wilcoxon")

        common = sorted(set(base.labels.index) & set(opt.labels.index))
        base_labels = base.labels.loc[common]
        opt_labels = opt.labels.loc[common]

        ari = float(adjusted_rand_score(base_labels, opt_labels))
        precision, recall, f1 = _pairwise_precision_recall_f1(base_labels, opt_labels)
        marker_f1_mean, marker_f1_by_cluster = _marker_f1_overlap(base.top_markers, opt.top_markers)

        perf_rows.append(
            {
                "repeat": i + 1,
                "seed": seed,
                "baseline_time_sec": base.wall_seconds,
                "optimized_time_sec": opt.wall_seconds,
                "speedup_pct": (base.wall_seconds - opt.wall_seconds) / max(base.wall_seconds, 1e-9) * 100.0,
            }
        )
        metric_rows.append(
            {
                "repeat": i + 1,
                "seed": seed,
                "cluster_ari": ari,
                "pairwise_precision": precision,
                "pairwise_recall": recall,
                "pairwise_f1": f1,
                "marker_f1": marker_f1_mean,
                "marker_f1_by_cluster": json.dumps(marker_f1_by_cluster),
            }
        )
        audit_rows.append(
            {
                "repeat": i + 1,
                "mode": "baseline",
                **base.qc_audit,
                "n_cells_final": base.n_cells,
                "n_genes_final": base.n_genes,
            }
        )
        audit_rows.append(
            {
                "repeat": i + 1,
                "mode": "optimized",
                **opt.qc_audit,
                "n_cells_final": opt.n_cells,
                "n_genes_final": opt.n_genes,
            }
        )

    perf_df = pd.DataFrame(perf_rows)
    metric_df = pd.DataFrame(metric_rows)
    audit_df = pd.DataFrame(audit_rows)
    perf_df.to_csv(out_dir / "perf_runs.csv", index=False)
    metric_df.to_csv(out_dir / "accuracy_runs.csv", index=False)
    audit_df.to_csv(out_dir / "data_quality_audit.csv", index=False)

    robustness_rows = _run_robustness_tests(args.lusc_path, seed=args.seed)
    pd.DataFrame(robustness_rows).to_csv(out_dir / "robustness_results.csv", index=False)

    runtime_p = float(
        stats.wilcoxon(
            perf_df["baseline_time_sec"],
            perf_df["optimized_time_sec"],
            alternative="greater",
        ).pvalue
    )
    runtime_p = max(0.0, min(1.0, runtime_p))

    acc_loss = (1.0 - metric_df[["cluster_ari", "pairwise_f1", "marker_f1"]]).clip(lower=0.0) * 100.0
    loss_per_repeat = acc_loss.max(axis=1)
    max_loss = float(loss_per_repeat.mean())
    worst_loss = float(loss_per_repeat.max())

    try:
        acc_p = float(
            stats.ttest_1samp(
                loss_per_repeat.to_numpy(),
                popmean=args.loss_threshold_pct,
                alternative="less",
            ).pvalue
        )
    except TypeError:
        acc_p = float(
            stats.ttest_1samp(loss_per_repeat.to_numpy(), popmean=args.loss_threshold_pct).pvalue
        )
    acc_p = max(0.0, min(1.0, acc_p))

    baseline_ci = _ci95(perf_df["baseline_time_sec"].to_numpy(dtype=float))
    optimized_ci = _ci95(perf_df["optimized_time_sec"].to_numpy(dtype=float))
    ari_ci = _ci95(metric_df["cluster_ari"].to_numpy(dtype=float))
    f1_ci = _ci95(metric_df["pairwise_f1"].to_numpy(dtype=float))
    marker_ci = _ci95(metric_df["marker_f1"].to_numpy(dtype=float))

    robust_pass = all(row["passed"] for row in robustness_rows)
    loss_pass = worst_loss <= args.loss_threshold_pct
    stat_pass = runtime_p < 0.05 and acc_p < 0.05

    # Separate algorithm-change assessment: same optimized path, only DE test changes.
    shift_scores = []
    for i in range(args.repeats):
        seed = args.seed + 100 + i
        opt_w = _run_optimized(_load_adata(args.lusc_path), seed=seed, de_method="wilcoxon")
        opt_t = _run_optimized(_load_adata(args.lusc_path), seed=seed, de_method="t-test_overestim_var")
        marker_f1_shift, _ = _marker_f1_overlap(opt_w.top_markers, opt_t.top_markers)
        shift_scores.append(marker_f1_shift)
    shift_mean = float(np.mean(shift_scores)) if shift_scores else 0.0

    ci_pass = (
        ari_ci[0] >= args.min_ari_lb
        and f1_ci[0] >= args.min_pairwise_f1_lb
        and marker_ci[0] >= args.min_marker_f1_lb
    )
    shift_pass = shift_mean >= args.min_shift_marker_f1

    if args.strict_mode:
        decision_ok = robust_pass and loss_pass and stat_pass and ci_pass and shift_pass
    else:
        decision_ok = robust_pass and loss_pass and stat_pass

    if decision_ok:
        decision_status = "PASS"
        decision_reason = "All thresholds met: speedup significant, accuracy loss <= threshold, robustness passed."
    else:
        decision_status = "ROLLBACK_REQUIRED"
        reasons = []
        if not robust_pass:
            reasons.append("robustness scenario failed")
        if not loss_pass:
            reasons.append("accuracy loss exceeded threshold")
        if not stat_pass:
            reasons.append("statistical significance criteria not satisfied")
        if args.strict_mode and not ci_pass:
            reasons.append("strict CI lower-bound criteria not satisfied")
        if args.strict_mode and not shift_pass:
            reasons.append("algorithm-shift marker consistency too low")
        decision_reason = "; ".join(reasons)

    report = {
        "dataset": "LUSC_3K_real_run",
        "repeats": int(args.repeats),
        "loss_threshold_pct": float(args.loss_threshold_pct),
        "strict_mode": bool(args.strict_mode),
        "strict_thresholds": {
            "min_ari_lb": float(args.min_ari_lb),
            "min_pairwise_f1_lb": float(args.min_pairwise_f1_lb),
            "min_marker_f1_lb": float(args.min_marker_f1_lb),
            "min_shift_marker_f1": float(args.min_shift_marker_f1),
        },
        "performance": {
            "baseline_mean_sec": float(perf_df["baseline_time_sec"].mean()),
            "optimized_mean_sec": float(perf_df["optimized_time_sec"].mean()),
            "baseline_ci95_sec": [baseline_ci[0], baseline_ci[1]],
            "optimized_ci95_sec": [optimized_ci[0], optimized_ci[1]],
            "speedup_pct_mean": float(perf_df["speedup_pct"].mean()),
            "runtime_p_value": runtime_p,
        },
        "accuracy": {
            "cluster_ari_mean": float(metric_df["cluster_ari"].mean()),
            "cluster_ari_ci95": [ari_ci[0], ari_ci[1]],
            "pairwise_precision_mean": float(metric_df["pairwise_precision"].mean()),
            "pairwise_recall_mean": float(metric_df["pairwise_recall"].mean()),
            "pairwise_f1_mean": float(metric_df["pairwise_f1"].mean()),
            "pairwise_f1_ci95": [f1_ci[0], f1_ci[1]],
            "marker_f1_mean": float(metric_df["marker_f1"].mean()),
            "marker_f1_ci95": [marker_ci[0], marker_ci[1]],
            "max_metric_loss_pct": max_loss,
            "worst_metric_loss_pct": worst_loss,
            "accuracy_p_value": acc_p,
            "ci_lower_bound_passed": bool(ci_pass),
        },
        "robustness": {
            "scenarios": robustness_rows,
        },
        "algorithm_shift": {
            "marker_f1_mean": shift_mean,
            "passed": bool(shift_pass),
            "interpretation": ">=0.85 means algorithm substitution is close; <0.85 suggests keeping wilcoxon in production.",
        },
        "decision": {
            "status": decision_status,
            "reason": decision_reason,
            "auto_rollback_steps": [
                "Set DE method back to wilcoxon in differential_expression module.",
                "Disable accelerated branch in clustering if loss persists.",
                "Rerun validation script with --repeats 5 and inspect robustness_results.csv.",
                "Compare accuracy_runs.csv per repeat to locate unstable metric.",
            ],
        },
    }

    _save_visuals(out_dir, perf_df, metric_df)
    _write_markdown_report(out_dir, report)
    (out_dir / "optimization_validation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    logger.info("Validation report generated under %s", out_dir)


if __name__ == "__main__":
    main()
