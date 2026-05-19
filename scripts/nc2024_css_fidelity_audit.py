from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import sys
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


DEFAULT_INPUT = Path(
    "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results/"
    "NC2024_NSCLC_FULL_COHORT_STAGE1_MASSIVE_FRESH_RERUN_AUTO_20260424_020329/"
    "final_adata.h5ad"
)
DEFAULT_OUTPUT_ROOT = Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/results")


@dataclass(frozen=True)
class AuditConfig:
    input_h5ad: Path
    output_dir: Path
    subset_size: int
    seed: int
    n_top_genes: int
    n_pcs: int
    n_neighbors: int
    leiden_resolution: float
    batch_key: str
    rare_min_cells: int
    rare_max_fraction: float
    min_per_stratum: int
    elf3_top_n: int


def _json_default(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default), encoding="utf-8")


def memory_math(input_h5ad: Path) -> dict:
    with h5py.File(input_h5ad, "r") as handle:
        x = handle["X"]
        shape = tuple(int(v) for v in x.attrs["shape"])
        parts = {
            name: {
                "shape": tuple(int(v) for v in x[name].shape),
                "dtype": str(x[name].dtype),
                "nbytes": int(x[name].size * x[name].dtype.itemsize),
            }
            for name in ("data", "indices", "indptr")
        }
        csr_core = sum(v["nbytes"] for v in parts.values())
    n_obs, n_vars = shape
    dense_float32 = n_obs * n_vars * 4
    dense_float64 = n_obs * n_vars * 8
    nnz = parts["data"]["shape"][0]
    return {
        "shape": [n_obs, n_vars],
        "x_encoding": "csr_matrix",
        "x_parts": parts,
        "x_csr_core_gb": csr_core / 1e9,
        "dense_float32_gb": dense_float32 / 1e9,
        "dense_float64_gb": dense_float64 / 1e9,
        "density": nnz / float(n_obs * n_vars),
    }


def _series_frame(obs: pd.DataFrame, cols: Iterable[str]) -> pd.Series:
    present = [c for c in cols if c in obs.columns]
    if not present:
        return pd.Series("all_cells", index=obs.index)
    frame = obs[present].astype(str).fillna("NA")
    return frame.agg("|".join, axis=1)


def stratified_indices(
    obs: pd.DataFrame,
    subset_size: int,
    seed: int,
    min_per_stratum: int,
    extra_indices: Iterable[int] = (),
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_obs = len(obs)
    target = min(int(subset_size), n_obs)
    selected: set[int] = {int(i) for i in extra_indices if 0 <= int(i) < n_obs}
    if len(selected) >= target:
        return np.array(sorted(rng.choice(np.array(sorted(selected)), size=target, replace=False)), dtype=np.int64)

    strata = _series_frame(obs, ("cell_type", "immune_subtype", "disease"))
    available = np.ones(n_obs, dtype=bool)
    if selected:
        available[np.array(sorted(selected), dtype=np.int64)] = False

    groups = defaultdict(list)
    for pos, label in enumerate(strata.to_numpy()):
        if available[pos]:
            groups[str(label)].append(pos)

    for _, positions in sorted(groups.items(), key=lambda item: len(item[1])):
        if len(selected) >= target:
            break
        quota = min(len(positions), max(1, min_per_stratum), target - len(selected))
        if quota <= 0:
            continue
        chosen = rng.choice(np.array(positions, dtype=np.int64), size=quota, replace=False)
        selected.update(int(i) for i in chosen)

    if len(selected) < target:
        remaining = np.setdiff1d(np.arange(n_obs, dtype=np.int64), np.array(sorted(selected), dtype=np.int64), assume_unique=True)
        fill_n = min(target - len(selected), len(remaining))
        if fill_n > 0:
            selected.update(int(i) for i in rng.choice(remaining, size=fill_n, replace=False))

    return np.array(sorted(selected), dtype=np.int64)


def _gene_index(var_names: pd.Index, gene: str) -> int | None:
    matches = np.where(var_names.astype(str).str.upper() == gene.upper())[0]
    if len(matches) == 0:
        return None
    return int(matches[0])


def _extract_column(adata: ad.AnnData, gene_idx: int) -> np.ndarray:
    x = adata[:, [gene_idx]].X
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel().astype(np.float32, copy=False)
    return np.asarray(x).ravel().astype(np.float32, copy=False)


def _top_expression_indices(values: np.ndarray, top_n: int) -> np.ndarray:
    positive = np.where(values > 0)[0]
    if positive.size == 0 or top_n <= 0:
        return np.array([], dtype=np.int64)
    n = min(int(top_n), positive.size)
    pos_values = values[positive]
    order = np.argpartition(pos_values, -n)[-n:]
    return np.array(sorted(positive[order]), dtype=np.int64)


def _select_hvg_mask(var: pd.DataFrame, n_top_genes: int) -> np.ndarray:
    if "highly_variable" in var.columns and var["highly_variable"].astype(bool).sum() > 0:
        mask = var["highly_variable"].astype(bool).to_numpy()
        if mask.sum() <= n_top_genes:
            return mask
        scores = var.get("dispersions_norm", pd.Series(np.arange(len(var)), index=var.index)).fillna(-np.inf).to_numpy()
        hvg_idx = np.where(mask)[0]
        keep = hvg_idx[np.argsort(scores[hvg_idx])[-n_top_genes:]]
        out = np.zeros(len(var), dtype=bool)
        out[keep] = True
        return out
    if "dispersions_norm" in var.columns:
        scores = var["dispersions_norm"].fillna(-np.inf).to_numpy()
        keep = np.argsort(scores)[-min(n_top_genes, len(scores)) :]
    else:
        keep = np.arange(min(n_top_genes, len(var)))
    out = np.zeros(len(var), dtype=bool)
    out[keep] = True
    return out


def compute_latent(adata: ad.AnnData, n_top_genes: int, n_pcs: int, seed: int) -> tuple[np.ndarray, int]:
    hvg_mask = _select_hvg_mask(adata.var, n_top_genes)
    x = adata[:, hvg_mask].X
    if sparse.issparse(x):
        x = x.tocsr().astype(np.float32)
    else:
        x = np.asarray(x, dtype=np.float32)
    n_components = min(max(2, n_pcs), 64, max(2, x.shape[1] - 1), max(2, x.shape[0] - 1))
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    latent = svd.fit_transform(x).astype(np.float32, copy=False)
    return latent, int(hvg_mask.sum())


def _leiden_from_rep(
    rep: np.ndarray,
    key: str,
    n_neighbors: int,
    resolution: float,
    seed: int,
) -> np.ndarray:
    import scanpy as sc

    holder = ad.AnnData(X=np.zeros((rep.shape[0], 1), dtype=np.float32))
    holder.obsm[key] = rep.astype(np.float32, copy=False)
    sc.pp.neighbors(holder, n_neighbors=min(n_neighbors, rep.shape[0] - 1), use_rep=key, method="umap")
    sc.tl.leiden(holder, resolution=resolution, flavor="igraph", directed=False, random_state=seed)
    return holder.obs["leiden"].astype(str).to_numpy()


def css_representation(
    latent: np.ndarray,
    obs: pd.DataFrame,
    batch_key: str,
    seed: int,
) -> tuple[np.ndarray, dict]:
    sample_labels = obs[batch_key].astype(str).to_numpy() if batch_key in obs.columns else np.array(["all"] * len(obs))
    sample_to_clusters: dict[str, list[int]] = defaultdict(list)
    centroid_blocks = []
    per_sample = []

    for sample in pd.unique(sample_labels):
        idx = np.where(sample_labels == sample)[0]
        if len(idx) == 0:
            continue
        n_clusters = max(2, min(12, len(idx) // 5000 if len(idx) >= 5000 else 2))
        if len(idx) < n_clusters:
            n_clusters = max(1, len(idx))
        km = MiniBatchKMeans(
            n_clusters=n_clusters,
            random_state=seed,
            batch_size=min(10000, max(1024, len(idx))),
            n_init=3,
        )
        km.fit(latent[idx])
        centroids = km.cluster_centers_.astype(np.float32, copy=False)
        start = sum(block.shape[0] for block in centroid_blocks)
        for c in range(centroids.shape[0]):
            sample_to_clusters[str(sample)].append(start + c)
        centroid_blocks.append(centroids)
        per_sample.append({"sample": str(sample), "n_cells": int(len(idx)), "n_css_clusters": int(centroids.shape[0])})

    centroids = np.vstack(centroid_blocks).astype(np.float32, copy=False)
    centroids = centroids / (np.linalg.norm(centroids, axis=1, keepdims=True) + 1e-8)

    css = np.empty((latent.shape[0], centroids.shape[0]), dtype=np.float32)
    chunk_size = 20000
    for start in range(0, latent.shape[0], chunk_size):
        stop = min(start + chunk_size, latent.shape[0])
        chunk = latent[start:stop].astype(np.float32, copy=False)
        chunk = chunk / (np.linalg.norm(chunk, axis=1, keepdims=True) + 1e-8)
        sims = (chunk @ centroids.T).astype(np.float32, copy=False)
        for cols in sample_to_clusters.values():
            block = sims[:, cols]
            if block.shape[1] == 1:
                block[:] = 0.0
            else:
                block[:] = (block - block.mean(axis=1, keepdims=True)) / (block.std(axis=1, keepdims=True) + 1e-6)
        css[start:stop] = sims

    return css, {
        "css_n_reference_clusters": int(css.shape[1]),
        "css_per_sample": per_sample,
    }


def cluster_majority_metrics(cluster_labels: np.ndarray, truth: pd.Series) -> dict:
    frame = pd.DataFrame({"cluster": cluster_labels.astype(str), "truth": truth.astype(str).to_numpy()})
    majority = frame.groupby("cluster")["truth"].agg(lambda values: values.value_counts().idxmax())
    predicted = frame["cluster"].map(majority)
    return {
        "majority_match_fraction": float((predicted.to_numpy() == frame["truth"].to_numpy()).mean()),
        "n_clusters": int(frame["cluster"].nunique()),
        "n_truth_labels": int(frame["truth"].nunique()),
    }


def rare_group_metrics(clean: np.ndarray, css: np.ndarray, obs: pd.DataFrame, cfg: AuditConfig) -> pd.DataFrame:
    rows = []
    for key in ("cell_type", "immune_subtype"):
        if key not in obs.columns:
            continue
        counts = obs[key].astype(str).value_counts()
        rare = counts[(counts >= cfg.rare_min_cells) & (counts <= max(cfg.rare_min_cells, cfg.rare_max_fraction * len(obs)))]
        for label, n_cells in rare.items():
            mask = obs[key].astype(str).to_numpy() == str(label)
            clean_sub = clean[mask]
            css_sub = css[mask]
            rows.append(
                {
                    "key": key,
                    "label": str(label),
                    "n_cells": int(n_cells),
                    "clean_n_clusters": int(pd.Series(clean_sub).nunique()),
                    "css_n_clusters": int(pd.Series(css_sub).nunique()),
                    "clean_top_cluster_fraction": float(pd.Series(clean_sub).value_counts(normalize=True).iloc[0]),
                    "css_top_cluster_fraction": float(pd.Series(css_sub).value_counts(normalize=True).iloc[0]),
                    "clean_css_ari_within_group": float(adjusted_rand_score(clean_sub, css_sub))
                    if len(np.unique(clean_sub)) > 1 or len(np.unique(css_sub)) > 1
                    else 1.0,
                    "clean_css_nmi_within_group": float(normalized_mutual_info_score(clean_sub, css_sub))
                    if len(np.unique(clean_sub)) > 1 or len(np.unique(css_sub)) > 1
                    else 1.0,
                }
            )
    return pd.DataFrame(rows)


def elf3_summary(
    expr: np.ndarray | None,
    clean: np.ndarray,
    css: np.ndarray,
    obs: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    if expr is None:
        return pd.DataFrame(), {"elf3_present": False}
    positive = expr > 0
    threshold = float(np.quantile(expr[positive], 0.95)) if positive.any() else float("nan")
    high = expr >= threshold if positive.any() else np.zeros(len(expr), dtype=bool)
    rows = []
    for key in ("cell_type", "immune_subtype", "leiden"):
        if key not in obs.columns:
            continue
        frame = pd.DataFrame({"label": obs[key].astype(str).to_numpy(), "expr": expr, "positive": positive, "high": high})
        grouped = frame.groupby("label", observed=False)
        for label, part in grouped:
            rows.append(
                {
                    "key": key,
                    "label": str(label),
                    "n_cells": int(len(part)),
                    "elf3_positive_n": int(part["positive"].sum()),
                    "elf3_positive_fraction": float(part["positive"].mean()),
                    "elf3_high_n": int(part["high"].sum()),
                    "elf3_mean": float(part["expr"].mean()),
                    "elf3_max": float(part["expr"].max()),
                }
            )
    metrics = {
        "elf3_present": True,
        "elf3_positive_n": int(positive.sum()),
        "elf3_positive_fraction": float(positive.mean()),
        "elf3_high_threshold": threshold,
        "elf3_high_n": int(high.sum()),
        "elf3_high_clean_css_ari": float(adjusted_rand_score(clean[high], css[high])) if high.sum() >= 2 else float("nan"),
        "elf3_high_clean_css_nmi": float(normalized_mutual_info_score(clean[high], css[high])) if high.sum() >= 2 else float("nan"),
        "elf3_high_clean_n_clusters": int(pd.Series(clean[high]).nunique()) if high.any() else 0,
        "elf3_high_css_n_clusters": int(pd.Series(css[high]).nunique()) if high.any() else 0,
    }
    return pd.DataFrame(rows), metrics


def expression_group_summary(
    expr: np.ndarray | None,
    obs: pd.DataFrame,
    gene: str,
    keys: Iterable[str] = ("cell_type", "immune_subtype", "leiden"),
) -> tuple[pd.DataFrame, dict]:
    if expr is None:
        return pd.DataFrame(), {f"{gene.lower()}_present_full": False}
    positive = expr > 0
    threshold = float(np.quantile(expr[positive], 0.95)) if positive.any() else float("nan")
    high = expr >= threshold if positive.any() else np.zeros(len(expr), dtype=bool)
    rows = []
    for key in keys:
        if key not in obs.columns:
            continue
        frame = pd.DataFrame({"label": obs[key].astype(str).to_numpy(), "expr": expr, "positive": positive, "high": high})
        for label, part in frame.groupby("label", observed=False):
            rows.append(
                {
                    "key": key,
                    "label": str(label),
                    "n_cells": int(len(part)),
                    f"{gene}_positive_n": int(part["positive"].sum()),
                    f"{gene}_positive_fraction": float(part["positive"].mean()),
                    f"{gene}_high_n": int(part["high"].sum()),
                    f"{gene}_mean": float(part["expr"].mean()),
                    f"{gene}_max": float(part["expr"].max()),
                }
            )
    metrics = {
        f"{gene.lower()}_present_full": True,
        f"{gene.lower()}_positive_n_full": int(positive.sum()),
        f"{gene.lower()}_positive_fraction_full": float(positive.mean()),
        f"{gene.lower()}_high_threshold_full": threshold,
        f"{gene.lower()}_high_n_full": int(high.sum()),
    }
    return pd.DataFrame(rows), metrics


def run_audit(cfg: AuditConfig) -> Path:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    mem = memory_math(cfg.input_h5ad)
    backed = ad.read_h5ad(cfg.input_h5ad, backed="r")
    try:
        obs = backed.obs.copy()
        var_names = backed.var_names.copy()
        elf3_idx = _gene_index(var_names, "ELF3")
        full_elf3 = _extract_column(backed, elf3_idx) if elf3_idx is not None else None
        full_elf3_table, full_elf3_metrics = expression_group_summary(full_elf3, obs, "ELF3")
        extra = _top_expression_indices(full_elf3, cfg.elf3_top_n) if full_elf3 is not None else np.array([], dtype=np.int64)
        chosen = stratified_indices(obs, cfg.subset_size, cfg.seed, cfg.min_per_stratum, extra)
        subset = backed[chosen, :].to_memory()
    finally:
        backed.file.close()

    subset.obs_names_make_unique()
    subset.var_names_make_unique()
    subset_obs = subset.obs.copy()
    latent, n_hvg = compute_latent(subset, cfg.n_top_genes, cfg.n_pcs, cfg.seed)

    clean_labels = _leiden_from_rep(latent, "X_clean_svd", cfg.n_neighbors, cfg.leiden_resolution, cfg.seed)
    css_rep, css_meta = css_representation(latent, subset_obs, cfg.batch_key, cfg.seed)
    css_labels = _leiden_from_rep(css_rep, "X_css", cfg.n_neighbors, cfg.leiden_resolution, cfg.seed)

    metrics = {
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "input_h5ad": str(cfg.input_h5ad),
        "subset_size_requested": int(cfg.subset_size),
        "subset_size_actual": int(subset.n_obs),
        "n_vars": int(subset.n_vars),
        "n_hvg_used": int(n_hvg),
        "seed": int(cfg.seed),
        "clean_n_clusters": int(pd.Series(clean_labels).nunique()),
        "css_n_clusters": int(pd.Series(css_labels).nunique()),
        "clean_css_ari": float(adjusted_rand_score(clean_labels, css_labels)),
        "clean_css_nmi": float(normalized_mutual_info_score(clean_labels, css_labels)),
        "peak_rss_mb": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024),
        **css_meta,
    }
    if "leiden" in subset_obs.columns:
        full_labels = subset_obs["leiden"].astype(str).to_numpy()
        metrics["clean_vs_full_massive_leiden_ari"] = float(adjusted_rand_score(clean_labels, full_labels))
        metrics["clean_vs_full_massive_leiden_nmi"] = float(normalized_mutual_info_score(clean_labels, full_labels))
        metrics["css_vs_full_massive_leiden_ari"] = float(adjusted_rand_score(css_labels, full_labels))
        metrics["css_vs_full_massive_leiden_nmi"] = float(normalized_mutual_info_score(css_labels, full_labels))

    label_rows = []
    for key in ("cell_type", "immune_subtype"):
        if key not in subset_obs.columns:
            continue
        truth = subset_obs[key].astype(str)
        clean_majority = cluster_majority_metrics(clean_labels, truth)
        css_majority = cluster_majority_metrics(css_labels, truth)
        label_rows.append({"label_key": key, "mode": "clean", **clean_majority})
        label_rows.append({"label_key": key, "mode": "css", **css_majority})
        metrics[f"clean_{key}_nmi"] = float(normalized_mutual_info_score(clean_labels, truth))
        metrics[f"css_{key}_nmi"] = float(normalized_mutual_info_score(css_labels, truth))
        metrics[f"clean_{key}_majority_match_fraction"] = clean_majority["majority_match_fraction"]
        metrics[f"css_{key}_majority_match_fraction"] = css_majority["majority_match_fraction"]

    subset_elf3_idx = _gene_index(subset.var_names, "ELF3")
    subset_elf3 = _extract_column(subset, subset_elf3_idx) if subset_elf3_idx is not None else None
    elf3_table, elf3_metrics = elf3_summary(subset_elf3, clean_labels, css_labels, subset_obs)
    metrics.update(elf3_metrics)
    metrics.update(full_elf3_metrics)
    metrics["memory_math"] = mem
    metrics["config"] = asdict(cfg)

    rare = rare_group_metrics(clean_labels, css_labels, subset_obs, cfg)
    labels = pd.DataFrame(
        {
            "cell_id": subset.obs_names.astype(str),
            "clean_leiden": clean_labels,
            "css_leiden": css_labels,
            "full_massive_leiden": subset_obs["leiden"].astype(str).to_numpy() if "leiden" in subset_obs.columns else "",
            "cell_type": subset_obs["cell_type"].astype(str).to_numpy() if "cell_type" in subset_obs.columns else "",
            "immune_subtype": subset_obs["immune_subtype"].astype(str).to_numpy() if "immune_subtype" in subset_obs.columns else "",
            "sample": subset_obs["sample"].astype(str).to_numpy() if "sample" in subset_obs.columns else "",
            "patient": subset_obs["patient"].astype(str).to_numpy() if "patient" in subset_obs.columns else "",
            "ELF3": subset_elf3 if subset_elf3 is not None else np.full(subset.n_obs, np.nan, dtype=np.float32),
        }
    )

    _write_json(cfg.output_dir / "summary.json", metrics)
    _write_json(cfg.output_dir / "memory_math.json", mem)
    pd.DataFrame([metrics]).drop(columns=["css_per_sample", "memory_math", "config"], errors="ignore").to_csv(
        cfg.output_dir / "metrics.csv", index=False
    )
    pd.DataFrame(label_rows).to_csv(cfg.output_dir / "label_majority_metrics.csv", index=False)
    rare.to_csv(cfg.output_dir / "rare_group_metrics.csv", index=False)
    if not elf3_table.empty:
        elf3_table = elf3_table.sort_values(["elf3_high_n", "elf3_positive_n", "elf3_mean"], ascending=False)
    elf3_table.to_csv(cfg.output_dir / "elf3_group_summary.csv", index=False)
    if not full_elf3_table.empty:
        full_elf3_table = full_elf3_table.sort_values(["ELF3_high_n", "ELF3_positive_n", "ELF3_mean"], ascending=False)
    full_elf3_table.to_csv(cfg.output_dir / "elf3_full_cohort_group_summary.csv", index=False)
    labels.to_csv(cfg.output_dir / "subset_cell_labels.csv", index=False)
    (cfg.output_dir / "README.md").write_text(
        "\n".join(
            [
                "# NC2024 CSS Fidelity Audit",
                "",
                f"- Input: `{cfg.input_h5ad}`",
                f"- Subset cells: `{subset.n_obs}`",
                f"- Clean-vs-CSS ARI: `{metrics['clean_css_ari']:.6f}`",
                f"- Clean-vs-CSS NMI: `{metrics['clean_css_nmi']:.6f}`",
                f"- ELF3 present: `{metrics.get('elf3_present')}`",
                f"- Full-cohort ELF3 positive cells: `{metrics.get('elf3_positive_n_full')}`",
                f"- Peak RSS MB: `{metrics['peak_rss_mb']:.1f}`",
                "",
                "`elf3_full_cohort_group_summary.csv` is the unbiased full-object ELF3 distribution.",
                "`elf3_group_summary.csv` is from the benchmark subset, which intentionally includes ELF3-high cells.",
                "See `summary.json`, `metrics.csv`, `rare_group_metrics.csv`, and ELF3 summaries.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return cfg.output_dir / "summary.json"


def parse_args(argv: list[str] | None = None) -> AuditConfig:
    parser = argparse.ArgumentParser(description="NC2024 clean-vs-CSS fidelity audit")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--subset-size", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260424)
    parser.add_argument("--n-top-genes", type=int, default=1000)
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--leiden-resolution", type=float, default=0.4)
    parser.add_argument("--batch-key", default="sample")
    parser.add_argument("--rare-min-cells", type=int, default=30)
    parser.add_argument("--rare-max-fraction", type=float, default=0.02)
    parser.add_argument("--min-per-stratum", type=int, default=100)
    parser.add_argument("--elf3-top-n", type=int, default=1000)
    args = parser.parse_args(argv)

    out = args.output_dir
    if out is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        out = DEFAULT_OUTPUT_ROOT / f"NC2024_NSCLC_CSS_FIDELITY_AUDIT_AUTO_{stamp}"

    return AuditConfig(
        input_h5ad=args.input_h5ad,
        output_dir=out,
        subset_size=args.subset_size,
        seed=args.seed,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        leiden_resolution=args.leiden_resolution,
        batch_key=args.batch_key,
        rare_min_cells=args.rare_min_cells,
        rare_max_fraction=args.rare_max_fraction,
        min_per_stratum=args.min_per_stratum,
        elf3_top_n=args.elf3_top_n,
    )


def main(argv: list[str] | None = None) -> int:
    cfg = parse_args(argv)
    if not cfg.input_h5ad.exists():
        raise FileNotFoundError(cfg.input_h5ad)
    summary = run_audit(cfg)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
