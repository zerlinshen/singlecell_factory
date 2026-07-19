#!/usr/bin/env python3
"""Produce immutable render-only artifacts for the public Figure 3A velocity lane."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.modular.velocity_render_contract import (  # noqa: E402
    choose_cells,
    write_render_bundle,
)

DEFAULT_MARKERS = (
    "SOX2",
    "PAX6",
    "EOMES",
    "NEUROD1",
    "NEUROD6",
    "SLC17A7",
    "SATB2",
    "GAD2",
)
TRUTH_BOUNDARY = (
    "Public processed-matrix RNA velocity reproduction from GSE162170; "
    "not FASTQ-level reproduction and not confirmatory biological evidence."
)


def _assert_output_outside_factory(output_dir: Path) -> None:
    factory_root = Path(__file__).resolve().parents[1]
    resolved_output = output_dir.expanduser().resolve()
    if resolved_output == factory_root or factory_root in resolved_output.parents:
        raise ValueError(
            "scientific render artifacts must be written under a governed project run, "
            "not inside the singlecell_factory tree"
        )


def _sparse_colsum(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).reshape(-1).astype(float, copy=False)
    return np.asarray(matrix).sum(axis=0).astype(float, copy=False)


def _software_versions(scv: Any) -> dict[str, str | None]:
    def package_version(name: str) -> str | None:
        try:
            return version(name)
        except PackageNotFoundError:
            return None

    return {
        "python": platform.python_version(),
        "anndata": package_version("anndata"),
        "numpy": package_version("numpy"),
        "pandas": package_version("pandas"),
        "scanpy": package_version("scanpy"),
        "scvelo": getattr(scv, "__version__", None) or package_version("scvelo"),
    }


def produce(
    velocity_h5ad: Path,
    pipeline_h5ad: Path,
    output_dir: Path,
    *,
    max_cells: int = 8000,
    top_genes: int = 2500,
    seed: int = 13,
    min_shared_counts: int = 20,
    n_pcs: int = 20,
    n_neighbors: int = 20,
    n_pseudotime_bins: int = 20,
    markers: Sequence[str] = DEFAULT_MARKERS,
    scv: Any | None = None,
) -> dict[str, Any]:
    """Run all Figure 3A analysis upstream and write strict render artifacts."""
    if scv is None:
        import scvelo as scv

    np.random.seed(seed)
    velocity = ad.read_h5ad(velocity_h5ad)
    pipeline = ad.read_h5ad(pipeline_h5ad)
    if not velocity.obs_names.is_unique or not pipeline.obs_names.is_unique:
        raise ValueError("both velocity and pipeline AnnData must have unique cell identifiers")

    common = sorted(
        set(velocity.obs_names.astype(str)).intersection(pipeline.obs_names.astype(str))
    )
    if not common:
        raise ValueError("no overlapping cells between velocity and pipeline AnnData")
    velocity = velocity[common, :].copy()
    pipeline = pipeline[common, :].copy()
    if "cell_type" not in pipeline.obs.columns:
        raise KeyError("pipeline AnnData missing required obs['cell_type']")
    if "X_umap" not in pipeline.obsm:
        raise KeyError("pipeline AnnData missing required obsm['X_umap']")
    velocity.obs["cell_type"] = pipeline.obs["cell_type"].astype(str).to_numpy()
    velocity.obsm["X_umap"] = np.asarray(pipeline.obsm["X_umap"], dtype=float)

    selected_cells = choose_cells(velocity.obs, max_cells=max_cells, seed=seed)
    velocity = velocity[selected_cells, :].copy()
    score = _sparse_colsum(velocity.X)
    for layer in ("spliced", "unspliced"):
        if layer not in velocity.layers:
            raise KeyError(f"velocity input missing required layer {layer!r}")
        score += _sparse_colsum(velocity.layers[layer])

    top_n = min(int(top_genes), velocity.n_vars)
    if top_n <= 0:
        raise ValueError("top_genes must be positive")
    top = np.argsort(score, kind="stable")[-top_n:]
    marker_indices = [
        velocity.var_names.get_loc(marker) for marker in markers if marker in velocity.var_names
    ]
    selected_genes = np.sort(
        np.unique(np.concatenate([top, np.asarray(marker_indices, dtype=int)]))
    )
    velocity = velocity[:, selected_genes].copy()
    retain_genes = [marker for marker in markers if marker in velocity.var_names]

    scv.pp.filter_and_normalize(
        velocity, min_shared_counts=min_shared_counts, retain_genes=retain_genes
    )
    scv.pp.moments(velocity, n_pcs=n_pcs, n_neighbors=n_neighbors)
    scv.tl.velocity(velocity, mode="stochastic")
    scv.tl.velocity_graph(velocity)
    scv.tl.velocity_embedding(velocity, basis="umap")
    scv.tl.velocity_confidence(velocity)
    scv.tl.velocity_pseudotime(velocity)

    if "velocity" not in velocity.layers:
        raise ValueError("scVelo did not produce layers['velocity']")
    velocity_layer = velocity.layers["velocity"]
    if sparse.issparse(velocity_layer):
        velocity.obs["velocity_length"] = np.sqrt(
            np.asarray(velocity_layer.multiply(velocity_layer).sum(axis=1)).reshape(-1)
        )
    else:
        velocity.obs["velocity_length"] = np.linalg.norm(np.asarray(velocity_layer), axis=1)

    parameters = {
        "max_cells": int(max_cells),
        "top_genes": int(top_genes),
        "min_shared_counts": int(min_shared_counts),
        "moments_n_pcs": int(n_pcs),
        "moments_n_neighbors": int(n_neighbors),
        "n_pseudotime_bins": int(n_pseudotime_bins),
        "velocity_mode": "stochastic",
        "basis": "pipeline_X_umap_aligned_by_cell_id",
        "input_overlap_cells": len(common),
        "computed_cells": int(velocity.n_obs),
        "computed_genes": int(velocity.n_vars),
    }
    return write_render_bundle(
        velocity,
        output_dir,
        markers=markers,
        source_files={"velocity_h5ad": velocity_h5ad, "pipeline_h5ad": pipeline_h5ad},
        parameters=parameters,
        seed=seed,
        software_versions=_software_versions(scv),
        producer_path=Path(__file__),
        truth_boundary=TRUTH_BOUNDARY,
        n_pseudotime_bins=n_pseudotime_bins,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--velocity-h5ad", required=True, type=Path)
    parser.add_argument("--pipeline-h5ad", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--max-cells", type=int, default=8000)
    parser.add_argument("--top-genes", type=int, default=2500)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--min-shared-counts", type=int, default=20)
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--n-pseudotime-bins", type=int, default=20)
    parser.add_argument("--marker", action="append", dest="markers")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    _assert_output_outside_factory(args.out_dir)
    manifest = produce(
        args.velocity_h5ad,
        args.pipeline_h5ad,
        args.out_dir,
        max_cells=args.max_cells,
        top_genes=args.top_genes,
        seed=args.seed,
        min_shared_counts=args.min_shared_counts,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        n_pseudotime_bins=args.n_pseudotime_bins,
        markers=tuple(args.markers) if args.markers else DEFAULT_MARKERS,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
