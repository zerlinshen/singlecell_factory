#!/usr/bin/env python3
"""C benchmark lane_manifest.json writer.

Walks a completed C-benchmark run directory, extracts provenance from
run_manifest.json + module_status.csv + the AnnData metadata, and emits a
lane_manifest.json conforming to G-C0 v2 contract sections A1-A8.

Usage:
    python scripts/c_benchmark_lane_manifest_writer.py \\
        --run-dir /home/zerlinshen/projects/nc-reproduction/runs/<run-id> \\
        --lane {1|2} \\
        --contract-sha b3ff5ba7110f45eacda786b11a169518c3396c9bb9a458c5c24b87324bb64240 \\
        --clustering-sha 200c31793be53eebcc9762beb67dd2c11e30ad4f625790fa37e885259deee4a2
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_zarr_dir(p: Path) -> str:
    """Recursive content hash of a zarr directory.

    Implements: find <dir> -type f | sort | xargs sha256sum | sha256sum
    """
    files = sorted([f for f in p.rglob("*") if f.is_file()])
    outer = hashlib.sha256()
    for fp in files:
        # Inner hash as "<sha>  <relpath>\n" — order-stable
        inner = sha256_file(fp)
        outer.update(f"{inner}  {fp.relative_to(p)}\n".encode())
    return outer.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="Path to completed run dir")
    ap.add_argument("--lane", required=True, choices=["1", "2"], help="Lane identifier")
    ap.add_argument(
        "--contract-sha",
        default="b3ff5ba7110f45eacda786b11a169518c3396c9bb9a458c5c24b87324bb64240",
        help="G-C0 v2 contract SHA256 (default: 2026-05-20 locked value)",
    )
    ap.add_argument(
        "--clustering-sha",
        default="200c31793be53eebcc9762beb67dd2c11e30ad4f625790fa37e885259deee4a2",
        help="clustering.py SHA256 at G-C0 lock-in (default: 2026-05-20 locked value)",
    )
    ap.add_argument(
        "--approval-timestamp",
        default="2026-05-20T00:00:00Z",
        help="G-C0 approval timestamp (default: 2026-05-20)",
    )
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.exists():
        raise SystemExit(f"run-dir does not exist: {run_dir}")

    # 1. Pipeline run_manifest.json
    run_manifest_path = run_dir / "run_manifest.json"
    if not run_manifest_path.exists():
        raise SystemExit(f"missing run_manifest.json: {run_manifest_path}")
    run_manifest = json.loads(run_manifest_path.read_text())

    # 2. module_status.csv
    module_status_path = run_dir / "module_status.csv"
    module_status_rows = []
    if module_status_path.exists():
        for line in module_status_path.read_text().splitlines()[1:]:  # skip header
            module_status_rows.append(line)

    # 3. Final adata SHA256 (best-effort — the file may live under a clustering subdir)
    final_adata_candidates = list(run_dir.rglob("final_adata.h5ad"))
    final_adata_sha = (
        sha256_file(final_adata_candidates[0]) if final_adata_candidates else None
    )

    # 4. Input zarr SHA256
    zarr_path = Path(
        "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/raw/nc2024_nsclc_emtab13526/full_cohort/prepared_input.zarr"
    )
    input_zarr_sha = sha256_zarr_dir(zarr_path) if zarr_path.exists() else None

    # 5. clustering.py SHA256 — recompute to detect post-approval drift
    clustering_py_path = Path(
        "/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/workflow/modular/modules/clustering.py"
    )
    clustering_py_sha_actual = sha256_file(clustering_py_path)
    drift_detected = clustering_py_sha_actual != args.clustering_sha

    # Extract clustering module metadata from run_manifest if present
    ctx_metadata = run_manifest.get("metadata", {}) or {}
    clustering_backend = ctx_metadata.get("clustering_backend", "unknown")

    # Compose the lane manifest
    manifest = {
        "g_c0_contract_sha256": args.contract_sha,
        "g_c0_approved_at": args.approval_timestamp,
        "clustering_py_sha256": args.clustering_sha,
        "clustering_py_sha256_actual_at_run": clustering_py_sha_actual,
        "clustering_py_drift_detected": drift_detected,
        "lane": int(args.lane),
        "lane_role": (
            "scanpy_cpu_sparse_exact_arpack" if args.lane == "1" else "rapids_singlecell_gpu_drop_in"
        ),
        # A1 — input
        "input_zarr_path": str(zarr_path),
        "input_zarr_sha256": input_zarr_sha,
        # A2 / A4b — numeric pins
        "qc_thresholds": {
            "min_genes": 200, "max_genes": 7000, "min_counts": 500, "max_counts": 50000,
            "max_mito_pct": 20.0, "max_ribo_pct": 50.0, "min_cells": 3,
        },
        "numeric_pins": {
            "n_pcs": 40, "n_neighbors": 15, "leiden_resolution": 0.8,
            "random_state": 42, "target_sum": 10000.0,
            "harmony_theta": 2.0, "harmony_sigma": 0.1, "harmony_max_iter": 10,
        },
        # A3 — HVG
        "hvg_flavor": "seurat",
        "hvg_n_top_genes": 3000,
        # A5 — Harmony
        "harmony_backend": "cpu_harmonypy",
        # A8 — ambient policy
        "ambient_correction_applied": False,
        "ambient_correction_skip_reason": "C benchmark interim policy; ambient module pending (see Plan G G-AMB-0..G-AMB-5)",
        # Provenance from run
        "run_id": run_dir.name,
        "run_dir": str(run_dir),
        "clustering_backend": clustering_backend,
        "final_adata_sha256": final_adata_sha,
        "module_status_rows": module_status_rows,
        "ctx_metadata_excerpt": {
            k: v for k, v in ctx_metadata.items()
            if k.startswith(("clustering_", "n_clusters", "gpu_", "css_", "hvg_", "harmony_", "pca_", "neighbors_", "leiden_"))
        },
        "manifest_written_at": datetime.now(timezone.utc).isoformat(),
        "manifest_writer_version": "v1_2026-05-20",
    }

    out_path = run_dir / "lane_manifest.json"
    out_path.write_text(json.dumps(manifest, indent=2, default=str))
    print(f"WROTE: {out_path}")
    print(f"  lane: {manifest['lane']} ({manifest['lane_role']})")
    print(f"  drift_detected: {drift_detected}")
    print(f"  clustering_backend: {clustering_backend}")
    print(f"  ambient_correction_applied: False")


if __name__ == "__main__":
    main()
