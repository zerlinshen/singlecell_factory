"""Python helper for nc_regression_gate.sh — compute the 4 sub-gates.

Usage:
    python scripts/ci/nc_regression_compare.py --baseline-run-dir <path> --candidate-run-dir <path>

Reads scripts/ci/nc_regression_tolerance_spec.yaml for thresholds.
Exit codes match the spec (0 pass; 10/11/12/13 per-gate fail; 2 missing input).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

SPEC_PATH = Path(__file__).resolve().parent / "nc_regression_tolerance_spec.yaml"


def _load_spec() -> dict:
    return yaml.safe_load(SPEC_PATH.read_text(encoding="utf-8"))


def _find_celltype_csv(run_dir: Path) -> Path | None:
    """Locate annotation/cell_type_annotation.csv for a given run dir."""
    candidates = list(run_dir.rglob("annotation/cell_type_annotation.csv"))
    return candidates[0] if candidates else None


def _find_de_markers_csv(run_dir: Path) -> Path | None:
    """Locate differential_expression/marker_genes.csv (top-N per cluster) for run dir."""
    candidates = list(run_dir.rglob("differential_expression/marker_genes.csv"))
    return candidates[0] if candidates else None


def gate_cluster_ari(baseline_dir: Path, candidate_dir: Path, threshold: float) -> tuple[bool, str]:
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score

    b_csv = _find_celltype_csv(baseline_dir)
    c_csv = _find_celltype_csv(candidate_dir)
    if not (b_csv and c_csv):
        return False, f"missing cell_type_annotation.csv (baseline={b_csv}, candidate={c_csv})"

    b = pd.read_csv(b_csv).set_index(pd.read_csv(b_csv).columns[0])
    c = pd.read_csv(c_csv).set_index(pd.read_csv(c_csv).columns[0])
    common = b.index.intersection(c.index)
    if len(common) < 100:
        return False, f"too few overlapping barcodes for ARI ({len(common)})"

    ari = float(adjusted_rand_score(b.loc[common, "leiden"], c.loc[common, "leiden"]))
    ok = ari >= threshold
    return ok, f"cluster_label_ari={ari:.4f} threshold={threshold}"


def gate_top_de_jaccard(baseline_dir: Path, candidate_dir: Path, top_k: int, per_cluster_min: float) -> tuple[bool, str]:
    import pandas as pd

    b_csv = _find_de_markers_csv(baseline_dir)
    c_csv = _find_de_markers_csv(candidate_dir)
    if not (b_csv and c_csv):
        return False, f"missing differential_expression/marker_genes.csv (baseline={b_csv}, candidate={c_csv})"

    b = pd.read_csv(b_csv)
    c = pd.read_csv(c_csv)
    # Expect columns: cluster, gene (or names), fdr/score — adapt to actual schema
    cluster_col = "cluster" if "cluster" in b.columns else b.columns[0]
    gene_col = "gene" if "gene" in b.columns else ("names" if "names" in b.columns else b.columns[1])
    sort_col = "fdr" if "fdr" in b.columns else ("pvals_adj" if "pvals_adj" in b.columns else None)

    per_cluster_jaccards: dict[str, float] = {}
    for cid in sorted(set(b[cluster_col]) | set(c[cluster_col])):
        b_top = b[b[cluster_col] == cid]
        c_top = c[c[cluster_col] == cid]
        if sort_col and sort_col in b_top.columns:
            b_top = b_top.sort_values(sort_col).head(top_k)
            c_top = c_top.sort_values(sort_col).head(top_k)
        else:
            b_top = b_top.head(top_k)
            c_top = c_top.head(top_k)
        b_set = set(b_top[gene_col])
        c_set = set(c_top[gene_col])
        if not (b_set | c_set):
            continue
        per_cluster_jaccards[str(cid)] = len(b_set & c_set) / len(b_set | c_set)

    if not per_cluster_jaccards:
        return False, "no overlapping clusters with DE outputs"

    failed = {cid: j for cid, j in per_cluster_jaccards.items() if j < per_cluster_min}
    ok = len(failed) == 0
    summary = f"min={min(per_cluster_jaccards.values()):.4f} median={sorted(per_cluster_jaccards.values())[len(per_cluster_jaccards)//2]:.4f} threshold={per_cluster_min}"
    if not ok:
        summary += f"; failed_clusters={list(failed.keys())[:5]}{'...' if len(failed) > 5 else ''}"
    return ok, summary


def gate_structural(baseline_dir: Path, candidate_dir: Path, fields: list[str]) -> tuple[bool, str]:
    """Check h5ad shape/dtype/obs schema parity. Uses backed-mode reads.

    Requires final_adata.h5ad in both runs (the canonical end-of-pipeline output).
    A missing or zero-byte final_adata.h5ad is itself a structural failure — the
    pipeline did not complete cleanly. No checkpoint substitution.
    """
    import anndata as ad

    b_h5ad = next(iter(baseline_dir.rglob("final_adata.h5ad")), None)
    c_h5ad = next(iter(candidate_dir.rglob("final_adata.h5ad")), None)
    if not (b_h5ad and c_h5ad and b_h5ad.stat().st_size > 0 and c_h5ad.stat().st_size > 0):
        return False, f"missing or empty final_adata.h5ad (baseline={b_h5ad}, candidate={c_h5ad})"

    b = ad.read_h5ad(b_h5ad, backed="r")
    c = ad.read_h5ad(c_h5ad, backed="r")
    diffs: list[str] = []
    if "shape" in fields and b.shape != c.shape:
        diffs.append(f"shape b={b.shape} c={c.shape}")
    if "X.dtype" in fields and b.X.dtype != c.X.dtype:
        diffs.append(f"X.dtype b={b.X.dtype} c={c.X.dtype}")
    if "obs.columns" in fields:
        b_cols = set(b.obs.columns)
        c_cols = set(c.obs.columns)
        if b_cols != c_cols:
            diffs.append(f"obs.columns symmetric_diff={sorted(b_cols ^ c_cols)}")
    if "var.columns" in fields:
        if set(b.var.columns) != set(c.var.columns):
            diffs.append("var.columns mismatch")
    ok = not diffs
    return ok, "structural=match" if ok else "; ".join(diffs)


def gate_schema_version(baseline_dir: Path, candidate_dir: Path, min_version: str) -> tuple[bool, str]:
    """Both runs' bundle_manifest.json must declare an accepted schema_version."""
    b_man = next(iter(baseline_dir.rglob("bundle_manifest.json")), None)
    c_man = next(iter(candidate_dir.rglob("bundle_manifest.json")), None)
    if not (b_man and c_man):
        # Bundle is optional — pass if both absent
        if not b_man and not c_man:
            return True, "no bundles produced; skipping schema_version check"
        return False, f"asymmetric bundle presence (baseline={bool(b_man)}, candidate={bool(c_man)})"
    b_doc = json.loads(b_man.read_text())
    c_doc = json.loads(c_man.read_text())
    accepted = {"singlecell_r_bundle_v2.1", "singlecell_r_bundle_v2.2"}
    b_ver = b_doc.get("schema_version", "")
    c_ver = c_doc.get("schema_version", "")
    ok = b_ver in accepted and c_ver in accepted
    return ok, f"baseline_schema={b_ver} candidate_schema={c_ver}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-run-dir", type=Path, required=True)
    parser.add_argument("--candidate-run-dir", type=Path, required=True)
    args = parser.parse_args()

    if not args.baseline_run_dir.exists():
        print(f"ERROR: baseline dir missing: {args.baseline_run_dir}", file=sys.stderr)
        return 2
    if not args.candidate_run_dir.exists():
        print(f"ERROR: candidate dir missing: {args.candidate_run_dir}", file=sys.stderr)
        return 2

    spec = _load_spec()

    # Gate 1: cluster ARI
    ok1, msg1 = gate_cluster_ari(
        args.baseline_run_dir, args.candidate_run_dir,
        threshold=float(spec["cluster_label_ari"]["threshold_min"]),
    )
    status1 = "PASS" if ok1 else "FAIL"
    print(f"[gate-1 cluster_label_ari] {status1}: {msg1}")
    if not ok1:
        return 10

    # Gate 2: top-50 DE marker Jaccard per cluster
    ok2, msg2 = gate_top_de_jaccard(
        args.baseline_run_dir, args.candidate_run_dir,
        top_k=int(spec["top_de_markers_jaccard"]["top_k"]),
        per_cluster_min=float(spec["top_de_markers_jaccard"]["per_cluster_threshold_min"]),
    )
    status2 = "PASS" if ok2 else "FAIL"
    print(f"[gate-2 top_de_markers_jaccard] {status2}: {msg2}")
    if not ok2:
        return 11

    # Gate 3: structural byte-identical
    ok3, msg3 = gate_structural(
        args.baseline_run_dir, args.candidate_run_dir,
        fields=list(spec["structural_byte_identical"]["fields"]),
    )
    status3 = "PASS" if ok3 else "FAIL"
    print(f"[gate-3 structural] {status3}: {msg3}")
    if not ok3:
        return 12

    # Gate 4: bundle schema version
    ok4, msg4 = gate_schema_version(
        args.baseline_run_dir, args.candidate_run_dir,
        min_version=spec["bundle_schema"]["schema_version_min"],
    )
    status4 = "PASS" if ok4 else "FAIL"
    print(f"[gate-4 bundle_schema] {status4}: {msg4}")
    if not ok4:
        return 13

    print("[nc_regression_gate] ALL 4 SUB-GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
