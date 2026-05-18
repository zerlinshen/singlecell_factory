#!/usr/bin/env python3
"""Audit Trevino Figure 3G/3I motif-synergy reproducibility resources.

Figures 3G and 3I depend on motif-cluster synergy scores. This audit records
which local/public resources are present and which method dependencies are
missing so the reproduction package does not fabricate synergy panels.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def file_info(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    }


def grep_count(root: Path, terms: list[str]) -> dict[str, int]:
    counts = {term: 0 for term in terms}
    if not root.exists():
        return counts
    for path in root.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        try:
            text = path.read_text(errors="ignore")
        except Exception:
            continue
        for term in terms:
            counts[term] += text.count(term)
    return counts


def run_r_probe(rscript: Path) -> dict[str, object]:
    expr = (
        "cat('chromVAR=', requireNamespace('chromVAR', quietly=TRUE), '\\n', sep='');"
        "if (requireNamespace('chromVAR', quietly=TRUE)) {"
        "cat('chromVAR_version=', as.character(packageVersion('chromVAR')), '\\n', sep='');"
        "cat('has_getAnnotationSynergy=', exists('getAnnotationSynergy', asNamespace('chromVAR')), '\\n', sep='');"
        "cat('has_getAnnotationCorrelation=', exists('getAnnotationCorrelation', asNamespace('chromVAR')), '\\n', sep='');"
        "};"
        "cat('motifmatchr=', requireNamespace('motifmatchr', quietly=TRUE), '\\n', sep='')"
    )
    if not rscript.exists():
        return {"rscript": str(rscript), "returncode": None, "stdout": "", "stderr": "Rscript not found"}
    proc = subprocess.run([str(rscript), "-e", expr], text=True, capture_output=True, check=False)
    parsed: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            key, val = line.split("=", 1)
            parsed[key] = val
    return {
        "rscript": str(rscript),
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "parsed": parsed,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--rscript", default=Path("/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript"), type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    inputs = args.project_root / "inputs"
    rds = inputs / "brainchromatin_s3/rds"
    external = inputs / "external_code"
    present_resources = {
        "tf_motif_expression_correlation": file_info(rds / "TF_MotifExpressionCorrelation.RDS"),
        "scatac_motif_match_matrix": file_info(rds / "scATAC_MotifMatchMatrix.RDS"),
        "glial_knn_atac_chromvar": file_info(rds / "Matrix_Glial_KNNpseudoBulk_ATAC_ChromVAR.RDS"),
        "glial_knn_atac_accessibility": file_info(rds / "Matrix_Glial_KNNpseudoBulk_ATAC_Accessibility.RDS"),
        "figure3f_selected_pairs": file_info(
            args.run_dir
            / "python/figures/public_matrix/figure_3F_tf_motif_heatmap/figure3f_selected_tf_motif_pairs.tsv"
        ),
    }
    search_terms = [
        "getAnnotationSynergy",
        "getAnnotationCorrelation",
        "Vierstra",
        "motif_clustering",
        "motif cluster",
        "synergy",
    ]
    code_counts = grep_count(external, search_terms)
    input_name_hits = [
        str(p.relative_to(inputs))
        for p in inputs.rglob("*")
        if p.is_file()
        and any(token in p.name.lower() for token in ["synergy", "motif_cluster", "motif-cluster", "vierstra"])
    ]
    r_probe = run_r_probe(args.rscript)
    r_parsed = r_probe.get("parsed", {})
    chromvar_available = r_parsed.get("chromVAR") == "TRUE"
    motifmatchr_available = r_parsed.get("motifmatchr") == "TRUE"
    has_synergy_function = r_parsed.get("has_getAnnotationSynergy") == "TRUE"
    has_correlation_function = r_parsed.get("has_getAnnotationCorrelation") == "TRUE"
    has_motif_cluster_resource = bool(input_name_hits) or code_counts.get("Vierstra", 0) > 0 or code_counts.get("motif_clustering", 0) > 0
    all_present = all(v["exists"] and v["size_bytes"] > 0 for v in present_resources.values())

    figure3g_status = (
        "METHOD_RESOURCE_GAP_SYNERGY_NOT_REPRODUCED"
        if not (all_present and chromvar_available and motifmatchr_available and has_synergy_function and has_motif_cluster_resource)
        else "READY_FOR_CHROMVAR_SYNERGY_RECOMPUTE"
    )
    figure3i_status = (
        "METHOD_RESOURCE_GAP_MEAN_SYNERGY_NOT_REPRODUCED"
        if figure3g_status.startswith("METHOD_RESOURCE_GAP")
        else "READY_FOR_MEAN_SYNERGY_SCATTER"
    )
    audit = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "figures": ["3G", "3I"],
        "statuses": {"3G": figure3g_status, "3I": figure3i_status},
        "truth_boundary": (
            "Figures 3G and 3I require motif-cluster synergy scores. Present public resources support motif-level "
            "correlation/proxy panels, but exact synergy reproduction requires chromVAR getAnnotationSynergy plus a "
            "24 motif-cluster mapping/synergy workflow that is not available in the local public resource set."
        ),
        "present_resources": present_resources,
        "r_environment_probe": r_probe,
        "external_code_term_counts": code_counts,
        "motif_cluster_or_synergy_named_files": input_name_hits,
        "decision_checks": {
            "all_basic_input_files_present": all_present,
            "chromvar_available": chromvar_available,
            "motifmatchr_available": motifmatchr_available,
            "has_getAnnotationSynergy": has_synergy_function,
            "has_getAnnotationCorrelation": has_correlation_function,
            "has_local_motif_cluster_resource": has_motif_cluster_resource,
        },
        "honest_decision": (
            "Do not render a synthetic synergy heatmap/scatter as Figure 3G/3I. Keep Figure 3H as correlation-only "
            "public-resource evidence and mark 3G/3I as method-resource gaps until chromVAR/motif-cluster resources "
            "and package versions are installed/pinned."
        ),
    }
    out = args.out_dir / "figure3g_3i_motif_synergy_resource_audit.json"
    out.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
