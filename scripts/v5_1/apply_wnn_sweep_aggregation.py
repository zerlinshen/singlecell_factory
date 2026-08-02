#!/usr/bin/env python3
"""Apply the SHA-pinned Wave 5 WNN sweep aggregation rule.

The executable is version controlled because tests, reruns, and fresh
worktrees must not depend on the gitignored ``.omc`` research scratch tree.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def passes(row: dict, constraints: list[dict]) -> tuple[bool, list[str]]:
    failed = []
    ops = {
        ">=": lambda a, b: a >= b,
        "<=": lambda a, b: a <= b,
        ">": lambda a, b: a > b,
        "<": lambda a, b: a < b,
        "==": lambda a, b: a == b,
    }
    for constraint in constraints:
        field = constraint["field"]
        value = row.get(field)
        if value is None or (isinstance(value, float) and value != value):
            failed.append(f"{field}=NA")
            continue
        operator = constraint["op"]
        if operator not in ops:
            raise ValueError(f"Unsupported aggregation operator: {operator}")
        threshold = constraint["value"]
        if not ops[operator](value, threshold):
            failed.append(f"{field}={value} !{operator} {threshold}")
    return not failed, failed


def _load_ok_rows(results_path: Path, search_sha: str, agg_sha: str) -> list[dict]:
    rows = []
    with results_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get("status") != "OK":
                continue
            if row.get("search_space_sha") != search_sha:
                raise ValueError(
                    "SHA mismatch: search_space "
                    f"({row.get('search_space_sha')} != {search_sha})"
                )
            if row.get("aggregation_rule_sha") != agg_sha:
                raise ValueError(
                    "SHA mismatch: aggregation_rule "
                    f"({row.get('aggregation_rule_sha')} != {agg_sha})"
                )
            rows.append(row)
    return rows


def _write_incomplete_report(
    out_dir: Path,
    rows: list[dict],
    constraints: list[dict],
    search_sha: str,
    agg_sha: str,
) -> None:
    best_violating = sorted(
        rows,
        key=lambda row: -row.get("ari_vs_seurat_clusters", -1),
    )[:10]
    report = out_dir / "hv2_incomplete_no_viable_fix.md"
    report_lines = [
        "# HV2 INCOMPLETE-NO-VIABLE-FIX",
        "",
        "Zero configurations satisfied the SHA-pinned aggregation guard.",
        "",
        f"- search_space_sha: `{search_sha}`",
        f"- aggregation_rule_sha: `{agg_sha}`",
        f"- n_results_total: {len(rows)}",
        "",
        "## Top-10 by ARI (each violates at least one guard)",
        "",
        "| config_id | ari | violations |",
        "|---|---|---|",
    ]
    for row in best_violating:
        _, failed = passes(row, constraints)
        report_lines.append(
            f"| {row['config_id']} | {row['ari_vs_seurat_clusters']:.4f} "
            f"| {', '.join(failed)} |"
        )
    report_lines.extend(
        [
            "",
            "The aggregation guard must not be silently relaxed.",
            "",
        ]
    )
    report.write_text("\n".join(report_lines), encoding="utf-8")

    marker = {
        "wnn_fix_status": "INCOMPLETE-NO-VIABLE-FIX",
        "search_space_sha": search_sha,
        "aggregation_rule_sha": agg_sha,
        "n_results_total": len(rows),
        "n_eligible": 0,
    }
    (out_dir / "selected_config.json").write_text(
        json.dumps(marker, indent=2),
        encoding="utf-8",
    )
    print(f"INCOMPLETE-NO-VIABLE-FIX written to {report}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--search-space", required=True, type=Path)
    parser.add_argument("--aggregation", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    rule = json.loads(args.aggregation.read_text(encoding="utf-8"))
    selection_rule = rule["selection_rule_machine_readable"]
    constraints = selection_rule["constraints"]
    objective = selection_rule["objective"]
    if objective != {"maximize": "ari_vs_seurat_clusters"}:
        raise ValueError(
            "Aggregation objective must maximize ari_vs_seurat_clusters"
        )

    search_sha = sha256_of_file(args.search_space)
    agg_sha = sha256_of_file(args.aggregation)
    rows = _load_ok_rows(args.results, search_sha, agg_sha)
    if not rows:
        raise SystemExit("No OK rows in sweep results")

    eligible = [row for row in rows if passes(row, constraints)[0]]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if not eligible:
        _write_incomplete_report(
            args.out_dir,
            rows,
            constraints,
            search_sha,
            agg_sha,
        )
        return

    def sort_key(row: dict) -> tuple:
        return (
            -row["ari_vs_seurat_clusters"],
            row["resolution"],
            row["n_neighbors"],
            row["wnn_weight_strategy"],
            row["knn_pruning"],
        )

    winner = sorted(eligible, key=sort_key)[0]
    selected = {
        "wnn_fix_status": "APPROVED",
        "resolution": winner["resolution"],
        "n_neighbors": winner["n_neighbors"],
        "wnn_weight_strategy": winner["wnn_weight_strategy"],
        "knn_pruning": winner["knn_pruning"],
        "ari_vs_seurat_clusters": winner["ari_vs_seurat_clusters"],
        "nmi_vs_seurat_clusters": winner["nmi_vs_seurat_clusters"],
        "silhouette_wnn": winner["silhouette_wnn"],
        "n_clusters": winner["n_clusters"],
        "max_cluster_fraction": winner["max_cluster_fraction"],
        "search_space_sha": search_sha,
        "aggregation_rule_sha": agg_sha,
        "n_results_total": len(rows),
        "n_eligible": len(eligible),
    }
    (args.out_dir / "selected_config.json").write_text(
        json.dumps(selected, indent=2),
        encoding="utf-8",
    )
    print(
        f"APPROVED: {winner['config_id']} "
        f"ARI={winner['ari_vs_seurat_clusters']:.4f}"
    )


if __name__ == "__main__":
    main(sys.argv[1:])
