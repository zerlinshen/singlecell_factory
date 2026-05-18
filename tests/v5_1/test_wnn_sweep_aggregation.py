from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
APPLY_AGGREGATION = REPO_ROOT / ".omc" / "research" / "wave5" / "wnn_sweep_draft" / "apply_aggregation.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("wave5_apply_aggregation", APPLY_AGGREGATION)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _base_search_space(tmp_path: Path) -> Path:
    path = tmp_path / "search_space.json"
    _write_json(path, {"grid": {}, "total_configurations": 0})
    return path


def _base_rule(tmp_path: Path) -> Path:
    path = tmp_path / "aggregation_rule.json"
    _write_json(
        path,
        {
            "selection_rule_machine_readable": {
                "objective": {"maximize": "ari_vs_seurat_clusters"},
                "constraints": [
                    {"field": "silhouette_wnn", "op": ">=", "value": 0.10},
                    {"field": "n_clusters", "op": ">=", "value": 10},
                    {"field": "n_clusters", "op": "<=", "value": 30},
                    {"field": "max_cluster_fraction", "op": "<=", "value": 0.4},
                ],
            }
        },
    )
    return path


def _row(search_sha: str, agg_sha: str, **overrides) -> dict:
    row = {
        "status": "OK",
        "config_id": "baseline",
        "resolution": 0.5,
        "n_neighbors": 20,
        "wnn_weight_strategy": "default",
        "knn_pruning": "none",
        "ari_vs_seurat_clusters": 0.4,
        "nmi_vs_seurat_clusters": 0.5,
        "silhouette_wnn": 0.2,
        "n_clusters": 14,
        "max_cluster_fraction": 0.25,
        "search_space_sha": search_sha,
        "aggregation_rule_sha": agg_sha,
    }
    row.update(overrides)
    return row


def test_passes_reports_failed_constraints():
    module = _load_module()
    constraints = [
        {"field": "silhouette_wnn", "op": ">=", "value": 0.10},
        {"field": "n_clusters", "op": "<=", "value": 30},
    ]

    ok, failed = module.passes({"silhouette_wnn": 0.05, "n_clusters": 31}, constraints)

    assert not ok
    assert "silhouette_wnn=0.05 !>= 0.1" in failed
    assert "n_clusters=31 !<= 30" in failed


def test_main_selects_highest_ari_eligible_config_with_sha_chain(tmp_path: Path):
    module = _load_module()
    search = _base_search_space(tmp_path)
    rule = _base_rule(tmp_path)
    search_sha = _sha(search)
    agg_sha = _sha(rule)
    results = tmp_path / "results.jsonl"
    rows = [
        _row(search_sha, agg_sha, config_id="eligible-low", ari_vs_seurat_clusters=0.40),
        _row(search_sha, agg_sha, config_id="ineligible-high", ari_vs_seurat_clusters=0.90, n_clusters=31),
        _row(search_sha, agg_sha, config_id="eligible-high", ari_vs_seurat_clusters=0.55, resolution=0.8),
    ]
    results.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    rc = module.main([
        "--search-space", str(search),
        "--aggregation", str(rule),
        "--results", str(results),
        "--out-dir", str(out_dir),
    ])

    assert rc is None
    selected = json.loads((out_dir / "selected_config.json").read_text(encoding="utf-8"))
    assert selected["wnn_fix_status"] == "APPROVED"
    assert selected["ari_vs_seurat_clusters"] == 0.55
    assert selected["resolution"] == 0.8
    assert selected["search_space_sha"] == search_sha
    assert selected["aggregation_rule_sha"] == agg_sha
    assert selected["n_results_total"] == 3
    assert selected["n_eligible"] == 2


def test_main_writes_incomplete_marker_when_no_config_passes(tmp_path: Path):
    module = _load_module()
    search = _base_search_space(tmp_path)
    rule = _base_rule(tmp_path)
    search_sha = _sha(search)
    agg_sha = _sha(rule)
    results = tmp_path / "results.jsonl"
    rows = [
        _row(search_sha, agg_sha, config_id="too-many-clusters", n_clusters=31),
        _row(search_sha, agg_sha, config_id="bad-silhouette", silhouette_wnn=0.01),
    ]
    results.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    rc = module.main([
        "--search-space", str(search),
        "--aggregation", str(rule),
        "--results", str(results),
        "--out-dir", str(out_dir),
    ])

    assert rc is None
    selected = json.loads((out_dir / "selected_config.json").read_text(encoding="utf-8"))
    assert selected["wnn_fix_status"] == "INCOMPLETE-NO-VIABLE-FIX"
    assert selected["n_results_total"] == 2
    assert selected["n_eligible"] == 0
    report = (out_dir / "hv2_incomplete_no_viable_fix.md").read_text(encoding="utf-8")
    assert "Zero configurations satisfied" in report
    assert search_sha in report
    assert agg_sha in report
