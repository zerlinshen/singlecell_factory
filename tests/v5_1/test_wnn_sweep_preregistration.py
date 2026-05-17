from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
WNN_SWEEP_DIR = REPO_ROOT / "ops" / "research" / "wave5" / "wnn_sweep"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_hv2_preregistration_files_parse_and_are_sha_addressable():
    search_space = WNN_SWEEP_DIR / "search_space.json"
    aggregation_rule = WNN_SWEEP_DIR / "aggregation_rule.json"

    assert search_space.exists()
    assert aggregation_rule.exists()
    assert len(_sha256(search_space)) == 64
    assert len(_sha256(aggregation_rule)) == 64

    assert _load_json(search_space)["plan_revision"] == "v5.1"
    assert _load_json(aggregation_rule)["plan_revision"] == "v5.1"


def test_hv2_search_space_total_matches_grid_product():
    spec = _load_json(WNN_SWEEP_DIR / "search_space.json")
    grid = spec["grid"]
    product = 1
    for values in grid.values():
        product *= len(values)

    assert product == spec["total_configurations"] == 288


def test_hv2_aggregation_rule_contains_required_guard_derivations():
    rule = _load_json(WNN_SWEEP_DIR / "aggregation_rule.json")
    derivation = rule["guard_threshold_derivation"]

    assert {"n_clusters_range", "silhouette_floor", "max_cluster_fraction"} <= set(derivation)
    assert derivation["n_clusters_range"]
    assert derivation["silhouette_floor"]
    assert derivation["max_cluster_fraction"]


def test_hv2_aggregation_guard_values_are_locked():
    rule = _load_json(WNN_SWEEP_DIR / "aggregation_rule.json")
    constraints = rule["selection_rule_machine_readable"]["constraints"]

    locked = {(c["field"], c["op"], c["value"]) for c in constraints}
    assert ("silhouette_wnn", ">=", 0.10) in locked
    assert ("n_clusters", ">=", 10) in locked
    assert ("n_clusters", "<=", 30) in locked
    assert ("max_cluster_fraction", "<=", 0.4) in locked

