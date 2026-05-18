#!/usr/bin/env python3
"""Validate NC2024 architecture contracts without rerunning the full cohort.

This is a controller-validation smoke: it checks current v2 run artifacts,
compact R bundle manifests, run ledger presence, and the singlecell -> multiomics
bridge layout. It deliberately avoids opening the 10s-of-GB H5AD objects.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path("/home/zerlinshen/singlecell_factory")
MULTIOMICS_ROOT = Path("/home/zerlinshen/r_multiomics_factory")
DEFAULT_RUNS = (
    "results/nc2024_tumor_20260426_v2",
    "results/nc2024_bh_20260426_v2",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(path: Path, label: str) -> Path:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    return path


def _load_json(path: Path) -> dict:
    return json.loads(_require(path, "json").read_text(encoding="utf-8"))


def _validate_bundle(run_dir: Path, verify_sha: bool) -> dict[str, object]:
    bundle_dir = _require(run_dir / "r_bundle", "r_bundle")
    manifest = _load_json(bundle_dir / "bundle_manifest.json")
    if manifest.get("schema_version") != "singlecell_r_bundle_v2":
        raise SystemExit(f"{bundle_dir}: expected singlecell_r_bundle_v2")

    source = manifest.get("source", {})
    bundle = manifest.get("bundle", {})
    files = manifest.get("files", {})
    exported = int(bundle.get("n_cells_exported", -1))
    source_n = int(source.get("n_obs", -2))
    if exported <= 0 or exported != source_n:
        raise SystemExit(f"{bundle_dir}: exported/source n_obs mismatch: {exported} vs {source_n}")

    required = set(bundle.get("required_files", []))
    required.update({"obs", "marker_expr", "X_umap", "X_pca", "uns_summary"})
    missing = sorted(required - set(files))
    if missing:
        raise SystemExit(f"{bundle_dir}: manifest missing file records: {missing}")

    for stem in sorted(required):
        rec = files[stem]
        path = _require(bundle_dir / rec["path"], f"bundle file {stem}")
        actual_bytes = path.stat().st_size
        if actual_bytes != int(rec["bytes"]):
            raise SystemExit(
                f"{bundle_dir}: byte mismatch for {stem}: "
                f"manifest={rec['bytes']} actual={actual_bytes}"
            )
        if verify_sha and rec.get("sha256") and _sha256(path) != str(rec["sha256"]).lower():
            raise SystemExit(f"{bundle_dir}: sha256 mismatch for {stem}")

    expression = manifest.get("expression", {})
    expected_guard = "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    if expression.get("claim_guard") != expected_guard:
        raise SystemExit(f"{bundle_dir}: unexpected claim_guard={expression.get('claim_guard')}")

    return {
        "bundle_dir": str(bundle_dir),
        "n_cells_exported": exported,
        "marker_format": bundle.get("marker_format"),
        "files": sorted(files),
    }


def _validate_run(run_dir: Path, verify_sha: bool) -> dict[str, object]:
    _require(run_dir, "NC2024 v2 run dir")
    final_h5ad = _require(run_dir / "final_adata.h5ad", "final_adata.h5ad")
    bundle_summary = _validate_bundle(run_dir, verify_sha)
    ledger_hits = sorted((REPO_ROOT / "ops" / "run_ledger").glob(f"{run_dir.name}_*.json"))
    if not ledger_hits:
        raise SystemExit(f"{run_dir}: no matching ops/run_ledger entry")
    return {
        "run_dir": str(run_dir),
        "final_h5ad_bytes": final_h5ad.stat().st_size,
        "latest_ledger": str(ledger_hits[-1]),
        "bundle": bundle_summary,
    }


def validate(verify_sha: bool) -> dict[str, object]:
    bridge = REPO_ROOT / "bridges" / "local_r_pipeline_macbook"
    for link_name in ("R", "R_bundle"):
        link = _require(bridge / link_name, f"bridge symlink {link_name}")
        if not link.is_symlink():
            raise SystemExit(f"{link} must remain a symlink into r_multiomics_factory")

    _require(MULTIOMICS_ROOT / "R_bundle" / "io_bundle.R", "multiomics bundle reader")
    _require(MULTIOMICS_ROOT / "scripts" / "plot_remote_bundle_large.R", "large bundle plotter")

    runs = [_validate_run(REPO_ROOT / rel, verify_sha) for rel in DEFAULT_RUNS]
    return {
        "status": "ok",
        "execution_mode": "controller_validation",
        "checked_runs": runs,
        "bridge": {
            "R": str((bridge / "R").resolve()),
            "R_bundle": str((bridge / "R_bundle").resolve()),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-sha", action="store_true", help="Hash every bundle file.")
    args = parser.parse_args()
    print(json.dumps(validate(verify_sha=args.verify_sha), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
