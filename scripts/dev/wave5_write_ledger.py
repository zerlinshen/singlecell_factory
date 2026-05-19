"""Wave-5 ledger writer (US-W5-10 / AC-LEDGER-1..3).

Assembles ops/run_ledger/wave5_trevino_<run_id>.json from:
  - <run-dir>/run_manifest.json           (driver-produced)
  - <run-dir>/metrics.json (or .omc/research/wave5/metrics.json)
  - environments/sc_gpu_stable.lock.txt   (env hash)
  - data/raw/trevino_2021_brain/SHA256SUMS (GEO sha256s)
  - data/external/trevino_2021_supp/SHA256SUMS (supplement sha256s)
  - module __references__ dicts at workflow/modular/modules/*.py

Validates the ledger against ops/run_ledger/schema/wave5.schema.json.

Wave-5 spec: .omc/specs/deep-interview-wave5-completion.md
Wave-5 plan: .omc/plans/wave5-completion-consensus-2026-05-16.md  (AC-LEDGER-1..3)
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory")
MODULES_DIR = REPO / "workflow" / "modular" / "modules"
SCHEMA_PATH = REPO / "ops" / "run_ledger" / "schema" / "wave5.schema.json"
LEDGER_DIR = REPO / "ops" / "run_ledger"


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def _load_sha_file(path: Path) -> dict[str, str]:
    """Parse 'sha256  filename' lines (gnu coreutils format)."""
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            sha, fname = parts
            out[fname.strip()] = sha
    return out


def _collect_module_dois() -> dict[str, list[str]]:
    """Import each module under its proper package path and harvest the
    top-level ``__references__`` dict for DOIs. Plain
    ``spec_from_file_location`` fails because the modules use relative
    imports (e.g. ``from .._contract_violation import ...``)."""
    import importlib
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    result: dict[str, list[str]] = {}
    for py in sorted(MODULES_DIR.glob("*.py")):
        if py.name.startswith("_"):
            continue
        qualname = f"workflow.modular.modules.{py.stem}"
        try:
            mod = importlib.import_module(qualname)
        except Exception:
            continue
        refs = getattr(mod, "__references__", None)
        if not isinstance(refs, dict):
            continue
        dois = sorted({str(v.get("doi")) for v in refs.values() if isinstance(v, dict) and v.get("doi")})
        if dois:
            rel = py.relative_to(REPO).as_posix()
            result[rel] = dois
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--regression-gate-exit-code", required=True, type=int)
    ap.add_argument("--metrics", default=None, type=Path,
                    help="Override metrics.json path (defaults to <run-dir>/metrics.json or .omc/research/wave5/metrics.json).")
    args = ap.parse_args()

    manifest_p = args.run_dir / "run_manifest.json"
    if not manifest_p.exists():
        print(f"FAIL: {manifest_p} missing", file=sys.stderr)
        return 2
    manifest = json.loads(manifest_p.read_text())

    metrics_p = args.metrics or (args.run_dir / "metrics.json")
    if not metrics_p.exists():
        alt = REPO / ".omc" / "research" / "wave5" / "metrics.json"
        if alt.exists():
            metrics_p = alt
        else:
            print(f"FAIL: metrics.json missing at {metrics_p}", file=sys.stderr)
            return 2
    metrics = json.loads(metrics_p.read_text())

    env_lock = REPO / "environments" / "sc_gpu_stable.lock.txt"
    env_hash = _sha256(env_lock) if env_lock.exists() else "missing"

    geo_sha = _load_sha_file(REPO / "data" / "raw" / "trevino_2021_brain" / "SHA256SUMS")
    supp_sha = _load_sha_file(REPO / "data" / "external" / "trevino_2021_supp" / "SHA256SUMS")
    geo_sha.update({f"trevino_2021_supp/{k}": v for k, v in supp_sha.items()})

    dois = _collect_module_dois()

    def metric_block(key: str, default_threshold: float) -> dict:
        meta = metrics.get(f"{key}_meta", {})
        value = metrics.get(key)
        threshold = meta.get("threshold", default_threshold)
        if value is None:
            status = "NOT-RUN"
        elif value >= threshold:
            status = "CLOSED"
        else:
            status = "CLOSED-PARTIAL"
        # Allow upstream meta to OVERRIDE only with a stricter status
        meta_status = meta.get("status")
        if meta_status in ("CLOSED-PARTIAL", "NOT-RUN"):
            status = meta_status
        return {"value": value, "threshold": threshold, "status": status}

    ledger = {
        "wave": "5",
        "run_id": manifest["run_id"],
        "project_root": manifest["project_root"],
        "run_dir": manifest["run_dir"],
        "run_manifest_path": str(manifest_p),
        "module_status_path": str(args.run_dir / "module_status.csv"),
        "conda_env_hash": env_hash,
        "geo_sha256s": geo_sha,
        "metrics": {
            "cell_type_ari":             metric_block("cell_type_ari", 0.70),
            "peak_gene_overlap_top1000": metric_block("peak_gene_overlap_top1000", 0.70),
            "pseudotime_spearman":       metric_block("pseudotime_spearman", 0.70),
        },
        "regression_gate_exit_code": args.regression_gate_exit_code,
        "dois_by_module": dois,
        "acquisition_status": {
            "trevino_S2F":   metrics.get("peak_gene_overlap_meta", {}).get("acquisition_status", "unknown"),
            "trevino_Fig4D": metrics.get("pseudotime_spearman_meta", {}).get("acquisition_status", "failed"),
        },
        "ari_evidence_path": metrics.get("ari_evidence_path", "SECONDARY"),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    if metrics.get("pseudotime_n_joined") is not None:
        ledger["metrics"]["pseudotime_n_joined"] = int(metrics["pseudotime_n_joined"])
    if metrics.get("pseudotime_status"):
        ledger["metrics"]["pseudotime_status"] = metrics["pseudotime_status"]

    # Validate against schema
    schema = json.loads(SCHEMA_PATH.read_text())
    try:
        import jsonschema
        jsonschema.Draft202012Validator(schema).validate(ledger)
    except Exception as exc:
        print(f"FAIL: ledger fails schema validation: {exc}", file=sys.stderr)
        return 2

    LEDGER_DIR.mkdir(parents=True, exist_ok=True)
    out_path = LEDGER_DIR / f"wave5_trevino_{manifest['run_id']}.json"
    out_path.write_text(json.dumps(ledger, indent=2, default=str), encoding="utf-8")
    print(f"VALID  wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
