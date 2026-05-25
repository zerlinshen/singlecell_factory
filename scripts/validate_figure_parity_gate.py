#!/usr/bin/env python3
"""Validate figure parity for current paper/process-data reproduction outputs.

The gate is manifest-driven. It always checks produced figure existence and file
headers. When a reference image exists, it compares dimensions plus exact hash
or quantitative drift metrics. Missing optional references are reported as a
conditional scientific boundary rather than silently passing.

Restored 2026-05-22 after the script was dropped in commit 04eeb39 alongside
the NC2024 abort. The NC2024-specific default-manifest construction has been
replaced with a registry-driven path: each entry in
`<suite_root>/governance/figure_parity_references_<date>.json` declares both
`produced` and `reference`, and the gate evaluates them directly without any
per-dataset Python coupling.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PNG_HEADER = b"\x89PNG\r\n\x1a\n"
PDF_HEADER = b"%PDF-"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "ops"
    / "governance_records"
    / "figure_parity_gate"
    / "latest.json"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _png_dimensions(path: Path) -> tuple[int, int]:
    data = path.read_bytes()[:24]
    if len(data) < 24 or data[:8] != PNG_HEADER or data[12:16] != b"IHDR":
        raise ValueError("invalid PNG header")
    width = int.from_bytes(data[16:20], "big")
    height = int.from_bytes(data[20:24], "big")
    return width, height


def _file_kind(path: Path) -> str:
    header = path.read_bytes()[:8]
    if path.suffix.lower() == ".png" and header == PNG_HEADER:
        return "png"
    if path.suffix.lower() == ".pdf" and header.startswith(PDF_HEADER):
        return "pdf"
    raise ValueError(f"unsupported or invalid figure header: {path}")


def _byte_mismatch_ratio(left: Path, right: Path) -> float:
    try:
        import numpy as np  # type: ignore
        a_arr = np.frombuffer(left.read_bytes(), dtype=np.uint8)
        b_arr = np.frombuffer(right.read_bytes(), dtype=np.uint8)
        max_len = max(a_arr.size, b_arr.size, 1)
        overlap = min(a_arr.size, b_arr.size)
        mismatches = int((a_arr[:overlap] != b_arr[:overlap]).sum()) + abs(a_arr.size - b_arr.size)
        return mismatches / max_len
    except Exception:
        a = left.read_bytes()
        b = right.read_bytes()
        max_len = max(len(a), len(b), 1)
        overlap = min(len(a), len(b))
        mismatches = sum(1 for i in range(overlap) if a[i] != b[i]) + abs(len(a) - len(b))
        return mismatches / max_len


def _pillow_rms(left: Path, right: Path) -> float | None:
    try:
        from PIL import Image, ImageChops  # type: ignore
    except Exception:
        return None
    try:
        with Image.open(left) as a_img, Image.open(right) as b_img:
            a = a_img.convert("RGB")
            b = b_img.convert("RGB")
            if a.size != b.size:
                return math.inf
            diff = ImageChops.difference(a, b)
            hist = diff.histogram()
            sq = sum(value * ((idx % 256) ** 2) for idx, value in enumerate(hist))
            rms = math.sqrt(sq / float(a.size[0] * a.size[1] * 3))
            return rms / 255.0
    except Exception:
        return None


def _resolve(path_value: str | None, *, base: Path) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value)
    return path if path.is_absolute() else (base / path)


def _evaluate_item(
    item: dict[str, Any],
    *,
    base: Path,
    allow_missing_produced: bool = False,
) -> dict[str, Any]:
    figure_id = str(item.get("id", "unnamed"))
    produced = _resolve(item.get("produced"), base=base)
    reference = _resolve(item.get("reference"), base=base)
    reference_required = bool(item.get("reference_required", False))
    max_rms = float(item.get("max_rms_normalized", 0.03))
    max_byte_ratio = float(item.get("max_byte_mismatch_ratio", 0.05))
    result: dict[str, Any] = {
        "id": figure_id,
        "produced": str(produced) if produced else None,
        "reference": str(reference) if reference else None,
    }
    if produced is None:
        result.update({"status": "conditional", "reason": "produced_path_unregistered"})
        return result
    if not produced.exists():
        status = "conditional" if allow_missing_produced else "fail"
        result.update({"status": status, "reason": "missing_produced_figure"})
        return result
    try:
        kind = _file_kind(produced)
    except ValueError as exc:
        result.update({"status": "fail", "reason": str(exc)})
        return result
    result["kind"] = kind
    result["produced_bytes"] = produced.stat().st_size
    result["produced_sha256"] = _sha256(produced)
    if kind == "png":
        result["produced_dimensions"] = list(_png_dimensions(produced))

    if reference is None or not reference.exists():
        status = "fail" if reference_required else "conditional"
        result.update({"status": status, "reason": "reference_missing"})
        return result
    try:
        ref_kind = _file_kind(reference)
    except ValueError as exc:
        result.update({"status": "fail", "reason": f"invalid_reference: {exc}"})
        return result
    if ref_kind != kind:
        result.update({"status": "fail", "reason": "kind_mismatch", "reference_kind": ref_kind})
        return result
    result["reference_bytes"] = reference.stat().st_size
    result["reference_sha256"] = _sha256(reference)
    if result["produced_sha256"] == result["reference_sha256"]:
        result.update({"status": "pass", "metric": "exact_sha256", "value": 0.0})
        return result
    if kind == "png":
        ref_dims = _png_dimensions(reference)
        result["reference_dimensions"] = list(ref_dims)
        if tuple(result["produced_dimensions"]) != ref_dims:
            result.update({"status": "fail", "reason": "dimension_mismatch"})
            return result
        rms = _pillow_rms(produced, reference)
        if rms is not None:
            result.update({"metric": "rms_normalized", "value": rms, "threshold": max_rms})
            result["status"] = "pass" if rms <= max_rms else "fail"
            if result["status"] == "fail":
                result["reason"] = "rms_above_threshold"
            return result
    ratio = _byte_mismatch_ratio(produced, reference)
    result.update({"metric": "byte_mismatch_ratio", "value": ratio, "threshold": max_byte_ratio})
    result["status"] = "pass" if ratio <= max_byte_ratio else "fail"
    if result["status"] == "fail":
        result["reason"] = "byte_mismatch_above_threshold"
    return result


class _RegistryError(SystemExit):
    pass


def _load_registry() -> tuple[Path, dict[str, Any]]:
    """Load the latest figure parity registry from the suite governance dir."""
    suite_root = Path(__file__).resolve().parents[2]
    governance = suite_root / "governance"
    if not governance.is_dir():
        raise _RegistryError(f"suite governance dir not found: {governance}")
    candidates = sorted(governance.glob("figure_parity_references_*.json"))
    if not candidates:
        raise _RegistryError(f"no figure_parity_references_*.json under {governance}")
    registry_path = candidates[-1]
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise _RegistryError(f"figure parity registry unreadable: {registry_path}: {exc}")
    except json.JSONDecodeError as exc:
        raise _RegistryError(f"figure parity registry not valid JSON: {registry_path}: {exc}")
    return registry_path, registry


def _registry_to_manifest(registry: dict[str, Any]) -> dict[str, Any]:
    """Convert a registry JSON to a gate-compatible manifest.

    Each registry `entries[*]` becomes one figure entry. The registry is the
    single source of truth for both `produced` and `reference` paths.
    """
    allowed_root = Path("/home/zerlinshen/data/external/published_figures").resolve()
    figures: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    defaults = registry.get("default_thresholds") or {}
    for item in registry.get("entries", []):
        item_id = item.get("id")
        if not item_id:
            raise _RegistryError("registry entry missing `id`")
        if item_id in seen_ids:
            raise _RegistryError(f"duplicate id {item_id!r} in registry")
        seen_ids.add(item_id)
        ref = item.get("reference")
        if ref:
            ref_resolved = Path(ref).resolve()
            try:
                ref_resolved.relative_to(allowed_root)
            except ValueError:
                raise _RegistryError(
                    f"reference {ref!r} for {item_id!r} must resolve under {allowed_root}"
                )
        figure: dict[str, Any] = {
            "id": item_id,
            "produced": item.get("produced"),
            "reference": ref,
            "reference_required": bool(item.get("reference_required", False)),
            "max_rms_normalized": float(
                item.get("max_rms_normalized", defaults.get("max_rms_normalized", 0.03))
            ),
            "max_byte_mismatch_ratio": float(
                item.get(
                    "max_byte_mismatch_ratio",
                    defaults.get("max_byte_mismatch_ratio", 0.05),
                )
            ),
        }
        figures.append(figure)
    return {
        "schema_version": "figure_parity_gate_v1",
        "scientific_boundary": registry.get("description"),
        "figures": figures,
    }


def validate_manifest(
    manifest: dict[str, Any],
    *,
    base: Path,
    allow_missing_produced: bool = False,
) -> dict[str, Any]:
    items = [
        _evaluate_item(item, base=base, allow_missing_produced=allow_missing_produced)
        for item in manifest.get("figures", [])
    ]
    counts: dict[str, int] = {}
    for item in items:
        counts[item["status"]] = counts.get(item["status"], 0) + 1
    if counts.get("fail", 0):
        verdict = "fail"
    elif counts.get("conditional", 0):
        verdict = "conditional"
    else:
        verdict = "pass"
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": manifest.get("schema_version", "figure_parity_gate_v1"),
        "verdict": verdict,
        "status_counts": counts,
        "scientific_boundary": manifest.get("scientific_boundary"),
        "figures": items,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--allow-conditional", action="store_true")
    parser.add_argument(
        "--allow-missing-produced",
        action="store_true",
        help="Downgrade missing produced figures to conditional for data-free CI.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.manifest:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        base = args.manifest.parent
    else:
        _, registry = _load_registry()
        manifest = _registry_to_manifest(registry)
        base = Path.cwd()
    report = validate_manifest(
        manifest,
        base=base,
        allow_missing_produced=args.allow_missing_produced,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["verdict"] == "pass" or (
        args.allow_conditional and report["verdict"] == "conditional"
    ):
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
