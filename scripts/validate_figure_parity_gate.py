#!/usr/bin/env python3
"""Validate registered figure artifacts or audit curated reference parity.

The enforced suite mode is ``artifact-integrity``: produced figures must be
registered, present, and valid PNG/PDF artifacts. It makes no reference-parity
claim. ``reference-parity`` is a separate manual/project audit; when zero
references are staged its aggregate verdict is ``not_run``, never pass.

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
import zlib
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
    data = path.read_bytes()
    if data[:8] != PNG_HEADER:
        raise ValueError("invalid PNG signature")
    offset = 8
    dimensions: tuple[int, int] | None = None
    first_chunk = True
    saw_image_data = False
    while offset < len(data):
        if offset + 12 > len(data):
            raise ValueError("truncated PNG chunk")
        length = int.from_bytes(data[offset : offset + 4], "big")
        chunk_type = data[offset + 4 : offset + 8]
        chunk_end = offset + 12 + length
        if chunk_end > len(data):
            raise ValueError("truncated PNG chunk payload")
        chunk_data = data[offset + 8 : offset + 8 + length]
        expected_crc = int.from_bytes(data[offset + 8 + length : chunk_end], "big")
        actual_crc = zlib.crc32(chunk_type)
        actual_crc = zlib.crc32(chunk_data, actual_crc) & 0xFFFFFFFF
        if actual_crc != expected_crc:
            raise ValueError("invalid PNG chunk CRC")
        if first_chunk:
            if chunk_type != b"IHDR" or length != 13:
                raise ValueError("PNG must start with a 13-byte IHDR chunk")
            width = int.from_bytes(chunk_data[0:4], "big")
            height = int.from_bytes(chunk_data[4:8], "big")
            if width <= 0 or height <= 0:
                raise ValueError("invalid PNG dimensions")
            dimensions = (width, height)
            first_chunk = False
        elif chunk_type == b"IHDR":
            raise ValueError("PNG contains multiple IHDR chunks")
        if chunk_type == b"IDAT":
            saw_image_data = True
        if chunk_type == b"IEND":
            if length != 0 or chunk_end != len(data):
                raise ValueError("invalid PNG IEND chunk")
            if dimensions is None:
                raise ValueError("PNG is missing IHDR")
            if not saw_image_data:
                raise ValueError("PNG is missing IDAT")
            return dimensions
        offset = chunk_end
    raise ValueError("PNG is missing IEND")


def _validate_pdf(path: Path) -> None:
    with path.open("rb") as handle:
        header = handle.read(8)
        if not header.startswith(PDF_HEADER):
            raise ValueError("invalid PDF header")
        handle.seek(0, 2)
        size = handle.tell()
        handle.seek(max(0, size - 2048))
        trailer = handle.read()
    if b"%%EOF" not in trailer:
        raise ValueError("PDF is missing EOF marker")


def _file_kind(path: Path) -> str:
    header = path.read_bytes()[:8]
    if path.suffix.lower() == ".png" and header == PNG_HEADER:
        _png_dimensions(path)
        return "png"
    if path.suffix.lower() == ".pdf" and header.startswith(PDF_HEADER):
        _validate_pdf(path)
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


def _evaluate_produced(
    item: dict[str, Any],
    *,
    base: Path,
    allow_missing_produced: bool = False,
    require_registration: bool = False,
) -> dict[str, Any]:
    """Validate only the registered produced artifact, without parity semantics."""
    figure_id = str(item.get("id", "unnamed"))
    produced = _resolve(item.get("produced"), base=base)
    result: dict[str, Any] = {
        "id": figure_id,
        "produced": str(produced) if produced else None,
    }
    if produced is None:
        result.update({
            "status": "fail" if require_registration else "conditional",
            "reason": "produced_path_unregistered",
        })
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
    result.update({
        "kind": kind,
        "produced_bytes": produced.stat().st_size,
        "produced_sha256": _sha256(produced),
    })
    if kind == "png":
        result["produced_dimensions"] = list(_png_dimensions(produced))
    return result


def _evaluate_item(
    item: dict[str, Any],
    *,
    base: Path,
    allow_missing_produced: bool = False,
) -> dict[str, Any]:
    result = _evaluate_produced(
        item,
        base=base,
        allow_missing_produced=allow_missing_produced,
    )
    if result.get("status") in {"fail", "conditional"}:
        result["reference"] = str(_resolve(item.get("reference"), base=base)) \
            if item.get("reference") else None
        return result

    produced = Path(str(result["produced"]))
    reference = _resolve(item.get("reference"), base=base)
    reference_required = bool(item.get("reference_required", False))
    max_rms = float(item.get("max_rms_normalized", 0.03))
    max_byte_ratio = float(item.get("max_byte_mismatch_ratio", 0.05))
    result["reference"] = str(reference) if reference else None
    kind = str(result["kind"])

    if reference is None or not reference.exists():
        result["reference_present"] = False
        status = "fail" if reference_required else "conditional"
        result.update({"status": status, "reason": "reference_missing"})
        return result
    result["reference_present"] = True
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
        historical_ref = item.get("reference")
        lifecycle = str(
            item.get(
                "reference_audit_lifecycle",
                registry.get("reference_audit_lifecycle", "active"),
            )
        )
        reference_is_active = lifecycle not in {"historical", "historical_not_staged", "retired"}
        ref = historical_ref if reference_is_active else None
        if historical_ref:
            ref_resolved = Path(historical_ref).resolve()
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
    declared = sum(bool(item.get("reference")) for item in items)
    staged = sum(bool(item.get("reference_present")) for item in items)
    evaluated = sum("metric" in item for item in items)
    if counts.get("fail", 0):
        verdict = "fail"
    elif declared == 0 or staged == 0:
        verdict = "not_run"
    elif counts.get("conditional", 0):
        verdict = "partial"
    else:
        verdict = "pass"
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": manifest.get("schema_version", "figure_parity_gate_v1"),
        "verdict": verdict,
        "status_counts": counts,
        "scientific_boundary": manifest.get("scientific_boundary"),
        "reference_coverage": {
            "declared": declared,
            "staged": staged,
            "evaluated": evaluated,
        },
        "figures": items,
    }


def validate_artifact_integrity(
    manifest: dict[str, Any],
    *,
    base: Path,
    allow_missing_produced: bool = False,
) -> dict[str, Any]:
    """Hard validation of produced artifacts with no paper-parity implication."""
    items = [
        _evaluate_produced(
            item,
            base=base,
            allow_missing_produced=allow_missing_produced,
            require_registration=True,
        )
        for item in manifest.get("figures", [])
    ]
    counts: dict[str, int] = {}
    for item in items:
        status = item.get("status", "pass")
        item["status"] = status
        counts[status] = counts.get(status, 0) + 1
    if not items or counts.get("fail", 0):
        verdict = "fail"
    elif counts.get("conditional", 0):
        verdict = "conditional"
    else:
        verdict = "pass"
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": "figure_artifact_integrity_v1",
        "verdict": verdict,
        "status_counts": counts,
        "registered_artifacts": len(items),
        "claim_boundary": (
            "Produced-artifact existence, format structure, dimensions, and observed "
            "size/hash recording only; "
            "no original-paper visual or scientific parity claim."
        ),
        "figures": items,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--mode",
        choices=("artifact-integrity", "reference-parity"),
        default="reference-parity",
    )
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
    if args.mode == "artifact-integrity":
        if args.allow_conditional:
            raise SystemExit("--allow-conditional is valid only with --mode reference-parity")
        report = validate_artifact_integrity(
            manifest,
            base=base,
            allow_missing_produced=args.allow_missing_produced,
        )
    else:
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
    allowed_nonpass = (
        args.mode == "reference-parity"
        and args.allow_conditional
        and report["verdict"] in {"not_run", "partial"}
    ) or (
        args.mode == "artifact-integrity"
        and args.allow_missing_produced
        and report["verdict"] == "conditional"
    )
    if report["verdict"] == "pass" or allowed_nonpass:
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
