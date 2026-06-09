#!/usr/bin/env python3
"""Dataset-agnostic real-data parity gate for the three-factory suite.

This gate asserts that a recorded source-of-truth (SOT) bounded real-data
validation run actually exists on disk as a complete, real evidence chain:

  1. the run manifest parses and carries the fields this gate reads;
  2. the run directory exists;
  3. the bundle directory (``outputs.bundle_dir`` under the run) exists;
  4. the recorded number of PNGs (``outputs.png_count``) and PDFs
     (``outputs.pdf_count``) exist under the run's R output dir
     (``outputs.r_dir``) as REAL, non-blank files (PNG/PDF magic bytes +
     nonzero size; if Pillow is importable, PNG min dimensions are also
     asserted against the manifest's recorded threshold, otherwise a byte-size
     heuristic is used);
  5. the bundle manifest (``bundle_manifest.json`` under the bundle dir)
     reports the expected exported cell count, an ``X_umap`` and ``X_pca``
     embedding, and the expected marker-gene count.

It is deliberately DATASET-AGNOSTIC: every expected quantity (cell count,
marker count, embedding names, image dimensions, output dir names) is read from
the manifest itself, NOT hard-coded to any particular dataset / paper / cohort.
Point ``--run-manifest`` at any manifest with the same shape and the gate
re-derives its expectations from that manifest.

This is REAL enforcement: if the SOT run is absent or incomplete the gate exits
nonzero with a clear message (no silent skip). ``--allow-missing-run`` is
provided ONLY for data-free CI / manual structural probes; it downgrades a
missing run directory to a structural-only pass (manifest parse + shape) and
logs the downgrade loudly. The default is strict.

A JSON report is written to ``--output`` (default under /tmp) so the gate never
dirties a working tree.

stdlib only (Pillow used opportunistically if present).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# --- path safety (reused allowlist + containment pattern) -------------------
# Run-id / output-dir come from the command line and are used to build
# filesystem paths. Validate them against a strict character allowlist +
# containment check before use so a typo or a pathological argument cannot
# escape the project tree. (Pattern carried over from the retired
# validate_two_realdata_final.py path-safety guard.)
_SAFE_RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9._:T/+-]{1,256}$")


def _is_relative_to(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def _safe_roots() -> tuple[Path, ...]:
    suite_root = Path(__file__).resolve().parents[2]
    return (
        Path("/home/zerlinshen/projects").resolve(),
        suite_root,
        Path("/tmp").resolve(),
    )


def _safe_path(value: str | Path, *, kind: str) -> Path:
    """Resolve ``value`` and require it under projects/, suite root, or /tmp."""
    path = Path(value).expanduser().resolve(strict=False)
    if not any(_is_relative_to(path, root) for root in _safe_roots()):
        raise argparse.ArgumentTypeError(
            f"{kind} {path} must live under projects/, suite root, or /tmp "
            f"(allowed roots: {[str(r) for r in _safe_roots()]})"
        )
    return path


def _safe_run_manifest(value: str) -> Path:
    return _safe_path(value, kind="run manifest")


def _safe_output_dir(value: str) -> Path:
    return _safe_path(value, kind="output path")


def _contained_child(base: Path, rel: str, *, label: str) -> Path:
    """Resolve ``base / rel`` and require the result to stay under ``base``.

    ``rel`` is read from the manifest (``outputs.*``); guard against a manifest
    that smuggles ``..`` or an absolute path that would escape the run dir.
    """
    if not _SAFE_RUN_ID_PATTERN.match(rel):
        raise SystemExit(
            f"{label} value {rel!r} is not in the allowed character class "
            f"[A-Za-z0-9._:T/+-]"
        )
    child = (base / rel).resolve(strict=False)
    if not _is_relative_to(child, base.resolve()):
        raise SystemExit(f"{label} {child} escapes the run directory {base}")
    return child


# --- magic-byte / image checks ----------------------------------------------
_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
_PDF_MAGIC = b"%PDF-"
# Byte-size floor for a "real" non-blank raster when Pillow is unavailable: a
# blank/placeholder PNG at publication dimensions is far smaller than this.
_PNG_MIN_BYTES = 10_000
_PDF_MIN_BYTES = 1_000


def _load_pillow():
    try:
        from PIL import Image  # noqa: WPS433 (runtime optional dep)

        return Image
    except ImportError:
        return None


def _parse_min_dims(verification: Any) -> tuple[int, int] | None:
    """Best-effort extract a recorded ``>=WxH`` PNG threshold from verification.

    The manifest records human-readable verification strings such as
    "PIL QA found 7/7 PNGs nonblank with dimensions >= 3200x1280". We parse the
    first ``>= WxH`` we find so the dimension assertion tracks the manifest's
    OWN recorded threshold rather than a hard-coded constant. Returns None when
    no threshold is recorded (then only the byte-size heuristic applies).
    """
    if isinstance(verification, list):
        text = "\n".join(str(item) for item in verification)
    else:
        text = str(verification or "")
    match = re.search(r">=\s*(\d{2,5})\s*[xX]\s*(\d{2,5})", text)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None


def _check_png(path: Path, min_dims: tuple[int, int] | None, image_mod) -> list[str]:
    errors: list[str] = []
    if not path.exists() or not path.is_file():
        return [f"missing PNG: {path}"]
    size = path.stat().st_size
    if size <= 0:
        return [f"empty PNG: {path}"]
    if path.read_bytes()[:8] != _PNG_MAGIC:
        errors.append(f"invalid PNG magic bytes: {path}")
    if image_mod is not None:
        try:
            with image_mod.open(path) as img:
                width, height = img.size
        except Exception as exc:  # pragma: no cover - corrupt file path
            return errors + [f"unreadable PNG ({exc}): {path}"]
        if min_dims is not None:
            min_w, min_h = min_dims
            if width < min_w or height < min_h:
                errors.append(
                    f"PNG below recorded min dimensions {min_w}x{min_h}: "
                    f"{path} is {width}x{height}"
                )
    elif size < _PNG_MIN_BYTES:
        errors.append(
            f"PNG suspiciously small ({size} bytes < {_PNG_MIN_BYTES} byte "
            f"non-blank floor; Pillow unavailable for a dimension check): {path}"
        )
    return errors


def _check_pdf(path: Path) -> list[str]:
    if not path.exists() or not path.is_file():
        return [f"missing PDF: {path}"]
    size = path.stat().st_size
    if size <= 0:
        return [f"empty PDF: {path}"]
    errors: list[str] = []
    if not path.read_bytes()[:8].startswith(_PDF_MAGIC):
        errors.append(f"invalid PDF magic bytes: {path}")
    if size < _PDF_MIN_BYTES:
        errors.append(
            f"PDF suspiciously small ({size} bytes < {_PDF_MIN_BYTES} byte floor): "
            f"{path}"
        )
    return errors


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


# --- bundle manifest assertions ---------------------------------------------

def _bundle_cell_count(manifest: dict[str, Any]) -> int | None:
    """Pull the exported cell count from a bundle manifest, schema-tolerantly."""
    bundle = manifest.get("bundle", {})
    for key in ("n_cells_exported", "n_cells", "cells_exported"):
        if isinstance(bundle.get(key), int):
            return int(bundle[key])
    return None


def _bundle_marker_count(manifest: dict[str, Any]) -> int | None:
    bundle = manifest.get("bundle", {})
    for key in ("markers_present", "markers_requested"):
        value = bundle.get(key)
        if isinstance(value, list):
            return len(value)
    return None


def _bundle_has_embeddings(manifest: dict[str, Any]) -> tuple[bool, bool]:
    """Return (has_X_umap, has_X_pca) across the bundle manifest's records."""
    files = manifest.get("files", {})
    keys = set(files) if isinstance(files, dict) else set()
    # Also accept embeddings recorded under an extensions/obsm listing.
    for extra_key in ("obsm", "embeddings", "extensions"):
        extra = manifest.get(extra_key)
        if isinstance(extra, dict):
            keys |= set(extra)
        elif isinstance(extra, list):
            keys |= {str(item) for item in extra}
    return ("X_umap" in keys, "X_pca" in keys)


def validate_run(
    manifest_path: Path,
    *,
    allow_missing_run: bool,
) -> dict[str, Any]:
    """Validate the real-data evidence chain; return a structured report.

    Raises SystemExit with a clear message on any hard failure (strict mode).
    """
    report: dict[str, Any] = {
        "gate": "realdata_parity",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_manifest": str(manifest_path),
        "strict": not allow_missing_run,
        "checks": [],
        "errors": [],
    }

    def record(name: str, ok: bool, detail: str = "") -> None:
        report["checks"].append({"check": name, "ok": ok, "detail": detail})
        if not ok:
            report["errors"].append(f"{name}: {detail}")

    # 1) manifest parses
    if not manifest_path.exists():
        if allow_missing_run:
            report["structural_only"] = True
            report["verdict"] = "structural_only_pass"
            print(
                "WARNING: --allow-missing-run: run manifest absent "
                f"({manifest_path}); downgrading to STRUCTURAL-ONLY pass (NO "
                "real-data evidence verified).",
                file=sys.stderr,
            )
            # Informational only — do NOT push to errors (a downgrade is a pass).
            report["checks"].append(
                {
                    "check": "manifest_parses",
                    "ok": False,
                    "detail": f"absent (downgraded, allow-missing-run): {manifest_path}",
                }
            )
            return report
        raise SystemExit(f"run manifest does not exist: {manifest_path}")
    try:
        manifest = _load_json(manifest_path)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"run manifest is not valid JSON: {manifest_path}: {exc}")
    record("manifest_parses", True, str(manifest_path))

    inputs = manifest.get("inputs", {})
    outputs = manifest.get("outputs", {})
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise SystemExit(
            f"run manifest missing inputs/outputs objects: {manifest_path}"
        )

    expected_cells = inputs.get("exported_cells")
    expected_markers = inputs.get("markers_requested")
    expected_marker_count = (
        len(expected_markers) if isinstance(expected_markers, list) else None
    )
    report["expected"] = {
        "exported_cells": expected_cells,
        "marker_count": expected_marker_count,
        "export_seed": inputs.get("export_seed"),
        "png_count": outputs.get("png_count"),
        "pdf_count": outputs.get("pdf_count"),
    }

    # 2) run directory exists (the manifest's own parent is the run dir)
    run_dir = manifest_path.parent
    run_dir_exists = run_dir.is_dir()
    if not run_dir_exists:
        if allow_missing_run:
            report["structural_only"] = True
            print(
                "WARNING: --allow-missing-run: run directory absent "
                f"({run_dir}); downgrading to STRUCTURAL-ONLY pass "
                "(manifest parse + shape only, NO real-data evidence verified).",
                file=sys.stderr,
            )
            # Informational only — do NOT push to errors (a downgrade is a pass).
            report["checks"].append(
                {
                    "check": "run_dir_exists",
                    "ok": False,
                    "detail": f"absent (downgraded, allow-missing-run): {run_dir}",
                }
            )
            report["verdict"] = "structural_only_pass"
            return report
        raise SystemExit(
            f"SOT run directory is ABSENT: {run_dir} (real-data enforcement; "
            "pass --allow-missing-run ONLY for data-free CI/manual probes)"
        )
    record("run_dir_exists", True, str(run_dir))

    image_mod = _load_pillow()
    report["pillow_available"] = image_mod is not None
    min_dims = _parse_min_dims(manifest.get("verification"))
    report["png_min_dims_from_manifest"] = list(min_dims) if min_dims else None

    # 3) bundle directory exists
    bundle_rel = outputs.get("bundle_dir")
    if not bundle_rel:
        raise SystemExit("run manifest outputs.bundle_dir is missing")
    bundle_dir = _contained_child(run_dir, str(bundle_rel), label="outputs.bundle_dir")
    record("bundle_dir_exists", bundle_dir.is_dir(), str(bundle_dir))

    # 4) png_count PNGs + pdf_count PDFs under r_dir are REAL non-blank files
    r_rel = outputs.get("r_dir")
    if not r_rel:
        raise SystemExit("run manifest outputs.r_dir is missing")
    r_dir = _contained_child(run_dir, str(r_rel), label="outputs.r_dir")
    record("r_dir_exists", r_dir.is_dir(), str(r_dir))

    png_count = int(outputs.get("png_count", 0) or 0)
    pdf_count = int(outputs.get("pdf_count", 0) or 0)

    pngs = sorted(r_dir.glob("*.png")) if r_dir.is_dir() else []
    pdfs = sorted(r_dir.glob("*.pdf")) if r_dir.is_dir() else []

    # The manifest records the EXPECTED count; require at least that many real
    # files (a renderer producing an extra composite PNG must not fail the gate,
    # but a missing one must).
    record(
        "png_count_satisfied",
        len(pngs) >= png_count,
        f"found {len(pngs)} PNG(s) under {r_dir}, manifest png_count={png_count}",
    )
    record(
        "pdf_count_satisfied",
        len(pdfs) >= pdf_count,
        f"found {len(pdfs)} PDF(s) under {r_dir}, manifest pdf_count={pdf_count}",
    )

    real_pngs = 0
    for png in pngs:
        errs = _check_png(png, min_dims, image_mod)
        if errs:
            for err in errs:
                record("png_real_nonblank", False, err)
        else:
            real_pngs += 1
    record(
        "png_real_nonblank_count",
        real_pngs >= png_count,
        f"{real_pngs} real non-blank PNG(s) >= png_count {png_count}",
    )

    real_pdfs = 0
    for pdf in pdfs:
        errs = _check_pdf(pdf)
        if errs:
            for err in errs:
                record("pdf_real_nonblank", False, err)
        else:
            real_pdfs += 1
    record(
        "pdf_real_nonblank_count",
        real_pdfs >= pdf_count,
        f"{real_pdfs} real non-blank PDF(s) >= pdf_count {pdf_count}",
    )

    # 5) bundle_manifest.json asserts cells + X_umap + X_pca + marker count
    bundle_manifest_path = bundle_dir / "bundle_manifest.json"
    if not bundle_manifest_path.exists():
        record(
            "bundle_manifest_exists", False, f"missing: {bundle_manifest_path}"
        )
    else:
        record("bundle_manifest_exists", True, str(bundle_manifest_path))
        try:
            bundle_manifest = _load_json(bundle_manifest_path)
        except json.JSONDecodeError as exc:
            record("bundle_manifest_parses", False, str(exc))
            bundle_manifest = {}
        else:
            record("bundle_manifest_parses", True, "")
            report["bundle_schema_version"] = bundle_manifest.get("schema_version")

            actual_cells = _bundle_cell_count(bundle_manifest)
            record(
                "bundle_cell_count",
                expected_cells is not None and actual_cells == expected_cells,
                f"bundle reports {actual_cells} cells, manifest expects "
                f"{expected_cells}",
            )

            has_umap, has_pca = _bundle_has_embeddings(bundle_manifest)
            record("bundle_has_X_umap", has_umap, "X_umap present in bundle manifest")
            record("bundle_has_X_pca", has_pca, "X_pca present in bundle manifest")

            actual_markers = _bundle_marker_count(bundle_manifest)
            record(
                "bundle_marker_count",
                expected_marker_count is not None
                and actual_markers == expected_marker_count,
                f"bundle reports {actual_markers} markers, manifest expects "
                f"{expected_marker_count}",
            )

    report["verdict"] = "pass" if not report["errors"] else "fail"
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-manifest",
        type=_safe_run_manifest,
        default=_safe_run_manifest(
            "/home/zerlinshen/projects/round9-singlecell-comparison/runs/"
            "2026-05-24T2347Z-0321773/manifest.json"
        ),
        help="canonical SOT run manifest (default: Round9 bounded real-data run)",
    )
    parser.add_argument(
        "--allow-missing-run",
        action="store_true",
        help="DATA-FREE probe only: downgrade an absent run dir to a "
        "structural-only pass (logged). Default is strict real enforcement.",
    )
    parser.add_argument(
        "--output",
        type=_safe_output_dir,
        default=None,
        help="path for the JSON report (default: a /tmp file). Must live under "
        "projects/, suite root, or /tmp.",
    )
    args = parser.parse_args()

    output_path = args.output
    if output_path is None:
        output_path = _safe_output_dir(
            f"/tmp/realdata_parity_gate.{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"
        )

    try:
        report = validate_run(
            args.run_manifest, allow_missing_run=args.allow_missing_run
        )
    except SystemExit as exc:
        # Persist a minimal failure report before propagating the nonzero exit.
        fail_report = {
            "gate": "realdata_parity",
            "run_manifest": str(args.run_manifest),
            "verdict": "fail",
            "fatal": str(exc),
        }
        try:
            output_path.write_text(json.dumps(fail_report, indent=2), encoding="utf-8")
        except OSError:
            pass
        print(f"FAIL realdata_parity_gate: {exc}", file=sys.stderr)
        return 1

    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    if report["errors"]:
        for error in report["errors"]:
            print(f"ERROR: {error}", file=sys.stderr)
        print(
            f"FAIL realdata_parity_gate verdict={report['verdict']} "
            f"errors={len(report['errors'])} report={output_path}",
            file=sys.stderr,
        )
        return 1

    print(
        f"OK realdata_parity_gate verdict={report['verdict']} "
        f"checks={len(report['checks'])} report={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
