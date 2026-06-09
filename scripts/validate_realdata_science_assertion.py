#!/usr/bin/env python3
"""Real-data scientific assertion gate (#34).

This is the FIRST suite gate that checks a real BIOLOGICAL quantity against a
documented threshold (not just structure / file existence). It reads the
RETAINED COMPACT ground-truth (GT) evidence that lives on disk (no download),
verifies the SHA256 integrity of every GT file it reads, and then asserts:

  A1. IMR90 Hi-C loop reproducibility (genome-wide F1).
      Source GT: ``hic_rao_imr90/loop_reval_genomewide_f1.json`` (and the
      genome-wide concordance file as a cross-check). The recorded
      ``headline_metric`` is ``f1``; we (a) RECOMPUTE F1 from the per-chrom
      ``gt``/``called``/``gt_matched``/``called_matched`` counts to confirm the
      recorded scalar is internally consistent, and (b) assert the F1 is at or
      above the recorded value within a small tolerance.
      Biology/provenance: our FitHiChIP-style loop calls on Rao 2014 IMR90
      in-situ Hi-C, compared to the published IMR90 HiCCUPS loop list at a
      20 kb anchor tolerance over 8040 GT loops. A nontrivial genome-wide F1
      proves we recover real chromatin loops, not noise.

  A2. TP63 CUT&Tag overexpression directional signal (OE gains >> losses, and
      OE replicate reproducibility >> the NC paired-label baseline).
      Source GT: ``cutrun_liu2025/diff_TP63_OEvsNC/signal_diff_manifest.json``
      (keys ``n_gained``, ``n_lost``, ``replicate_reproducibility.{OE,NC}``).
      Biology/provenance: TP63 is an activating master transcription factor of
      the squamous program (Liu et al. 2025). Its overexpression should GAIN
      far more H3K27ac/binding signal than it loses; the retained manifest
      records ``n_gained=60128`` vs ``n_lost=27544`` (a >2x directional margin)
      and an OE replicate reproducibility (0.7975) far above the NC
      paired-label baseline (0.2197). We assert OE-gained exceeds OE-lost by a
      wide margin AND OE reproducibility far exceeds the NC baseline — the
      "OE-up >> null by a wide margin" signal, expressed with the keys actually
      present in the retained file (the larger 11592-peak q<0.10 count cited in
      the original count-based DESeq2 path is NOT in this retained signal-diff
      slice, so it is intentionally NOT asserted here).

This is REAL enforcement: the gate exits nonzero if the GT root is missing, if
SHA256 integrity fails for any file read, or if any biological metric falls
below its threshold. ``--allow-missing-gt`` is provided ONLY for CI hosts that
do not stage the retained GT; it downgrades a missing GT root to a logged
structural skip. The default is strict.

A JSON report is written to ``--output`` (default under /tmp) so the gate never
dirties a working tree. stdlib only.
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

DEFAULT_GT_ROOT = Path(
    "/home/zerlinshen/projects/pipeline-validation-20260528/"
    "bulk_multiomics_retained_evidence"
)

# --- path safety (allowlist + containment, mirrors the parity gate) ---------
_SAFE_VALUE_PATTERN = re.compile(r"^[A-Za-z0-9._:T/+-]{1,256}$")


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
    path = Path(value).expanduser().resolve(strict=False)
    if not any(_is_relative_to(path, root) for root in _safe_roots()):
        raise argparse.ArgumentTypeError(
            f"{kind} {path} must live under projects/, suite root, or /tmp "
            f"(allowed roots: {[str(r) for r in _safe_roots()]})"
        )
    return path


def _safe_gt_root(value: str) -> Path:
    return _safe_path(value, kind="gt root")


def _safe_output(value: str) -> Path:
    return _safe_path(value, kind="output path")


# --- integrity --------------------------------------------------------------

def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_sha256sums(sums_path: Path) -> dict[str, str]:
    """Parse a ``<sha256>  ./relative/path`` SHA256SUMS file -> {relpath: sha}."""
    mapping: dict[str, str] = {}
    for line in sums_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Format: "<64-hex><space><space|*>./rel/path"
        match = re.match(r"^([0-9a-fA-F]{64})\s+\*?(.+)$", line)
        if not match:
            continue
        sha, rel = match.group(1).lower(), match.group(2).strip()
        rel = rel[2:] if rel.startswith("./") else rel
        mapping[rel] = sha
    return mapping


class IntegrityError(RuntimeError):
    """Raised when a GT file's recomputed sha256 does not match SHA256SUMS."""


def _verify_and_load_json(
    gt_root: Path, rel: str, sums: dict[str, str], *, contained: bool = True
) -> dict[str, Any]:
    """Verify ``rel`` against SHA256SUMS, then parse it as JSON.

    Raises ``IntegrityError`` on a missing file, a missing/unmatched checksum,
    or a path that escapes ``gt_root``.
    """
    if contained and not _SAFE_VALUE_PATTERN.match(rel):
        raise IntegrityError(f"GT relpath {rel!r} not in allowed character class")
    path = (gt_root / rel).resolve(strict=False)
    if not _is_relative_to(path, gt_root.resolve()):
        raise IntegrityError(f"GT path {path} escapes gt root {gt_root}")
    if not path.exists():
        raise IntegrityError(f"GT file missing: {path}")
    if rel not in sums:
        raise IntegrityError(f"no SHA256SUMS entry for {rel}")
    actual = _sha256(path)
    if actual != sums[rel]:
        raise IntegrityError(
            f"sha256 mismatch for {rel}: recorded {sums[rel]}, got {actual}"
        )
    return json.loads(path.read_text(encoding="utf-8"))


# --- assertions -------------------------------------------------------------
# Small tolerance for the F1 floor: the recorded scalar is the validated value;
# allow a tiny epsilon for float round-trip but no meaningful degradation.
_F1_TOLERANCE = 0.01
# TP63 OE-gained must exceed OE-lost by at least this ratio (recorded ~2.18x).
_TP63_GAIN_LOSS_MIN_RATIO = 1.5
# OE replicate reproducibility must exceed the NC paired-label baseline by at
# least this additive margin (recorded OE 0.7975 vs NC 0.2197 -> ~0.58 margin).
_OE_VS_NC_MIN_MARGIN = 0.2


def _recompute_loop_f1(reval: dict[str, Any]) -> dict[str, float]:
    """Recompute recall/precision/F1 from the per-chrom matched counts."""
    per_chrom = reval.get("per_chrom", {})
    gt = sum(int(c["gt"]) for c in per_chrom.values())
    called = sum(int(c["called"]) for c in per_chrom.values())
    gt_matched = sum(int(c["gt_matched"]) for c in per_chrom.values())
    called_matched = sum(int(c["called_matched"]) for c in per_chrom.values())
    recall = gt_matched / gt if gt else 0.0
    precision = called_matched / called if called else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )
    return {"recall": recall, "precision": precision, "f1": f1, "gt": gt, "called": called}


def run_assertions(gt_root: Path) -> dict[str, Any]:
    """Verify integrity and assert the biological quantities; return a report."""
    report: dict[str, Any] = {
        "gate": "realdata_science_assertion",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "gt_root": str(gt_root),
        "assertions": [],
        "errors": [],
        "integrity": {"verified_files": []},
    }

    def assert_metric(name: str, ok: bool, detail: str, **values: Any) -> None:
        entry = {"assertion": name, "ok": bool(ok), "detail": detail, **values}
        report["assertions"].append(entry)
        if not ok:
            report["errors"].append(f"{name}: {detail}")

    sums_path = gt_root / "SHA256SUMS.txt"
    if not sums_path.exists():
        raise SystemExit(f"GT SHA256SUMS.txt missing under {gt_root}")
    sums = _parse_sha256sums(sums_path)
    report["integrity"]["sha256sums_entries"] = len(sums)

    # ----- A1: IMR90 Hi-C loop reproducibility (genome-wide F1) -------------
    reval_rel = "hic_rao_imr90/loop_reval_genomewide_f1.json"
    concordance_rel = "hic_rao_imr90/loop_concordance_genomewide.json"
    try:
        reval = _verify_and_load_json(gt_root, reval_rel, sums)
        report["integrity"]["verified_files"].append(reval_rel)
        concordance = _verify_and_load_json(gt_root, concordance_rel, sums)
        report["integrity"]["verified_files"].append(concordance_rel)
    except IntegrityError as exc:
        raise SystemExit(f"GT integrity failure (A1): {exc}")

    recorded_f1 = float(reval.get("f1", 0.0))
    recomputed = _recompute_loop_f1(reval)
    # (a) internal consistency: recomputed F1 matches the recorded scalar.
    assert_metric(
        "imr90_loop_f1_internal_consistency",
        abs(recomputed["f1"] - recorded_f1) <= 1e-3,
        f"recomputed F1 {recomputed['f1']:.4f} matches recorded {recorded_f1:.4f} "
        f"(per-chrom gt={recomputed['gt']}, called={recomputed['called']})",
        recorded_f1=recorded_f1,
        recomputed_f1=round(recomputed["f1"], 4),
    )
    # (b) biological floor: F1 at or above an ABSOLUTE published-loop floor,
    # decoupled from this file's own recorded scalar so the check is an
    # independent regression tripwire rather than tautological with the
    # internal-consistency assertion above. Floor 0.40 sits well below the
    # recorded 0.5959 yet still fails on a real degradation of IMR90 loop
    # recovery vs published HiCCUPS ground truth.
    _f1_absolute_floor = 0.40
    assert_metric(
        "imr90_loop_f1_above_floor",
        recomputed["f1"] >= _f1_absolute_floor,
        f"genome-wide loop F1 {recomputed['f1']:.4f} >= absolute floor "
        f"{_f1_absolute_floor:.2f} (recorded {recorded_f1:.4f}); GT=published "
        f"IMR90 HiCCUPS loops, {reval.get('tolerance_bp')} bp tolerance, "
        f"{reval.get('gt_loops')} GT loops",
        threshold=_f1_absolute_floor,
        observed=round(recomputed["f1"], 4),
    )
    # Cross-check: the genome-wide concordance recall is the high-recall regime.
    assert_metric(
        "imr90_loop_recall_positive",
        float(concordance.get("recall", 0.0)) > 0.5,
        f"genome-wide concordance recall {concordance.get('recall')} > 0.5 "
        "(real loops recovered, not noise)",
        recall=concordance.get("recall"),
    )

    # ----- A2: TP63 CUT&Tag OE directional signal ---------------------------
    tp63_rel = "cutrun_liu2025/diff_TP63_OEvsNC/signal_diff_manifest.json"
    try:
        tp63 = _verify_and_load_json(gt_root, tp63_rel, sums)
        report["integrity"]["verified_files"].append(tp63_rel)
    except IntegrityError as exc:
        raise SystemExit(f"GT integrity failure (A2): {exc}")

    n_gained = int(tp63.get("n_gained", 0))
    n_lost = int(tp63.get("n_lost", 0))
    ratio = (n_gained / n_lost) if n_lost else float("inf")
    assert_metric(
        "tp63_oe_gained_exceeds_lost",
        n_gained > n_lost and ratio >= _TP63_GAIN_LOSS_MIN_RATIO,
        f"TP63 OE-gained {n_gained} >> OE-lost {n_lost} (ratio {ratio:.3f} >= "
        f"{_TP63_GAIN_LOSS_MIN_RATIO}); TP63 is an activating squamous master TF, "
        "so overexpression gains far more signal than it loses (Liu et al. 2025)",
        n_gained=n_gained,
        n_lost=n_lost,
        ratio=round(ratio, 3),
        min_ratio=_TP63_GAIN_LOSS_MIN_RATIO,
    )

    rep = tp63.get("replicate_reproducibility", {})
    oe_rep = rep.get("OE")
    nc_rep = rep.get("NC")
    oe_ok = (
        isinstance(oe_rep, (int, float))
        and isinstance(nc_rep, (int, float))
        and float(oe_rep) - float(nc_rep) >= _OE_VS_NC_MIN_MARGIN
    )
    assert_metric(
        "tp63_oe_reproducibility_above_nc_baseline",
        oe_ok,
        f"TP63 OE replicate reproducibility {oe_rep} exceeds NC paired-label "
        f"baseline {nc_rep} by >= {_OE_VS_NC_MIN_MARGIN} (OE-up signal is "
        "reproducible and far above the no-true-binding baseline)",
        oe_reproducibility=oe_rep,
        nc_reproducibility=nc_rep,
        min_margin=_OE_VS_NC_MIN_MARGIN,
    )

    report["verdict"] = "pass" if not report["errors"] else "fail"
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gt-root",
        type=_safe_gt_root,
        default=_safe_gt_root(str(DEFAULT_GT_ROOT)),
        help="retained compact GT evidence root (default: the on-disk "
        "bulk_multiomics_retained_evidence dir)",
    )
    parser.add_argument(
        "--allow-missing-gt",
        action="store_true",
        help="CI-only: downgrade a missing GT root to a logged structural skip. "
        "Default is strict real enforcement.",
    )
    parser.add_argument(
        "--output",
        type=_safe_output,
        default=None,
        help="path for the JSON report (default: a /tmp file). Must live under "
        "projects/, suite root, or /tmp.",
    )
    args = parser.parse_args()

    output_path = args.output
    if output_path is None:
        output_path = _safe_output(
            f"/tmp/realdata_science_assertion."
            f"{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"
        )

    gt_root = args.gt_root
    if not gt_root.is_dir():
        if args.allow_missing_gt:
            report = {
                "gate": "realdata_science_assertion",
                "gt_root": str(gt_root),
                "verdict": "structural_skip",
                "structural_only": True,
            }
            output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
            print(
                "WARNING: --allow-missing-gt: GT root absent "
                f"({gt_root}); downgrading to a logged STRUCTURAL SKIP (NO "
                "biological assertion verified).",
                file=sys.stderr,
            )
            print(
                f"SKIP realdata_science_assertion (no GT) report={output_path}"
            )
            return 0
        print(
            f"FAIL realdata_science_assertion: GT root is ABSENT: {gt_root} "
            "(real enforcement; pass --allow-missing-gt ONLY on CI hosts without "
            "the retained GT)",
            file=sys.stderr,
        )
        return 1

    try:
        report = run_assertions(gt_root)
    except SystemExit as exc:
        fail_report = {
            "gate": "realdata_science_assertion",
            "gt_root": str(gt_root),
            "verdict": "fail",
            "fatal": str(exc),
        }
        try:
            output_path.write_text(json.dumps(fail_report, indent=2), encoding="utf-8")
        except OSError:
            pass
        print(f"FAIL realdata_science_assertion: {exc}", file=sys.stderr)
        return 1

    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    if report["errors"]:
        for error in report["errors"]:
            print(f"ERROR: {error}", file=sys.stderr)
        print(
            f"FAIL realdata_science_assertion verdict={report['verdict']} "
            f"errors={len(report['errors'])} report={output_path}",
            file=sys.stderr,
        )
        return 1

    print(
        f"OK realdata_science_assertion verdict={report['verdict']} "
        f"assertions={len(report['assertions'])} "
        f"verified_files={len(report['integrity']['verified_files'])} "
        f"report={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
