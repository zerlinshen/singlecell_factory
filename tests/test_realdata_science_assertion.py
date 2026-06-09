"""Hermetic tests for scripts/validate_realdata_science_assertion.py (#34).

Build a synthetic GT directory in tmp_path that mirrors the retained-evidence
key shapes (IMR90 loop reval per-chrom counts + TP63 OE/NC signal-diff
manifest), write a real SHA256SUMS.txt over those files, then:
  * assert the biological assertions PASS on integrity-clean GT;
  * tamper a metric below threshold (TP63 OE gains < losses) -> FAIL;
  * corrupt a file after the checksum is written (sha256 mismatch) -> hard
    SystemExit (integrity failure).

No real GT path is touched; everything lives under tmp_path.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "validate_realdata_science_assertion.py"
)
_spec = importlib.util.spec_from_file_location(
    "validate_realdata_science_assertion", MODULE_PATH
)
gate = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(gate)


def _per_chrom(recall: float, precision: float) -> dict[str, dict[str, int]]:
    """Build a 2-chrom per-chrom block whose totals yield ~(recall, precision).

    gt=1000, called=1000 across two chroms; gt_matched=recall*gt,
    called_matched=precision*called.
    """
    return {
        "1": {
            "gt": 600,
            "called": 600,
            "gt_matched": round(recall * 600),
            "called_matched": round(precision * 600),
        },
        "2": {
            "gt": 400,
            "called": 400,
            "gt_matched": round(recall * 400),
            "called_matched": round(precision * 400),
        },
    }


def _f1(recall: float, precision: float) -> float:
    return 2 * precision * recall / (precision + recall)


def _build_gt(
    tmp_path: Path,
    *,
    tp63_gained: int = 60128,
    tp63_lost: int = 27544,
    oe_rep: float = 0.7975,
    nc_rep: float = 0.2197,
    recall: float = 0.6642,
    precision: float = 0.5404,
) -> Path:
    """Create a synthetic GT root with a valid SHA256SUMS; return the root."""
    gt_root = tmp_path / "gt"
    (gt_root / "hic_rao_imr90").mkdir(parents=True)
    (gt_root / "cutrun_liu2025" / "diff_TP63_OEvsNC").mkdir(parents=True)

    per_chrom = _per_chrom(recall, precision)
    gt_total = sum(c["gt"] for c in per_chrom.values())
    called_total = sum(c["called"] for c in per_chrom.values())
    gt_matched = sum(c["gt_matched"] for c in per_chrom.values())
    called_matched = sum(c["called_matched"] for c in per_chrom.values())
    eff_recall = gt_matched / gt_total
    eff_precision = called_matched / called_total
    eff_f1 = _f1(eff_recall, eff_precision)

    reval = {
        "headline_metric": "f1",
        "gt_loops": gt_total,
        "called_loops": called_total,
        "recall": round(eff_recall, 4),
        "precision": round(eff_precision, 4),
        "f1": round(eff_f1, 4),
        "tolerance_bp": 20000,
        "per_chrom": per_chrom,
    }
    concordance = {
        "gt_loops": gt_total,
        "called_loops": 2000,
        "recall": 0.7405,
        "precision": 0.2383,
        "f1": 0.3606,
        "tolerance_bp": 20000,
    }
    tp63 = {
        "lane": "chip_atac_signal_diff",
        "test": "OE",
        "ref": "NC",
        "n_peaks": 164175,
        "n_gained": tp63_gained,
        "n_lost": tp63_lost,
        "replicate_reproducibility": {"NC": nc_rep, "OE": oe_rep},
        "median_abs_log2FC": 1.0968,
    }

    files = {
        "hic_rao_imr90/loop_reval_genomewide_f1.json": reval,
        "hic_rao_imr90/loop_concordance_genomewide.json": concordance,
        "cutrun_liu2025/diff_TP63_OEvsNC/signal_diff_manifest.json": tp63,
    }
    for rel, payload in files.items():
        (gt_root / rel).write_text(json.dumps(payload, indent=2), encoding="utf-8")

    _write_sha256sums(gt_root, list(files))
    return gt_root


def _write_sha256sums(gt_root: Path, rels: list[str]) -> None:
    lines = []
    for rel in rels:
        sha = hashlib.sha256((gt_root / rel).read_bytes()).hexdigest()
        lines.append(f"{sha}  ./{rel}")
    (gt_root / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_assertions_pass_on_clean_gt(tmp_path: Path) -> None:
    gt_root = _build_gt(tmp_path)
    assert gate._is_relative_to(gt_root, Path("/tmp").resolve()), (
        f"tmp_path {gt_root} not under /tmp; gate allowlist would reject it"
    )
    report = gate.run_assertions(gt_root)
    assert report["verdict"] == "pass", report["errors"]
    assert report["errors"] == []
    assert len(report["integrity"]["verified_files"]) == 3
    names = {a["assertion"] for a in report["assertions"]}
    assert "imr90_loop_f1_above_floor" in names
    assert "tp63_oe_gained_exceeds_lost" in names


def test_tp63_gains_below_threshold_fails(tmp_path: Path) -> None:
    # OE gains < losses violates the directional biological assertion.
    gt_root = _build_gt(tmp_path, tp63_gained=10000, tp63_lost=27544)
    report = gate.run_assertions(gt_root)
    assert report["verdict"] == "fail"
    assert any("tp63_oe_gained_exceeds_lost" in err for err in report["errors"]), (
        report["errors"]
    )


def test_low_loop_f1_fails(tmp_path: Path) -> None:
    # Drive recall/precision down so the recomputed F1 sits below the recorded
    # floor it is compared against (the recorded scalar is derived from the same
    # low per-chrom counts, so internal consistency still holds, but we then
    # tamper the recorded f1 UP to create an above-floor violation).
    gt_root = _build_gt(tmp_path, recall=0.30, precision=0.20)
    reval_path = gt_root / "hic_rao_imr90" / "loop_reval_genomewide_f1.json"
    reval = json.loads(reval_path.read_text())
    reval["f1"] = 0.90  # claim a high F1 the per-chrom counts cannot support
    reval_path.write_text(json.dumps(reval, indent=2))
    _write_sha256sums(
        gt_root,
        [
            "hic_rao_imr90/loop_reval_genomewide_f1.json",
            "hic_rao_imr90/loop_concordance_genomewide.json",
            "cutrun_liu2025/diff_TP63_OEvsNC/signal_diff_manifest.json",
        ],
    )
    report = gate.run_assertions(gt_root)
    assert report["verdict"] == "fail"
    # Either the internal-consistency check or the above-floor check must fire.
    assert any("imr90_loop_f1" in err for err in report["errors"]), report["errors"]


def test_sha256_mismatch_is_hard_failure(tmp_path: Path) -> None:
    # Corrupt a GT file AFTER the checksum is recorded -> integrity SystemExit.
    gt_root = _build_gt(tmp_path)
    tp63_path = gt_root / "cutrun_liu2025" / "diff_TP63_OEvsNC" / "signal_diff_manifest.json"
    tampered = json.loads(tp63_path.read_text())
    tampered["n_gained"] = 999999  # change content without updating SHA256SUMS
    tp63_path.write_text(json.dumps(tampered, indent=2))
    with pytest.raises(SystemExit):
        gate.run_assertions(gt_root)


def test_missing_sha256sums_is_hard_failure(tmp_path: Path) -> None:
    gt_root = _build_gt(tmp_path)
    (gt_root / "SHA256SUMS.txt").unlink()
    with pytest.raises(SystemExit):
        gate.run_assertions(gt_root)
