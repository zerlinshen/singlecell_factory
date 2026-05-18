from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_project_governance.py"
_SPEC = importlib.util.spec_from_file_location("validate_project_governance", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
validator = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(validator)


class ProjectGovernanceValidatorTests(unittest.TestCase):
    def test_validate_project_governance_accepts_legacy_missing_manifests(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp) / "legacy-project"
            run = project / "runs" / "2026-05-18T0101Z-abcdef0"
            (run / "python" / "figures").mkdir(parents=True)
            (run / "r").mkdir()
            (run / "logs").mkdir()
            (project / "inputs").mkdir()
            (project / "configs").mkdir()
            (project / "notebooks").mkdir()
            (project / "project.yaml").write_text("project_id: legacy-project\n", encoding="utf-8")
            (run / "python" / "figures" / "umap.png").write_bytes(b"fake")

            report = validator.validate_project(project)

        self.assertEqual(report["overall_status"], "pass")
        self.assertGreaterEqual(report["severity_counts"].get("warning", 0), 1)
        self.assertIsNone(report["runs"][0]["root_manifest"])
        codes = {f["code"] for f in report["runs"][0]["findings"]}
        self.assertIn("missing_root_manifest", codes)

    def test_validate_project_governance_separates_root_and_producer_manifests(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp) / "new-project"
            run = project / "runs" / "2026-05-18T0101Z-abcdef0"
            (run / "python" / "bundle").mkdir(parents=True)
            (run / "r").mkdir()
            (run / "logs").mkdir()
            (project / "project.yaml").write_text("project_id: new-project\n", encoding="utf-8")
            (run / "manifest.json").write_text("{}\n", encoding="utf-8")
            (run / "python" / "run_manifest.json").write_text("{}\n", encoding="utf-8")
            (run / "python" / "module_status.csv").write_text("module,status\nqc,ok\n", encoding="utf-8")
            (run / "python" / "bundle" / "provenance.json").write_text("{}\n", encoding="utf-8")

            report = validator.validate_project(project, run_id=run.name)

        run_report = report["runs"][0]
        self.assertEqual(report["overall_status"], "pass")
        self.assertTrue(run_report["root_manifest"].endswith("manifest.json"))
        self.assertEqual(run_report["producer_native_provenance"]["python_run_manifests"], ["run_manifest.json"])
        self.assertEqual(run_report["producer_native_provenance"]["python_module_status"], ["module_status.csv"])
        self.assertEqual(run_report["producer_native_provenance"]["bundle_provenance"], ["bundle/provenance.json"])

    def test_validate_project_governance_missing_project_root_is_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = validator.validate_project(Path(tmp) / "missing")

        self.assertEqual(report["overall_status"], "fail")
        self.assertEqual(report["severity_counts"]["error"], 1)


class C2ReproductionFieldsValidatorTests(unittest.TestCase):
    def test_missing_all_c2_fields_warns_but_does_not_error(self) -> None:
        data: dict = {}
        warnings = validator._validate_c2_reproduction_fields(data, Path("project.yaml"))
        codes = {w["code"] for w in warnings}
        for field in validator._C2_TOP_LEVEL_FIELDS:
            self.assertIn(f"c2_missing_{field}", codes)
        for w in warnings:
            self.assertEqual(w["severity"], "warning")

    def test_invalid_parity_class_warns(self) -> None:
        data = {
            "data_object_reproduction": {"obs": "exact", "obsm": "bogus_class"},
        }
        warnings = validator._validate_c2_reproduction_fields(data, Path("project.yaml"))
        codes = {w["code"] for w in warnings}
        self.assertIn("c2_data_object_reproduction_invalid_parity_class", codes)

    def test_invalid_gap_target_warns(self) -> None:
        data = {
            "module_gap_decisions": [
                {"module": "atac", "gap_class": "approximate", "target": "made_up_factory"},
            ],
        }
        warnings = validator._validate_c2_reproduction_fields(data, Path("project.yaml"))
        codes = {w["code"] for w in warnings}
        self.assertIn("c2_module_gap_decisions_invalid_target", codes)

    def test_fully_populated_c2_fields_silent(self) -> None:
        data = {
            "upstream_repository": {
                "url": "https://example.com/repo.git",
                "commit_or_tag": "abc1234",
                "doi": "10.0/example",
                "license": "MIT",
                "forked_at": "2026-05-18",
                "fork_path": "projects/<id>/upstream/repo/",
            },
            "raw_data_reproduction": {"attempted": True, "success": True},
            "data_object_reproduction": {"obs": "exact", "obsm": "approximate"},
            "figure_reproduction": {"figure1a": "proxy"},
            "module_gap_decisions": [
                {"module": "atac", "gap_class": "approximate", "target": "r_multiomics_factory", "rationale": "x"},
            ],
            "context_optimization_decisions": [
                {"parameter": "n_pcs", "upstream_value": 30, "our_value": 50, "rationale": "y", "parity_class": "approximate"},
            ],
        }
        warnings = validator._validate_c2_reproduction_fields(data, Path("project.yaml"))
        # upstream_repository fully populated → no warnings; other fields valid → no warnings.
        # Should produce empty warnings list.
        self.assertEqual(warnings, [])


if __name__ == "__main__":
    unittest.main()
