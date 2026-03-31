from pathlib import Path

from workflow.velocity import VelocityWorkflowConfig, run_velocity_workflow


def test_run_velocity_workflow(monkeypatch, tmp_path):
    import workflow.velocity as mod

    root = Path("/home/zerlinshen/singlecell_factory")
    scripts = root / "rna_velocity_pseudotime_analysis" / "scripts"
    scripts.mkdir(parents=True, exist_ok=True)
    (scripts / "generate_velocity_loom.sh").write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")
    (scripts / "run_rna_velocity_pseudotime_analysis.sh").write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")

    calls = []

    def fake_run(cmd, check=True):
        calls.append(cmd)
        return 0

    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    cfg = VelocityWorkflowConfig(
        samplefolder=tmp_path / "count",
        gtf_path=tmp_path / "ref.gtf",
        tenx_dir=tmp_path / "tenx",
        run_id="r1",
        config_path=tmp_path / "cfg.yaml",
        python_bin=Path("/usr/bin/python3"),
        threads=8,
    )
    out = run_velocity_workflow(cfg)
    assert len(calls) == 2
    assert out.name == "r1"


def test_run_velocity_workflow_missing_scripts(tmp_path):
    import workflow.velocity as mod

    cfg = VelocityWorkflowConfig(
        samplefolder=tmp_path / "count",
        gtf_path=tmp_path / "ref.gtf",
        tenx_dir=tmp_path / "tenx",
        run_id="r1",
        config_path=tmp_path / "cfg.yaml",
        python_bin=Path("/usr/bin/python3"),
    )
    orig = Path.exists

    def fake_exists(self):
        if str(self).endswith("generate_velocity_loom.sh"):
            return False
        return orig(self)

    import pathlib

    pathlib.Path.exists = fake_exists
    try:
        try:
            run_velocity_workflow(cfg)
        except FileNotFoundError:
            pass
        else:
            raise AssertionError("Expected FileNotFoundError")
    finally:
        pathlib.Path.exists = orig


def test_run_velocity_workflow_missing_run_script(tmp_path):
    import workflow.velocity as mod
    import pathlib

    cfg = VelocityWorkflowConfig(
        samplefolder=tmp_path / "count",
        gtf_path=tmp_path / "ref.gtf",
        tenx_dir=tmp_path / "tenx",
        run_id="r1",
        config_path=tmp_path / "cfg.yaml",
        python_bin=Path("/usr/bin/python3"),
    )
    orig = pathlib.Path.exists

    def fake_exists(self):
        s = str(self)
        if s.endswith("generate_velocity_loom.sh"):
            return True
        if s.endswith("run_rna_velocity_pseudotime_analysis.sh"):
            return False
        return orig(self)

    pathlib.Path.exists = fake_exists
    try:
        try:
            run_velocity_workflow(cfg)
        except FileNotFoundError:
            pass
        else:
            raise AssertionError("Expected FileNotFoundError")
    finally:
        pathlib.Path.exists = orig


def test_velocity_main(monkeypatch, tmp_path):
    import workflow.velocity as mod

    monkeypatch.setattr(
        mod,
        "parse_args",
        lambda: type(
            "Args",
            (),
            {
                "samplefolder": str(tmp_path / "count"),
                "gtf_path": str(tmp_path / "a.gtf"),
                "tenx_dir": str(tmp_path / "tenx"),
                "run_id": "r2",
                "config_path": str(tmp_path / "cfg.yaml"),
                "python_bin": "/usr/bin/python3",
                "threads": 4,
            },
        )(),
    )
    monkeypatch.setattr(mod, "run_velocity_workflow", lambda cfg: tmp_path / "out")
    mod.main()


def test_velocity_parse_args(monkeypatch):
    import workflow.velocity as mod

    monkeypatch.setattr(
        "sys.argv",
        [
            "prog",
            "--samplefolder",
            "/tmp/s",
            "--gtf-path",
            "/tmp/g.gtf",
            "--tenx-dir",
            "/tmp/t",
            "--run-id",
            "r9",
            "--config-path",
            "/tmp/c.yaml",
            "--python-bin",
            "/usr/bin/python3",
            "--threads",
            "6",
        ],
    )
    args = mod.parse_args()
    assert args.run_id == "r9"
    assert args.threads == 6
