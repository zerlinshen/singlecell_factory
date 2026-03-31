from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import argparse


@dataclass(frozen=True)
class VelocityWorkflowConfig:
    """RNA velocity workflow configuration.

    Parameters
    ----------
    samplefolder
        Cell Ranger count folder containing `outs/`.
    gtf_path
        Path to annotation GTF or GTF.GZ file.
    tenx_dir
        Matrix directory used for velocity integration.
    run_id
        Unique run identifier.
    config_path
        YAML config for velocity runner.
    python_bin
        Python interpreter path for isolated execution.
    """

    samplefolder: Path
    gtf_path: Path
    tenx_dir: Path
    run_id: str
    config_path: Path
    python_bin: Path
    threads: int = 4


def run_velocity_workflow(cfg: VelocityWorkflowConfig) -> Path:
    """Run loom generation and velocity estimation workflow.

    Returns
    -------
    Path
        Path to the generated velocity run result directory.
    """

    root = Path("/home/zerlinshen/singlecell_factory")
    loom_script = root / "rna_velocity_pseudotime_analysis" / "scripts" / "generate_velocity_loom.sh"
    run_script = root / "rna_velocity_pseudotime_analysis" / "scripts" / "run_rna_velocity_pseudotime_analysis.sh"
    if not loom_script.exists():
        raise FileNotFoundError(f"Missing loom generation script: {loom_script}")
    if not run_script.exists():
        raise FileNotFoundError(f"Missing velocity run script: {run_script}")

    subprocess.run(
        [
            "bash",
            str(loom_script),
            "--samplefolder",
            str(cfg.samplefolder),
            "--gtf",
            str(cfg.gtf_path),
            "--threads",
            str(cfg.threads),
            "--run-tag",
            cfg.run_id,
        ],
        check=True,
    )

    subprocess.run(
        [
            "bash",
            str(run_script),
            "--dataset",
            cfg.samplefolder.name,
            "--run-id",
            cfg.run_id,
            "--config",
            str(cfg.config_path),
            "--python-bin",
            str(cfg.python_bin),
        ],
        check=True,
    )

    return root / "rna_velocity_pseudotime_analysis" / "runtime" / "results" / cfg.run_id


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for velocity workflow."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--samplefolder", required=True)
    parser.add_argument("--gtf-path", required=True)
    parser.add_argument("--tenx-dir", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--config-path", required=True)
    parser.add_argument("--python-bin", required=True)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    """CLI entrypoint for velocity workflow."""

    args = parse_args()
    cfg = VelocityWorkflowConfig(
        samplefolder=Path(args.samplefolder),
        gtf_path=Path(args.gtf_path),
        tenx_dir=Path(args.tenx_dir),
        run_id=args.run_id,
        config_path=Path(args.config_path),
        python_bin=Path(args.python_bin),
        threads=args.threads,
    )
    out = run_velocity_workflow(cfg)
    print(out)


if __name__ == "__main__":
    main()
