"""Execute the governed scfactory quick-start documented in README.md."""

from __future__ import annotations

import re
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
QUICKSTART_MARKER = "<!-- scfactory-governed-h5ad-quickstart -->"


def _documented_command() -> list[str]:
    text = README.read_text(encoding="utf-8")
    match = re.search(
        re.escape(QUICKSTART_MARKER) + r"\s*```bash\s*(.*?)\s*```",
        text,
        flags=re.DOTALL,
    )
    assert match is not None, "README must expose the governed h5ad quick-start marker"
    command = match.group(1).replace("\\\n", " ")
    lines = [line for line in command.splitlines() if not line.lstrip().startswith("#")]
    return shlex.split(" ".join(lines))


def test_readme_governed_h5ad_quickstart_is_executable(tmp_path):
    anndata = pytest.importorskip("anndata")
    input_h5ad = tmp_path / "arbitrary-source-name.h5ad"
    anndata.AnnData(X=np.zeros((4, 3), dtype=np.float32)).write_h5ad(input_h5ad)
    project_root = tmp_path / "project"

    command = _documented_command()
    assert command[:3] == ["python", "scripts/scfactory.py", "run"]
    assert "--project-root" in command
    assert "--run-id" in command
    assert "--dry-run" in command

    command[0] = sys.executable
    command = [
        str(input_h5ad) if token == "INPUT_H5AD" else
        str(project_root) if token == "PROJECT_ROOT" else token
        for token in command
    ]
    proc = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert proc.returncode == 0, proc.stderr
    assert f"--input-h5ad {input_h5ad}" in proc.stdout
    assert f"--project-root {project_root}" in proc.stdout
