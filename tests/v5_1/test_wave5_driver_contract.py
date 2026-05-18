from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
DRIVER = REPO_ROOT / "scripts" / "dev" / "run_wave5_trevino_pcw21.py"


def _load_driver():
    spec = importlib.util.spec_from_file_location("wave5_trevino_driver", DRIVER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ctx():
    return SimpleNamespace(
        adata=SimpleNamespace(obsm={}, uns={}),
        set_module_dir=lambda name: None,
    )


def test_driver_rejects_silent_multimodal_skip(monkeypatch):
    driver = _load_driver()

    class FakeSkippedModule:
        name = "multimodal_integration"

        def __init__(self, config):
            self.config = config

        def run(self, ctx):
            ctx.adata.uns["multimodal_status"] = {
                "engine": "wnn",
                "status": "skipped",
                "reason": "Rscript not found",
            }

    monkeypatch.setattr(driver, "MultimodalIntegrationModule", FakeSkippedModule)

    try:
        driver._run_multimodal(_ctx())
    except RuntimeError as exc:
        assert "completed without 'X_wnn'" in str(exc)
        assert "Rscript not found" in str(exc)
    else:
        raise AssertionError("silent WNN skip was not rejected")


def test_driver_accepts_multimodal_when_x_wnn_is_written(monkeypatch):
    driver = _load_driver()

    class FakeOkModule:
        name = "multimodal_integration"

        def __init__(self, config):
            self.config = config

        def run(self, ctx):
            ctx.adata.obsm["X_wnn"] = np.zeros((3, 2), dtype=np.float32)
            ctx.adata.uns["multimodal_status"] = {"engine": "wnn", "status": "ok"}

    monkeypatch.setattr(driver, "MultimodalIntegrationModule", FakeOkModule)
    ctx = _ctx()

    driver._run_multimodal(ctx)

    assert ctx.adata.obsm["X_wnn"].shape == (3, 2)
