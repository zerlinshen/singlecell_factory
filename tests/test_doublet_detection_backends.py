from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.modular.config import CellRangerConfig, DoubletConfig, PipelineConfig, QCConfig
from workflow.modular.context import PipelineContext
from workflow.modular.modules.doublet_detection import DoubletDetectionModule


def _adata(n_obs: int = 30, n_vars: int = 60) -> ad.AnnData:
    rng = np.random.default_rng(7)
    x = sparse.csr_matrix(rng.poisson(1.2, size=(n_obs, n_vars)).astype("float32"))
    obs = pd.DataFrame(index=[f"cell{i:03d}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(n_vars)])
    return ad.AnnData(X=x, obs=obs, var=var)


def _ctx(tmp_path: Path, doublet: DoubletConfig) -> PipelineContext:
    cfg = PipelineConfig(
        project="doublet-test",
        output_dir=tmp_path,
        cellranger=CellRangerConfig(sample_root=tmp_path, outs_dir=tmp_path),
        qc=QCConfig(min_genes=0, min_counts=0),
        doublet=doublet,
        random_state=11,
    )
    fig_dir = tmp_path / "figures"
    table_dir = tmp_path / "tables"
    fig_dir.mkdir()
    table_dir.mkdir()
    return PipelineContext(cfg=cfg, run_dir=tmp_path, figure_dir=fig_dir, table_dir=table_dir)


def test_scdblfinder_backend_records_effective_metadata(monkeypatch, tmp_path):
    adata = _adata()
    calls = np.zeros(adata.n_obs, dtype=bool)
    calls[[1, 3, 5]] = True
    scores = np.linspace(0.0, 1.0, adata.n_obs, dtype=np.float32)

    def fake_scdbl(self, adata_arg, cfg, ctx):
        assert adata_arg is adata
        return scores, calls, None

    monkeypatch.setattr(DoubletDetectionModule, "_run_scdblfinder_via_r", fake_scdbl)
    ctx = _ctx(tmp_path, DoubletConfig(backend="scdblfinder", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_backend_requested"] == "scdblfinder"
    assert ctx.metadata["doublet_backend"] == "scdblfinder"
    assert ctx.metadata["doublet_method"] == "scdblfinder"
    assert ctx.metadata["doublets_detected"] == 3
    assert ctx.adata.obs["predicted_doublet"].sum() == 3


def test_scdbl_alias_maps_to_scdblfinder(monkeypatch, tmp_path):
    adata = _adata()
    calls = np.zeros(adata.n_obs, dtype=bool)
    scores = np.linspace(0.0, 1.0, adata.n_obs, dtype=np.float32)

    def fake_scdbl(self, adata_arg, cfg, ctx):
        return scores, calls, None

    monkeypatch.setattr(DoubletDetectionModule, "_run_scdblfinder_via_r", fake_scdbl)
    ctx = _ctx(tmp_path, DoubletConfig(backend="scdbl", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_backend_requested"] == "scdbl"
    assert ctx.metadata["doublet_backend"] == "scdblfinder"
    assert ctx.metadata["doublet_method"] == "scdblfinder"


def test_unknown_backend_records_fallback_reason_on_small_dataset(tmp_path):
    adata = _adata(n_obs=10, n_vars=20)
    ctx = _ctx(tmp_path, DoubletConfig(backend="bogus", remove_doublets=False))
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_backend_requested"] == "bogus"
    assert ctx.metadata["doublet_backend"] == "scrublet"
    assert ctx.metadata["doublet_backend_fallback_reason"] == "unknown_backend:bogus"
    assert ctx.metadata["doublet_method"] == "fallback_all_singlets"


def test_consensus_can_use_scrublet_scdblfinder_pair(monkeypatch, tmp_path):
    adata = _adata()
    scrub_calls = np.zeros(adata.n_obs, dtype=bool)
    scrub_calls[[0, 1]] = True
    scdbl_calls = np.zeros(adata.n_obs, dtype=bool)
    scdbl_calls[[1, 2, 3]] = True

    def fake_scrub(self, adata_arg, cfg, ctx):
        return np.linspace(0.0, 1.0, adata_arg.n_obs, dtype=np.float32), scrub_calls, 0.5

    def fake_scdbl(self, adata_arg, cfg, ctx):
        return np.linspace(1.0, 0.0, adata_arg.n_obs, dtype=np.float32), scdbl_calls, None

    monkeypatch.setattr(DoubletDetectionModule, "_run_scrublet_only", fake_scrub)
    monkeypatch.setattr(DoubletDetectionModule, "_run_scdblfinder_via_r", fake_scdbl)
    ctx = _ctx(
        tmp_path,
        DoubletConfig(
            backend="consensus",
            consensus_pair="scrublet_scdblfinder",
            consensus_logic="or",
            remove_doublets=False,
        ),
    )
    ctx.adata = adata

    DoubletDetectionModule().run(ctx)

    assert ctx.metadata["doublet_method"] == "consensus_or_scrublet_scdblfinder"
    assert ctx.metadata["doublet_consensus_pair"] == "scrublet_scdblfinder"
    assert ctx.metadata["doublet_consensus_backends"] == "scrublet,scdblfinder"
    assert ctx.metadata["doublet_consensus_per_backend_rates"]["final"] == 4 / adata.n_obs
    assert "scrublet_score" in ctx.adata.obs
    assert "scdblfinder_score" in ctx.adata.obs
    assert ctx.adata.obs["predicted_doublet"].sum() == 4


def test_resolve_scdblfinder_driver_from_script_dir(monkeypatch, tmp_path):
    driver = tmp_path / "scdblfinder_run.R"
    driver.write_text("# dummy\n", encoding="utf-8")
    monkeypatch.setenv("R_MULTIOMICS_SCRIPT_DIR", str(tmp_path))

    assert DoubletDetectionModule._resolve_scdblfinder_driver() == driver
