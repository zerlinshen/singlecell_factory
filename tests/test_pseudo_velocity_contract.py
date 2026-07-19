"""Contract tests for the exploratory pseudo-velocity proxy."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from workflow.modular.modules.pseudo_velocity import PseudoVelocityModule


class _Ctx:
    def __init__(self, adata, tmp_path: Path):
        self.adata = adata
        self.metadata = {}
        self.module_status = []
        self.table_dir = tmp_path / "tables"
        self.figure_dir = tmp_path / "figures"
        self.table_dir.mkdir()
        self.figure_dir.mkdir()

    def status(self, module, ok, message):
        self.module_status.append({"module": module, "status": ok, "message": message})


def _make_adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    adata = ad.AnnData(np.zeros((60, 3), dtype=np.float32))
    adata.obsm["X_umap"] = rng.normal(size=(60, 2)).astype(np.float32)
    adata.obs["dpt_pseudotime"] = np.linspace(0, 1, 60, dtype=np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 3) for i in range(60)])
    return adata


def test_every_status_and_table_marks_exploratory_proxy(monkeypatch, tmp_path):
    adata = _make_adata()
    ctx = _Ctx(adata, tmp_path)
    for method_name in (
        "_plot_speed_umap",
        "_plot_arrows",
        "_plot_stream",
        "_plot_arrows_temporal",
        "_plot_celltype_arrows",
        "_plot_celltype_stream",
        "_plot_speed_boxplot",
    ):
        monkeypatch.setattr(PseudoVelocityModule, method_name, staticmethod(lambda *args: None))

    PseudoVelocityModule().run(ctx)

    provenance = adata.uns["pseudo_velocity"]
    assert provenance["inference_class"] == "exploratory_proxy"
    assert provenance["claim_status"] == "non_claimable_as_canonical_rna_velocity"
    assert provenance["canonical_rna_velocity"] is False
    assert ctx.metadata["pseudo_velocity_inference_class"] == "exploratory_proxy"
    assert ctx.metadata["pseudo_velocity_canonical_rna_velocity"] is False
    assert ctx.module_status == [{
        "module": "pseudo_velocity",
        "status": "ok",
        "message": "exploratory_proxy; non-claimable as canonical RNA velocity",
    }]

    for table_name in ("pseudo_velocity_speed.csv", "pseudo_velocity_per_cluster.csv"):
        table = pd.read_csv(ctx.table_dir / table_name)
        assert set(table["inference_class"]) == {"exploratory_proxy"}
        assert set(table["claim_status"]) == {
            "non_claimable_as_canonical_rna_velocity"
        }
        assert not table["canonical_rna_velocity"].any()


def test_plot_title_visibly_rejects_canonical_rna_velocity(monkeypatch, tmp_path):
    adata = _make_adata()
    adata.obs["pseudo_velocity_speed"] = np.linspace(0, 1, adata.n_obs)
    ctx = _Ctx(adata, tmp_path)
    saved_titles = []

    def _capture_titles(*args, **kwargs):
        saved_titles.extend(axis.get_title() for axis in plt.gcf().axes)

    monkeypatch.setattr(plt, "savefig", _capture_titles)

    PseudoVelocityModule._plot_speed_umap(adata, ctx)

    joined = " ".join(saved_titles)
    assert "Exploratory proxy" in joined
    assert "NOT canonical RNA velocity" in joined
