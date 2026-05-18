"""Wave 2B / US-013 — Ribo-seq vertical slice integration smoke."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


FACTORY_ROOT = Path(__file__).resolve().parent.parent
R_FACTORY_ROOT = FACTORY_ROOT.parent / "multiomics_r_factory"


@pytest.fixture
def synthetic_adata():
    rng = np.random.default_rng(42)
    n_cells, n_genes = 40, 30
    X = sp.random(n_cells, n_genes, density=0.3, dtype=np.float32, random_state=42).tocsr()
    barcodes = [f"BCRIB{i:04d}-1" for i in range(n_cells)]
    obs = pd.DataFrame(
        {
            "leiden": rng.integers(0, 3, size=n_cells).astype(str),
            "sample": [f"S{i % 3}" for i in range(n_cells)],
        },
        index=barcodes,
    )
    var = pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def synthetic_ribo_footprints(tmp_path):
    """TSV with gene/sample/footprint_count for 3 samples × 30 genes."""
    rng = np.random.default_rng(0)
    rows = []
    for s in ("S0", "S1", "S2"):
        for i in range(30):
            rows.append({"gene": f"GENE_{i:03d}", "sample": s,
                         "footprint_count": int(rng.integers(5, 100))})
    df = pd.DataFrame(rows)
    path = tmp_path / "footprints.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_ribo_ingest_skips_without_footprints(synthetic_adata, tmp_path):
    from workflow.modular.modules.ribo_ingest import RiboIngestModule
    cfg = SimpleNamespace()
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    RiboIngestModule().run(ctx)
    assert ctx.metadata["ribo_ingest_status"] == "skipped_no_input"


def test_ribo_ingest_computes_te_with_rna_from_adata(synthetic_adata, synthetic_ribo_footprints, tmp_path):
    from workflow.modular.modules.ribo_ingest import RiboIngestModule
    cfg = SimpleNamespace(
        ribo_footprints_path=str(synthetic_ribo_footprints),
        ribo_min_footprint_count=1,
    )
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    RiboIngestModule().run(ctx)
    te_df = synthetic_adata.uns["ribo_translation_efficiency"]
    assert {"gene", "sample", "footprint_count", "rna_count", "te"}.issubset(te_df.columns)
    assert te_df.shape[0] > 0
    assert (te_df["te"] >= 0).all()
    meta = synthetic_adata.uns["ribo_ingest_metadata"]
    assert meta["n_samples"] >= 3
    assert meta["te_mean"] >= 0
    # Outputs written
    assert (tmp_path / "ribo_ingest" / "translation_efficiency.csv").exists()
    assert (tmp_path / "ribo_ingest" / "ribo_summary.json").exists()


def test_schema_parity_after_ribo_contract_check():
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    rf = R_FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    assert sf.read_bytes() == rf.read_bytes()


def test_ribo_slot_reserved_until_bundle_vertical_slice():
    import yaml
    doc = yaml.safe_load((FACTORY_ROOT / "contracts" / "bundle_schema.yaml").read_text())
    ribo = doc["extensions_v22"]["ribo"]
    assert ribo["status"] == "reserved"
    assert "table" in ribo
    assert "te_path" in ribo["table"]


def test_r_extension_registry_includes_ribo():
    src = (R_FACTORY_ROOT / "R_bundle" / "io_bundle.R").read_text()
    for line in src.splitlines():
        if line.strip().startswith("KNOWN_BUNDLE_EXTENSIONS"):
            assert '"ribo"' in line
            return
    raise AssertionError("KNOWN_BUNDLE_EXTENSIONS line not found")


def test_r_ribo_module_present():
    rmod = R_FACTORY_ROOT / "R" / "ribo_module.R"
    assert rmod.exists()
    src = rmod.read_text()
    for fn in ("load_ribo_extension", "plot_translation_efficiency", "plot_te_volcano"):
        assert fn in src


def test_ribo_module_in_registry():
    from workflow.modular.pipeline import _build_registry, MUTATING_MODULES, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES
    reg = _build_registry()
    assert "ribo_ingest" in reg
    assert "ribo_ingest" not in MUTATING_MODULES
    order = _resolve_execution_order(list(MANDATORY_MODULES), ["ribo_ingest"])
    assert "ribo_ingest" in order
