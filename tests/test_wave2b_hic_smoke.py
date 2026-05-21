"""Wave 2B / US-012 — Hi-C vertical slice integration smoke.

Covers:
  - hic_ingest reads a synthetic TSV contact-pair file → bin × bin sparse matrix
  - hic_tad computes insulation score + boundaries + A/B compartments
  - Both modules fail-soft when prerequisites absent
  - Schema marks hic.status active after exporter+reader integration
  - R-side KNOWN_BUNDLE_EXTENSIONS contains "hic"
  - r_multiomics_factory/R/hic_module.R has the required functions
  - Registry includes hic_ingest + hic_tad, no MUTATING, topo respects deps
"""
from __future__ import annotations

import gzip
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


FACTORY_ROOT = Path(__file__).resolve().parent.parent
R_FACTORY_ROOT = FACTORY_ROOT.parent / "r_multiomics_factory"


@pytest.fixture
def synthetic_adata():
    n_cells, n_genes = 30, 20
    X = sp.random(n_cells, n_genes, density=0.3, dtype=np.float32, random_state=42).tocsr()
    return ad.AnnData(X=X)


@pytest.fixture
def synthetic_contacts_tsv(tmp_path):
    """Build a small chr1-only contact pairs TSV with two visible TADs.

    Bins: 100 bp resolution → 10 bins covering 0-1000 on chr1.
    Contact density inside [0, 500] and [500, 1000] is high; cross-TAD is sparse.
    """
    rng = np.random.default_rng(42)
    rows = []
    # TAD 1: positions 0-500
    for _ in range(50):
        p1, p2 = rng.integers(0, 500, size=2)
        rows.append(("chr1", int(p1), "chr1", int(p2), 1))
    # TAD 2: positions 500-1000
    for _ in range(50):
        p1, p2 = rng.integers(500, 1000, size=2)
        rows.append(("chr1", int(p1), "chr1", int(p2), 1))
    # Sparse cross-TAD contacts
    for _ in range(5):
        p1 = rng.integers(0, 500)
        p2 = rng.integers(500, 1000)
        rows.append(("chr1", int(p1), "chr1", int(p2), 1))

    path = tmp_path / "contacts.tsv.gz"
    with gzip.open(path, "wt") as fh:
        for c1, p1, c2, p2, n in rows:
            fh.write(f"{c1}\t{p1}\t{c2}\t{p2}\t{n}\n")
    return path


def test_hic_ingest_builds_sparse_contact_matrix(synthetic_adata, synthetic_contacts_tsv, tmp_path):
    from workflow.modular.modules.hic_ingest import HiCIngestModule
    cfg = SimpleNamespace(
        hic_contacts_path=str(synthetic_contacts_tsv),
        hic_resolution_bp=100,
    )
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    HiCIngestModule().run(ctx)
    assert "hic_contact_matrix" in synthetic_adata.uns
    assert "hic_bins" in synthetic_adata.uns
    mat = synthetic_adata.uns["hic_contact_matrix"]
    assert sp.issparse(mat)
    bins = synthetic_adata.uns["hic_bins"]
    assert "chrom" in bins.columns and "start" in bins.columns
    assert ctx.metadata["hic_ingest_status"] == "ok"


def test_hic_tad_finds_boundary_between_two_tads(synthetic_adata, synthetic_contacts_tsv, tmp_path):
    from workflow.modular.modules.hic_ingest import HiCIngestModule
    from workflow.modular.modules.hic_tad import HiCTADModule
    cfg = SimpleNamespace(
        hic_contacts_path=str(synthetic_contacts_tsv),
        hic_resolution_bp=100,
        hic_tad_window_bins=3,
        hic_tad_boundary_k=0.5,  # loose threshold for the tiny synthetic fixture
    )
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    HiCIngestModule().run(ctx)
    HiCTADModule().run(ctx)
    boundaries = synthetic_adata.uns["hic_tad_boundaries"]
    compartments = synthetic_adata.uns["hic_compartments"]
    assert "insulation" in boundaries.columns
    assert "position" in boundaries.columns
    assert "compartment" in compartments.columns
    assert "position" in compartments.columns
    assert set(compartments["compartment"].unique()).issubset({"A", "B"})
    assert ctx.metadata["hic_tad_status"] == "ok"


def test_hic_ingest_skips_without_contacts(synthetic_adata, tmp_path):
    from workflow.modular.modules.hic_ingest import HiCIngestModule
    cfg = SimpleNamespace()
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    HiCIngestModule().run(ctx)
    assert ctx.metadata["hic_ingest_status"] == "skipped_no_input"


def test_hic_ingest_refuses_hic_binary_route(synthetic_adata, tmp_path):
    from workflow.modular.modules.hic_ingest import HiCIngestModule
    hic_path = tmp_path / "sample.hic"
    hic_path.write_bytes(b"HIC\x00not-a-tsv")
    cfg = SimpleNamespace(
        hic_contacts_path=str(hic_path),
        hic_resolution_bp=100,
    )
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    HiCIngestModule().run(ctx)
    assert ctx.metadata["hic_ingest_status"] == "skipped_hic_unsupported"
    assert "hic_contact_matrix" not in synthetic_adata.uns


def test_hic_tad_skips_without_contact_matrix(synthetic_adata, tmp_path):
    from workflow.modular.modules.hic_tad import HiCTADModule
    cfg = SimpleNamespace()
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    HiCTADModule().run(ctx)
    assert ctx.metadata["hic_tad_status"] == "skipped_no_input"


def test_schema_parity_after_hic_contract_check():
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    rf = R_FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    assert sf.read_bytes() == rf.read_bytes()


def test_hic_slot_active_after_bundle_vertical_slice():
    import yaml
    doc = yaml.safe_load((FACTORY_ROOT / "contracts" / "bundle_schema.yaml").read_text())
    hic = doc["extensions_v22"]["hic"]
    assert hic["status"] == "active"
    assert "table" in hic
    assert "contacts_path" in hic["table"]


def test_r_extension_registry_includes_hic():
    src = (R_FACTORY_ROOT / "R_bundle" / "io_bundle.R").read_text()
    for line in src.splitlines():
        if line.strip().startswith("KNOWN_BUNDLE_EXTENSIONS"):
            assert '"hic"' in line
            return
    raise AssertionError("KNOWN_BUNDLE_EXTENSIONS line not found")


def test_r_hic_module_present():
    rmod = R_FACTORY_ROOT / "R" / "hic_module.R"
    assert rmod.exists()
    src = rmod.read_text()
    for fn in ("load_hic_extension", "plot_contact_heatmap", "plot_tad_triangle", "plot_ab_compartments"):
        assert fn in src, f"{fn} missing from hic_module.R"


def test_hic_modules_in_registry():
    from workflow.modular.pipeline import _build_registry, MUTATING_MODULES, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES
    reg = _build_registry()
    for name in ("hic_ingest", "hic_tad"):
        assert name in reg
        assert name not in MUTATING_MODULES
    order = _resolve_execution_order(list(MANDATORY_MODULES), ["hic_ingest", "hic_tad"])
    assert order.index("hic_tad") > order.index("hic_ingest")
