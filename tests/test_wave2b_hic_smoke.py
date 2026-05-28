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
import json
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
    assert set(compartments["compartment"].unique()).issubset({"A", "B", "low_information"})
    assert ctx.metadata["hic_tad_status"] == "ok"
    assert "hic_tad_metadata" in synthetic_adata.uns
    assert synthetic_adata.uns["hic_tad_metadata"]["compartment_status"] in {
        "confident",
        "partial_low_information",
        "low_information",
    }


def test_hic_tad_marks_low_information_chromosomes_without_warnings(synthetic_adata, tmp_path):
    from workflow.modular.modules.hic_tad import HiCTADModule

    bins = pd.DataFrame(
        {
            "bin_id": np.arange(10, dtype=np.int64),
            "chrom": ["chr1"] * 5 + ["chr2"] * 5,
            "start": np.tile(np.arange(0, 5000, 1000, dtype=np.int64), 2),
            "end": np.tile(np.arange(1000, 6000, 1000, dtype=np.int64), 2),
        }
    )
    informative = np.array(
        [
            [8, 2, 1, 0, 0],
            [2, 7, 2, 1, 0],
            [1, 2, 9, 2, 1],
            [0, 1, 2, 6, 3],
            [0, 0, 1, 3, 5],
        ],
        dtype=np.float32,
    )
    zero_chrom = np.zeros((5, 5), dtype=np.float32)
    synthetic_adata.uns["hic_bins"] = bins
    synthetic_adata.uns["hic_contact_matrix"] = sp.block_diag((informative, zero_chrom), format="csr")

    cfg = SimpleNamespace(hic_tad_window_bins=1, hic_tad_boundary_k=0.5)
    ctx = SimpleNamespace(
        adata=synthetic_adata,
        cfg=cfg,
        random_state=42,
        run_dir=tmp_path,
        metadata={},
        status=lambda *a, **k: None,
    )
    HiCTADModule().run(ctx)

    meta = synthetic_adata.uns["hic_tad_metadata"]
    assert meta["compartment_status"] == "partial_low_information"
    assert meta["low_information_chromosomes"] == ["chr2"]
    assert meta["compartment_status_by_chrom"] == {"chr1": "confident", "chr2": "low_information"}
    assert ctx.metadata["hic_tad_metadata"] == meta

    compartments = synthetic_adata.uns["hic_compartments"]
    chr2_scores = compartments.loc[compartments["chrom"] == "chr2", "eigenvector_1"].to_numpy()
    assert np.allclose(chr2_scores, 0.0)
    assert set(compartments.loc[compartments["chrom"] == "chr2", "compartment"]) == {"low_information"}

    summary = json.loads((tmp_path / "hic_tad" / "tad_summary.json").read_text())
    assert summary["hic_tad_metadata"]["compartment_status"] == "partial_low_information"
    assert ctx.metadata["hic_tad_status"] == "ok"


def test_hic_tad_keeps_informative_chromosome_with_sparse_tail_bins(synthetic_adata, tmp_path):
    from workflow.modular.modules.hic_tad import HiCTADModule

    bins = pd.DataFrame(
        {
            "bin_id": np.arange(7, dtype=np.int64),
            "chrom": ["chr1"] * 7,
            "start": np.arange(0, 7000, 1000, dtype=np.int64),
            "end": np.arange(1000, 8000, 1000, dtype=np.int64),
        }
    )
    contact_matrix = np.array(
        [
            [9, 3, 1, 0, 0, 0, 0],
            [3, 8, 3, 1, 0, 0, 0],
            [1, 3, 8, 3, 1, 0, 0],
            [0, 1, 3, 9, 2, 0, 0],
            [0, 0, 1, 2, 7, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
        ],
        dtype=np.float32,
    )
    synthetic_adata.uns["hic_bins"] = bins
    synthetic_adata.uns["hic_contact_matrix"] = sp.csr_matrix(contact_matrix)

    cfg = SimpleNamespace(hic_tad_window_bins=1, hic_tad_boundary_k=0.5)
    ctx = SimpleNamespace(
        adata=synthetic_adata,
        cfg=cfg,
        random_state=42,
        run_dir=tmp_path,
        metadata={},
        status=lambda *a, **k: None,
    )
    HiCTADModule().run(ctx)

    meta = synthetic_adata.uns["hic_tad_metadata"]
    assert meta["compartment_status"] == "confident"
    assert meta["low_information_chromosomes"] == []

    compartments = synthetic_adata.uns["hic_compartments"]
    assert set(compartments.loc[compartments["bin_id"] < 5, "compartment"]) <= {"A", "B"}
    assert set(compartments.loc[compartments["bin_id"] >= 5, "compartment"]) == {"low_information"}
    assert compartments.loc[compartments["bin_id"] < 5, "eigenvector_1"].notna().all()
    assert compartments.loc[compartments["bin_id"] >= 5, "eigenvector_1"].isna().all()


def test_insulation_score_marks_empty_windows_nan():
    from workflow.modular.modules.hic_tad import _insulation_score

    mat = sp.csr_matrix((5, 5), dtype=np.float32)
    scores = _insulation_score(mat, window_bins=1)
    assert np.isnan(scores).all()


def test_hic_tad_marks_eigendecomposition_failure_low_information(synthetic_adata, tmp_path, monkeypatch):
    from workflow.modular.modules.hic_tad import HiCTADModule

    bins = pd.DataFrame(
        {
            "bin_id": np.arange(5, dtype=np.int64),
            "chrom": ["chr1"] * 5,
            "start": np.arange(0, 5000, 1000, dtype=np.int64),
            "end": np.arange(1000, 6000, 1000, dtype=np.int64),
        }
    )
    contact_matrix = np.array(
        [
            [8, 2, 1, 0, 0],
            [2, 7, 2, 1, 0],
            [1, 2, 9, 2, 1],
            [0, 1, 2, 6, 3],
            [0, 0, 1, 3, 5],
        ],
        dtype=np.float32,
    )
    synthetic_adata.uns["hic_bins"] = bins
    synthetic_adata.uns["hic_contact_matrix"] = sp.csr_matrix(contact_matrix)

    def raise_linalg_error(_corr):
        raise np.linalg.LinAlgError("forced eigendecomposition failure")

    monkeypatch.setattr(np.linalg, "eigh", raise_linalg_error)

    cfg = SimpleNamespace(hic_tad_window_bins=1, hic_tad_boundary_k=0.5)
    ctx = SimpleNamespace(
        adata=synthetic_adata,
        cfg=cfg,
        random_state=42,
        run_dir=tmp_path,
        metadata={},
        status=lambda *a, **k: None,
    )
    HiCTADModule().run(ctx)

    meta = synthetic_adata.uns["hic_tad_metadata"]
    assert meta["compartment_status"] == "low_information"
    assert meta["low_information_chromosomes"] == ["chr1"]
    assert meta["compartment_status_by_chrom"] == {"chr1": "low_information"}
    assert np.allclose(synthetic_adata.uns["hic_compartments"]["eigenvector_1"].to_numpy(), 0.0)
    assert set(synthetic_adata.uns["hic_compartments"]["compartment"]) == {"low_information"}


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
    assert hic["metadata"]["hic_tad_metadata"]["fields"]["compartment_status"]["values"] == [
        "confident",
        "partial_low_information",
        "low_information",
        "unknown_or_unvalidated",
    ]


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
    from workflow.modular.modules.hic_tad import HiCTADModule
    reg = _build_registry()
    for name in ("hic_ingest", "hic_tad"):
        assert name in reg
        assert name not in MUTATING_MODULES
    assert "hic_tad_metadata" in HiCTADModule.provides_keys["uns"]
    order = _resolve_execution_order(list(MANDATORY_MODULES), ["hic_ingest", "hic_tad"])
    assert order.index("hic_tad") > order.index("hic_ingest")
