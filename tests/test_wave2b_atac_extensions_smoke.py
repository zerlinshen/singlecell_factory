"""Wave 2B / US-014 — ATAC extensions smoke (atac_qc + peak_to_gene).

Covers:
  - atac_qc reports TSS enrichment + FRiP per cell with ENCODE-style pass/warn/fail
  - atac_qc fail-soft skips when --atac-fragments-path absent
  - peak_to_gene maps peaks to nearest TSS within window
  - peak_to_gene fail-soft skips when uns['atac_peaks'] absent
  - Registry includes both modules, no cycle, no MUTATING entry
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


@pytest.fixture
def synthetic_adata():
    rng = np.random.default_rng(42)
    n_cells, n_genes = 50, 20
    X = sp.random(n_cells, n_genes, density=0.3, dtype=np.float32, random_state=42).tocsr()
    barcodes = [f"BCQC{i:04d}-1" for i in range(n_cells)]
    obs = pd.DataFrame({"leiden": rng.integers(0, 3, size=n_cells).astype(str)}, index=barcodes)
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))


@pytest.fixture
def synthetic_fragments(tmp_path, synthetic_adata):
    """fragments.tsv.gz — 80% of cells have ATAC fragments concentrated near TSSs."""
    barcodes = list(synthetic_adata.obs_names)
    rng = np.random.default_rng(42)
    lines = []
    for bc in barcodes[:40]:  # 80% with fragments
        for _ in range(rng.integers(30, 60)):
            chrom = "chr1"
            # Cluster fragments near TSS at 1000 (within ±1000 bp)
            tss_mid = rng.choice([1000, 5000, 10000])
            offset = rng.integers(-1500, 1500)
            start = max(0, tss_mid + offset)
            end = start + rng.integers(50, 200)
            lines.append(f"{chrom}\t{start}\t{end}\t{bc}")
    frag_path = tmp_path / "fragments.tsv.gz"
    with gzip.open(frag_path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return frag_path


@pytest.fixture
def synthetic_tss_bed(tmp_path):
    bed = tmp_path / "tss.bed"
    bed.write_text("chr1\t900\t1100\nchr1\t4900\t5100\nchr1\t9900\t10100\n")
    return bed


@pytest.fixture
def synthetic_peaks_bed(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t950\t1050\nchr1\t4950\t5050\nchr1\t9950\t10050\nchr2\t1000\t1100\n")
    return bed


@pytest.fixture
def synthetic_gene_tss_bed(tmp_path):
    bed = tmp_path / "genes.bed"
    bed.write_text(
        "chr1\t1000\t1001\tGENE_A\nchr1\t5000\t5001\tGENE_B\nchr1\t10000\t10001\tGENE_C\nchr2\t1000\t1001\tGENE_D\n"
    )
    return bed


def test_atac_qc_skips_without_fragments(synthetic_adata, tmp_path):
    from workflow.modular.modules.atac_qc import ATACQCModule
    cfg = SimpleNamespace()
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    ATACQCModule().run(ctx)
    assert ctx.metadata["atac_qc_status_summary"] == "skipped_no_fragments"
    assert "atac_qc_status" not in synthetic_adata.obs.columns


def test_atac_qc_computes_frip_tss_status(synthetic_adata, synthetic_fragments, synthetic_tss_bed, synthetic_peaks_bed, tmp_path):
    from workflow.modular.modules.atac_qc import ATACQCModule
    cfg = SimpleNamespace(
        atac_fragments_path=str(synthetic_fragments),
        atac_tss_bed_path=str(synthetic_tss_bed),
        atac_peaks_bed_path=str(synthetic_peaks_bed),
    )
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    ATACQCModule().run(ctx)
    adata = synthetic_adata
    for col in ("frip", "tss_enrichment", "atac_qc_status"):
        assert col in adata.obs.columns, f"missing column: {col}"
    assert ctx.metadata["atac_qc_status_summary"] == "ok"
    # First 40 cells have fragments; remaining 10 should be no_data
    n_no_data = int((adata.obs["atac_qc_status"] == "no_data").sum())
    assert n_no_data == 10
    summary = adata.uns["atac_qc_summary"]
    assert "thresholds" in summary
    # JSON written
    assert (tmp_path / "atac_qc" / "atac_qc.json").exists()


def test_peak_to_gene_skips_without_peaks(synthetic_adata, synthetic_gene_tss_bed, tmp_path):
    from workflow.modular.modules.peak_to_gene import PeakToGeneModule
    cfg = SimpleNamespace(gene_tss_bed_path=str(synthetic_gene_tss_bed))
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    PeakToGeneModule().run(ctx)
    assert ctx.metadata["peak_to_gene_status"] == "skipped_no_peaks"


def test_peak_to_gene_links_peaks_to_genes(synthetic_adata, synthetic_gene_tss_bed, tmp_path):
    from workflow.modular.modules.peak_to_gene import PeakToGeneModule

    # Seed adata.uns["atac_peaks"] with peaks adjacent to GENE_A and GENE_B
    synthetic_adata.uns["atac_peaks"] = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr2"],
        "start": [950, 4950, 9950, 1000],
        "end":   [1050, 5050, 10050, 1100],
    })
    synthetic_adata.obsm["atac_peaks"] = sp.random(
        synthetic_adata.n_obs,
        4,
        density=0.2,
        format="csr",
        dtype=np.float32,
        random_state=7,
    )
    cfg = SimpleNamespace(
        gene_tss_bed_path=str(synthetic_gene_tss_bed),
        peak_to_gene_window=5_000,
    )
    ctx = SimpleNamespace(adata=synthetic_adata, cfg=cfg, random_state=42,
                          run_dir=tmp_path, metadata={}, status=lambda *a, **k: None)
    PeakToGeneModule().run(ctx)
    link = synthetic_adata.uns["peak_to_gene"]
    # All 4 peaks should link to their respective genes within the window
    assert link.shape[0] == 4
    assert set(link["gene"]) == {"GENE_A", "GENE_B", "GENE_C", "GENE_D"}
    assert (link["dist_bp"] <= 100).all()
    # CSV written
    assert (tmp_path / "peak_to_gene" / "peak_to_gene.csv").exists()


def test_atac_extensions_in_registry():
    from workflow.modular.pipeline import _build_registry, MUTATING_MODULES, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES
    reg = _build_registry()
    for name in ("atac_qc", "peak_to_gene"):
        assert name in reg
        assert name not in MUTATING_MODULES
    # peak_to_gene must come after atac_ingest in topo
    optional = ["atac_ingest", "atac_qc", "peak_to_gene"]
    order = _resolve_execution_order(list(MANDATORY_MODULES), optional)
    assert order.index("peak_to_gene") > order.index("atac_ingest")


def test_atac_ingest_lsi_peak_to_gene_chain(synthetic_adata, synthetic_gene_tss_bed, tmp_path):
    """Synthetic ATAC chain uses one canonical sparse peak matrix contract."""
    import scipy.io as sio

    from workflow.modular.modules.atac_ingest import ATACIngestModule
    from workflow.modular.modules.atac_lsi import ATACLSIModule
    from workflow.modular.modules.peak_to_gene import PeakToGeneModule

    n_cells = synthetic_adata.n_obs
    n_peaks = 8
    peak_mat = sp.random(
        n_cells,
        n_peaks,
        density=0.25,
        format="csr",
        dtype=np.float32,
        random_state=42,
    )
    peak_mtx_path = tmp_path / "peaks.mtx"
    sio.mmwrite(str(peak_mtx_path), peak_mat)

    peaks_bed = tmp_path / "peaks.bed"
    peaks_bed.write_text(
        "\n".join([
            "chr1\t950\t1050",
            "chr1\t4950\t5050",
            "chr1\t9950\t10050",
            "chr2\t1000\t1100",
            "chr1\t20000\t20100",
            "chr1\t25000\t25100",
            "chr2\t30000\t30100",
            "chr2\t35000\t35100",
        ])
        + "\n"
    )

    cfg = SimpleNamespace(
        atac_peak_matrix_path=str(peak_mtx_path),
        atac_peaks_bed_path=str(peaks_bed),
        atac_n_components=4,
        gene_tss_bed_path=str(synthetic_gene_tss_bed),
        peak_to_gene_window=5_000,
    )
    ctx = SimpleNamespace(
        adata=synthetic_adata,
        cfg=cfg,
        random_state=42,
        run_dir=tmp_path,
        figure_dir=tmp_path,
        metadata={},
        status=lambda *a, **k: None,
    )

    ATACIngestModule().run(ctx)
    assert sp.issparse(synthetic_adata.obsm["atac_peaks"])
    assert synthetic_adata.obsm["atac_peaks"].shape == (n_cells, n_peaks)
    assert "atac_var" in synthetic_adata.uns

    ATACLSIModule().run(ctx)
    assert "X_lsi" in synthetic_adata.obsm
    assert synthetic_adata.obsm["X_lsi"].shape[0] == n_cells

    PeakToGeneModule().run(ctx)
    assert "peak_to_gene" in synthetic_adata.uns
    assert "peak_to_gene_linkages" in synthetic_adata.uns
