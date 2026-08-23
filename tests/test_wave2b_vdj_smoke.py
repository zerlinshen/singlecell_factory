"""Wave 2B / US-011 — VDJ vertical slice integration smoke.

Covers:
  - vdj_ingest module on a synthetic cellranger-vdj contigs CSV
  - vdj_metrics computes Shannon + Gini + clonal_expansion classes
  - Registry includes vdj_ingest + vdj_metrics, no MUTATING entry, no cycle
  - Bundle schema keeps vdj.status reserved until exporter+reader are complete
  - R-side KNOWN_BUNDLE_EXTENSIONS contains "vdj"
  - r_multiomics_factory/R/vdj_module.R has the required functions
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


FACTORY_ROOT = Path(__file__).resolve().parent.parent


def _resolve_r_factory_root() -> Path:
    candidates = [
        FACTORY_ROOT.parent / "r_multiomics_factory",
        FACTORY_ROOT.parents[2] / "r_multiomics_factory",
        Path("/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory"),
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


R_FACTORY_ROOT = _resolve_r_factory_root()


@pytest.fixture
def synthetic_adata_with_barcodes():
    """Small adata with realistic 10x-style barcodes (so vdj_ingest can map onto them)."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 80, 30
    X = sp.random(n_cells, n_genes, density=0.3, dtype=np.float32, random_state=42).tocsr()
    barcodes = [f"AAAC{i:04d}-1" for i in range(n_cells)]
    obs = pd.DataFrame(
        {
            "leiden": rng.integers(0, 4, size=n_cells).astype(str),
            "sample": [f"S{i % 4}" for i in range(n_cells)],
        },
        index=barcodes,
    )
    var = pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def synthetic_vdj_contigs(tmp_path, synthetic_adata_with_barcodes):
    """Build a small cellranger-vdj-like filtered_contig_annotations.csv.

    Cells 0-30 get TRA+TRB pairs (split into 5 distinct clonotypes by shared CDR3).
    Cells 40-50 get only TRA (orphan).
    Cells 60-70 get IGH+IGK (B cell pattern).
    """
    barcodes = list(synthetic_adata_with_barcodes.obs_names)
    rows = []

    # Clonotype 1: cells 0-9, TRA+TRB
    for bc in barcodes[0:10]:
        rows.append({"barcode": bc, "chain": "TRA", "v_gene": "TRAV1-1", "j_gene": "TRAJ12",
                     "cdr3_nt": "ACGTACGT", "cdr3": "AAB", "is_cell": "true", "productive": "true"})
        rows.append({"barcode": bc, "chain": "TRB", "v_gene": "TRBV1-1", "j_gene": "TRBJ1-1",
                     "cdr3_nt": "GTCAGTCA", "cdr3": "CCD", "is_cell": "true", "productive": "true"})
    # Clonotype 2: cells 10-20, TRA+TRB (different CDR3)
    for bc in barcodes[10:21]:
        rows.append({"barcode": bc, "chain": "TRA", "v_gene": "TRAV2", "j_gene": "TRAJ20",
                     "cdr3_nt": "TTTAGGAA", "cdr3": "EEF", "is_cell": "true", "productive": "true"})
        rows.append({"barcode": bc, "chain": "TRB", "v_gene": "TRBV2", "j_gene": "TRBJ2-1",
                     "cdr3_nt": "CCAATTGG", "cdr3": "GGH", "is_cell": "true", "productive": "true"})
    # Clonotype 3: cells 21-30, TRA+TRB (third CDR3)
    for bc in barcodes[21:31]:
        rows.append({"barcode": bc, "chain": "TRA", "v_gene": "TRAV3", "j_gene": "TRAJ30",
                     "cdr3_nt": "AAACCCGG", "cdr3": "III", "is_cell": "true", "productive": "true"})
        rows.append({"barcode": bc, "chain": "TRB", "v_gene": "TRBV3", "j_gene": "TRBJ3",
                     "cdr3_nt": "TTGGAACC", "cdr3": "JJK", "is_cell": "true", "productive": "true"})
    # Orphan TRA only: cells 40-49
    for bc in barcodes[40:50]:
        rows.append({"barcode": bc, "chain": "TRA", "v_gene": "TRAV4", "j_gene": "TRAJ40",
                     "cdr3_nt": "ORPHANTA", "cdr3": "LLM", "is_cell": "true", "productive": "true"})
    # B cell pattern IGH+IGK: cells 60-69
    for bc in barcodes[60:70]:
        rows.append({"barcode": bc, "chain": "IGH", "v_gene": "IGHV1", "j_gene": "IGHJ1",
                     "cdr3_nt": "BCELLIGH", "cdr3": "NNN", "is_cell": "true", "productive": "true"})
        rows.append({"barcode": bc, "chain": "IGK", "v_gene": "IGKV1", "j_gene": "IGKJ1",
                     "cdr3_nt": "BCELLIGK", "cdr3": "OOO", "is_cell": "true", "productive": "true"})

    df = pd.DataFrame(rows)
    contigs_path = tmp_path / "filtered_contig_annotations.csv"
    df.to_csv(contigs_path, index=False)
    return contigs_path


def test_vdj_ingest_writes_obs_and_uns(synthetic_adata_with_barcodes, synthetic_vdj_contigs, tmp_path):
    from workflow.modular.modules.vdj_ingest import VDJIngestModule

    cfg = SimpleNamespace(vdj_contigs_path=str(synthetic_vdj_contigs))
    ctx = SimpleNamespace(
        adata=synthetic_adata_with_barcodes, cfg=cfg, random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    VDJIngestModule().run(ctx)
    adata = synthetic_adata_with_barcodes
    assert "clonotype_id" in adata.obs.columns
    assert "chain_pairing" in adata.obs.columns
    assert "vdj_clonotypes" in adata.uns
    # Cells 0-30 are 3 TRA+TRB clonotypes; cells 40-49 are orphan TRA; cells 60-69 are IGH+IGK
    pairing = adata.obs["chain_pairing"].value_counts().to_dict()
    assert pairing.get("TRA+TRB", 0) == 31, f"expected 31 TRA+TRB cells, got {pairing}"
    assert pairing.get("orphan_TRA", 0) == 10, f"expected 10 orphan_TRA, got {pairing}"
    assert pairing.get("IGH+IGK", 0) == 10, f"expected 10 IGH+IGK, got {pairing}"
    # No-VDJ cells fill the rest
    assert pairing.get("no_chains", 0) == 80 - 31 - 10 - 10
    # Clonotype summary has at least 5 distinct clonotypes
    clonotypes_df = adata.uns["vdj_clonotypes"]
    assert clonotypes_df.shape[0] >= 4
    # Summary JSON written
    assert (tmp_path / "vdj_ingest" / "vdj_summary.json").exists()


def test_vdj_metrics_computes_shannon_gini_expansion(synthetic_adata_with_barcodes, synthetic_vdj_contigs, tmp_path):
    from workflow.modular.modules.vdj_ingest import VDJIngestModule
    from workflow.modular.modules.vdj_metrics import VDJMetricsModule

    cfg = SimpleNamespace(vdj_contigs_path=str(synthetic_vdj_contigs))
    ctx = SimpleNamespace(
        adata=synthetic_adata_with_barcodes, cfg=cfg, random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    VDJIngestModule().run(ctx)
    VDJMetricsModule().run(ctx)
    adata = synthetic_adata_with_barcodes
    assert "clonal_expansion" in adata.obs.columns
    assert "vdj_diversity" in adata.uns
    div = adata.uns["vdj_diversity"]
    assert list(div.columns) == ["sample", "n_cells_with_vdj", "n_clonotypes", "shannon", "gini"]
    # With 4 samples (S0..S3), 10-11-cell clonotypes are split sample-locally into 2-3 cells per sample
    expansion_counts = adata.obs["clonal_expansion"].value_counts().to_dict()
    assert "small_2-5" in expansion_counts
    assert "no_vdj" in expansion_counts
    assert "medium_6-20" not in expansion_counts  # proves sample-local counting (global count was 10-11)
    # Shannon entropy must be positive on a multi-clonotype sample
    assert (div["shannon"] > 0).any()
    # CSV and summary JSON written
    assert (tmp_path / "vdj_metrics" / "vdj_diversity.csv").exists()
    summary_path = tmp_path / "vdj_metrics" / "vdj_diversity.summary.json"
    assert summary_path.exists()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["provenance"]["expansion_unit"] == "obs.sample"
    assert summary["provenance"]["sample_local_clonotypes"] is True
    assert ctx.metadata["vdj_expansion_unit"] == "obs.sample"
    assert ctx.metadata["vdj_sample_local_clonotypes"] is True


def test_vdj_metrics_two_samples_same_clonotype_are_singletons(tmp_path):
    """Counterexample proof: cellranger clonotype IDs with same string across 2 samples are singletons."""
    from workflow.modular.modules.vdj_metrics import VDJMetricsModule

    obs = pd.DataFrame(
        {
            "sample": ["S1", "S2"],
            "clonotype_id": ["clonotype1", "clonotype1"],
        },
        index=["cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((2, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    ctx = SimpleNamespace(
        adata=adata, cfg=SimpleNamespace(), random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    VDJMetricsModule().run(ctx)
    # Both cells must be singleton (not small_2-5 from global pooling)
    assert adata.obs.loc["cell_1", "clonal_expansion"] == "singleton"
    assert adata.obs.loc["cell_2", "clonal_expansion"] == "singleton"

    div = adata.uns["vdj_diversity"]
    assert div.shape[0] == 2
    assert sorted(div["sample"].tolist()) == ["S1", "S2"]
    assert div["n_cells_with_vdj"].tolist() == [1, 1]
    assert div["n_clonotypes"].tolist() == [1, 1]
    assert ctx.metadata["vdj_expansion_unit"] == "obs.sample"
    assert ctx.metadata["vdj_sample_local_clonotypes"] is True


def test_vdj_metrics_sample_absent_fallback_all(tmp_path):
    """Fallback proof: when obs.sample is absent, fallback to _ALL_ and count globally."""
    from workflow.modular.modules.vdj_metrics import VDJMetricsModule

    obs = pd.DataFrame(
        {
            "clonotype_id": ["clonotype1", "clonotype1"],
        },
        index=["cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((2, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    ctx = SimpleNamespace(
        adata=adata, cfg=SimpleNamespace(), random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    VDJMetricsModule().run(ctx)
    # Both cells pool into small_2-5
    assert adata.obs.loc["cell_1", "clonal_expansion"] == "small_2-5"
    assert adata.obs.loc["cell_2", "clonal_expansion"] == "small_2-5"

    div = adata.uns["vdj_diversity"]
    assert div.shape[0] == 1
    assert div.loc[0, "sample"] == "_ALL_"
    assert div.loc[0, "n_cells_with_vdj"] == 2
    assert div.loc[0, "n_clonotypes"] == 1
    assert ctx.metadata["vdj_expansion_unit"] == "_ALL_"
    assert ctx.metadata["vdj_fallback_all_repertoire"] is True


def test_vdj_metrics_fails_loud_on_blank_sample_for_vdj_cells(tmp_path):
    """Fail-closed proof: VDJ-bearing cell with blank/null sample raises ValueError."""
    from workflow.modular.modules.vdj_metrics import VDJMetricsModule

    obs = pd.DataFrame(
        {
            "sample": ["S1", ""],
            "clonotype_id": ["clonotype1", "clonotype1"],
        },
        index=["cell_1", "cell_2"],
    )
    adata = ad.AnnData(X=np.zeros((2, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    ctx = SimpleNamespace(
        adata=adata, cfg=SimpleNamespace(), random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    with pytest.raises(ValueError, match="VDJ-bearing cells contain null, NaN, or blank sample identifiers"):
        VDJMetricsModule().run(ctx)


def test_vdj_metrics_skips_without_clonotype_column(synthetic_adata_with_barcodes, tmp_path):
    from workflow.modular.modules.vdj_metrics import VDJMetricsModule
    cfg = SimpleNamespace()
    ctx = SimpleNamespace(
        adata=synthetic_adata_with_barcodes, cfg=cfg, random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    # Don't run vdj_ingest first — column absent
    VDJMetricsModule().run(ctx)
    assert ctx.metadata.get("vdj_metrics_status") == "skipped_no_clonotypes"


def test_vdj_ingest_skips_without_contigs(synthetic_adata_with_barcodes, tmp_path):
    from workflow.modular.modules.vdj_ingest import VDJIngestModule
    cfg = SimpleNamespace()
    ctx = SimpleNamespace(
        adata=synthetic_adata_with_barcodes, cfg=cfg, random_state=42,
        run_dir=tmp_path, metadata={}, status=lambda *a, **k: None,
    )
    VDJIngestModule().run(ctx)
    assert ctx.metadata.get("vdj_ingest_status") == "skipped_no_input"
    assert "clonotype_id" not in synthetic_adata_with_barcodes.obs.columns


def test_schema_parity_after_vdj_contract_check():
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    rf = R_FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    assert sf.read_bytes() == rf.read_bytes(), "schema parity broken after vdj flip"


def test_vdj_slot_reserved_until_bundle_vertical_slice():
    import yaml
    sf = FACTORY_ROOT / "contracts" / "bundle_schema.yaml"
    doc = yaml.safe_load(sf.read_text())
    vdj = doc["extensions_v22"]["vdj"]
    assert vdj["status"] == "reserved"
    assert "table" in vdj
    assert "clonotypes_path" in vdj["table"]


def test_r_extension_registry_includes_vdj():
    io_bundle = R_FACTORY_ROOT / "R_bundle" / "io_bundle.R"
    src = io_bundle.read_text()
    for line in src.splitlines():
        if line.strip().startswith("KNOWN_BUNDLE_EXTENSIONS"):
            assert '"vdj"' in line, f"vdj not in registry line: {line}"
            return
    raise AssertionError("KNOWN_BUNDLE_EXTENSIONS line not found")


def test_r_vdj_module_present():
    rmod = R_FACTORY_ROOT / "R" / "vdj_module.R"
    assert rmod.exists()
    src = rmod.read_text()
    for fn in ("load_vdj_extension", "plot_clonotype_expansion", "plot_clonotype_chord", "plot_vdj_diversity"):
        assert fn in src, f"{fn} missing from vdj_module.R"


def test_vdj_modules_in_registry():
    from workflow.modular.pipeline import _build_registry, MUTATING_MODULES, _resolve_execution_order
    from workflow.modular.module_catalog import MANDATORY_MODULES
    reg = _build_registry()
    for name in ("vdj_ingest", "vdj_metrics"):
        assert name in reg
        assert name not in MUTATING_MODULES
    # Topo: vdj_metrics depends on vdj_ingest
    optional = ["vdj_ingest", "vdj_metrics"]
    order = _resolve_execution_order(list(MANDATORY_MODULES), optional)
    assert order.index("vdj_metrics") > order.index("vdj_ingest")
