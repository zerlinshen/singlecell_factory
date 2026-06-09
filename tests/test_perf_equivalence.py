"""Performance-fix equivalence tests.

Each test pins a surgical speedup against a reference implementation of the
*original* (pre-optimization) code path and asserts numerically identical
output. These are speedups, not behavior changes — the bar is exact equality
(or float32-rounding-level equality where a float64-vs-float32 reduction order
legitimately differs, documented per test).

Covers:
  #22 hic_ingest._binize_contacts / _build_bins  — vectorized vs iterrows
  #23 cell_communication                          — sparse vs dense mean engine
  #25 spatial_neighborhoods._select_morans_genes  — sparse-safe vs dense mean
  #26 export_singlecell_r_bundle.write_marker_expr_parquet — single-slice vs chunk+concat
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

ad = pytest.importorskip("anndata")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


# ---------------------------------------------------------------------------
# #22 — Hi-C contact binning: vectorized COO must equal iterrows accumulation
# ---------------------------------------------------------------------------

def _ref_build_bins(chromsizes: pd.DataFrame, resolution: int) -> pd.DataFrame:
    """Original row-wise bin builder (pre-#22)."""
    rows = []
    bin_id = 0
    for _, r in chromsizes.iterrows():
        starts = np.arange(0, int(r["length"]), resolution, dtype=int)
        for s in starts:
            rows.append({
                "bin_id": bin_id,
                "chrom": r["chrom"],
                "start": int(s),
                "end": int(min(s + resolution, r["length"])),
            })
            bin_id += 1
    return pd.DataFrame(rows)


def _ref_binize_contacts(contacts: pd.DataFrame, bins: pd.DataFrame, resolution: int) -> sp.csr_matrix:
    """Original iterrows accumulation (pre-#22)."""
    bin_lookup: dict[tuple[str, int], int] = {}
    for _, r in bins.iterrows():
        bin_lookup[(r["chrom"], int(r["start"]) // resolution)] = int(r["bin_id"])

    n_bins = bins.shape[0]
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for _, c in contacts.iterrows():
        b1 = bin_lookup.get((c["chrom1"], int(c["pos1"]) // resolution))
        b2 = bin_lookup.get((c["chrom2"], int(c["pos2"]) // resolution))
        if b1 is None or b2 is None:
            continue
        rows.append(b1)
        cols.append(b2)
        data.append(float(c["count"]))
        if b1 != b2:
            rows.append(b2)
            cols.append(b1)
            data.append(float(c["count"]))
    return sp.csr_matrix((data, (rows, cols)), shape=(n_bins, n_bins), dtype=np.float32)


def _make_contacts(seed: int = 7) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Multi-chrom contacts with duplicate pairs, off-bin positions, and an
    out-of-range contact that must be dropped by the bin lookup."""
    resolution = 100
    chromsizes = pd.DataFrame({"chrom": ["chr1", "chr2"], "length": [500, 300]})
    rng = np.random.default_rng(seed)
    recs = []
    for _ in range(200):
        c1 = "chr1" if rng.random() < 0.6 else "chr2"
        c2 = "chr1" if rng.random() < 0.6 else "chr2"
        lim1 = 500 if c1 == "chr1" else 300
        lim2 = 500 if c2 == "chr1" else 300
        p1 = int(rng.integers(0, lim1))
        p2 = int(rng.integers(0, lim2))
        cnt = float(rng.integers(1, 5))
        recs.append((c1, p1, c2, p2, cnt))
    # Force duplicate (bin1,bin2) pairs to exercise duplicate-summing.
    recs.append(("chr1", 10, "chr1", 250, 3.0))
    recs.append(("chr1", 40, "chr1", 230, 2.0))   # same bins (0, 2) as above
    recs.append(("chr1", 10, "chr1", 10, 5.0))    # diagonal, no mirror
    recs.append(("chr1", 60, "chr1", 60, 1.0))    # diagonal same bin again
    # Out-of-range contact (no matching bin) -> must be skipped by both impls.
    recs.append(("chrX", 10, "chr1", 10, 9.0))
    contacts = pd.DataFrame(recs, columns=["chrom1", "pos1", "chrom2", "pos2", "count"]).astype(
        {"chrom1": str, "pos1": int, "chrom2": str, "pos2": int, "count": float}
    )
    bins = _ref_build_bins(chromsizes, resolution)
    return contacts, bins, resolution


def test_hic_build_bins_vectorized_equals_iterrows():
    from workflow.modular.modules.hic_ingest import _build_bins

    chromsizes = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3"], "length": [550, 300, 99]})
    resolution = 100
    ref = _ref_build_bins(chromsizes, resolution).reset_index(drop=True)
    opt = _build_bins(chromsizes, resolution).reset_index(drop=True)
    # Same columns, order, dtypes, values.
    pd.testing.assert_frame_equal(opt, ref, check_dtype=True)


def test_hic_binize_vectorized_equals_iterrows():
    from workflow.modular.modules.hic_ingest import _binize_contacts

    contacts, bins, resolution = _make_contacts()
    ref = _ref_binize_contacts(contacts, bins, resolution)
    opt = _binize_contacts(contacts, bins, resolution)

    assert opt.dtype == ref.dtype == np.float32
    assert opt.shape == ref.shape
    # Exact structural + value equality (symmetry, duplicate-summing, diagonal).
    diff = (opt - ref)
    assert diff.nnz == 0, f"matrices differ in {diff.nnz} entries"
    a = opt.toarray()
    b = ref.toarray()
    max_abs = float(np.max(np.abs(a - b))) if a.size else 0.0
    assert max_abs == 0.0, f"max abs diff {max_abs}"
    # Sanity: matrix is symmetric exactly as the original built it.
    assert (opt - opt.T).nnz == 0


def test_hic_binize_empty_contacts():
    from workflow.modular.modules.hic_ingest import _binize_contacts

    chromsizes = pd.DataFrame({"chrom": ["chr1"], "length": [300]})
    bins = _ref_build_bins(chromsizes, 100)
    empty = pd.DataFrame(
        {"chrom1": [], "pos1": [], "chrom2": [], "pos2": [], "count": []}
    ).astype({"chrom1": str, "pos1": int, "chrom2": str, "pos2": int, "count": float})
    ref = _ref_binize_contacts(empty, bins, 100)
    opt = _binize_contacts(empty, bins, 100)
    assert opt.shape == ref.shape == (bins.shape[0], bins.shape[0])
    assert opt.nnz == ref.nnz == 0
    assert opt.dtype == np.float32


# ---------------------------------------------------------------------------
# #23 — cell_communication: sparse engine (new default) == dense engine
# ---------------------------------------------------------------------------

def _make_cc_adata(seed: int = 3):
    rng = np.random.default_rng(seed)
    genes = [
        "CD274", "PDCD1", "CXCL12", "CXCR4", "VEGFA", "KDR", "TGFB1", "TGFBR2",
        "IL6", "IL6R", "CCL5", "CCR5", "FILLER1", "FILLER2",
    ]
    n_cells = 60
    dense = rng.random((n_cells, len(genes))).astype(np.float32)
    dense[dense < 0.2] = 0.0  # inject sparsity
    X = sp.csr_matrix(dense)
    # 3 cell types, each with >= 5 cells, in a deterministic order.
    cts = (["T"] * 20) + (["Myeloid"] * 20) + (["Tumor"] * 20)
    obs = pd.DataFrame({"cell_type": cts}, index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=genes)
    return ad.AnnData(X=X, obs=obs, var=var)


def _run_cc(adata, engine: str, tmp_path: Path) -> pd.DataFrame:
    import os
    from workflow.modular.modules.cell_communication import CellCommunicationModule

    table_dir = tmp_path / engine
    table_dir.mkdir(parents=True, exist_ok=True)
    ctx = SimpleNamespace(
        adata=adata.copy(),
        cfg=SimpleNamespace(),
        metadata={},
        table_dir=table_dir,
        figure_dir=table_dir,
        status=lambda *a, **k: None,
    )
    prev = os.environ.get("SC_CELLCOMM_ENGINE")
    os.environ["SC_CELLCOMM_ENGINE"] = engine
    try:
        # Force the manual L-R fallback (no LIANA in the equivalence harness).
        CellCommunicationModule()._run_manual_lr(ctx.adata, ctx)
    finally:
        if prev is None:
            os.environ.pop("SC_CELLCOMM_ENGINE", None)
        else:
            os.environ["SC_CELLCOMM_ENGINE"] = prev
    return pd.read_csv(table_dir / "cell_communication_lr.csv")


def test_cellcomm_sparse_equals_dense(tmp_path):
    adata = _make_cc_adata()
    dense_df = _run_cc(adata, "dense", tmp_path)
    sparse_df = _run_cc(adata, "sparse", tmp_path)

    # Same rows/columns and (after deterministic sort) identical scored output.
    assert list(dense_df.columns) == list(sparse_df.columns)
    sort_cols = ["ligand", "receptor", "source", "target"]
    d = dense_df.sort_values(sort_cols).reset_index(drop=True)
    s = sparse_df.sort_values(sort_cols).reset_index(drop=True)
    assert d.shape == s.shape
    # Label columns must match exactly.
    for col in sort_cols:
        assert (d[col] == s[col]).all(), f"label mismatch in {col}"
    # Numeric columns: the CSV is rounded to 4 decimals; float64-mean (sparse)
    # vs float32-mean (dense) of the same data agree to that precision.
    for col in ["ligand_expr", "receptor_expr", "lr_score"]:
        max_abs = float(np.max(np.abs(d[col].to_numpy() - s[col].to_numpy())))
        assert max_abs <= 1e-4, f"{col} max abs diff {max_abs} exceeds rounding tolerance"


def test_cellcomm_auto_default_uses_sparse_for_sparse_X(tmp_path, monkeypatch):
    """With no env override (auto), sparse X must take the single-pass sparse
    path — i.e. produce output identical to an explicit sparse run."""
    import os
    monkeypatch.delenv("SC_CELLCOMM_ENGINE", raising=False)
    adata = _make_cc_adata()

    from workflow.modular.modules.cell_communication import CellCommunicationModule

    table_dir = tmp_path / "auto"
    table_dir.mkdir(parents=True, exist_ok=True)
    ctx = SimpleNamespace(
        adata=adata.copy(), cfg=SimpleNamespace(), metadata={},
        table_dir=table_dir, figure_dir=table_dir, status=lambda *a, **k: None,
    )
    assert os.environ.get("SC_CELLCOMM_ENGINE") is None
    CellCommunicationModule()._run_manual_lr(ctx.adata, ctx)
    auto_df = pd.read_csv(table_dir / "cell_communication_lr.csv")
    sparse_df = _run_cc(adata, "sparse", tmp_path)
    sort_cols = ["ligand", "receptor", "source", "target"]
    a = auto_df.sort_values(sort_cols).reset_index(drop=True)
    s = sparse_df.sort_values(sort_cols).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, s)


# ---------------------------------------------------------------------------
# #25 — spatial_neighborhoods._select_morans_genes sparse-safe mean
# ---------------------------------------------------------------------------

def _ref_select_morans_genes_dense(adata, n_top_genes: int) -> list[str]:
    """Original implementation with the dense else-branch (pre-#25)."""
    if "highly_variable" in adata.var.columns:
        hv = adata.var.index[adata.var["highly_variable"].astype(bool)].tolist()
        if hv:
            return hv[:n_top_genes]
    try:
        from scipy import sparse

        if sparse.issparse(adata.X):
            means = np.asarray(adata.X.mean(axis=0)).ravel()
        else:
            means = np.asarray(adata.X).mean(axis=0)
    except Exception:
        return list(adata.var_names[:n_top_genes])
    order = np.argsort(-means)[:n_top_genes]
    return [str(adata.var_names[i]) for i in order]


def test_spatial_select_genes_sparse_equals_dense_reference():
    from workflow.modular.modules.spatial_neighborhoods import SpatialNeighborhoodsModule

    rng = np.random.default_rng(11)
    n_cells, n_genes = 40, 25
    dense = rng.random((n_cells, n_genes)).astype(np.float32)
    dense[dense < 0.3] = 0.0
    genes = [f"g{i}" for i in range(n_genes)]
    var = pd.DataFrame(index=genes)

    # Sparse X: optimized path vs original sparse branch -> identical.
    adata_sp = ad.AnnData(X=sp.csr_matrix(dense), var=var.copy())
    ref_sp = _ref_select_morans_genes_dense(adata_sp, 10)
    opt_sp = SpatialNeighborhoodsModule._select_morans_genes(adata_sp, 10)
    assert opt_sp == ref_sp

    # Dense X: the dropped else-branch used np.asarray(X).mean(axis=0); the new
    # X.mean(axis=0) on an ndarray is the identical reduction -> same selection.
    adata_de = ad.AnnData(X=dense.copy(), var=var.copy())
    ref_de = _ref_select_morans_genes_dense(adata_de, 10)
    opt_de = SpatialNeighborhoodsModule._select_morans_genes(adata_de, 10)
    assert opt_de == ref_de
    # And dense vs sparse of the same data agree on the ordering.
    assert opt_de == opt_sp


# ---------------------------------------------------------------------------
# #26 — bundle marker export: single-slice == chunk+concat (byte-identical parquet)
# ---------------------------------------------------------------------------

def _ref_write_marker_expr_parquet_chunked(adata, cell_idx, markers_present, config):
    """Original chunk+concat materialization (pre-#26)."""
    from scripts.export_singlecell_r_bundle import (
        materialize_small_matrix, _write_parquet, file_record,
    )
    cells = adata.obs_names[cell_idx].astype(str)
    chunks = []
    for start in range(0, len(cell_idx), config.marker_chunk_size):
        stop = min(start + config.marker_chunk_size, len(cell_idx))
        rows = cell_idx[start:stop]
        block = materialize_small_matrix(adata[rows, markers_present].X)
        block = block.astype(np.float32, copy=False)
        chunks.append(pd.DataFrame(block, index=cells[start:stop], columns=markers_present))
    df = pd.concat(chunks, axis=0)
    df.index.name = "cell"
    path = config.output_dir / "marker_expr.parquet"
    _write_parquet(df, path)
    rec = file_record(path, config.output_dir, df.shape[0], df.shape[1], {"format": "parquet"})
    return path, rec


def _make_marker_h5ad(path: Path, n_cells: int = 120, n_genes: int = 8) -> list[str]:
    rng = np.random.default_rng(5)
    dense = rng.random((n_cells, n_genes)).astype(np.float32)
    dense[dense < 0.4] = 0.0
    X = sp.csr_matrix(dense)
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(n_cells)])
    genes = [f"GENE{i}" for i in range(n_genes)]
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.write_h5ad(path)
    return genes


def test_marker_expr_single_slice_equals_chunked_backed(tmp_path):
    from scripts.export_singlecell_r_bundle import (
        ExportConfig, write_marker_expr_parquet,
    )

    h5_path = tmp_path / "marker.h5ad"
    genes = _make_marker_h5ad(h5_path, n_cells=120, n_genes=8)
    markers_present = genes[:5]
    # Subset that crosses several chunk boundaries (chunk_size < n_cells).
    cell_idx = np.sort(np.random.default_rng(2).choice(120, size=90, replace=False))

    # --- reference: chunked, on a backed='r' adata (matches export_bundle) ---
    ref_dir = tmp_path / "ref"
    ref_dir.mkdir()
    ref_cfg = ExportConfig(
        input_h5ad=h5_path, output_dir=ref_dir, marker_chunk_size=25,
    )
    backed_ref = ad.read_h5ad(h5_path, backed="r")
    try:
        ref_path, ref_rec = _ref_write_marker_expr_parquet_chunked(
            backed_ref, cell_idx, markers_present, ref_cfg
        )
    finally:
        backed_ref.file.close()

    # --- optimized: single slice, on its own backed='r' adata ---
    opt_dir = tmp_path / "opt"
    opt_dir.mkdir()
    opt_cfg = ExportConfig(
        input_h5ad=h5_path, output_dir=opt_dir, marker_chunk_size=25,
    )
    backed_opt = ad.read_h5ad(h5_path, backed="r")
    try:
        opt_path, opt_rec = write_marker_expr_parquet(
            backed_opt, cell_idx, markers_present, opt_cfg
        )
    finally:
        backed_opt.file.close()

    # Byte-identical parquet payload (same schema, order, dtypes, values).
    ref_bytes = ref_path.read_bytes()
    opt_bytes = opt_path.read_bytes()
    assert ref_rec["sha256"] == opt_rec["sha256"], "marker_expr.parquet sha256 mismatch"
    assert ref_bytes == opt_bytes, "marker_expr.parquet bytes differ"

    # Content cross-check: same frame including the literal 'cell' index column.
    ref_df = pd.read_parquet(ref_path)
    opt_df = pd.read_parquet(opt_path)
    pd.testing.assert_frame_equal(opt_df, ref_df, check_dtype=True)
    assert list(opt_df.columns) == ["cell"] + markers_present
    assert opt_df["cell"].tolist() == [f"cell{i}" for i in cell_idx]
    max_abs = float(
        np.max(np.abs(opt_df[markers_present].to_numpy() - ref_df[markers_present].to_numpy()))
    )
    assert max_abs == 0.0, f"marker values differ by {max_abs}"
