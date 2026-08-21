from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip("anndata")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.export_singlecell_r_bundle import ExportConfig, export_bundle, maybe_export_ribo


def read_bundle_csv(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return pd.read_csv(handle, index_col=0)


def make_tiny_h5ad(
    path: Path,
    include_pca: bool = True,
    include_hic: bool = False,
    include_ribo: bool = False,
) -> None:
    x = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 2, 0],
                [3, 0, 0, 4],
                [0, 5, 0, 0],
                [6, 0, 7, 0],
            ],
            dtype=np.float32,
        )
    )
    obs = pd.DataFrame(
        {
            "cell_type": ["T", "Myeloid", "Tumor", "Tumor"],
            "leiden": ["0", "1", "1", "2"],
            "unused": ["a", "b", "c", "d"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    var = pd.DataFrame(index=["CD3E", "LYZ", "ELF3", "KRT8"])
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_umap"] = np.arange(8, dtype=np.float32).reshape(4, 2)
    if include_pca:
        adata.obsm["X_pca"] = np.arange(12, dtype=np.float32).reshape(4, 3)
    if include_hic:
        bins = pd.DataFrame(
            {
                "bin_id": [0, 1, 2, 3],
                "chrom": ["chr1", "chr1", "chr1", "chr1"],
                "start": [0, 1000, 2000, 3000],
                "end": [1000, 2000, 3000, 4000],
            }
        )
        adata.uns["hic_bins"] = bins
        adata.uns["hic_contact_matrix"] = sparse.csr_matrix(
            np.array(
                [
                    [10, 2, 0, 0],
                    [2, 9, 1, 0],
                    [0, 1, 8, 3],
                    [0, 0, 3, 7],
                ],
                dtype=np.float32,
            )
        )
        adata.uns["hic_tad_boundaries"] = pd.DataFrame(
            {
                "bin_id": [0, 1, 2, 3],
                "chrom": ["chr1"] * 4,
                "start": [0, 1000, 2000, 3000],
                "end": [1000, 2000, 3000, 4000],
                "insulation": [0.2, -1.4, -0.8, 0.1],
                "is_boundary": [False, True, False, False],
            }
        )
        adata.uns["hic_compartments"] = pd.DataFrame(
            {
                "bin_id": [0, 1, 2, 3],
                "chrom": ["chr1"] * 4,
                "start": [0, 1000, 2000, 3000],
                "end": [1000, 2000, 3000, 4000],
                "eigenvector_1": [0.3, 0.2, -0.4, -0.2],
                "compartment": ["A", "A", "B", "B"],
            }
        )
        adata.uns["hic_ingest_metadata"] = {
            "resolution_bp": 1000,
            "matrix_format": "csr_sparse",
            "input_format": "tsv_contact_pairs",
            "normalization": "raw_counts_or_unbalanced",
        }
        adata.uns["hic_tad_metadata"] = {
            "compartment_status": "confident",
            "low_information_chromosomes": [],
            "compartment_status_by_chrom": {"chr1": "confident"},
        }
    if include_ribo:
        te = pd.DataFrame(
            {
                "gene": ["CD3E", "LYZ", "ELF3", "KRT8"],
                "sample": ["S0", "S0", "S1", "S1"],
                "footprint_count": [10.0, 20.0, 15.0, 5.0],
                "rna_count": [100.0, 50.0, 30.0, 10.0],
                "te": [10.0 / 101.0, 20.0 / 51.0, 15.0 / 31.0, 5.0 / 11.0],
            }
        )
        adata.uns["ribo_translation_efficiency"] = te
        adata.uns["ribo_ingest_metadata"] = {
            "n_genes_with_footprints": 4,
            "n_samples": 2,
            "te_mean": float(te["te"].mean()),
            "te_median": float(te["te"].median()),
            "min_footprint_threshold": 1,
        }
    adata.write_h5ad(path)


def test_export_bundle_writes_sparse_safe_compact_contract(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path)

    # Pin to v1 explicitly: this test exercises v1-specific schema invariants
    # (CSV files, bundle_manifest.tsv, schema_version == "singlecell_r_bundle_v1").
    # The default writer now emits v2.1; pinning preserves v1 code-path coverage.
    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            source_run_dir=tmp_path,
            obs_columns=("cell_type", "leiden", "missing_col"),
            markers=("ELF3", "CD3E", "MISSING"),
            obsm_keys=("X_umap", "X_pca"),
            marker_chunk_size=2,
            schema_version="v1",
        )
    )

    assert manifest["schema_version"] == "singlecell_r_bundle_v1"
    assert manifest["source"]["n_obs"] == 4
    assert manifest["source"]["n_vars"] == 4
    assert manifest["bundle"]["markers_present"] == ["ELF3", "CD3E"]
    assert manifest["bundle"]["required_files"] == ["X_pca", "X_umap", "marker_expr", "obs"]
    assert manifest["expression"]["source_slot"] == "X"
    assert manifest["expression"]["value_scale"] == "source_X_as_stored"
    assert manifest["expression"]["export_dtype"] == "float32"
    assert manifest["expression"]["intended_use"] == "plotting_and_visual_summary_only"
    assert (
        manifest["expression"]["claim_guard"]
        == "not_for_de_or_new_quantitative_claims_without_full_object_validation"
    )
    assert "full expression matrix is not converted to dense" in manifest["precision_policy"]

    obs = read_bundle_csv(out_dir / "obs.csv.gz")
    marker = read_bundle_csv(out_dir / "marker_expr.csv.gz")
    umap = read_bundle_csv(out_dir / "X_umap.csv.gz")
    pca = read_bundle_csv(out_dir / "X_pca.csv.gz")

    assert obs.columns.tolist() == ["cell_type", "leiden"]
    assert marker.loc["c2", "CD3E"] == 3
    assert marker.loc["c3", "ELF3"] == 0
    assert umap.shape == (4, 2)
    assert pca.shape == (4, 3)

    manifest_disk = json.loads((out_dir / "bundle_manifest.json").read_text())
    assert manifest_disk["files"]["marker_expr"]["sha256"]
    tsv = pd.read_csv(out_dir / "bundle_manifest.tsv", sep="\t")
    assert "precision_policy" in set(tsv["key"])
    tsv_map = dict(zip(tsv["key"], tsv["value"]))
    assert tsv_map["expression_value_scale"] == "source_X_as_stored"
    assert tsv_map["file_marker_expr_n_rows"] == "4"
    assert tsv_map["file_marker_expr_n_cols"] == "2"


def test_export_bundle_can_deterministically_subset_cells(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path)

    export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            markers=("CD3E",),
            max_cells=2,
            seed=7,
        )
    )

    # Default schema is v2.1: outputs are parquet, not CSV.
    # v2.1 marker_expr.parquet has a "cell" column + one column per marker gene.
    obs = pd.read_parquet(out_dir / "obs.parquet")
    marker = pd.read_parquet(out_dir / "marker_expr.parquet")
    assert obs.shape[0] == 2
    assert marker.shape[0] == 2
    assert "CD3E" in marker.columns


def test_export_bundle_fails_when_required_embedding_missing(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_missing_pca.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path, include_pca=False)

    with pytest.raises(ValueError, match="X_pca"):
        export_bundle(ExportConfig(input_h5ad=input_path, output_dir=out_dir))


def test_export_bundle_rejects_invalid_marker_chunk_size(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle"
    make_tiny_h5ad(input_path)

    with pytest.raises(ValueError, match="marker-chunk-size"):
        export_bundle(ExportConfig(input_h5ad=input_path, output_dir=out_dir, marker_chunk_size=0))


def test_export_bundle_v2_mtx_records_sidecar_integrity(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny.h5ad"
    out_dir = tmp_path / "bundle_v2_mtx"
    make_tiny_h5ad(input_path)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            schema_version="v2",
            format="mtx",
            markers=("CD3E", "LYZ", "ELF3"),
        )
    )

    assert manifest["bundle"]["marker_format"] == "mtx_gz"
    assert set(manifest["files"]).issuperset(
        {"marker_expr", "marker_expr_barcodes", "marker_expr_genes"}
    )
    assert set(manifest["bundle"]["required_files"]).issuperset(
        {"marker_expr", "marker_expr_barcodes", "marker_expr_genes"}
    )
    assert manifest["files"]["marker_expr_barcodes"]["sha256"]
    assert manifest["files"]["marker_expr_genes"]["n_rows"] == 3
    assert (out_dir / "marker_expr.barcodes.tsv.gz").exists()
    assert (out_dir / "marker_expr.genes.tsv.gz").exists()


def test_export_bundle_v22_ribo_extension(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_ribo.h5ad"
    out_dir = tmp_path / "bundle_v22_ribo"
    make_tiny_h5ad(input_path, include_ribo=True)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            schema_version="v2.2",
            format="parquet",
            markers=("CD3E", "LYZ"),
            include_ribo=True,
        )
    )

    ribo_ext = manifest["extensions"]["ribo"]
    assert ribo_ext["status"] == "active"
    assert ribo_ext["files"] == ["ribo_te"]
    assert ribo_ext["n_genes"] == 4
    assert ribo_ext["n_samples"] == 2
    assert ribo_ext["table"]["te_path"] == "extensions/ribo/translation_efficiency.parquet"
    assert "ribo_te" in manifest["files"]
    assert "ribo_te" in manifest["bundle"]["required_files"]

    te = pd.read_parquet(out_dir / "extensions" / "ribo" / "translation_efficiency.parquet")
    assert list(te.columns) == ["gene", "sample", "footprint_count", "rna_count", "te"]
    assert te.shape[0] == 4
    assert (te["te"] >= 0).all()


@pytest.mark.parametrize(
    ("bad_case", "message"),
    [
        ("blank_id", "non-empty gene and sample"),
        ("duplicate_key", "duplicate gene/sample"),
        ("non_numeric", "finite numeric"),
        ("non_finite", "finite numeric"),
        ("negative_count", "non-negative"),
        ("fractional_count", "integer-valued"),
        ("te_mismatch", "does not match footprint_count /"),
        ("empty", "at least one row"),
    ],
)
def test_ribo_extension_rejects_malformed_payload_before_writing(
    tmp_path: Path,
    bad_case: str,
    message: str,
) -> None:
    te = pd.DataFrame(
        {
            "gene": ["CD3E", "LYZ"],
            "sample": ["S0", "S0"],
            "footprint_count": [10.0, 20.0],
            "rna_count": [100.0, 50.0],
            "te": [10.0 / 101.0, 20.0 / 51.0],
        }
    )
    if bad_case == "blank_id":
        te.loc[0, "gene"] = "  "
    elif bad_case == "duplicate_key":
        te.loc[1, ["gene", "sample"]] = te.loc[0, ["gene", "sample"]]
    elif bad_case == "non_numeric":
        te["footprint_count"] = te["footprint_count"].astype(object)
        te.loc[0, "footprint_count"] = "not-a-number"
    elif bad_case == "non_finite":
        te.loc[0, "te"] = np.inf
    elif bad_case == "negative_count":
        te.loc[0, "rna_count"] = -1
    elif bad_case == "fractional_count":
        te.loc[0, "footprint_count"] = 1.5
    elif bad_case == "te_mismatch":
        te.loc[0, "te"] = 999.0
    elif bad_case == "empty":
        te = te.iloc[0:0]

    output_dir = tmp_path / "bundle"
    adata = SimpleNamespace(uns={"ribo_translation_efficiency": te})
    with pytest.raises(ValueError, match=message):
        maybe_export_ribo({}, output_dir, adata)
    assert not (output_dir / "extensions" / "ribo").exists()


def test_export_bundle_ribo_requires_v22(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_ribo_v21.h5ad"
    out_dir = tmp_path / "bundle_ribo_v21"
    make_tiny_h5ad(input_path, include_ribo=True)
    with pytest.raises(ValueError, match="include-ribo requires"):
        export_bundle(
            ExportConfig(
                input_h5ad=input_path,
                output_dir=out_dir,
                schema_version="v2.1",
                format="parquet",
                include_ribo=True,
            )
        )


def test_export_bundle_v22_hic_extension(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_hic.h5ad"
    out_dir = tmp_path / "bundle_v22_hic"
    make_tiny_h5ad(input_path, include_hic=True)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            schema_version="v2.2",
            format="parquet",
            markers=("CD3E", "LYZ"),
            include_hic=True,
            hic_max_contacts=100,
        )
    )

    hic_ext = manifest["extensions"]["hic"]
    assert hic_ext["status"] == "active"
    assert hic_ext["files"] == ["hic_bins", "hic_contacts", "hic_boundaries", "hic_compartments"]
    assert hic_ext["n_bins"] == 4
    assert hic_ext["n_boundary_bins"] == 1
    assert hic_ext["resolution_bp"] == 1000
    assert hic_ext["table"]["contacts_path"] == "extensions/hic/contacts.parquet"
    assert hic_ext["hic_tad_metadata"]["compartment_status"] == "confident"
    assert hic_ext["hic_tad_metadata"]["low_information_chromosomes"] == []
    assert hic_ext["hic_tad_metadata"]["compartment_status_by_chrom"] == {"chr1": "confident"}
    for stem in ("hic_bins", "hic_contacts", "hic_boundaries", "hic_compartments"):
        assert stem in manifest["files"]
        assert stem in manifest["bundle"]["required_files"]

    contacts = pd.read_parquet(out_dir / "extensions" / "hic" / "contacts.parquet")
    assert list(contacts.columns) == ["row", "col", "count"]
    assert "cell" not in contacts.columns
    assert contacts.shape[0] == hic_ext["n_contacts"]
    boundaries = pd.read_parquet(out_dir / "extensions" / "hic" / "boundaries.parquet")
    assert list(boundaries.columns) == ["bin_id", "chrom", "position", "insulation", "is_boundary"]
    assert boundaries.loc[boundaries["bin_id"] == 1, "position"].item() == 1500
    compartments = pd.read_parquet(out_dir / "extensions" / "hic" / "compartments.parquet")
    assert set(compartments["compartment"]) == {"A", "B"}


def test_export_bundle_v22_hic_extension_marks_missing_tad_metadata_unknown(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_hic_without_metadata.h5ad"
    out_dir = tmp_path / "bundle_v22_hic_unknown"
    make_tiny_h5ad(input_path, include_hic=True)
    adata = ad.read_h5ad(input_path)
    del adata.uns["hic_tad_metadata"]
    adata.write_h5ad(input_path)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            schema_version="v2.2",
            format="parquet",
            markers=("CD3E", "LYZ"),
            include_hic=True,
            hic_max_contacts=100,
        )
    )

    assert manifest["extensions"]["hic"]["hic_tad_metadata"] == {
        "compartment_status": "unknown_or_unvalidated",
        "low_information_chromosomes": [],
        "compartment_status_by_chrom": {},
    }


def test_export_bundle_v22_hic_extension_allows_low_information_compartment_label(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_hic_low_information.h5ad"
    out_dir = tmp_path / "bundle_v22_hic_low_information"
    make_tiny_h5ad(input_path, include_hic=True)
    adata = ad.read_h5ad(input_path)
    compartments = adata.uns["hic_compartments"].copy()
    compartments.loc[compartments["bin_id"].isin([2, 3]), "eigenvector_1"] = 0.0
    compartments.loc[compartments["bin_id"].isin([2, 3]), "compartment"] = "low_information"
    adata.uns["hic_compartments"] = compartments
    adata.uns["hic_tad_metadata"] = {
        "compartment_status": "partial_low_information",
        "low_information_chromosomes": ["chr1"],
        "compartment_status_by_chrom": {"chr1": "low_information"},
    }
    adata.write_h5ad(input_path)

    manifest = export_bundle(
        ExportConfig(
            input_h5ad=input_path,
            output_dir=out_dir,
            schema_version="v2.2",
            format="parquet",
            markers=("CD3E", "LYZ"),
            include_hic=True,
            hic_max_contacts=100,
        )
    )

    hic_ext = manifest["extensions"]["hic"]
    exported = pd.read_parquet(out_dir / "extensions" / "hic" / "compartments.parquet")
    assert "low_information" in set(exported["compartment"])
    assert hic_ext["hic_tad_metadata"]["compartment_status"] == "partial_low_information"
    assert hic_ext["hic_tad_metadata"]["low_information_chromosomes"] == ["chr1"]


def test_export_bundle_hic_requires_v22_and_contact_cap(tmp_path: Path) -> None:
    input_path = tmp_path / "tiny_hic.h5ad"
    make_tiny_h5ad(input_path, include_hic=True)

    with pytest.raises(ValueError, match="include-hic requires --schema-version v2.2"):
        export_bundle(
            ExportConfig(
                input_h5ad=input_path,
                output_dir=tmp_path / "bad_schema",
                schema_version="v2.1",
                markers=("CD3E",),
                include_hic=True,
            )
        )

    with pytest.raises(ValueError, match="above --hic-max-contacts"):
        export_bundle(
            ExportConfig(
                input_h5ad=input_path,
                output_dir=tmp_path / "too_many_contacts",
                schema_version="v2.2",
                markers=("CD3E",),
                include_hic=True,
                hic_max_contacts=1,
            )
        )


def test_export_cli_output_not_required_with_project_root() -> None:
    from scripts.export_singlecell_r_bundle import build_parser

    args = build_parser().parse_args([
        "--input", "/tmp/final_adata.h5ad",
        "--project-root", "/tmp/project",
        "--run-id", "2026-05-18T0900Z-13c2c88",
    ])
    assert args.output is None
    assert str(args.project_root) == "/tmp/project"


# ---------------------------------------------------------------------------
# Cross-factory manifest: bundle_sha256 back-fill
# ---------------------------------------------------------------------------
# The pipeline writes <run>/manifest.json BEFORE any bundle exists, and emits
# ctx.metadata.get("bundle_sha256", "") — a key nothing ever sets. The field was
# therefore structurally empty in every run, including runs whose overall_status
# was "complete": the cross-factory manifest advertised a bundle checksum it
# never carried, so a consumer could not detect a truncated or substituted
# bundle. The authoritative digest is computed at export time into
# bundle/provenance.json and is now copied into the manifest.

from scripts.export_singlecell_r_bundle import _backfill_project_manifest_bundle_sha256

DIGEST = "5b2277592860747a6092ca79b2a835103d975ef265eeed9cb11f6fa86e474d52"


def _make_run(tmp_path: Path, *, manifest: dict | None, provenance: dict | None):
    run_dir = tmp_path / "runs" / "2026-01-01T0000Z-abcdef0"
    bundle = run_dir / "python" / "bundle"
    bundle.mkdir(parents=True)
    if manifest is not None:
        (run_dir / "manifest.json").write_text(json.dumps(manifest))
    if provenance is not None:
        (bundle / "provenance.json").write_text(json.dumps(provenance))
    return run_dir, bundle


def test_backfill_writes_the_real_digest_into_the_project_manifest(tmp_path):
    run_dir, bundle = _make_run(
        tmp_path,
        manifest={"project_id": "p", "run_id": "r", "bundle_sha256": ""},
        provenance={"bundle_sha256": DIGEST},
    )

    returned = _backfill_project_manifest_bundle_sha256(run_dir, bundle)

    assert returned == DIGEST
    written = json.loads((run_dir / "manifest.json").read_text())
    assert written["bundle_sha256"] == DIGEST
    assert written["bundle_sha256_source"].endswith("provenance.json")
    # Pre-existing manifest content must survive the update.
    assert written["project_id"] == "p" and written["run_id"] == "r"


def test_backfill_is_a_noop_without_a_project_manifest(tmp_path):
    run_dir, bundle = _make_run(tmp_path, manifest=None, provenance={"bundle_sha256": DIGEST})
    assert _backfill_project_manifest_bundle_sha256(run_dir, bundle) == ""


def test_backfill_is_a_noop_without_provenance_or_digest(tmp_path):
    run_dir, bundle = _make_run(tmp_path, manifest={"bundle_sha256": ""}, provenance=None)
    assert _backfill_project_manifest_bundle_sha256(run_dir, bundle) == ""

    run_dir2, bundle2 = _make_run(
        tmp_path / "second", manifest={"bundle_sha256": ""}, provenance={"bundle_sha256": ""}
    )
    assert _backfill_project_manifest_bundle_sha256(run_dir2, bundle2) == ""
    # An empty digest must never overwrite the field with a misleading value.
    assert json.loads((run_dir2 / "manifest.json").read_text())["bundle_sha256"] == ""
