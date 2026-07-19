from __future__ import annotations

import json
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from scripts.produce_public_rna_velocity_figure3a_artifacts import (
    _assert_governed_project_output,
    produce,
)
from workflow.modular.velocity_render_contract import (
    CELL_ARTIFACT_NAME,
    CELL_COLUMNS,
    MANIFEST_NAME,
    MARKER_ARTIFACT_NAME,
    MARKER_COLUMNS,
    SCHEMA_NAME,
    SCHEMA_VERSION,
    choose_cells,
    validate_render_bundle,
    validate_render_tables,
)


def _write_inputs(tmp_path):
    rng = np.random.default_rng(7)
    n_cells, n_genes = 36, 12
    cell_ids = [f"cell-{i:03d}" for i in range(n_cells)]
    gene_ids = ["SOX2", "PAX6", "EOMES"] + [f"gene-{i}" for i in range(n_genes - 3)]
    counts = sparse.csr_matrix(rng.poisson(3, size=(n_cells, n_genes)).astype(np.float32))

    velocity = AnnData(counts.copy())
    velocity.obs_names = cell_ids
    velocity.var_names = gene_ids
    velocity.layers["spliced"] = counts.copy()
    velocity.layers["unspliced"] = sparse.csr_matrix(
        rng.poisson(1, size=(n_cells, n_genes)).astype(np.float32)
    )

    pipeline = AnnData(counts.copy())
    pipeline.obs_names = cell_ids
    pipeline.var_names = gene_ids
    pipeline.obs["cell_type"] = pd.Categorical(
        ["radial_glia"] * 12 + ["intermediate_progenitor"] * 12 + ["neuron"] * 12
    )
    pipeline.obsm["X_umap"] = rng.normal(size=(n_cells, 2))

    velocity_path = tmp_path / "velocity.h5ad"
    pipeline_path = tmp_path / "pipeline.h5ad"
    velocity.write_h5ad(velocity_path)
    pipeline.write_h5ad(pipeline_path)
    return velocity_path, pipeline_path


def _fake_scvelo():
    fake = types.SimpleNamespace(__version__="0.test")
    fake.pp = types.SimpleNamespace()
    fake.tl = types.SimpleNamespace()
    fake.pp.filter_and_normalize = lambda adata, **kwargs: None
    fake.pp.moments = lambda adata, **kwargs: None

    def velocity(adata, mode):
        adata.layers["velocity"] = sparse.csr_matrix(
            np.full((adata.n_obs, adata.n_vars), 0.25, dtype=np.float32)
        )

    fake.tl.velocity = velocity
    fake.tl.velocity_graph = lambda adata: adata.uns.__setitem__("velocity_graph", "computed")
    fake.tl.velocity_embedding = lambda adata, basis: adata.obsm.__setitem__(
        "velocity_umap", np.column_stack([np.full(adata.n_obs, 0.1), np.full(adata.n_obs, -0.2)])
    )
    fake.tl.velocity_confidence = lambda adata: adata.obs.__setitem__(
        "velocity_confidence", np.linspace(0.2, 0.9, adata.n_obs)
    )
    fake.tl.velocity_pseudotime = lambda adata: adata.obs.__setitem__(
        "velocity_pseudotime", np.linspace(0.0, 1.0, adata.n_obs)
    )
    return fake


def test_choose_cells_is_seed_reproducible_and_sorted():
    obs = pd.DataFrame(
        {"cell_type": ["common"] * 90 + ["rare"] * 10},
        index=[f"cell-{i:03d}" for i in range(100)],
    )
    first = choose_cells(obs, max_cells=40, seed=13)
    second = choose_cells(obs.sample(frac=1, random_state=4), max_cells=40, seed=13)
    assert first == sorted(first)
    assert len(first) == len(set(first)) == 40
    assert sum(cell in set(obs.index[-10:]) for cell in first) == 10
    assert second == first


def test_producer_writes_strict_named_hash_linked_bundle(tmp_path):
    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    manifest = produce(
        velocity_path,
        pipeline_path,
        output_dir,
        max_cells=30,
        top_genes=8,
        seed=13,
        min_shared_counts=1,
        n_pcs=5,
        n_neighbors=5,
        n_pseudotime_bins=6,
        markers=("SOX2", "PAX6", "EOMES", "MISSING"),
        scv=_fake_scvelo(),
    )

    cells = pd.read_csv(output_dir / CELL_ARTIFACT_NAME)
    trends = pd.read_csv(output_dir / MARKER_ARTIFACT_NAME)
    assert tuple(cells.columns) == CELL_COLUMNS
    assert tuple(trends.columns) == MARKER_COLUMNS
    assert cells["cell_id"].is_unique
    assert len(cells) == 30
    assert set(trends["marker"]) == {"SOX2", "PAX6", "EOMES"}
    assert manifest["schema_version"] == SCHEMA_VERSION
    assert manifest["schema_name"] == SCHEMA_NAME
    assert manifest["status"] == "complete"
    assert manifest["claim_class"] == "exploratory"
    assert manifest["seed"] == 13
    assert manifest["markers"]["missing"] == ["MISSING"]
    assert validate_render_bundle(output_dir) == manifest


def test_provenance_and_artifacts_reproduce_byte_for_byte(tmp_path):
    velocity_path, pipeline_path = _write_inputs(tmp_path)
    kwargs = dict(
        max_cells=24,
        top_genes=8,
        seed=9,
        min_shared_counts=1,
        n_pcs=5,
        n_neighbors=5,
        n_pseudotime_bins=5,
        markers=("SOX2", "PAX6"),
    )
    first = tmp_path / "first"
    second = tmp_path / "second"
    produce(velocity_path, pipeline_path, first, scv=_fake_scvelo(), **kwargs)
    produce(velocity_path, pipeline_path, second, scv=_fake_scvelo(), **kwargs)

    for filename in (CELL_ARTIFACT_NAME, MARKER_ARTIFACT_NAME, MANIFEST_NAME):
        assert (first / filename).read_bytes() == (second / filename).read_bytes()


def test_validator_rejects_value_tampering_by_hash(tmp_path):
    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    produce(
        velocity_path,
        pipeline_path,
        output_dir,
        max_cells=20,
        top_genes=8,
        min_shared_counts=1,
        markers=("SOX2",),
        scv=_fake_scvelo(),
    )
    cells = pd.read_csv(output_dir / CELL_ARTIFACT_NAME)
    cells.loc[0, "umap1"] += 1.0
    cells.to_csv(output_dir / CELL_ARTIFACT_NAME, index=False)
    with pytest.raises(ValueError, match="sha256 does not match"):
        validate_render_bundle(output_dir)

    manifest = json.loads((output_dir / MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["artifacts"]["cells"]["sha256"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda manifest: manifest.pop("truth_boundary"), "truth_boundary"),
        (
            lambda manifest: manifest["sources"]["velocity_h5ad"].update({"sha256": "not-a-sha"}),
            "sources.velocity_h5ad.sha256",
        ),
        (
            lambda manifest: manifest.update({"reproducibility_key_sha256": "0" * 64}),
            "reproducibility_key_sha256 does not match",
        ),
        (
            lambda manifest: manifest.update({"schema_name": "wrong_schema"}),
            "schema name",
        ),
        (
            lambda manifest: manifest["parameters"].pop("velocity_mode"),
            "missing parameters",
        ),
        (
            lambda manifest: manifest["software_versions"].pop("scvelo"),
            "missing software versions",
        ),
        (
            lambda manifest: manifest.update(
                {"truth_boundary": "different but still nonblank scientific boundary"}
            ),
            "reproducibility_key_sha256 does not match",
        ),
        (
            lambda manifest: manifest["markers"]["missing"].append("NEW_MISSING"),
            "reproducibility_key_sha256 does not match",
        ),
    ],
)
def test_validator_rejects_incomplete_or_tampered_manifest(tmp_path, mutation, message):
    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    produce(
        velocity_path,
        pipeline_path,
        output_dir,
        max_cells=20,
        top_genes=8,
        min_shared_counts=1,
        markers=("SOX2",),
        scv=_fake_scvelo(),
    )
    manifest_path = output_dir / MANIFEST_NAME
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    mutation(manifest)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match=message):
        validate_render_bundle(output_dir)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("pseudotime_bin_index", 0.5, "must contain integers"),
        ("pseudotime_bin_midpoint", 1.1, r"within \[0, 1\]"),
    ],
)
def test_table_validator_rejects_invalid_marker_bins(tmp_path, column, value, message):
    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    produce(
        velocity_path,
        pipeline_path,
        output_dir,
        max_cells=20,
        top_genes=8,
        min_shared_counts=1,
        markers=("SOX2",),
        scv=_fake_scvelo(),
    )
    cells = pd.read_csv(output_dir / CELL_ARTIFACT_NAME)
    trends = pd.read_csv(output_dir / MARKER_ARTIFACT_NAME)
    trends[column] = trends[column].astype(float)
    trends.loc[0, column] = value
    with pytest.raises(ValueError, match=message):
        validate_render_tables(cells, trends)


def test_producer_refuses_to_overwrite_immutable_bundle(tmp_path):
    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    kwargs = dict(
        max_cells=20,
        top_genes=8,
        min_shared_counts=1,
        markers=("SOX2",),
        scv=_fake_scvelo(),
    )
    produce(velocity_path, pipeline_path, output_dir, **kwargs)
    with pytest.raises(FileExistsError, match="immutable render artifact"):
        produce(velocity_path, pipeline_path, output_dir, **kwargs)


def test_atomic_publication_cleans_partial_staging_and_allows_retry(tmp_path, monkeypatch):
    import workflow.modular.velocity_render_contract as contract

    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    kwargs = dict(
        max_cells=20,
        top_genes=8,
        min_shared_counts=1,
        markers=("SOX2",),
        scv=_fake_scvelo(),
    )
    original_validate = contract.validate_render_bundle
    monkeypatch.setattr(
        contract,
        "validate_render_bundle",
        lambda bundle_dir: (_ for _ in ()).throw(ValueError("injected validation failure")),
    )
    with pytest.raises(ValueError, match="injected validation failure"):
        produce(velocity_path, pipeline_path, output_dir, **kwargs)
    assert not output_dir.exists()
    assert list(tmp_path.glob(".bundle.*.tmp")) == []

    monkeypatch.setattr(contract, "validate_render_bundle", original_validate)
    produce(velocity_path, pipeline_path, output_dir, **kwargs)
    assert (output_dir / CELL_ARTIFACT_NAME).is_file()
    assert (output_dir / MARKER_ARTIFACT_NAME).is_file()
    assert (output_dir / MANIFEST_NAME).is_file()


def test_atomic_publication_cleans_staging_on_source_record_failure(tmp_path, monkeypatch):
    import workflow.modular.velocity_render_contract as contract

    velocity_path, pipeline_path = _write_inputs(tmp_path)
    output_dir = tmp_path / "bundle"
    monkeypatch.setattr(
        contract,
        "_file_record",
        lambda path: (_ for _ in ()).throw(OSError("injected source hash failure")),
    )
    with pytest.raises(OSError, match="injected source hash failure"):
        produce(
            velocity_path,
            pipeline_path,
            output_dir,
            max_cells=20,
            top_genes=8,
            min_shared_counts=1,
            markers=("SOX2",),
            scv=_fake_scvelo(),
        )
    assert not output_dir.exists()
    assert list(tmp_path.glob(".bundle.*.tmp")) == []


def test_cli_guard_requires_governed_project_run_path():
    factory_output = Path(__file__).resolve().parents[1] / "forbidden-scientific-output"
    with pytest.raises(ValueError, match="governed project run"):
        _assert_governed_project_output(factory_output)
    with pytest.raises(ValueError, match="governed project run"):
        _assert_governed_project_output(Path("/tmp/figure3a-output"))
    with pytest.raises(ValueError, match="governed project run"):
        _assert_governed_project_output(Path("/home/zerlinshen/projects/project-only"))

    governed = Path(
        "/home/zerlinshen/projects/wave5-trevino/runs/run-id/python/figure3a_velocity_render"
    )
    assert _assert_governed_project_output(governed) is None


def test_direct_script_help_bootstraps_repository_imports():
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "scripts" / "produce_public_rna_velocity_figure3a_artifacts.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "--velocity-h5ad" in result.stdout
    assert "--pipeline-h5ad" in result.stdout
    assert "--out-dir" in result.stdout
