from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil

import anndata as ad
from anndata.experimental import concat_on_disk
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


def load_prepare_module():
    repo_root = Path(__file__).resolve().parents[1]
    module_path = repo_root / "scripts" / "prepare_emtab13526_full_cohort_zarr.py"
    spec = importlib.util.spec_from_file_location("prepare_emtab13526_full_cohort_zarr", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_build_sample_meta_yields_required_columns(tmp_path):
    module = load_prepare_module()
    sdrf = tmp_path / "E-MTAB-13526.sdrf.txt"
    pd.DataFrame(
        {
            "Source Name": ["P1_T1", "P1_B1"],
            "Characteristics[individual]": ["P1", "P1"],
            "Characteristics[disease]": ["Tumor", "normal"],
            "Characteristics[FACS]": ["CD45+", "CD45+"],
            "Characteristics[sampling site]": ["site1", "site2"],
            "Characteristics[sex]": ["F", "F"],
            "Characteristics[original source name]": ["orig1", "orig2"],
        }
    ).to_csv(sdrf, sep="\t", index=False)

    meta = module.build_sample_meta(sdrf)
    assert set(module.REQUIRED_OBS_COLUMNS).issubset(meta.columns)
    condition_map = dict(zip(meta["sample"], meta["condition"]))
    assert condition_map["P1_T1"] == "tumor"
    assert condition_map["P1_B1"] == "healthy_background"


def make_valid_merged_zarr(tmp_path: Path, include_condition: bool = True) -> Path:
    module = load_prepare_module()
    parts_root = tmp_path / "parts"
    parts_root.mkdir(parents=True, exist_ok=True)

    var = pd.DataFrame(
        {"gene_symbol": ["g1", "g2", "g3"]},
        index=pd.Index(["g1", "g2", "g3"], name=None),
    )
    for idx, sample in enumerate(["S1", "S2"], start=1):
        x = sparse.csr_matrix(np.array([[1, 0, 2], [0, 3, 0]], dtype=np.float32) + idx)
        obs = pd.DataFrame(
            {
                "sample": [sample, sample],
                "patient": [f"P{idx}", f"P{idx}"],
                "batch": [sample, sample],
                "disease": ["Tumor", "Tumor"],
                "sorting": ["CD45+", "CD45+"],
                "sampling_site": ["site1", "site1"],
                "sex": ["F", "F"],
                "original_source_name": ["orig", "orig"],
                "tumor_type": ["NSCLC", "NSCLC"],
            },
            index=[f"{sample}_1", f"{sample}_2"],
        )
        if include_condition:
            obs["condition"] = ["tumor", "tumor"]
        obs = obs[
            [
                "sample",
                "patient",
                "batch",
                "disease",
                *([] if not include_condition else ["condition"]),
                "sorting",
                "sampling_site",
                "sex",
                "original_source_name",
                "tumor_type",
            ]
        ]
        adata = ad.AnnData(X=x, obs=obs, var=var)
        adata.write_zarr(parts_root / f"{sample}.zarr")

    merged = tmp_path / "merged.zarr"
    concat_on_disk(
        [parts_root / "S1.zarr", parts_root / "S2.zarr"],
        merged,
        axis=0,
        join="outer",
    )
    return merged


def write_summary_json(tmp_path: Path, prepared_zarr: Path, *, n_samples: int = 2, retained_barcodes_total: int = 4) -> Path:
    summary_path = tmp_path / "prepared_input.summary.json"
    payload = {
        "prepared_zarr": str(prepared_zarr),
        "n_samples": n_samples,
        "raw_barcodes_total": retained_barcodes_total,
        "retained_barcodes_total": retained_barcodes_total,
        "nnz_total": retained_barcodes_total,
        "samples": [],
    }
    summary_path.write_text(json.dumps(payload), encoding="utf-8")
    return summary_path


def test_validate_prepared_input_passes_on_toy_csr_zarr_parts(tmp_path):
    module = load_prepare_module()
    merged = make_valid_merged_zarr(tmp_path)
    summary = write_summary_json(tmp_path, merged, n_samples=2, retained_barcodes_total=4)
    payload = module.validate_prepared_input(merged, expected_n_samples=2, summary_json=summary)
    assert payload["n_samples"] == 2
    assert payload["retained_barcodes_total"] == 4
    assert payload["x_encoding"] == "csr_matrix"


def test_validate_prepared_input_fails_when_x_storage_is_incomplete(tmp_path):
    module = load_prepare_module()
    merged = make_valid_merged_zarr(tmp_path)
    summary = write_summary_json(tmp_path, merged, n_samples=2, retained_barcodes_total=4)
    shutil.rmtree(merged / "X" / "indptr")
    with pytest.raises(ValueError, match="missing X/indptr"):
        module.validate_prepared_input(merged, expected_n_samples=2, summary_json=summary)


def test_validate_prepared_input_fails_when_required_obs_columns_are_missing(tmp_path):
    module = load_prepare_module()
    merged = make_valid_merged_zarr(tmp_path, include_condition=False)
    summary = write_summary_json(tmp_path, merged, n_samples=2, retained_barcodes_total=4)
    with pytest.raises(ValueError, match="required obs columns"):
        module.validate_prepared_input(merged, expected_n_samples=2, summary_json=summary)


def test_validate_prepared_input_requires_summary_json(tmp_path):
    module = load_prepare_module()
    merged = make_valid_merged_zarr(tmp_path)
    with pytest.raises(FileNotFoundError, match="summary missing"):
        module.validate_prepared_input(
            merged,
            expected_n_samples=2,
            summary_json=tmp_path / "missing.summary.json",
        )


def test_validate_prepared_input_fails_when_summary_obs_count_mismatches(tmp_path):
    module = load_prepare_module()
    merged = make_valid_merged_zarr(tmp_path)
    summary = write_summary_json(tmp_path, merged, n_samples=2, retained_barcodes_total=999)
    with pytest.raises(ValueError, match="obs count mismatch"):
        module.validate_prepared_input(merged, expected_n_samples=2, summary_json=summary)


def load_fallback_script_text() -> str:
    repo_root = Path(__file__).resolve().parents[1]
    return (repo_root / "scripts" / "run_emtab13526_full_cohort_with_fallback.sh").read_text()


def test_full_cohort_fallback_uses_directory_checks_for_prepared_zarr():
    script = load_fallback_script_text()
    assert '[[ ! -d "${PREPARED_ZARR}" ]]' in script
    assert '[[ ! -f "${PREPARED_ZARR}" ]]' not in script


def test_full_cohort_cleanup_preserves_existing_prepared_input():
    script = load_fallback_script_text()
    cleanup_block = script.split("cleanup_derived_artifacts() {", 1)[1].split("\n}\n", 1)[0]
    assert '"${FULL_ROOT}/prepared_input.zarr"' not in cleanup_block
    assert '"${FULL_ROOT}/prepared_input.ready"' not in cleanup_block
