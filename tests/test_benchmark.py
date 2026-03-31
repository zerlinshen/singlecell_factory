from pathlib import Path
import json

import numpy as np
from anndata import AnnData

import workflow.benchmark as bm


def test_cluster_ari(tmp_path):
    a = AnnData(np.ones((3, 2)))
    a.obs_names = ["c1", "c2", "c3"]
    a.obs["leiden"] = ["0", "1", "1"]
    b = AnnData(np.ones((3, 2)))
    b.obs_names = ["c1", "c2", "c3"]
    b.obs["leiden"] = ["0", "1", "1"]
    p1 = tmp_path / "a.h5ad"
    p2 = tmp_path / "b.h5ad"
    a.write(p1)
    b.write(p2)
    ari = bm.cluster_ari(p1, p2)
    assert ari == 1.0


def test_velocity_direction_consistency(tmp_path):
    from scipy.sparse import csr_matrix

    a1 = AnnData(np.ones((3, 2)))
    a1.uns["velocity_graph"] = csr_matrix(
        np.array(
            [
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 0.0],
            ]
        )
    )
    a2 = a1.copy()
    p1 = tmp_path / "v1.h5ad"
    p2 = tmp_path / "v2.h5ad"
    a1.write(p1)
    a2.write(p2)
    score = bm.velocity_direction_consistency(p1, p2)
    assert 0.9 <= score <= 1.0


def test_cluster_ari_no_overlap(tmp_path):
    a = AnnData(np.ones((1, 1)))
    a.obs_names = ["a1"]
    a.obs["leiden"] = ["0"]
    b = AnnData(np.ones((1, 1)))
    b.obs_names = ["b1"]
    b.obs["leiden"] = ["1"]
    p1 = tmp_path / "a1.h5ad"
    p2 = tmp_path / "b1.h5ad"
    a.write(p1)
    b.write(p2)
    assert np.isnan(bm.cluster_ari(p1, p2))


def test_velocity_direction_consistency_nan(tmp_path):
    a = AnnData(np.ones((2, 2)))
    b = AnnData(np.ones((2, 2)))
    p1 = tmp_path / "none1.h5ad"
    p2 = tmp_path / "none2.h5ad"
    a.write(p1)
    b.write(p2)
    assert np.isnan(bm.velocity_direction_consistency(p1, p2))


def test_run_benchmark(monkeypatch, tmp_path):
    fake_root = tmp_path
    out_root = fake_root / "output" / "workflow_benchmark"
    out_root.mkdir(parents=True, exist_ok=True)

    def fake_run_standard(cfg):
        p = out_root / cfg.project / "processed.h5ad"
        p.parent.mkdir(parents=True, exist_ok=True)
        a = AnnData(np.ones((2, 2)))
        a.obs_names = ["c1", "c2"]
        a.obs["leiden"] = ["0", "1"]
        a.write(p)
        return p

    def fake_cluster_ari(new_h5ad, baseline_h5ad):
        return 0.99

    def fake_velocity_consistency(path1, path2):
        return 0.95

    monkeypatch.setattr(bm, "run_standard_workflow", fake_run_standard)
    monkeypatch.setattr(bm, "cluster_ari", fake_cluster_ari)
    monkeypatch.setattr(bm, "velocity_direction_consistency", fake_velocity_consistency)

    out = bm.run_benchmark(fake_root)
    data = json.loads(Path(out).read_text(encoding="utf-8"))
    assert len(data["datasets"]) == 3
    assert data["velocity_direction_consistency"] == 0.95
