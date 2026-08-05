#!/usr/bin/env python3
"""Joint spatial × communication method smoke (Wave-8).

Loads the 10x Visium mouse-brain sample, builds spatial neighbors, clusters
spots, then runs LIANA consensus treating Leiden clusters as "cell types".

Claim honesty
-------------
- Methodology smoke only: proves the **joint software path**
  (spatial graph + LIANA) is executable on a real Visium object.
- Tissue is **mouse brain**, not NSCLC. Do **not** upgrade F4-02 to
  NSCLC localization or paper multi-condition claims from this smoke.

Acceptance
----------
- spatial_neighbors succeeds
- LIANA returns a non-empty table
- At least one VEGF-family or checkpoint-family panel row after filter
  (if resource has human symbols; mouse may miss — then empty panels are
  recorded but LIANA non-empty still passes engine claim)

Run::
  python scripts/test_joint_spatial_liana_smoke.py [--output-dir DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

VISIUM = Path("/home/zerlinshen/data/raw/visium_v1_mouse_brain")
_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional directory for liana CSV + summary JSON",
    )
    args = p.parse_args(argv)

    if not (VISIUM / "filtered_feature_bc_matrix.h5").exists():
        print(f"FAIL: Visium missing at {VISIUM}", file=sys.stderr)
        return 2

    import scanpy as sc
    import squidpy as sq
    import liana as li

    from scripts.filter_liana_claim_panels import filter_panels

    adata = sc.read_visium(str(VISIUM))
    adata.var_names_make_unique()
    # Human-symbol upper for LIANA resource matching (many genes shared orthologs)
    adata.var_names = pd.Index([str(g).upper() for g in adata.var_names])
    adata.var_names_make_unique()

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sq.gr.spatial_neighbors(adata, coord_type="generic")
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.4)
    adata.obs["cell_type"] = "cluster_" + adata.obs["leiden"].astype(str)

    # Cap spots for runtime
    if adata.n_obs > 2500:
        rng = np.random.default_rng(41)
        idx = np.sort(rng.choice(adata.n_obs, size=2500, replace=False))
        adata = adata[idx].copy()

    li.mt.rank_aggregate(
        adata,
        groupby="cell_type",
        resource_name="consensus",
        use_raw=False,
        verbose=False,
    )
    res = adata.uns.get("liana_res")
    if res is None or len(res) == 0:
        print("FAIL: empty LIANA results", file=sys.stderr)
        return 3

    panels = filter_panels(res)
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "tissue": "mouse_brain_visium_demo",
        "claimable_as_nsclc": False,
        "n_spots": int(adata.n_obs),
        "n_clusters": int(adata.obs["cell_type"].nunique()),
        "spatial_connectivities": "spatial_connectivities" in adata.obsp,
        "n_liana": int(len(res)),
        "n_vegf_egfr": int(len(panels["vegf_egfr"])),
        "n_checkpoint_any": int(len(panels["checkpoint_any"])),
        "engine": "liana_consensus",
        "scientific_boundary": (
            "Joint spatial_neighbors + LIANA path smoke on non-NSCLC Visium; "
            "not tissue localization of NSCLC claims."
        ),
    }

    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        res.to_csv(args.output_dir / "joint_spatial_liana.csv", index=False)
        (args.output_dir / "joint_spatial_liana_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )

    print(json.dumps(summary, indent=2))
    if not summary["spatial_connectivities"]:
        print("FAIL: missing spatial_connectivities", file=sys.stderr)
        return 4
    if summary["n_liana"] < 1:
        print("FAIL: n_liana < 1", file=sys.stderr)
        return 5
    print("JOINT SPATIAL×LIANA VISIUM SMOKE PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
