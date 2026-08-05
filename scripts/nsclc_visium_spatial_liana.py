#!/usr/bin/env python3
"""NSCLC Visium paper-tissue spatial × communication lane (Wave-9).

Loads De Zuani et al. 2024 (Nat Commun 15:4388) 10x Visium sections staged
from E-MTAB-13530, and per section runs the validated joint path:
spatial_neighbors -> Leiden spot clusters -> LIANA consensus -> claim-panel
filters, plus De Zuani / macrophage panel ``score_genes`` and Moran's I
spatial autocorrelation. Finally emits a paired tumour-vs-background
descriptive contrast within patients.

Science / scope
---------------
- First **NSCLC-tissue** spatial lane: tissue is the paper's own Visium
  cohort (human LUAD/LUSC tumour + adjacent background), so panel-level
  spatial summaries are claimable as NSCLC-tissue method evidence.
- NOT paper figure parity: spot groups are Leiden clusters (the paper used
  pathologist annotation + cell2location deconvolution); LIANA consensus is
  our engine, not the paper's CellphoneDB multi-condition design. All
  outputs are descriptive and stamped ``figure_parity=False``.
- Do not upgrade F4-01/F4-02 beyond ``partial`` from this lane alone.

Usage
-----
  python scripts/nsclc_visium_spatial_liana.py \\
    --sections-root /home/zerlinshen/data/raw/nc2024_nsclc_visium_emtab13530/sections \\
    --section-map  /home/zerlinshen/data/raw/nc2024_nsclc_visium_emtab13530/section_map.csv \\
    --output-dir   <run_dir>/python/nsclc_spatial_liana
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.filter_liana_claim_panels import filter_panels  # noqa: E402
from workflow.modular.modules.gene_signature_scoring import (  # noqa: E402
    BUILTIN_SIGNATURES,
)

SCORE_PANELS = [
    "de_zuani_stab1_foetal",
    "de_zuani_stab1_signature_reconstructed",
    "de_zuani_cholesterol_export",
    "foetal_like_mac",
    "spp1_mac",
]

BOUNDARY = (
    "De Zuani 2024 paper-cohort NSCLC Visium (E-MTAB-13530); Leiden spot "
    "clusters as proxy groups (not cell2location deconvolution, no "
    "pathologist annotation); LIANA consensus engine (not paper CellphoneDB "
    "multi-condition design). Descriptive spatial method lane — not paper "
    "figure parity."
)


def _moran_i(values: np.ndarray, W) -> float:
    """Moran's I against a spatial weights matrix (esda 'r' convention).

    Used for per-spot panel scores, which squidpy's ``spatial_autocorr``
    cannot take (var_names only in squidpy 1.8). Weights are row-standardised
    before applying I = (n/S0) * z'Wz / z'z.
    """
    v = np.asarray(values, dtype=float)
    mask = ~np.isnan(v)
    if int(mask.sum()) < 10:
        return float("nan")
    idx = np.where(mask)[0]
    v = v[idx]
    Ws = W[idx][:, idx]
    row_sums = np.asarray(Ws.sum(axis=1)).ravel()
    nz = row_sums > 0
    v = v[nz]
    Ws = Ws[nz][:, nz]
    row_sums = row_sums[nz]
    Ws = Ws.multiply(1.0 / row_sums[:, None]).tocsr()
    z = v - v.mean()
    denom = float(z @ z)
    s0 = float(Ws.sum())
    if denom == 0.0 or s0 == 0.0:
        return float("nan")
    num = float(z @ (Ws @ z))
    return (len(v) / s0) * (num / denom)


def analyze_section(
    section_dir: Path,
    *,
    max_spots: int = 2500,
    seed: int = 41,
    leiden_resolution: float = 0.4,
) -> dict:
    """Ingest one staged Visium section and run the joint path on it."""
    import scanpy as sc

    h5 = section_dir / "filtered_feature_bc_matrix.h5"
    if not h5.exists():
        raise FileNotFoundError(f"missing {h5}")
    if not (section_dir / "spatial").exists():
        raise FileNotFoundError(f"missing {section_dir / 'spatial'}")

    adata = sc.read_visium(str(section_dir))
    adata.var_names_make_unique()
    return analyze_adata(
        adata, max_spots=max_spots, seed=seed, leiden_resolution=leiden_resolution
    )


def analyze_adata(
    adata,
    *,
    max_spots: int = 2500,
    seed: int = 41,
    leiden_resolution: float = 0.4,
):
    """Joint spatial x LIANA path on a Visium AnnData (count space).

    Separated from ingest so tests can drive it with synthetic objects.
    """
    import scanpy as sc
    import squidpy as sq
    import liana as li

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sq.gr.spatial_neighbors(adata, coord_type="generic")
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution)
    adata.obs["cell_type"] = "cluster_" + adata.obs["leiden"].astype(str)

    if adata.n_obs > max_spots:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(adata.n_obs, size=max_spots, replace=False))
        adata = adata[idx].copy()

    # --- panel scores (log-normalised space) ---
    panel_coverage: dict[str, str] = {}
    score_cols: list[str] = []
    for panel in SCORE_PANELS:
        genes = BUILTIN_SIGNATURES[panel]
        present = [g for g in genes if g in adata.var_names]
        panel_coverage[panel] = f"{len(present)}/{len(genes)}"
        col = f"score_{panel}"
        if len(present) >= 3:
            sc.tl.score_genes(adata, present, score_name=col, use_raw=False)
        else:
            adata.obs[col] = np.nan
        score_cols.append(col)

    # --- Moran's I: genes via squidpy; panel scores via weights helper ---
    # squidpy 1.8 spatial_autocorr only accepts var_names, so per-spot panel
    # scores are computed against the same spatial_connectivities graph with
    # row-standardised weights (esda 'r' convention) in _moran_i. The squidpy
    # gene path needs numba JIT; where JIT is unavailable (e.g. pytest runs
    # with NUMBA_DISABLE_JIT=1) genes fall back to the same weights helper and
    # the engine is recorded in moran_engine.
    panel_genes = sorted(
        {g for p in SCORE_PANELS for g in BUILTIN_SIGNATURES[p] if g in adata.var_names}
    )
    moran_targets = [c for c in score_cols if not adata.obs[c].isna().all()]
    W = adata.obsp["spatial_connectivities"]
    moran_engine = "squidpy_spatial_autocorr(genes)+weights_moran(panel_scores)"
    try:
        sq.gr.spatial_autocorr(adata, genes=panel_genes, mode="moran")
        moran = adata.uns["moranI"]
        if not isinstance(moran, pd.DataFrame):
            moran = pd.DataFrame(moran)
        moran = moran.reset_index().rename(
            columns={moran.reset_index().columns[0]: "feature"}
        )
        moran["feature_kind"] = "gene"
    except Exception as exc:
        moran_engine = f"weights_moran(all) [squidpy unavailable: {type(exc).__name__}]"
        gene_rows = []
        for g in panel_genes:
            x = adata[:, g].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            gene_rows.append(
                {
                    "feature": g,
                    "I": _moran_i(np.asarray(x).ravel().astype(float), W),
                    "feature_kind": "gene",
                }
            )
        moran = pd.DataFrame(gene_rows)

    score_rows = []
    for col in moran_targets:
        score_rows.append(
            {
                "feature": col,
                "I": _moran_i(adata.obs[col].to_numpy(dtype=float), W),
                "feature_kind": "panel_score",
            }
        )
    moran = pd.concat(
        [moran, pd.DataFrame(score_rows)], ignore_index=True
    )

    # --- LIANA consensus on Leiden spot clusters ---
    li.mt.rank_aggregate(
        adata,
        groupby="cell_type",
        resource_name="consensus",
        use_raw=False,
        verbose=False,
    )
    res = adata.uns.get("liana_res")
    if res is None or len(res) == 0:
        raise RuntimeError("empty LIANA results")
    panels = filter_panels(res)

    score_means = {
        c: (None if adata.obs[c].isna().all() else float(adata.obs[c].mean()))
        for c in score_cols
    }
    moran_scores = (
        moran[moran["feature"].isin(moran_targets)]
        .set_index("feature")["I"]
        .to_dict()
    )

    return {
        "adata_n_spots": int(adata.n_obs),
        "n_clusters": int(adata.obs["cell_type"].nunique()),
        "liana_res": res,
        "moran_table": moran,
        "summary": {
            "n_spots_analyzed": int(adata.n_obs),
            "n_clusters": int(adata.obs["cell_type"].nunique()),
            "spatial_connectivities": "spatial_connectivities" in adata.obsp,
            "n_liana": int(len(res)),
            "n_vegf_egfr": int(len(panels["vegf_egfr"])),
            "n_classic_vegf_vegfr": int(len(panels["classic_vegf_vegfr"])),
            "n_vegfa_egfr": int(len(panels["vegfa_egfr"])),
            "n_checkpoint_any": int(len(panels["checkpoint_any"])),
            "n_checkpoint_immune": int(len(panels["checkpoint_immune"])),
            "checkpoint_immune_note": (
                "checkpoint_immune requires curated immune cell-type labels; "
                "Leiden cluster_X labels never match, so it is 0 by "
                "construction in this lane — use checkpoint_any."
            ),
            "panel_gene_coverage": panel_coverage,
            "panel_score_means": score_means,
            "panel_score_moranI": {k: float(v) for k, v in moran_scores.items()},
            "moran_engine": moran_engine,
            "engine": "liana_consensus",
            "spot_groups": "leiden_clusters",
        },
    }


def build_contrast(
    section_rows: list[dict],
    section_map: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    """Paired tumour-vs-background descriptive contrast within patients."""
    rows = []
    for r in section_rows:
        sec = r["section"]
        meta = section_map.set_index("section").loc[sec]
        s = r["summary"]
        row = {
            "section": sec,
            "patient": meta["patient"],
            "histology": meta["histology"],
            "condition": meta["condition"],
            "n_spots": s["n_spots_analyzed"],
            "n_liana": s["n_liana"],
            "n_vegf_egfr": s["n_vegf_egfr"],
            "n_classic_vegf_vegfr": s["n_classic_vegf_vegfr"],
            "n_checkpoint_any": s["n_checkpoint_any"],
            "n_checkpoint_immune": s["n_checkpoint_immune"],
        }
        for panel in SCORE_PANELS:
            row[f"mean_{panel}"] = s["panel_score_means"].get(f"score_{panel}")
            row[f"moranI_{panel}"] = s["panel_score_moranI"].get(f"score_{panel}")
        rows.append(row)
    table = pd.DataFrame(rows)

    contrasts: list[dict] = []
    for (patient, hist), grp in table.groupby(["patient", "histology"]):
        t = grp[grp["condition"] == "tumor"]
        b = grp[grp["condition"] == "background"]
        if t.empty or b.empty:
            continue
        entry: dict = {
            "patient": patient,
            "histology": hist,
            "n_tumor_sections": int(len(t)),
            "n_background_sections": int(len(b)),
            "delta_n_checkpoint_any": float(
                t["n_checkpoint_any"].mean() - b["n_checkpoint_any"].mean()
            ),
            # checkpoint_immune is 0 by construction for Leiden cluster labels
            # (see per-section checkpoint_immune_note); retained for schema
            # parity with scRNA lanes, not for interpretation.
            "delta_n_checkpoint_immune": float(
                t["n_checkpoint_immune"].mean() - b["n_checkpoint_immune"].mean()
            ),
            "delta_n_classic_vegf_vegfr": float(
                t["n_classic_vegf_vegfr"].mean() - b["n_classic_vegf_vegfr"].mean()
            ),
            "panel_score_deltas": {},
            "panel_moranI_deltas": {},
        }
        for panel in SCORE_PANELS:
            mt, mb = t[f"mean_{panel}"].mean(), b[f"mean_{panel}"].mean()
            it, ib = t[f"moranI_{panel}"].mean(), b[f"moranI_{panel}"].mean()
            entry["panel_score_deltas"][panel] = (
                None if np.isnan(mt) or np.isnan(mb) else float(mt - mb)
            )
            entry["panel_moranI_deltas"][panel] = (
                None if np.isnan(it) or np.isnan(ib) else float(it - ib)
            )
        contrasts.append(entry)

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "contrast_unit": "within-patient tumour-minus-background section means",
        "patients": contrasts,
        "claimable_as_nsclc_tissue": True,
        "figure_parity": False,
        "scientific_boundary": BOUNDARY,
    }
    return table, summary


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sections-root", required=True, type=Path)
    p.add_argument("--section-map", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--max-spots", type=int, default=2500)
    p.add_argument("--seed", type=int, default=41)
    args = p.parse_args(argv)

    section_map = pd.read_csv(args.section_map)
    sections = section_map["section"].tolist()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    section_rows: list[dict] = []
    errors: list[dict] = []
    for sec in sections:
        sec_dir = args.sections_root / sec
        out = args.output_dir / sec
        out.mkdir(parents=True, exist_ok=True)
        print(f"=== {sec} ===", flush=True)
        try:
            result = analyze_section(
                sec_dir, max_spots=args.max_spots, seed=args.seed
            )
        except Exception as exc:  # record, do not silently skip
            errors.append({"section": sec, "error": str(exc)})
            print(f"ERROR {sec}: {exc}", file=sys.stderr, flush=True)
            continue
        result["liana_res"].to_csv(out / "liana_consensus.csv", index=False)
        result["moran_table"].to_csv(out / "moranI.csv", index=False)
        summary = {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "section": sec,
            **result["summary"],
            "claimable_as_nsclc_tissue": True,
            "figure_parity": False,
            "scientific_boundary": BOUNDARY,
        }
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        section_rows.append({"section": sec, "summary": result["summary"]})
        print(json.dumps({k: summary[k] for k in ("n_liana", "n_clusters")}, indent=None), flush=True)

    table, contrast = build_contrast(section_rows, section_map)
    table.to_csv(args.output_dir / "spatial_section_table.csv", index=False)
    (args.output_dir / "spatial_contrast_summary.json").write_text(
        json.dumps(contrast, indent=2) + "\n", encoding="utf-8"
    )

    run_summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "lane": "W9_nsclc_visium_spatial_liana",
        "dataset": "E-MTAB-13530 (De Zuani et al. 2024 Nat Commun 15:4388)",
        "sections_requested": len(sections),
        "sections_ok": len(section_rows),
        "errors": errors,
        "claimable_as_nsclc_tissue": True,
        "figure_parity": False,
        "scientific_boundary": BOUNDARY,
    }
    (args.output_dir / "run_summary.json").write_text(
        json.dumps(run_summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(run_summary, indent=2))
    if not section_rows:
        print("FAIL: no section completed", file=sys.stderr)
        return 3
    if errors:
        print(f"WARN: {len(errors)} section(s) failed (recorded)", file=sys.stderr)
    print("NSCLC VISIUM SPATIAL×LIANA LANE COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
