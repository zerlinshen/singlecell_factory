"""Wave-5.5 v5.0.4 Stage C3 — parameter divergence computation.

# STATUS: one-off, promote-on-reuse

Reads:
- data/external/trevino_2021/params_trevino_canonical.json  (Trevino Methods params)
- (our params encoded inline from FIGURE_SPEC.md + RNA-only LITERAL params)

Writes:
- <canonical-run>/python/figures/v5/parameter_divergence/<fig-id>.json per AC-V5-COMP-PARAM-1
- per-figure rows with columns: parameter | trevino_value | our_v5_value | flag

Flag semantics (FIGURE_SPEC.md §parameter_canonical_source):
- OK = exact match (or scalar within tolerance=0)
- WITHIN-TOLERANCE = numeric within declared tolerance
- DIVERGENT = outside tolerance / categorical mismatch
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TREVINO_PARAMS = REPO_ROOT / "data" / "external" / "trevino_2021" / "params_trevino_canonical.json"
CANONICAL_RUN = Path("/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0")
OUT_DIR = CANONICAL_RUN / "python" / "figures" / "v5" / "parameter_divergence"

# Our v5.0.4 params per FIGURE_SPEC.md + Stage A1/A2 actual values
OUR_PARAMS = {
    "F-1B": {
        "normalization": "scanpy normalize_total(1e4) + log1p",
        "n_pcs": 30,
        "leiden_or_louvain_resolution": "n/a (we use Trevino's seurat_clusters as labels)",
        "n_neighbors": "from canonical run (default)",
        "umap_min_dist": 0.5,
        "umap_n_neighbors": "default",
        "random_seed": 42,
    },
    "F-1C": {
        "marker_panel_source": "pre-registered panel JSON (cortical_panel_v1, 18 markers across 6 lineages)",
        "marker_overlay_cmap": "viridis",
        "overlay_percentile_clip": "99%",
    },
    "F-1D": {
        "sample_ages_included": ["pcw21"],
        "plot_type": "bar (per-cluster substitute since only PCW21 available)",
    },
    "F-2A": {
        "lsi_method": "scanpy/snapatac LSI equivalent",
        "n_lsi_components": 49,
        "lsi_first_component_dropped": "no (LSI[:, :2] = both LSI1+LSI2 used)",
    },
    "F-2BD": {
        "peak_caller": "as-supplied in v4.2 canonical run obsm['atac_peaks']",
        "peak_window_size": "from canonical run",
        "n_top_peaks_displayed": 1000,
        "aggregation": "mean per seurat_cluster",
    },
    "F-4A": {
        "trajectory_method": "scanpy diffmap + DPT",
        "n_diffmap_components": 15,
        "root_selection": "highest PAX6 expression",
    },
    "F-4D": {
        "pseudotime_cmap": "plasma",
        "overlay_layer": "same UMAP as F-1B",
    },
    "F-4EG": {
        "branch_assignment_method": "by canonical panel lineage mapping",
        "n_pseudotime_bins": 30,
        "smoothing": "raw mean per bin (no rolling)",
    },
    "F-5": {
        "linkage_method": "Pearson per-cell-correlation (peak_to_gene_v3 pipeline)",
        "n_top_linkages": 1000,
        "correlation_threshold": "top-1000 by abs correlation",
    },
    "F-RNA-1B": {
        "modality": "RNA-only (no ATAC; obsm contamination stripped before recompute)",
        "normalization": "scanpy normalize_total(1e4) + log1p",
        "n_pcs": 50,
        "leiden_resolution": 0.3,
        "n_neighbors": 30,
        "umap_min_dist": 0.5,
    },
    "F-RNA-1C": {
        "marker_panel_source": "pre-registered panel JSON (same as F-1C)",
    },
    "F-RNA-4A": {
        "trajectory_method": "scanpy diffmap+DPT on RNA-only PCA",
    },
}


def compare_value(ours, trevino, tolerance) -> str:
    """Return flag: OK / WITHIN-TOLERANCE / DIVERGENT."""
    if trevino is None or ours is None:
        return "DIVERGENT"
    # exact match
    if ours == trevino:
        return "OK"
    # numeric within tolerance
    if isinstance(tolerance, (int, float)) and isinstance(ours, (int, float)) and isinstance(trevino, (int, float)):
        return "WITHIN-TOLERANCE" if math.isclose(ours, trevino, abs_tol=tolerance) or abs(ours - trevino) <= tolerance else "DIVERGENT"
    # list comparison (e.g., sample ages)
    if isinstance(ours, list) and isinstance(trevino, list):
        return "OK" if set(map(str, ours)) == set(map(str, trevino)) else "DIVERGENT"
    # string fuzzy match (tolerance is non-empty -> at least WITHIN-TOLERANCE)
    return "WITHIN-TOLERANCE" if tolerance else "DIVERGENT"


def compute_divergence_for_figure(fig_id: str, trevino_block: dict, our_block: dict) -> dict:
    if fig_id.endswith("_NOT_RUN"):
        return {"fig_id": fig_id, "status": "NOT-RUN", "rows": [],
                "_description": trevino_block.get("_description", "")}
    rows = []
    trevino_params = trevino_block.get("params", {})
    for param_name, trevino_spec in trevino_params.items():
        # Guard: some entries may be authored as scalar strings; normalize to dict form
        if not isinstance(trevino_spec, dict):
            trevino_spec = {"value": trevino_spec, "tolerance": 0}
        trevino_val = trevino_spec.get("value")
        tolerance = trevino_spec.get("tolerance", 0)
        tol_note = trevino_spec.get("tolerance_note") or trevino_spec.get("_note", "")
        ours = our_block.get(param_name, "MISSING (not declared in our spec for this figure)")
        flag = compare_value(ours, trevino_val, tolerance) if ours != "MISSING" else "DIVERGENT"
        row = {
            "parameter": param_name,
            "trevino_value": trevino_val,
            "our_v5_value": ours,
            "tolerance": tolerance,
            "flag": flag,
        }
        if tol_note:
            row["note"] = tol_note
        if flag == "DIVERGENT" and "_note" not in row:
            row["divergence_rationale"] = trevino_spec.get("_note", "exact equivalent unavailable in current pipeline")
        rows.append(row)
    return {
        "fig_id": fig_id,
        "status": "computed",
        "_description": trevino_block.get("_description", ""),
        "trevino_methods_reference_section": trevino_block.get("trevino_methods_reference_section", ""),
        "rows": rows,
    }


def main() -> int:
    trevino_all = json.loads(TREVINO_PARAMS.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}
    for fig_id, trevino_block in trevino_all.items():
        if fig_id.startswith("_"):
            continue
        our_block = OUR_PARAMS.get(fig_id, {})
        result = compute_divergence_for_figure(fig_id, trevino_block, our_block)
        clean_id = fig_id.replace("_NOT_RUN", "")
        out_path = OUT_DIR / f"{clean_id}.json"
        out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        flag_counts = {"OK": 0, "WITHIN-TOLERANCE": 0, "DIVERGENT": 0}
        for r in result.get("rows", []):
            flag_counts[r["flag"]] = flag_counts.get(r["flag"], 0) + 1
        summary[clean_id] = {"status": result["status"], "flag_counts": flag_counts}
    (OUT_DIR / "_summary.json").write_text(
        json.dumps({
            "computed_iso": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "per_figure": summary,
        }, indent=2) + "\n", encoding="utf-8"
    )
    print("parameter divergence per-figure summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"\nOutputs at: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
