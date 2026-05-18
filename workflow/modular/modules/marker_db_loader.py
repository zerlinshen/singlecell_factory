"""Marker DB loader: vendors offline marker DBs into adata.uns["marker_db_index"].

Resolves precedence: project override > curated YAML > DB consensus.
No online I/O. All data read from ref/markers/<db>/<version>/markers.tsv.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)

__references__ = {
    "CellMarker2": {
        "title": "CellMarker 2.0: an updated database of manually curated cell markers in human/mouse",
        "authors": "Hu et al.",
        "journal": "Nucleic Acids Research",
        "year": "2023",
        "doi": "10.1093/nar/gkac947",
        "description": "Curated cell type marker database used for tissue/condition-aware annotation",
    },
    "PanglaoDB": {
        "title": "PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data",
        "authors": "Franzen et al.",
        "journal": "Database",
        "year": "2019",
        "doi": "10.1093/database/baz046",
        "description": "Broad cell type marker resource used as consensus fallback",
    },
    "CellTypist": {
        "title": "Integrated single-cell atlases reveal an oral SARS-CoV-2 infection and a persistently activated immune axis in long COVID",
        "authors": "Dominguez Conde et al.",
        "journal": "Science",
        "year": "2022",
        "doi": "10.1126/science.abq1006",
        "description": "CellTypist model-based marker references",
    },
    "scTypeDB": {
        "title": "Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data",
        "authors": "Ianevski et al.",
        "journal": "Nature Communications",
        "year": "2022",
        "doi": "10.1038/s41467-022-28803-w",
        "description": "scType marker DB used for condition-specific sub-state detection",
    },
}

# Factory-relative path to the offline marker asset root.
_REF_MARKERS_ROOT = Path(__file__).resolve().parent.parent.parent.parent / "ref" / "markers"

# Curated-overlay directory (precedence above all DBs; matches plan's
# "project override > curated YAML > DB consensus" rule).
_CURATED_DIR = _REF_MARKERS_ROOT / "curated"

# Ordered list of DB names; first match wins in consensus resolution.
# "curated" is loaded first (highest priority) when matching YAMLs exist.
_DB_PRIORITY = ["curated", "cellmarker2", "panglaodb", "celltypist", "sctypedb"]

# Latest pinned version per DB (updated when a new snapshot is vendored).
_DB_VERSIONS: dict[str, str] = {
    "curated": "v1",
    "cellmarker2": "v2.0",
    "panglaodb": "v2021-03-18",
    "celltypist": "v2.6.0",
    "sctypedb": "v1.0",
}


def _load_curated_yaml(tissue: str, condition: str) -> pd.DataFrame | None:
    """Load curated marker overlays from ref/markers/curated/*.yaml.

    YAML schema:
        tissue: lung
        condition: NSCLC
        markers:
          <cell_type>:
            positive: [GENE1, GENE2, ...]
            negative: [GENE3, ...]   # optional

    Returns a DataFrame matching the DB-TSV schema (tissue, condition,
    cell_type, marker_gene, marker_type) or None if no matching curated
    file exists for the requested (tissue, condition).
    """
    if not _CURATED_DIR.exists():
        return None
    try:
        import yaml
    except ImportError:
        logger.warning("MarkerDbLoader: pyyaml not installed — skipping curated overlays")
        return None

    rows: list[dict] = []
    for yml in sorted(_CURATED_DIR.glob("*.yaml")):
        try:
            doc = yaml.safe_load(yml.read_text(encoding="utf-8"))
        except Exception as exc:
            logger.warning("MarkerDbLoader: failed to parse %s: %s", yml, exc)
            continue
        if not isinstance(doc, dict):
            continue
        if str(doc.get("tissue", "")).lower() != tissue.lower():
            continue
        if str(doc.get("condition", "")).lower() != condition.lower():
            continue
        markers = doc.get("markers") or {}
        if not isinstance(markers, dict):
            continue
        for cell_type, lists in markers.items():
            if not isinstance(lists, dict):
                continue
            for marker_type in ("positive", "negative"):
                for gene in lists.get(marker_type, []) or []:
                    rows.append({
                        "tissue": tissue,
                        "condition": condition,
                        "cell_type": cell_type,
                        "marker_gene": gene,
                        "marker_type": marker_type,
                        "source_file": yml.name,
                    })
    if not rows:
        return None
    return pd.DataFrame(rows)


class MarkerDbLoaderModule:
    """Load offline-vendored marker DBs and index them by (tissue, condition).

    Writes adata.uns["marker_db_index"]: dict keyed by db_name, each value a
    DataFrame with columns [tissue, condition, cell_type, marker_gene, marker_type].

    Does NOT depend on clustering; clustering dependency belongs exclusively on
    context_aware_annotation (plan P1A.S2 improvement #1).
    """

    name = "marker_db_loader"
    required = False
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"uns": ["marker_db_index"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ValueError(f"{self.name} requires loaded AnnData.")

        tissue = getattr(ctx.cfg, "tissue", None) or "lung"
        condition = getattr(ctx.cfg, "condition", None) or "NSCLC"
        logger.info("MarkerDbLoader: resolving markers for tissue=%s condition=%s", tissue, condition)

        ref_root = _REF_MARKERS_ROOT
        index: dict[str, pd.DataFrame] = {}
        loaded_dbs: list[str] = []
        missing_dbs: list[str] = []

        for db_name in _DB_PRIORITY:
            version = _DB_VERSIONS[db_name]

            # Curated-overlay branch: read YAMLs from ref/markers/curated/
            if db_name == "curated":
                filtered = _load_curated_yaml(tissue, condition)
                if filtered is None or filtered.empty:
                    logger.info(
                        "MarkerDbLoader: no curated YAML for (%s, %s) — skipping curated overlay",
                        tissue, condition,
                    )
                    continue
                filtered["source_db"] = db_name
                filtered["source_version"] = version
                filtered["snapshot_date"] = ""
                index[db_name] = filtered.reset_index(drop=True)
                loaded_dbs.append(f"{db_name}/{version}")
                logger.info(
                    "MarkerDbLoader: loaded %d entries from curated overlay for (%s, %s)",
                    len(filtered), tissue, condition,
                )
                continue

            tsv_path = ref_root / db_name / version / "markers.tsv"
            manifest_path = ref_root / db_name / version / "MANIFEST.json"

            if not tsv_path.exists():
                logger.warning("MarkerDbLoader: missing TSV for %s/%s at %s", db_name, version, tsv_path)
                missing_dbs.append(db_name)
                continue

            try:
                df = pd.read_csv(tsv_path, sep="\t", dtype=str)
            except Exception as exc:
                logger.warning("MarkerDbLoader: failed to read %s: %s", tsv_path, exc)
                missing_dbs.append(db_name)
                continue

            # Filter to requested (tissue, condition) — case-insensitive match
            mask = (
                df["tissue"].str.lower() == tissue.lower()
            ) & (
                df["condition"].str.lower() == condition.lower()
            )
            filtered = df[mask].copy()
            if filtered.empty:
                logger.info(
                    "MarkerDbLoader: %s/%s has no entries for (%s, %s) — skipping",
                    db_name, version, tissue, condition,
                )
                continue

            # Attach provenance columns
            filtered["source_db"] = db_name
            filtered["source_version"] = version

            # Read manifest for version provenance
            if manifest_path.exists():
                try:
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    filtered["snapshot_date"] = manifest.get("snapshot_date", "")
                except Exception:
                    filtered["snapshot_date"] = ""
            else:
                filtered["snapshot_date"] = ""

            index[db_name] = filtered.reset_index(drop=True)
            loaded_dbs.append(f"{db_name}/{version}")
            logger.info(
                "MarkerDbLoader: loaded %d entries from %s/%s for (%s, %s)",
                len(filtered), db_name, version, tissue, condition,
            )

        if not index:
            logger.warning(
                "MarkerDbLoader: no marker entries found for (%s, %s) across all DBs. "
                "Check ref/markers/COVERAGE.md for supported (tissue, condition) pairs.",
                tissue, condition,
            )
            ctx.status(
                self.name,
                "skipped",
                f"No marker entries for ({tissue}, {condition}) in any vendored DB",
            )
            ctx.adata.uns["marker_db_index"] = {}
            ctx.metadata["marker_db_loader_status"] = "skipped_no_coverage"
            ctx.metadata["marker_db_tissue"] = tissue
            ctx.metadata["marker_db_condition"] = condition
            return

        ctx.adata.uns["marker_db_index"] = index
        ctx.metadata["marker_db_loaded"] = loaded_dbs
        ctx.metadata["marker_db_missing"] = missing_dbs
        ctx.metadata["marker_db_tissue"] = tissue
        ctx.metadata["marker_db_condition"] = condition
        ctx.metadata["marker_db_loader_status"] = "ok"
        ctx.metadata["marker_db_total_entries"] = sum(len(df) for df in index.values())

        logger.info(
            "MarkerDbLoader: loaded %d DBs, %d total marker entries for (%s, %s)",
            len(loaded_dbs), ctx.metadata["marker_db_total_entries"], tissue, condition,
        )
