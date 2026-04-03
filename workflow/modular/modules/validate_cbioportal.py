from __future__ import annotations

import json
from collections import Counter
from typing import Any

import pandas as pd
import requests

from ..context import PipelineContext

BASE_URL = "https://www.cbioportal.org/api"


def _get_entrez_id(symbol: str, timeout: int) -> int | None:
    try:
        resp = requests.get(f"{BASE_URL}/genes/{symbol}", timeout=timeout)
        if resp.status_code == 200:
            return resp.json().get("entrezGeneId")
    except requests.RequestException:
        pass
    return None


def _get_study_sample_count(study_id: str, timeout: int) -> int:
    try:
        resp = requests.get(f"{BASE_URL}/studies/{study_id}", timeout=timeout)
        if resp.status_code == 200:
            return resp.json().get("sequencedSampleCount", 0)
    except requests.RequestException:
        pass
    return 0


def _fetch_mutations(profile_id: str, entrez_ids: list[int], study_id: str, timeout: int) -> list[dict]:
    try:
        resp = requests.post(
            f"{BASE_URL}/molecular-profiles/{profile_id}/mutations/fetch",
            params={"projection": "SUMMARY"},
            json={"entrezGeneIds": entrez_ids, "sampleListId": f"{study_id}_sequenced"},
            timeout=timeout,
        )
        if resp.status_code == 200:
            return resp.json()
    except requests.RequestException:
        pass
    return []


def query_gene_mutations(
    gene_symbol: str,
    study_id: str,
    total_samples: int,
    timeout: int,
) -> dict[str, Any]:
    """Query mutation frequency for one gene from a cBioPortal study."""
    entrez_id = _get_entrez_id(gene_symbol, timeout)
    if entrez_id is None:
        return {
            "gene": gene_symbol,
            "study_id": study_id,
            "entrez_id": None,
            "total_samples": total_samples,
            "mutated_samples": 0,
            "mutation_rate_pct": 0.0,
            "mutation_types": {},
            "error": "gene symbol not found",
        }

    profile_id = f"{study_id}_mutations"
    mutations = _fetch_mutations(profile_id, [entrez_id], study_id, timeout)
    mutated_samples = len({m["sampleId"] for m in mutations})
    mutation_types = dict(
        Counter(m.get("mutationType", "Unknown") for m in mutations).most_common()
    )

    return {
        "gene": gene_symbol,
        "study_id": study_id,
        "entrez_id": entrez_id,
        "total_samples": total_samples,
        "mutated_samples": mutated_samples,
        "mutation_rate_pct": round(mutated_samples / total_samples * 100, 2) if total_samples else 0.0,
        "mutation_types": mutation_types,
        "error": None,
    }


class ValidateCbioPortalModule:
    """Optional module: validate gene mutation status via cBioPortal API.

    Gene list is resolved from two sources (merged and deduplicated):
      1. Explicit genes provided via cfg.cbioportal.genes
      2. Top-scoring DE genes from marker_genes.csv (if cfg.cbioportal.use_de_genes=True)

    Outputs:
      tables/cbioportal_mutation_summary.csv   — per-gene mutation frequencies
      tables/cbioportal_validation_report.json — summary report
    """

    name = "validate_cbioportal"

    def run(self, ctx: PipelineContext) -> None:
        cfg = ctx.cfg.cbioportal
        study_id = cfg.study_id
        timeout = cfg.timeout

        genes: list[str] = list(cfg.genes)
        if cfg.use_de_genes:
            de_dir = ctx.module_output_dir("differential_expression")
            de_csv = de_dir / "marker_genes.csv" if de_dir else None
            if de_csv is not None and de_csv.exists():
                de_df = pd.read_csv(de_csv)
                top_de = (
                    de_df.sort_values("scores", ascending=False)["names"]
                    .drop_duplicates()
                    .head(cfg.top_n_de_genes)
                    .tolist()
                )
                # Preserve explicit genes first, append DE genes not already listed.
                seen = set(genes)
                for g in top_de:
                    if g not in seen:
                        genes.append(g)
                        seen.add(g)

        if not genes:
            raise ValueError(
                "No genes to validate. Provide --cbioportal-genes or run "
                "differential_expression first."
            )

        total_samples = _get_study_sample_count(study_id, timeout)
        if total_samples == 0:
            raise ValueError(
                f"Study '{study_id}' returned 0 samples. "
                "Check the study ID with: python -m workflow.query_cbioportal --list-studies"
            )

        results = [
            query_gene_mutations(gene, study_id, total_samples, timeout)
            for gene in genes
        ]

        summary_rows = [
            {
                "gene": r["gene"],
                "study_id": r["study_id"],
                "total_samples": r["total_samples"],
                "mutated_samples": r["mutated_samples"],
                "mutation_rate_pct": r["mutation_rate_pct"],
                "mutation_types": json.dumps(r["mutation_types"]),
                "error": r.get("error"),
            }
            for r in results
        ]
        pd.DataFrame(summary_rows).to_csv(
            ctx.table_dir / "cbioportal_mutation_summary.csv", index=False
        )

        report = {
            "study_id": study_id,
            "total_samples": total_samples,
            "genes_queried": len(results),
            "genes_with_mutations": sum(1 for r in results if r["mutated_samples"] > 0),
            "results": [{k: v for k, v in r.items()} for r in results],
        }
        (ctx.table_dir / "cbioportal_validation_report.json").write_text(
            json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        ctx.metadata["cbioportal_validation"] = {
            "study_id": study_id,
            "total_samples": total_samples,
            "genes_queried": len(results),
            "genes_with_mutations": report["genes_with_mutations"],
        }
