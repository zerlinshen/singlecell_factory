#!/usr/bin/env python3
"""Standalone CLI for quick cBioPortal gene mutation queries.

Usage examples:
  python -m workflow.query_cbioportal --genes ELF3
  python -m workflow.query_cbioportal --genes ELF3 TP53 CDKN2A
  python -m workflow.query_cbioportal --genes ELF3 --study lusc_tcga
  python -m workflow.query_cbioportal --genes ELF3 --output results.csv
  python -m workflow.query_cbioportal --list-studies
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

import requests

BASE_URL = "https://www.cbioportal.org/api"
DEFAULT_STUDY = "lusc_tcga_pan_can_atlas_2018"


def _get_entrez_id(symbol: str, timeout: int) -> int | None:
    try:
        resp = requests.get(f"{BASE_URL}/genes/{symbol}", timeout=timeout)
        if resp.status_code == 200:
            return resp.json().get("entrezGeneId")
    except requests.RequestException as exc:
        print(f"[warn] Could not resolve gene {symbol}: {exc}", file=sys.stderr)
    return None


def _get_study_info(study_id: str, timeout: int) -> dict:
    try:
        resp = requests.get(f"{BASE_URL}/studies/{study_id}", timeout=timeout)
        if resp.status_code == 200:
            return resp.json()
    except requests.RequestException as exc:
        print(f"[warn] Could not fetch study {study_id}: {exc}", file=sys.stderr)
    return {}


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
    except requests.RequestException as exc:
        print(f"[warn] Mutation fetch failed: {exc}", file=sys.stderr)
    return []


def list_lusc_studies(timeout: int = 30) -> None:
    """Print LUSC-related studies available on cBioPortal."""
    try:
        resp = requests.get(
            f"{BASE_URL}/studies",
            params={"keyword": "squamous", "pageSize": 50},
            timeout=timeout,
        )
        studies = resp.json() if resp.status_code == 200 else []
    except requests.RequestException as exc:
        print(f"[error] {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{'Study ID':<48} {'Samples':>8}  Name")
    print("-" * 95)
    for s in studies:
        print(
            f"{s['studyId']:<48} {s.get('allSampleCount', '?'):>8}  "
            f"{s.get('name', '')[:35]}"
        )
    print(f"\n{len(studies)} studies found.")


def query_genes(
    genes: list[str],
    study_id: str = DEFAULT_STUDY,
    timeout: int = 30,
    output: str | None = None,
) -> None:
    """Query mutation status for a list of genes and print a formatted table."""
    study_info = _get_study_info(study_id, timeout)
    total_samples = study_info.get("sequencedSampleCount", 0)
    study_name = study_info.get("name", study_id)
    profile_id = f"{study_id}_mutations"

    if total_samples == 0:
        print(
            f"[error] Study '{study_id}' not found or returned 0 samples.",
            file=sys.stderr,
        )
        print(
            "Tip: run with --list-studies to see available LUSC studies.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"\n=== cBioPortal Mutation Query ===")
    print(f"Study : {study_id}")
    print(f"Name  : {study_name}")
    print(f"Total : {total_samples} samples\n")

    rows = []
    for symbol in genes:
        entrez_id = _get_entrez_id(symbol, timeout)
        if entrez_id is None:
            rows.append(
                {
                    "gene": symbol,
                    "mutated_samples": "N/A",
                    "total_samples": total_samples,
                    "mutation_rate_pct": "N/A",
                    "top_mutation_types": "gene not found",
                }
            )
            continue

        mutations = _fetch_mutations(profile_id, [entrez_id], study_id, timeout)
        mutated = len({m["sampleId"] for m in mutations})
        rate = round(mutated / total_samples * 100, 2) if total_samples else 0.0
        types = Counter(m.get("mutationType", "Unknown") for m in mutations)
        type_str = ", ".join(f"{t}({n})" for t, n in types.most_common(3)) or "—"

        rows.append(
            {
                "gene": symbol,
                "mutated_samples": mutated,
                "total_samples": total_samples,
                "mutation_rate_pct": rate,
                "top_mutation_types": type_str,
            }
        )

    _print_table(rows)

    if output:
        import pandas as pd

        pd.DataFrame(rows).to_csv(output, index=False)
        print(f"\nSaved to {output}")


def _print_table(rows: list[dict]) -> None:
    header = (
        f"{'Gene':<12}  {'Mutated':>8}  {'Total':>8}  {'Rate%':>7}  "
        f"Top Mutation Types"
    )
    print(header)
    print("-" * 90)
    for r in rows:
        print(
            f"{str(r['gene']):<12}  "
            f"{str(r['mutated_samples']):>8}  "
            f"{str(r['total_samples']):>8}  "
            f"{str(r['mutation_rate_pct']):>7}  "
            f"{r['top_mutation_types']}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quick cBioPortal gene mutation query tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m workflow.query_cbioportal --genes ELF3
  python -m workflow.query_cbioportal --genes ELF3 TP53 CDKN2A
  python -m workflow.query_cbioportal --genes ELF3 --study lusc_tcga
  python -m workflow.query_cbioportal --genes ELF3 --output results.csv
  python -m workflow.query_cbioportal --list-studies
""",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        metavar="GENE",
        help="Gene symbols to query, e.g. --genes ELF3 TP53",
    )
    parser.add_argument(
        "--study",
        default=DEFAULT_STUDY,
        metavar="STUDY_ID",
        help=f"cBioPortal study ID (default: {DEFAULT_STUDY})",
    )
    parser.add_argument(
        "--output",
        metavar="FILE",
        help="Save results to CSV file",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        metavar="SEC",
        help="Request timeout in seconds (default: 30)",
    )
    parser.add_argument(
        "--list-studies",
        action="store_true",
        help="List available LUSC-related studies and exit",
    )
    args = parser.parse_args()

    if args.list_studies:
        list_lusc_studies(args.timeout)
        return

    if not args.genes:
        parser.error("--genes is required unless --list-studies is used")

    query_genes(
        genes=args.genes,
        study_id=args.study,
        timeout=args.timeout,
        output=args.output,
    )


if __name__ == "__main__":
    main()
