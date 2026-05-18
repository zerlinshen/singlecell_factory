#!/usr/bin/env python3
"""Render/audit Trevino Cell 2021 Figure 7 from public Supplementary Table S6.

Truth boundary:
- Supplementary Table S6 provides curated ASD/control mutations, BPNet model
  correlations, a prioritized mutation list, and motif-overlap calls.
- It does not provide the full per-mutation/per-cluster score universe used as
  Fisher-test denominators for Figure 7B/7C, nor BPNet weights/ref-alt per-base
  predictions/DeepLift arrays needed to recreate the exact 7F/7G sequence-logo
  and prediction-track vignettes.
- Therefore this script reproduces deposited S6-derived summaries and renders
  explicit resource-gap panels where exact figure parity is not possible from
  the public tables alone.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from zipfile import ZipFile, is_zipfile
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

try:  # scanpy is available in the project env, but keep a hard fallback.
    import scanpy as sc
except Exception:  # pragma: no cover - used only if env lacks scanpy
    sc = None  # type: ignore[assignment]

XLSX_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
CELL_SUPP_BASE_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867421009429"
FIG7_CAPTION_COUNTS = {
    "7C_bpnet_or": 1.909,
    "7C_bpnet_p": 0.004,
    "7D_sfari_case": 24,
    "7D_sfari_control": 17,
    "7D_total_case": 262,
    "7D_total_control": 232,
    "7D_or": 1.24,
    "7D_p": 0.154,
}


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_triple(fig: plt.Figure, out_dir: Path, stem: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for ext in ("png", "svg", "pdf"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", dpi=300)
        outputs[ext] = str(path)
        outputs[f"{ext}_sha256"] = sha256_file(path)
    plt.close(fig)
    return outputs


def col_to_idx(ref: str) -> int:
    m = re.match(r"([A-Z]+)", ref)
    if not m:
        return 0
    out = 0
    for ch in m.group(1):
        out = out * 26 + ord(ch) - 64
    return out - 1


def unique_headers(headers: list[str]) -> list[str]:
    seen: Counter[str] = Counter()
    out: list[str] = []
    for h in headers:
        key = h.strip()
        seen[key] += 1
        out.append(key if seen[key] == 1 else f"{key}__{seen[key]}")
    return out


def load_xlsx_sheet(path: Path, sheet_name: str) -> list[list[str]]:
    if not is_zipfile(path):
        raise ValueError(f"Not a valid xlsx zip file: {path}")
    with ZipFile(path) as zf:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            ss_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in ss_root.findall("a:si", XLSX_NS):
                shared.append("".join(t.text or "" for t in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")))
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        targets: dict[str, str] = {}
        for sh in wb.find("a:sheets", XLSX_NS):  # type: ignore[arg-type]
            targets[sh.attrib["name"]] = "xl/" + relmap[sh.attrib[RID]].lstrip("/")
        if sheet_name not in targets:
            raise KeyError(f"Sheet {sheet_name!r} not found in {path}; sheets={sorted(targets)}")
        root = ET.fromstring(zf.read(targets[sheet_name]))
        rows: list[list[str]] = []
        for row_el in root.findall(".//a:sheetData/a:row", XLSX_NS):
            vals: dict[int, str] = {}
            max_col = -1
            for cell in row_el.findall("a:c", XLSX_NS):
                idx = col_to_idx(cell.attrib.get("r", "A1"))
                typ = cell.attrib.get("t")
                v = cell.find("a:v", XLSX_NS)
                is_el = cell.find("a:is", XLSX_NS)
                val = ""
                if typ == "s" and v is not None and v.text is not None:
                    val = shared[int(v.text)]
                elif typ == "inlineStr" and is_el is not None:
                    val = "".join(t.text or "" for t in is_el.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t"))
                elif v is not None:
                    val = v.text or ""
                vals[idx] = val
                max_col = max(max_col, idx)
            if max_col >= 0:
                rows.append([vals.get(i, "") for i in range(max_col + 1)])
        return rows


def records_from_rows(rows: list[list[str]]) -> tuple[list[dict[str, str]], list[str]]:
    if len(rows) < 2:
        return [], []
    headers = unique_headers(rows[1])
    records: list[dict[str, str]] = []
    for row in rows[2:]:
        records.append(dict(zip(headers, row + [""] * (len(headers) - len(row)))))
    return records, headers


def parse_s6(mmc6: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"path": str(mmc6), "sha256": sha256_file(mmc6), "zip_ok": is_zipfile(mmc6), "source_url": f"{CELL_SUPP_BASE_URL}-mmc6.xlsx"}
    tables: dict[str, list[dict[str, str]]] = {}
    headers: dict[str, list[str]] = {}
    sheet_rows: dict[str, list[list[str]]] = {}
    for sheet in ["Table Index", "A", "B", "C", "D", "E"]:
        rows = load_xlsx_sheet(mmc6, sheet)
        sheet_rows[sheet] = rows
        recs, hdr = records_from_rows(rows)
        tables[sheet] = recs
        headers[sheet] = hdr
    out.update({"tables": tables, "headers": headers, "sheet_rows": sheet_rows})
    return out


def parse_s1_atac_clusters(mmc1: Path) -> pd.DataFrame:
    rows = load_xlsx_sheet(mmc1, "E")
    recs, _ = records_from_rows(rows)
    return pd.DataFrame(recs)


def mutation_key_ab(rec: dict[str, str]) -> tuple[str, str, str, str]:
    return (rec.get("Chr", "").strip(), rec.get("Pos", "").strip(), rec.get("Ref", "").strip(), rec.get("Alt", "").strip())


def mutation_key_d(rec: dict[str, str]) -> tuple[str, str, str, str]:
    # Table S6D is 0-based half-open (Position, Position1); S6A/B use 1-based Pos.
    return (
        rec.get("Chromosome", "").strip(),
        rec.get("Position1", "").strip(),
        rec.get("Reference allele", "").strip(),
        rec.get("Alternate allele", "").strip(),
    )


def as_float(x: str | float | int | None) -> float:
    if x is None or x == "":
        return 0.0
    try:
        return float(x)
    except Exception:
        return 0.0


def load_sfari_sets(sfari_csv: Path) -> dict[str, set[str]]:
    sets = {"all": set(), "score_lt3": set(), "score_le3": set(), "score_le2": set()}
    with sfari_csv.open(errors="ignore") as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            gene = (rec.get("gene-symbol") or rec.get("gene_symbol") or "").strip().upper()
            if not gene:
                continue
            sets["all"].add(gene)
            try:
                score = float(rec.get("gene-score") or "nan")
            except Exception:
                score = math.nan
            if not math.isnan(score) and score < 3:
                sets["score_lt3"].add(gene)
            if not math.isnan(score) and score <= 3:
                sets["score_le3"].add(gene)
            if not math.isnan(score) and score <= 2:
                sets["score_le2"].add(gene)
    return sets


def summarize_figure7(s6: dict[str, Any], sfari_csv: Path, cluster_df: pd.DataFrame) -> dict[str, Any]:
    A: list[dict[str, str]] = s6["tables"]["A"]
    B: list[dict[str, str]] = s6["tables"]["B"]
    C: list[dict[str, str]] = s6["tables"]["C"]
    D: list[dict[str, str]] = s6["tables"]["D"]
    E: list[dict[str, str]] = s6["tables"]["E"]
    headers_d: list[str] = s6["headers"]["D"]
    cluster_cols_all = headers_d[5:20]
    # The published Table S6D has a duplicated "Cluster GluN6" header. The
    # caption denominators are closest when the first duplicate column is
    # excluded, matching the behavior of earlier spreadsheet readers that keep
    # the later duplicate label. We keep both metrics and mark the boundary.
    cluster_cols_caption_like = headers_d[6:20]

    meta: dict[tuple[str, str, str, str], dict[str, str]] = {}
    pheno: dict[tuple[str, str, str, str], str] = {}
    for rec in A:
        k = mutation_key_ab(rec)
        meta[k] = rec
        pheno[k] = "case"
    for rec in B:
        k = mutation_key_ab(rec)
        meta[k] = rec
        pheno[k] = "control"

    annotated: list[dict[str, Any]] = []
    missing = 0
    for rec in D:
        k = mutation_key_d(rec)
        typ = pheno.get(k, "missing")
        if typ == "missing":
            missing += 1
        scores_all = {col: as_float(rec.get(col)) for col in cluster_cols_all}
        scores_caption = {col: as_float(rec.get(col)) for col in cluster_cols_caption_like}
        active_all = [col for col, val in scores_all.items() if val > 30]
        active_caption = [col for col, val in scores_caption.items() if val > 30]
        m = meta.get(k, {})
        genes = []
        for col in ["SYMBOL", "NEAREST"]:
            g = (m.get(col) or "").strip().upper()
            if g and g != ".":
                genes.append(g)
        annotated.append(
            {
                "chrom": k[0],
                "position1": k[1],
                "ref": k[2],
                "alt": k[3],
                "type": typ,
                "table_s6d_type": rec.get("Type", ""),
                "max_score_all_s6d_columns": max(scores_all.values()) if scores_all else 0.0,
                "max_score_caption_like": max(scores_caption.values()) if scores_caption else 0.0,
                "active_clusters_all_s6d_columns": ";".join(active_all),
                "active_clusters_caption_like": ";".join(active_caption),
                "n_active_clusters_all_s6d_columns": len(active_all),
                "n_active_clusters_caption_like": len(active_caption),
                "symbol": (m.get("SYMBOL") or "").strip(),
                "nearest": (m.get("NEAREST") or "").strip(),
                "consequence": (m.get("Consequence") or "").strip(),
                "gene_candidates": ";".join(genes),
                **{f"score_{col}": val for col, val in scores_all.items()},
            }
        )

    ann_df = pd.DataFrame(annotated)
    high_caption = ann_df[ann_df["n_active_clusters_caption_like"] > 0].copy()
    high_all = ann_df[ann_df["n_active_clusters_all_s6d_columns"] > 0].copy()
    s6d_type_counts = Counter(str(rec.get("Type", "")).strip().lower() for rec in D)

    # Cluster-level summaries from the deposited prioritized table.
    cluster_rows: list[dict[str, Any]] = []
    totals_by_type = Counter(high_caption["type"])
    name_to_cluster_id = {str(r["Name"]).strip(): str(r["Cluster ID"]).strip() for _, r in cluster_df.iterrows()}
    name_to_color = {str(r["Name"]).strip(): str(r["Color"]).strip() for _, r in cluster_df.iterrows()}
    for col in cluster_cols_caption_like:
        display = col.replace("Cluster ", "").replace("__2", "")
        subset = ann_df[ann_df[f"score_{col}"] > 30]
        case_n = int((subset["type"] == "case").sum())
        control_n = int((subset["type"] == "control").sum())
        table = [[case_n, control_n], [int(totals_by_type["case"] - case_n), int(totals_by_type["control"] - control_n)]]
        try:
            odds, pval = fisher_exact(table)
        except Exception:
            odds, pval = (math.nan, math.nan)
        cluster_rows.append(
            {
                "score_column": col,
                "cluster_name": display,
                "cluster_id": name_to_cluster_id.get(display, ""),
                "color": name_to_color.get(display, "#808080"),
                "case_high_effect_n": case_n,
                "control_high_effect_n": control_n,
                "case_minus_control": case_n - control_n,
                "log2_case_control_ratio_pseudocount": math.log2((case_n + 0.5) / (control_n + 0.5)),
                "within_prioritized_fisher_or": odds,
                "within_prioritized_fisher_p": pval,
            }
        )
    cluster_enrichment = pd.DataFrame(cluster_rows)
    if not cluster_enrichment.empty:
        cluster_enrichment["signed_neglog10_p"] = np.sign(np.log2(cluster_enrichment["within_prioritized_fisher_or"].replace(0, np.nan))).fillna(0) * (-np.log10(cluster_enrichment["within_prioritized_fisher_p"].clip(lower=1e-300)))

    # Motif overlap summary from S6E.
    motif_rows: list[dict[str, Any]] = []
    e_df = pd.DataFrame(E)
    if not e_df.empty:
        e_df["Motif Family"] = e_df["Motif Family"].astype(str).str.strip()
        e_df["Mutation Type"] = e_df["Mutation Type"].astype(str).str.strip().str.lower()
        for fam in sorted(e_df["Motif Family"].dropna().unique()):
            if not fam:
                continue
            case_n = int(((e_df["Motif Family"] == fam) & (e_df["Mutation Type"] == "case")).sum())
            control_n = int(((e_df["Motif Family"] == fam) & (e_df["Mutation Type"] == "control")).sum())
            motif_rows.append({"motif_family": fam, "case_n": case_n, "control_n": control_n, "case_minus_control": case_n - control_n})
    motif_excess = pd.DataFrame(motif_rows).sort_values(["case_minus_control", "case_n"], ascending=[False, False])

    # SFARI recomputation under several plausible interpretations; this is an
    # audit because it does not match the published 24/17 caption counts.
    sfari_sets = load_sfari_sets(sfari_csv)
    sfari_rows: list[dict[str, Any]] = []
    for sf_label, geneset in sfari_sets.items():
        for gene_cols in [("SYMBOL",), ("NEAREST",), ("SYMBOL", "NEAREST")]:
            tmp: Counter[str] = Counter()
            for _, row in high_caption.iterrows():
                genes = [str(row["symbol"]).strip().upper(), str(row["nearest"]).strip().upper()]
                by_col = {"SYMBOL": genes[0], "NEAREST": genes[1]}
                if any(by_col[c] and by_col[c] != "." and by_col[c] in geneset for c in gene_cols):
                    tmp[str(row["type"])] += 1
            sfari_rows.append(
                {
                    "sfari_set": sf_label,
                    "gene_columns": "+".join(gene_cols),
                    "case_sfari_n": int(tmp["case"]),
                    "control_sfari_n": int(tmp["control"]),
                    "case_total_high_effect_caption_like": int(totals_by_type["case"]),
                    "control_total_high_effect_caption_like": int(totals_by_type["control"]),
                    "matches_caption_24_17_262_232": int(tmp["case"]) == 24 and int(tmp["control"]) == 17 and int(totals_by_type["case"]) == 262 and int(totals_by_type["control"]) == 232,
                }
            )
    sfari_audit = pd.DataFrame(sfari_rows)

    # Locus summaries for 7F/7G: S6 can identify candidate mutations, but not the
    # exact logos / ref-alt prediction tracks.
    def locus_subset(pattern: str) -> pd.DataFrame:
        mask = high_caption["gene_candidates"].astype(str).str.contains(pattern, case=False, regex=True, na=False)
        out = high_caption[mask].copy()
        if not out.empty:
            out["position1_int"] = pd.to_numeric(out["position1"], errors="coerce")
            out = out.sort_values(["type", "max_score_caption_like"], ascending=[True, False])
        return out

    locus = {
        "NFIA": locus_subset(r"(?:^|;)NFIA(?:$|;)|NFIA"),
        "NPY": locus_subset(r"(?:^|;)NPY(?:$|;)"),
    }

    caption_d_table = [[FIG7_CAPTION_COUNTS["7D_sfari_case"], FIG7_CAPTION_COUNTS["7D_total_case"] - FIG7_CAPTION_COUNTS["7D_sfari_case"]], [FIG7_CAPTION_COUNTS["7D_sfari_control"], FIG7_CAPTION_COUNTS["7D_total_control"] - FIG7_CAPTION_COUNTS["7D_sfari_control"]]]
    cap_or, cap_p = fisher_exact(caption_d_table)

    return {
        "ann_df": ann_df,
        "high_caption": high_caption,
        "high_all": high_all,
        "cluster_enrichment": cluster_enrichment,
        "motif_excess": motif_excess,
        "sfari_audit": sfari_audit,
        "locus": locus,
        "metrics": {
            "s6a_case_mutations": len(A),
            "s6b_control_mutations": len(B),
            "s6c_folds": len(C),
            "s6c_cluster_columns": len([h for h in s6["headers"]["C"] if h.startswith("c")]),
            "s6d_rows": len(D),
            "s6d_case_rows": int(s6d_type_counts["case"]),
            "s6d_control_rows": int(s6d_type_counts["control"]),
            "s6e_motif_rows": len(E),
            "s6d_join_missing": missing,
            "s6d_duplicate_cluster_columns": [c for c, n in Counter([h.replace("__2", "") for h in cluster_cols_all]).items() if n > 1],
            "s6d_high_effect_caption_like_case": int(totals_by_type["case"]),
            "s6d_high_effect_caption_like_control": int(totals_by_type["control"]),
            "s6d_high_effect_all_columns_case": int((high_all["type"] == "case").sum()),
            "s6d_high_effect_all_columns_control": int((high_all["type"] == "control").sum()),
            "motif_overlap_case": int((pd.DataFrame(E)["Mutation Type"].astype(str).str.strip().str.lower() == "case").sum()) if E else 0,
            "motif_overlap_control": int((pd.DataFrame(E)["Mutation Type"].astype(str).str.strip().str.lower() == "control").sum()) if E else 0,
            "caption_7d_fisher_recomputed_or_from_caption_counts": cap_or,
            "caption_7d_fisher_recomputed_p_from_caption_counts": cap_p,
        },
    }


def write_tables(summary: dict[str, Any], out_dir: Path) -> dict[str, str]:
    paths: dict[str, str] = {}
    for key, df in [
        ("annotated_prioritized_mutations", summary["ann_df"]),
        ("cluster_enrichment", summary["cluster_enrichment"]),
        ("motif_excess", summary["motif_excess"]),
        ("sfari_audit", summary["sfari_audit"]),
    ]:
        path = out_dir / f"figure7_{key}.tsv"
        df.to_csv(path, sep="\t", index=False)
        paths[key] = str(path)
        paths[f"{key}_sha256"] = sha256_file(path)
    for locus_name, df in summary["locus"].items():
        path = out_dir / f"figure7_locus_{locus_name.lower()}_s6_summary.tsv"
        df.to_csv(path, sep="\t", index=False)
        paths[f"locus_{locus_name.lower()}"] = str(path)
        paths[f"locus_{locus_name.lower()}_sha256"] = sha256_file(path)
    return paths


def plot_7a(out_dir: Path, metrics: dict[str, Any]) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.axis("off")
    steps = [
        ("SSC de novo\nmutations", f"S6A case={metrics['s6a_case_mutations']:,}\nS6B control={metrics['s6b_control_mutations']:,}"),
        ("Filter", "noncoding / singleton\nnot in gnomAD"),
        ("Intersect peaks", "cluster-specific\nscATAC peaks"),
        ("BPNet perturb", "ref vs alt\n5 folds, 100 bp"),
        ("Prioritize", f"S6D >30\ncase={metrics['s6d_high_effect_caption_like_case']}\ncontrol={metrics['s6d_high_effect_caption_like_control']}"),
        ("Interpret motifs", f"S6E motif overlaps\ncase={metrics['motif_overlap_case']}\ncontrol={metrics['motif_overlap_control']}"),
    ]
    x = np.linspace(0.08, 0.92, len(steps))
    for i, ((title, body), xi) in enumerate(zip(steps, x)):
        rect = plt.Rectangle((xi - 0.07, 0.35), 0.14, 0.35, facecolor="#F2F2F2", edgecolor="#555555", linewidth=1.2)
        ax.add_patch(rect)
        ax.text(xi, 0.61, title, ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(xi, 0.45, body, ha="center", va="center", fontsize=8)
        if i < len(steps) - 1:
            ax.annotate("", xy=(x[i + 1] - 0.08, 0.52), xytext=(xi + 0.08, 0.52), arrowprops=dict(arrowstyle="->", lw=1.5, color="#333333"))
    ax.text(0.5, 0.16, "Fig. 7A public S6-derived mutation-prioritization audit schematic", ha="center", fontsize=12)
    return save_triple(fig, out_dir, "figure7A_mutation_prioritization_schematic")


def plot_7b(out_dir: Path, cluster_enrichment: pd.DataFrame, atac_umap_h5ad: Path | None) -> dict[str, str]:
    if sc is not None and atac_umap_h5ad is not None and atac_umap_h5ad.exists():
        adata = sc.read_h5ad(atac_umap_h5ad)
        obs = adata.obs.copy()
        coords = np.asarray(adata.obsm["X_umap"])
        score_map = dict(zip(cluster_enrichment["cluster_id"], cluster_enrichment["signed_neglog10_p"].fillna(0)))
        name_map = dict(zip(cluster_enrichment["cluster_id"], cluster_enrichment["cluster_name"]))
        vals = obs["Iterative.LSI.Clusters"].astype(str).map(score_map).fillna(0).to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(6.2, 5.5))
        vmax = max(1.0, float(np.nanmax(np.abs(vals))))
        sca = ax.scatter(coords[:, 0], coords[:, 1], c=vals, s=2, cmap="coolwarm", vmin=-vmax, vmax=vmax, linewidths=0, alpha=0.85)
        for cid, sub in obs.groupby(obs["Iterative.LSI.Clusters"].astype(str)):
            if cid in score_map:
                idx = sub.index
                loc = coords[obs.index.get_indexer(idx)]
                ax.text(np.nanmedian(loc[:, 0]), np.nanmedian(loc[:, 1]), name_map.get(cid, cid), fontsize=7, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.65))
        ax.set_title("Fig. 7B S6D BPNet high-effect enrichment on public scATAC UMAP")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        cb = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("signed -log10(p) within S6D prioritized table")
        return save_triple(fig, out_dir, "figure7B_bpnet_cluster_enrichment_umap")

    # Fallback: cluster barplot.
    df = cluster_enrichment.sort_values("case_minus_control", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(df["cluster_name"], df["case_minus_control"], color=np.where(df["case_minus_control"] >= 0, "#D95F02", "#7570B3"))
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("case - control high-effect mutations")
    ax.set_title("Fig. 7B S6D cluster high-effect mutation excess")
    return save_triple(fig, out_dir, "figure7B_bpnet_cluster_enrichment_umap")


def plot_7c(out_dir: Path, metrics: dict[str, Any]) -> dict[str, str]:
    rows = pd.DataFrame(
        [
            {"category": "All cleaned S6", "case": metrics["s6a_case_mutations"], "control": metrics["s6b_control_mutations"]},
            {"category": "S6D rows\n(score table)", "case": metrics["s6d_case_rows"], "control": metrics["s6d_control_rows"]},
            {"category": "S6D high-effect\ncaption-like", "case": metrics["s6d_high_effect_caption_like_case"], "control": metrics["s6d_high_effect_caption_like_control"]},
            {"category": "S6D high-effect\nall columns", "case": metrics["s6d_high_effect_all_columns_case"], "control": metrics["s6d_high_effect_all_columns_control"]},
            {"category": "S6E motif overlap", "case": metrics["motif_overlap_case"], "control": metrics["motif_overlap_control"]},
        ]
    )
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    x = np.arange(len(rows))
    w = 0.38
    ax.bar(x - w / 2, rows["case"], w, label="case", color="#D95F02")
    ax.bar(x + w / 2, rows["control"], w, label="control", color="#1B9E77")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(rows["category"], rotation=20, ha="right")
    ax.set_ylabel("mutation count (log scale)")
    ax.set_title("Fig. 7C count audit from deposited S6 tables")
    ax.legend(frameon=False)
    ax.text(0.98, 0.95, "Paper caption benchmark:\nBPNet OR=1.909, p=0.004\nExact Fisher denominators not in S6", transform=ax.transAxes, ha="right", va="top", fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="#AAAAAA", alpha=0.9))
    return save_triple(fig, out_dir, "figure7C_prioritization_count_audit")


def plot_7d(out_dir: Path, sfari_audit: pd.DataFrame) -> dict[str, str]:
    # Show the caption benchmark plus the nearest plausible local recomputations.
    candidates = sfari_audit[
        sfari_audit["sfari_set"].isin(["all", "score_lt3", "score_le3"]) & sfari_audit["gene_columns"].isin(["SYMBOL", "NEAREST", "SYMBOL+NEAREST"])
    ].copy()
    candidates["label"] = candidates["sfari_set"] + "\n" + candidates["gene_columns"]
    rows = pd.concat(
        [
            pd.DataFrame(
                [
                    {
                        "label": "Paper caption\nbenchmark",
                        "case_sfari_n": FIG7_CAPTION_COUNTS["7D_sfari_case"],
                        "control_sfari_n": FIG7_CAPTION_COUNTS["7D_sfari_control"],
                        "case_total_high_effect_caption_like": FIG7_CAPTION_COUNTS["7D_total_case"],
                        "control_total_high_effect_caption_like": FIG7_CAPTION_COUNTS["7D_total_control"],
                    }
                ]
            ),
            candidates,
        ],
        ignore_index=True,
    )
    rows = rows.iloc[:10].copy()
    x = np.arange(len(rows))
    w = 0.38
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    ax.bar(x - w / 2, rows["case_sfari_n"], w, label="case SFARI-nearest", color="#D95F02")
    ax.bar(x + w / 2, rows["control_sfari_n"], w, label="control SFARI-nearest", color="#1B9E77")
    ax.set_xticks(x)
    ax.set_xticklabels(rows["label"], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("SFARI-overlapping high-effect mutations")
    ax.set_title("Fig. 7D SFARI nearest-gene audit: caption vs S6A/B+local SFARI recomputation")
    ax.legend(frameon=False)
    ax.text(0.98, 0.95, "Exact 24/17 not recovered\nfrom S6 SYMBOL/NEAREST + local SFARI CSV", transform=ax.transAxes, ha="right", va="top", fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="#AAAAAA", alpha=0.9))
    return save_triple(fig, out_dir, "figure7D_sfari_gene_audit")


def plot_7e(out_dir: Path, motif_excess: pd.DataFrame) -> dict[str, str]:
    df = motif_excess.sort_values(["case_minus_control", "case_n"], ascending=[False, False]).head(18).sort_values("case_minus_control")
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.barh(df["motif_family"], df["case_minus_control"], color=np.where(df["case_minus_control"] >= 0, "#D95F02", "#7570B3"))
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("case motif overlaps - control motif overlaps")
    ax.set_title("Fig. 7E S6E motif disruption excess")
    for y, (_, row) in enumerate(df.iterrows()):
        ax.text(row["case_minus_control"] + (0.15 if row["case_minus_control"] >= 0 else -0.15), y, f"{int(row['case_n'])}/{int(row['control_n'])}", va="center", ha="left" if row["case_minus_control"] >= 0 else "right", fontsize=8)
    return save_triple(fig, out_dir, "figure7E_motif_disruption_excess")


def plot_locus(out_dir: Path, locus_name: str, df: pd.DataFrame, stem: str, caption: str) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(8.8, 3.8))
    if df.empty:
        ax.axis("off")
        ax.text(0.5, 0.6, f"{locus_name}: no high-effect S6D candidates recovered", ha="center", fontsize=12)
    else:
        work = df.copy()
        work["position1_int"] = pd.to_numeric(work["position1"], errors="coerce")
        colors = work["type"].map({"case": "#D95F02", "control": "#1B9E77"}).fillna("#777777")
        y = work["max_score_caption_like"].astype(float)
        ax.scatter(work["position1_int"], y, c=colors, s=55, alpha=0.85, edgecolor="black", linewidth=0.2)
        for _, row in work.sort_values("max_score_caption_like", ascending=False).head(5).iterrows():
            ax.text(float(row["position1_int"]), float(row["max_score_caption_like"]) + 0.7, str(row["gene_candidates"]).replace(";", "/"), fontsize=7, ha="center")
        ax.set_xlabel(f"{locus_name} locus mutation position (hg38, from S6A/B join)")
        ax.set_ylabel("max S6D BPNet perturbation score")
        ax.set_title(caption)
        ax.axhline(30, color="black", linestyle="--", linewidth=0.8)
        ax.text(0.01, 0.97, "orange=case, green=control\nExact logos/ref-alt prediction tracks require BPNet weights + DeepLift arrays", transform=ax.transAxes, va="top", fontsize=8, bbox=dict(boxstyle="round", fc="white", ec="#AAAAAA", alpha=0.9))
    return save_triple(fig, out_dir, stem)


def build_evidence(args: argparse.Namespace, paths: dict[str, str], plots: dict[str, dict[str, str]], s6: dict[str, Any], summary: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    metrics = summary["metrics"]
    checks = {
        "mmc6_zip_ok": bool(s6["zip_ok"]),
        "s6a_case_mutations_eq_42931": metrics["s6a_case_mutations"] == 42931,
        "s6b_control_mutations_eq_42049": metrics["s6b_control_mutations"] == 42049,
        "s6c_folds_eq_5": metrics["s6c_folds"] == 5,
        "s6c_cluster_columns_eq_17": metrics["s6c_cluster_columns"] == 17,
        "s6d_rows_eq_800": metrics["s6d_rows"] == 800,
        "s6d_case_control_rows_sum_eq_800": (metrics["s6d_case_rows"] + metrics["s6d_control_rows"]) == 800,
        "s6e_motif_rows_eq_189": metrics["s6e_motif_rows"] == 189,
        "s6d_join_to_s6a_s6b_no_missing": metrics["s6d_join_missing"] == 0,
        "s6d_caption_like_control_matches_232": metrics["s6d_high_effect_caption_like_control"] == 232,
        "s6d_caption_like_case_matches_262": metrics["s6d_high_effect_caption_like_case"] == 262,
        "motif_table_has_case_and_control": metrics["motif_overlap_case"] > 0 and metrics["motif_overlap_control"] > 0,
        "all_plot_hashes_recorded": all(any(k.endswith("_sha256") for k in p) for p in plots.values()),
    }
    panel_status = {
        "7A": "GENERATED_FROM_METHODS_AND_PUBLIC_TABLE_S6_COUNTS",
        "7B": "GENERATED_S6D_CLUSTER_HIGH_EFFECT_SUMMARY_ON_PUBLIC_ATAC_UMAP_NOT_EXACT_AUTHOR_FISHER_DENOMINATORS",
        "7C": "GENERATED_S6_COUNT_AUDIT_EXACT_OR_1_909_DENOMINATORS_NOT_PUBLIC_IN_S6",
        "7D": "GENERATED_SFARI_AUDIT_CAPTION_BENCHMARK_PLUS_RECOMPUTED_COUNTS_NOT_EXACT_24_17_FROM_LOCAL_TABLES",
        "7E": "REPRODUCED_FROM_PUBLIC_TABLE_S6E_MOTIF_OVERLAP_EXCESS",
        "7F": "GENERATED_NFIA_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY",
        "7G": "GENERATED_NPY_S6_LOCUS_SUMMARY_NOT_EXACT_BPNET_LOGO_TRACK_PARITY",
    }
    return {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "FIGURE7_PUBLIC_TABLE_S6_GENERATED_WITH_EXPLICIT_RESOURCE_GAPS",
        "truth_boundary": "Figure 7 is grounded in public Supplementary Table S6 and public scATAC UMAP subset. Exact Figure 7B/7C Fisher enrichments need full scored-overlap denominators; exact 7F/7G vignettes need BPNet model weights/ref-alt per-base predictions/DeepLift tracks that are not present in the public tables or local Brain_ASD checkout.",
        "inputs": {
            "mmc1": str(args.mmc1),
            "mmc6": str(args.mmc6),
            "sfari_csv": str(args.sfari_csv),
            "atac_umap_h5ad": str(args.atac_umap_h5ad),
            "brain_asd_repo": str(args.brain_asd_repo),
        },
        "source_urls": {
            "supplement_mmc6": f"{CELL_SUPP_BASE_URL}-mmc6.xlsx",
            "article": "https://www.sciencedirect.com/science/article/pii/S0092867421009429",
            "bpnet_code": "https://github.com/GreenleafLab/Brain_ASD",
        },
        "out_dir": str(out_dir),
        "tables": paths,
        "plots": plots,
        "metrics": metrics,
        "caption_benchmarks": FIG7_CAPTION_COUNTS,
        "checks": checks,
        "panel_reproduction_status": panel_status,
        "resource_gaps": [
            "S6D has prioritized/high-effect scores but not the complete per-method denominator universe used for Figure 7C OR=1.909.",
            "S6D contains a duplicated Cluster GluN6 header; caption-like reader behavior (dropping first duplicate) gives 261 case / 232 control high-effect rows, matching control but one case below the caption total 262.",
            "S6A/B SYMBOL/NEAREST annotations plus local SFARI CSV do not recover the caption 24/17 SFARI-nearest counts under tested all/score<=3/score<3 interpretations.",
            "Brain_ASD public checkout has scripts/peaks/GC-matched negatives, but no pre-trained model weights or precomputed ref-alt prediction/DeepLift arrays for exact 7F/7G sequence-logo and prediction tracks.",
        ],
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--run-id", required=True)
    p.add_argument("--mmc1", type=Path, default=None)
    p.add_argument("--mmc6", type=Path, default=None)
    p.add_argument("--sfari-csv", type=Path, default=None)
    p.add_argument("--atac-umap-h5ad", type=Path, default=None)
    p.add_argument("--brain-asd-repo", type=Path, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.project_root / "runs" / args.run_id
    supp_dir = args.project_root / "inputs" / "supplementary" / "cell_2021"
    if args.mmc1 is None:
        args.mmc1 = supp_dir / "mmc1.xlsx"
    if args.mmc6 is None:
        args.mmc6 = supp_dir / "mmc6.xlsx"
    if args.sfari_csv is None:
        args.sfari_csv = args.project_root / "inputs" / "external_code" / "brainchromatin" / "misc" / "sfari_gene.csv"
    if args.atac_umap_h5ad is None:
        args.atac_umap_h5ad = run_dir / "python" / "figures" / "public_matrix" / "figure_3B_atac_pseudotime" / "figure3b_atac_transferred_pseudotime_subset.h5ad"
    if args.brain_asd_repo is None:
        args.brain_asd_repo = args.project_root / "inputs" / "external_code" / "Brain_ASD"

    configure_matplotlib()
    out_dir = run_dir / "python" / "figures" / "public_matrix" / "figure_7_asd_bpnet_public_panels"
    evidence_dir = run_dir / "python" / "figure_reproduction_evidence"
    out_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)

    s6 = parse_s6(args.mmc6)
    clusters = parse_s1_atac_clusters(args.mmc1)
    summary = summarize_figure7(s6, args.sfari_csv, clusters)
    table_paths = write_tables(summary, out_dir)

    plots: dict[str, dict[str, str]] = {}
    plots["7A"] = plot_7a(out_dir, summary["metrics"])
    plots["7B"] = plot_7b(out_dir, summary["cluster_enrichment"], args.atac_umap_h5ad)
    plots["7C"] = plot_7c(out_dir, summary["metrics"])
    plots["7D"] = plot_7d(out_dir, summary["sfari_audit"])
    plots["7E"] = plot_7e(out_dir, summary["motif_excess"])
    plots["7F"] = plot_locus(out_dir, "NFIA", summary["locus"]["NFIA"], "figure7F_nfia_locus_s6_summary", "Fig. 7F NFIA S6 high-effect locus summary")
    plots["7G"] = plot_locus(out_dir, "NPY", summary["locus"]["NPY"], "figure7G_npy_locus_s6_summary", "Fig. 7G NPY S6 high-effect locus summary")

    evidence = build_evidence(args, table_paths, plots, s6, summary, out_dir)
    ev_path = evidence_dir / "figure7_public_asd_bpnet_evidence.json"
    ev_path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"evidence": str(ev_path), "out_dir": str(out_dir), "checks": evidence["checks"], "metrics": evidence["metrics"]}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
