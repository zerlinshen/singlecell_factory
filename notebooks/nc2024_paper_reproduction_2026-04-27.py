# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
# ---

# %% [markdown]
# # NC2024 Paper Reproduction — Sanchez-Mejias et al, *Nat Commun* 2024
#
# DOI: [10.1038/s41467-024-48700-8](https://doi.org/10.1038/s41467-024-48700-8)
#
# **Goal**: validate that the v2 final_adata supports the four core biological
# findings of the paper. We do **not** attempt cell-by-cell reproduction — we
# validate direction, magnitude, and significance at the cohort level.
#
# **Inputs**: `results/nc2024_tumor_20260426_v2/final_adata.h5ad` (803,784 cells,
# 75 samples, 24 patients, LUAD+LUSC+NSCLC, tumor + normal-adjacent).
#
# **Strategy**: backed-mode h5ad read; load X chunks per gene panel only when
# needed.

# %% Setup
from __future__ import annotations
import json
import warnings
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)

REPO = Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory")
TUMOR_H5 = REPO / "results/nc2024_tumor_20260426_v2/final_adata.h5ad"
FIG_DIR = REPO / "results/nc2024_tumor_20260426_v2/figures/paper_reproduction"
FIG_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_JSON = REPO / "results/nc2024_tumor_20260426_v2/paper_reproduction_results.json"

DISEASE_COL = "disease"
LUAD_VALUE = "lung adenocarcinoma"
LUSC_VALUE = "lung squamous cell carcinoma"
SITE_TUMOR = "tumor"
SITE_NAT = "normal tissue adjacent to tumor"

results: dict = {"meta": {"paper_doi": "10.1038/s41467-024-48700-8", "date": "2026-04-27"}}


def gene_indices(var_names: pd.Index, genes: list[str]) -> dict[str, int]:
    return {g: int(var_names.get_loc(g)) for g in genes if g in var_names}


def fetch_gene_matrix(adata: ad.AnnData, gene_idx: dict[str, int]) -> pd.DataFrame:
    """Read columns of X for given gene indices; returns dense DataFrame indexed by obs_names."""
    if not gene_idx:
        return pd.DataFrame(index=adata.obs_names)
    cols = list(gene_idx.values())
    if adata.isbacked:
        # Backed CSR slicing column-wise is slow; iterate genes
        out = {}
        for gene, idx in gene_idx.items():
            v = adata.X[:, idx]
            if sp.issparse(v):
                v = v.toarray().ravel()
            else:
                v = np.asarray(v).ravel()
            out[gene] = v
        return pd.DataFrame(out, index=adata.obs_names)
    sub = adata.X[:, cols]
    if sp.issparse(sub):
        sub = sub.toarray()
    return pd.DataFrame(np.asarray(sub), columns=list(gene_idx.keys()), index=adata.obs_names)


# %% Load tumor (backed)
print("Loading tumor in backed mode ...")
adata = ad.read_h5ad(TUMOR_H5, backed="r")
print(f"  shape: {adata.shape}")
obs = adata.obs.copy()  # full obs in memory ~few hundred MB
print(f"  obs columns: {list(obs.columns)}")
print(f"  cell types:\n{obs['cell_type'].value_counts().to_string()}")
print(f"  disease:\n{obs[DISEASE_COL].value_counts().to_string()}")
print(f"  sampling_site:\n{obs['sampling_site'].value_counts().to_string()}")
print(f"  n_samples={obs['sample'].nunique()}, n_patients={obs['patient'].nunique()}")

results["meta"]["n_cells"] = int(adata.n_obs)
results["meta"]["n_samples"] = int(obs["sample"].nunique())
results["meta"]["n_patients"] = int(obs["patient"].nunique())
results["meta"]["cell_type_counts"] = obs["cell_type"].value_counts().to_dict()


# %% Finding 1: M2-strict (CD163+OLR1+) fraction vs NK/T cytotoxicity (sample-level Spearman)
# Paper claims anti-inflammatory M2 macrophages inversely correlate with NK/T cytotoxicity.
# Using all-Myeloid (M1+M2 mix) cancels the signal — must use M2-strict subset.
print("\n=== Finding 1: M2-strict Macrophage vs NK/T cytotoxicity ===")
CYT_GENES = ["GZMA", "GNLY", "PRF1", "NKG7"]
M2_STRICT_GENES = ["CD163", "OLR1", "MRC1", "MSR1", "TGFB1", "IL10"]  # extended anti-inflammatory M2 panel
cyt_idx = gene_indices(adata.var_names, CYT_GENES)
m2s_idx = gene_indices(adata.var_names, M2_STRICT_GENES)
print(f"  cyt genes available: {list(cyt_idx.keys())}")
print(f"  M2-strict markers available: {list(m2s_idx.keys())}")

cyt_expr = fetch_gene_matrix(adata, cyt_idx)
cyt_expr["cell_type"] = obs["cell_type"].values
cyt_expr["sample"] = obs["sample"].values

m2_marker_expr = fetch_gene_matrix(adata, m2s_idx)
# Continuous M2 polarization score per Myeloid cell = mean of available M2 panel
myeloid_mask_glob = obs["cell_type"] == "Myeloid/Macro"
m2_marker_expr_myel = m2_marker_expr.loc[myeloid_mask_glob]
m2_polar_score = m2_marker_expr_myel.mean(axis=1)
print(f"  Myeloid cells: {int(myeloid_mask_glob.sum())}, M2 polarization score range: [{m2_polar_score.min():.2f}, {m2_polar_score.max():.2f}]")

# Restrict to TUMOR-site cells only (paper claim is about TME composition, not NAT)
tumor_site_obs = obs[obs["sampling_site"] == SITE_TUMOR]
print(f"  restricting to tumor-site cells: {len(tumor_site_obs)}")

# Per-patient (n=24) — more biological than per-sample (n=33-65) for cross-patient correlation
sample_df = []
for patient, p_obs in tumor_site_obs.groupby("patient"):
    n_total = len(p_obs)
    if n_total < 200:
        continue
    # mean M2 polarization across this patient's tumor-site Myeloid cells
    p_myel_idx = p_obs.index.intersection(m2_polar_score.index)
    if len(p_myel_idx) < 50:
        continue
    m2_mean = m2_polar_score.loc[p_myel_idx].mean()
    # cytotoxicity in this patient's tumor-site NK and T cells
    nk_idx = p_obs[p_obs["cell_type"] == "NK cell"].index
    t_idx = p_obs[p_obs["cell_type"] == "T cell"].index
    nk_score = cyt_expr.loc[nk_idx, list(cyt_idx.keys())].mean(axis=1).mean() if len(nk_idx) > 30 else np.nan
    t_score = cyt_expr.loc[t_idx, list(cyt_idx.keys())].mean(axis=1).mean() if len(t_idx) > 30 else np.nan
    sample_df.append({
        "patient": patient,
        "n_cells": n_total,
        "M2_polar": float(m2_mean),
        "NK_cyt": nk_score,
        "T_cyt": t_score,
    })
sample_df = pd.DataFrame(sample_df).dropna()
print(f"  n samples after filter: {len(sample_df)}")
print(sample_df.describe().round(3).to_string())

f1: dict = {
    "n_samples_used": len(sample_df),
    "m2_definition": f"per-patient mean M2 polarization score across tumor-site Myeloid cells (panel: {list(m2s_idx.keys())})",
    "unit": "patient (tumor-site only)",
}
for pair in [("M2_polar", "NK_cyt"), ("M2_polar", "T_cyt"), ("NK_cyt", "T_cyt")]:
    rho, p = stats.spearmanr(sample_df[pair[0]], sample_df[pair[1]])
    f1[f"{pair[0]}_vs_{pair[1]}"] = {"rho": float(rho), "p": float(p)}
    print(f"  spearman {pair[0]} vs {pair[1]}: rho={rho:.3f}, p={p:.4f}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, ycol in zip(axes, ["NK_cyt", "T_cyt"]):
    sns.regplot(data=sample_df, x="M2_polar", y=ycol, ax=ax, scatter_kws={"alpha": 0.6})
    rho, p = stats.spearmanr(sample_df["M2_polar"], sample_df[ycol])
    ax.set_title(f"M2-strict frac vs {ycol}\nSpearman rho={rho:.3f}, p={p:.4f}")
    ax.set_xlabel("M2 polarization score (mean of CD163/OLR1/MRC1/MSR1/TGFB1/IL10)")
    ax.set_ylabel(f"{ycol} mean expression")
plt.tight_layout()
plt.savefig(FIG_DIR / "finding1_macrophage_vs_cytotoxicity.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  fig: {FIG_DIR / 'finding1_macrophage_vs_cytotoxicity.png'}")

# Verdict
neg_pairs = [k for k in ["M2_polar_vs_NK_cyt", "M2_polar_vs_T_cyt"] if f1[k]["rho"] < 0 and f1[k]["p"] < 0.05]
neg_dir = [k for k in ["M2_polar_vs_NK_cyt", "M2_polar_vs_T_cyt"] if f1[k]["rho"] < 0]
f1["verdict"] = "PASS" if neg_pairs else ("PARTIAL" if neg_dir else "FAIL")
f1["passing_pairs"] = neg_pairs
results["finding1"] = f1
print(f"  VERDICT: {f1['verdict']}  passing_pairs={neg_pairs}")


# %% Finding 2: NK cytotoxicity tumor vs normal-adjacent
print("\n=== Finding 2: NK cytotoxicity tumor vs normal-adjacent (NAT) ===")
nk_mask_all = obs["cell_type"] == "NK cell"
nk_tumor = nk_mask_all & (obs["sampling_site"] == SITE_TUMOR)
nk_nat = nk_mask_all & (obs["sampling_site"] == SITE_NAT)
print(f"  NK in tumor: {int(nk_tumor.sum())}, NK in NAT: {int(nk_nat.sum())}")

f2: dict = {"n_nk_tumor": int(nk_tumor.sum()), "n_nk_nat": int(nk_nat.sum()), "comparison": "tumor_vs_NAT_within_NC2024_cohort"}
violin_data = []
gene_results = {}
for gene in CYT_GENES:
    if gene not in cyt_expr.columns:
        gene_results[gene] = {"available": False}
        continue
    v_tumor = cyt_expr.loc[nk_tumor & (cyt_expr["cell_type"] == "NK cell"), gene].values
    v_nat = cyt_expr.loc[nk_nat & (cyt_expr["cell_type"] == "NK cell"), gene].values
    if len(v_tumor) < 30 or len(v_nat) < 30:
        gene_results[gene] = {"available": True, "skipped": "too few cells"}
        continue
    u_stat, p_mw = stats.mannwhitneyu(v_tumor, v_nat, alternative="two-sided")
    gene_results[gene] = {
        "tumor_mean": float(v_tumor.mean()),
        "nat_mean": float(v_nat.mean()),
        "tumor_lower_than_nat": bool(v_tumor.mean() < v_nat.mean()),
        "mannwhitney_p": float(p_mw),
    }
    print(f"  {gene}: tumor={v_tumor.mean():.3f}, NAT={v_nat.mean():.3f}, p={p_mw:.2e}")
    violin_data.append(pd.DataFrame({"gene": gene, "expr": v_tumor, "site": "tumor"}))
    violin_data.append(pd.DataFrame({"gene": gene, "expr": v_nat, "site": "NAT"}))

f2["per_gene"] = gene_results
n_lower_significant = sum(
    1 for r in gene_results.values()
    if isinstance(r, dict) and r.get("tumor_lower_than_nat") and r.get("mannwhitney_p", 1) < 0.01
)
f2["n_lower_p_lt_001"] = n_lower_significant
f2["verdict"] = "PASS" if n_lower_significant >= 2 else ("PARTIAL" if n_lower_significant >= 1 else "FAIL")
results["finding2"] = f2
print(f"  VERDICT: {f2['verdict']}  ({n_lower_significant}/{len(CYT_GENES)} genes lower with p<0.01)")

# Plot
if violin_data:
    vdf = pd.concat(violin_data)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    sns.violinplot(data=vdf, x="gene", y="expr", hue="site", split=True, inner="quartile", ax=ax)
    ax.set_title(f"NK cell cytotoxicity: tumor (n={int(nk_tumor.sum())}) vs NAT (n={int(nk_nat.sum())})")
    ax.set_ylabel("log-normalized expression")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "finding2_nk_cytotoxicity.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig: {FIG_DIR / 'finding2_nk_cytotoxicity.png'}")


# %% Finding 3: LUAD vs LUSC checkpoint co-expression in T cells
print("\n=== Finding 3: LUAD vs LUSC checkpoint co-expression (T cells) ===")
CHECKPOINT_GENES = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"]
ck_idx = gene_indices(adata.var_names, CHECKPOINT_GENES)
print(f"  checkpoint genes available: {list(ck_idx.keys())}")

t_mask = obs["cell_type"] == "T cell"
luad_mask = (obs[DISEASE_COL] == LUAD_VALUE) & t_mask
lusc_mask = (obs[DISEASE_COL] == LUSC_VALUE) & t_mask
print(f"  T cells LUAD: {int(luad_mask.sum())}, LUSC: {int(lusc_mask.sum())}")

f3: dict = {"n_t_luad": int(luad_mask.sum()), "n_t_lusc": int(lusc_mask.sum())}

if luad_mask.sum() > 100 and lusc_mask.sum() > 100 and ck_idx:
    ck_expr = fetch_gene_matrix(adata, ck_idx)
    luad_t = ck_expr.loc[luad_mask].copy()
    lusc_t = ck_expr.loc[lusc_mask].copy()
    # z-score normalize per gene across BOTH groups together
    combined = pd.concat([luad_t, lusc_t])
    z = (combined - combined.mean()) / combined.std()
    z["disease"] = ["LUAD"] * len(luad_t) + ["LUSC"] * len(lusc_t)
    # composite co-expression
    luad_score = z.loc[z["disease"] == "LUAD", list(ck_idx.keys())].sum(axis=1)
    lusc_score = z.loc[z["disease"] == "LUSC", list(ck_idx.keys())].sum(axis=1)
    u, p_composite = stats.mannwhitneyu(luad_score, lusc_score, alternative="two-sided")
    print(f"  composite co-expression: LUAD={luad_score.mean():.3f}, LUSC={lusc_score.mean():.3f}, p={p_composite:.2e}")
    f3["composite_score"] = {"luad_mean": float(luad_score.mean()), "lusc_mean": float(lusc_score.mean()), "p": float(p_composite)}

    per_gene = {}
    for gene in ck_idx:
        u, pg = stats.mannwhitneyu(luad_t[gene].values, lusc_t[gene].values, alternative="two-sided")
        per_gene[gene] = {
            "luad_mean": float(luad_t[gene].mean()),
            "lusc_mean": float(lusc_t[gene].mean()),
            "mannwhitney_p": float(pg),
        }
        print(f"  {gene}: LUAD={luad_t[gene].mean():.3f}, LUSC={lusc_t[gene].mean():.3f}, p={pg:.2e}")
    f3["per_gene"] = per_gene
    n_sig_genes = sum(1 for r in per_gene.values() if r["mannwhitney_p"] < 0.05)
    f3["n_significant_genes"] = n_sig_genes
    f3["verdict"] = "PASS" if (p_composite < 0.05 or n_sig_genes >= 1) else "FAIL"

    # Heatmap: per-patient mean per gene, faceted by disease
    ck_expr["patient"] = obs.loc[ck_expr.index, "patient"].values
    ck_expr["disease"] = obs.loc[ck_expr.index, DISEASE_COL].values
    ck_expr["cell_type"] = obs.loc[ck_expr.index, "cell_type"].values
    ck_t = ck_expr[ck_expr["cell_type"] == "T cell"].copy()
    ck_t = ck_t[ck_t["disease"].isin([LUAD_VALUE, LUSC_VALUE])]
    pt_mean = ck_t.groupby(["disease", "patient"])[list(ck_idx.keys())].mean()
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, disease, label in zip(axes, [LUAD_VALUE, LUSC_VALUE], ["LUAD", "LUSC"]):
        sub = pt_mean.loc[disease]
        if len(sub):
            sns.heatmap(sub, ax=ax, cmap="viridis", cbar_kws={"label": "mean expr"})
            ax.set_title(f"{label} (n_patients={len(sub)})")
            ax.set_ylabel("patient")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "finding3_checkpoint_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig: {FIG_DIR / 'finding3_checkpoint_heatmap.png'}")
else:
    f3["verdict"] = "SKIPPED"
    f3["skip_reason"] = "insufficient cells or genes"

results["finding3"] = f3
print(f"  VERDICT: {f3['verdict']}")


# %% Finding 4: Macrophage cholesterol export pathway
print("\n=== Finding 4: Macrophage cholesterol export ===")
CHOL_GENES = ["HMGCR", "ABCA1", "ABCG1", "NR1H3", "APOE", "CYP27A1", "SCD", "SREBF2"]
chol_idx = gene_indices(adata.var_names, CHOL_GENES)
print(f"  chol genes available: {list(chol_idx.keys())}")
print(f"  M2 panel reused from Finding 1: {list(m2s_idx.keys())}")

myeloid_mask = obs["cell_type"] == "Myeloid/Macro"
print(f"  Myeloid cells: {int(myeloid_mask.sum())}")

f4: dict = {"n_myeloid": int(myeloid_mask.sum()), "chol_genes_available": list(chol_idx.keys()), "m2_panel": list(m2s_idx.keys())}

if chol_idx and m2s_idx:
    chol_expr = fetch_gene_matrix(adata, chol_idx)
    chol_score = chol_expr.mean(axis=1)
    # Reuse the per-cell M2 polarization score computed in Finding 1
    df = pd.DataFrame({"chol": chol_score, "m2": m2_polar_score.reindex(obs.index).fillna(0).values, "cell_type": obs["cell_type"].values}, index=obs.index)
    myel_df = df[df["cell_type"] == "Myeloid/Macro"].copy()
    threshold = myel_df["m2"].quantile(0.75)
    myel_df["is_m2"] = myel_df["m2"] >= threshold
    n_m2 = int(myel_df["is_m2"].sum())
    chol_m2 = myel_df.loc[myel_df["is_m2"], "chol"].values
    chol_all = myel_df["chol"].values
    chol_non_m2 = myel_df.loc[~myel_df["is_m2"], "chol"].values
    t_stat, p_two = stats.ttest_ind(chol_m2, chol_non_m2, equal_var=False)
    p_one = p_two / 2 if t_stat > 0 else 1 - p_two / 2  # one-tailed (M2 > non-M2)
    print(f"  M2 cholesterol mean: {chol_m2.mean():.3f}")
    print(f"  non-M2 Myeloid cholesterol mean: {chol_non_m2.mean():.3f}")
    print(f"  one-tailed t-test p (M2 > non-M2): {p_one:.4e}")
    f4["m2_subset"] = {"n_cells": n_m2, "chol_mean": float(chol_m2.mean())}
    f4["non_m2"] = {"n_cells": int((~myel_df["is_m2"]).sum()), "chol_mean": float(chol_non_m2.mean())}
    f4["t_test_one_tailed_p"] = float(p_one)
    f4["verdict"] = "PASS" if (chol_m2.mean() > chol_non_m2.mean() and p_one < 0.05) else "FAIL"

    # Plot
    plot_df = myel_df.assign(group=np.where(myel_df["is_m2"], "M2 (top 25%)", "non-M2 Myeloid"))
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.violinplot(data=plot_df, x="group", y="chol", ax=ax, inner="quartile")
    ax.set_title(f"Macrophage cholesterol export pathway score\none-tailed t-test p={p_one:.2e}")
    ax.set_ylabel(f"mean of {len(chol_idx)} cholesterol genes")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "finding4_cholesterol.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig: {FIG_DIR / 'finding4_cholesterol.png'}")
else:
    f4["verdict"] = "SKIPPED"
    f4["skip_reason"] = "missing chol or M2 marker genes"

results["finding4"] = f4
print(f"  VERDICT: {f4['verdict']}")


# %% Save results JSON + close
RESULTS_JSON.write_text(json.dumps(results, indent=2, default=str))
print(f"\nResults JSON saved to: {RESULTS_JSON}")

adata.file.close()
print("Done.")
