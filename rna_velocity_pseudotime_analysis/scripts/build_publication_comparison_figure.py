import math
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from scipy.stats import kruskal, mannwhitneyu


def bh_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adj, 1.0)
    return out


def p_to_star(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def build():
    root = Path("/home/zerlinshen/singlecell_factory")
    result_dir = root / "rna_velocity_pseudotime_analysis" / "runtime" / "results" / "lung_with_loom_v25"
    fig_dir = result_dir / "figures"
    table_dir = result_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(result_dir / "integrated.h5ad")
    obs = adata.obs.copy()
    obs["leiden"] = obs["leiden"].astype(str)

    metrics = ["dpt_pseudotime", "velocity_speed"]
    count_by_cluster = obs["leiden"].value_counts()
    valid_clusters = sorted(count_by_cluster[count_by_cluster >= 50].index, key=lambda x: int(x))
    obs = obs[obs["leiden"].isin(valid_clusters)].copy()

    summary_rows = []
    for m in metrics:
        g = obs.groupby("leiden", observed=False)[m]
        s = g.agg(["count", "mean", "median", "std"]).reset_index()
        s["sem"] = s["std"] / np.sqrt(s["count"].clip(lower=1))
        s["metric"] = m
        summary_rows.append(s)
    summary = pd.concat(summary_rows, ignore_index=True)
    summary.to_csv(table_dir / "publication_comparison_summary.csv", index=False)

    ref_cluster = (
        summary[summary["metric"] == "dpt_pseudotime"]
        .sort_values("mean")
        .iloc[0]["leiden"]
    )

    stats_rows = []
    for m in metrics:
        groups = [obs.loc[obs["leiden"] == c, m].values for c in valid_clusters]
        h, p_global = kruskal(*groups)
        pvals = []
        clusters_to_test = []
        ref = obs.loc[obs["leiden"] == ref_cluster, m].values
        for c in valid_clusters:
            if c == ref_cluster:
                continue
            x = obs.loc[obs["leiden"] == c, m].values
            _, p = mannwhitneyu(ref, x, alternative="two-sided")
            pvals.append(p)
            clusters_to_test.append(c)
        qvals = bh_adjust(pvals) if pvals else []
        stats_rows.append(
            {
                "metric": m,
                "comparison": "global_kruskal",
                "cluster": "all",
                "p_value": float(p_global),
                "q_value": float(p_global),
                "star": p_to_star(float(p_global)),
            }
        )
        for c, p, q in zip(clusters_to_test, pvals, qvals):
            stats_rows.append(
                {
                    "metric": m,
                    "comparison": f"{ref_cluster}_vs_{c}",
                    "cluster": c,
                    "p_value": float(p),
                    "q_value": float(q),
                    "star": p_to_star(float(q)),
                }
            )
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(table_dir / "publication_comparison_stats.csv", index=False)

    sns.set_theme(style="whitegrid")
    palette = sns.color_palette("colorblind")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=320, constrained_layout=True)

    config = [
        ("dpt_pseudotime", "DPT Pseudotime (a.u.)", palette[0]),
        ("velocity_speed", "Velocity Speed (a.u.)", palette[1]),
    ]

    for ax, (metric, ylab, color) in zip(axes, config):
        sub = summary[summary["metric"] == metric].set_index("leiden").loc[valid_clusters].reset_index()
        x = np.arange(len(valid_clusters))
        ax.bar(
            x,
            sub["mean"].values,
            yerr=sub["sem"].values,
            color=color,
            edgecolor="black",
            linewidth=0.6,
            capsize=3,
            alpha=0.92,
            label=f"{metric} mean±SEM",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(valid_clusters)
        ax.set_xlabel("Leiden Cluster")
        ax.set_ylabel(ylab)
        ax.set_title("Cluster-level Comparison")
        ax.legend(loc="upper right", frameon=True, fontsize=8)
        max_y = float((sub["mean"] + sub["sem"]).max())
        gap = max(0.02, max_y * 0.06)
        y_top = max_y + gap
        ax.set_ylim(0, y_top + gap * 3.2)

        sig = stats_df[(stats_df["metric"] == metric) & (stats_df["cluster"] != "all")]
        sig_map = {r["cluster"]: r["star"] for _, r in sig.iterrows()}
        for i, c in enumerate(valid_clusters):
            if c == ref_cluster:
                ax.text(i, y_top + gap * 0.35, "ref", ha="center", va="bottom", fontsize=8)
                continue
            star = sig_map.get(c, "ns")
            ax.text(i, y_top + gap * 0.35, star, ha="center", va="bottom", fontsize=8, fontweight="bold")
        global_p = float(stats_df[(stats_df["metric"] == metric) & (stats_df["comparison"] == "global_kruskal")]["q_value"].iloc[0])
        ax.text(
            0.01,
            0.98,
            f"Kruskal-Wallis q={global_p:.2e}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.8),
        )

    fig.suptitle(
        "Lung Carcinoma scRNA-seq: Cluster-wise Pseudotime and Velocity Comparison\n"
        f"Reference Cluster for pairwise tests: {ref_cluster}",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.01,
        0.01,
        "Data source: lung_with_loom_v25/integrated.h5ad; Statistics: Kruskal-Wallis + Mann-Whitney U (BH corrected).",
        fontsize=8,
    )

    out_png = fig_dir / "lung_cluster_comparison_pubgrade.png"
    out_tiff = fig_dir / "lung_cluster_comparison_pubgrade.tiff"
    fig.savefig(out_png, dpi=320, bbox_inches="tight")
    fig.savefig(out_tiff, dpi=320, bbox_inches="tight")
    plt.close(fig)

    with Image.open(out_png) as img:
        dpi = img.info.get("dpi", (0, 0))

    caption = (
        "图X｜肺癌样本各Leiden亚群的DPT伪时序与RNA velocity速度对比。柱形为均值±SEM；"
        "以伪时序最低亚群作为参考，进行Mann-Whitney U两侧检验并作BH校正，星号表示显著性。"
        "整体组间差异采用Kruskal-Wallis检验。数据来源为lung_with_loom_v25整合对象。"
    )
    (result_dir / "caption_pubgrade.txt").write_text(caption, encoding="utf-8")

    checklist = f"""# 发表级图像修改清单

1. 标题与信息完整性  
- 已添加总标题与子标题，明确样本、对比对象与参考组。  
- 专业依据：图题应独立传达研究对象与比较维度，满足期刊图表可读性标准。

2. 坐标轴与单位  
- 已补全X轴“Leiden Cluster”、Y轴“(a.u.)”单位说明。  
- 专业依据：定量图必须给出测量量纲或标准化单位，避免歧义。

3. 图例与数据来源  
- 已添加每面板图例（mean±SEM）及图内数据来源注释。  
- 专业依据：学术与商业报告均要求可追溯来源与统计定义。

4. 分组逻辑与对比组数量  
- 使用n≥50的13个亚群，按Leiden编号统一分组。  
- 专业依据：保证组内样本量稳定，降低小样本波动风险。

5. 统计显著性标记  
- 已加入Kruskal-Wallis全局检验与参考组配对Mann-Whitney U检验，BH校正后以星号标注。  
- 专业依据：多组比较需全局检验+多重校正，控制假阳性率。

6. 色彩无障碍  
- 采用seaborn colorblind调色方案并保持高对比边框。  
- 专业依据：符合常用色觉缺陷可访问性建议（避免红绿混淆）。

7. 分辨率  
- 导出PNG/TIFF均为320 dpi（≥300 dpi）。  
- 实测PNG DPI: {dpi}
- 专业依据：满足多数期刊与商业印刷最小分辨率要求。

8. 输出文件  
- 图像：figures/lung_cluster_comparison_pubgrade.png  
- 图像：figures/lung_cluster_comparison_pubgrade.tiff  
- 图注：caption_pubgrade.txt  
- 统计明细：tables/publication_comparison_stats.csv  
- 汇总表：tables/publication_comparison_summary.csv
"""
    (result_dir / "modification_checklist_pubgrade.md").write_text(checklist, encoding="utf-8")


if __name__ == "__main__":
    build()
