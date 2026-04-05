from __future__ import annotations

import ast
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

def _load_dependencies(project_root: Path) -> dict[str, set[str]]:
    pipeline_py = project_root / "workflow" / "modular" / "pipeline.py"
    tree = ast.parse(pipeline_py.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if node.target.id == "MODULE_DEPENDENCIES":
                raw = ast.literal_eval(node.value)
                return {k: set(v) for k, v in raw.items()}
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "MODULE_DEPENDENCIES":
                    raw = ast.literal_eval(node.value)
                    return {k: set(v) for k, v in raw.items()}
    raise RuntimeError("Failed to load MODULE_DEPENDENCIES from pipeline.py")


def _resolve_execution_order(
    module_dependencies: dict[str, set[str]],
    mandatory: list[str],
    optional: list[str],
) -> list[str]:
    all_requested = set(mandatory) | set(optional)
    to_process = list(all_requested)
    while to_process:
        mod = to_process.pop()
        for dep in module_dependencies.get(mod, set()):
            if dep not in all_requested:
                all_requested.add(dep)
                to_process.append(dep)

    in_degree: dict[str, int] = {m: 0 for m in all_requested}
    for mod in all_requested:
        for dep in module_dependencies.get(mod, set()):
            if dep in all_requested:
                in_degree[mod] += 1

    queue = sorted([m for m, d in in_degree.items() if d == 0])
    order: list[str] = []
    while queue:
        node = queue.pop(0)
        order.append(node)
        for mod in sorted(all_requested):
            if node in module_dependencies.get(mod, set()):
                in_degree[mod] -= 1
                if in_degree[mod] == 0:
                    queue.append(mod)
    return order


def _build_tiers(module_dependencies: dict[str, set[str]], order: list[str]) -> list[list[str]]:
    remaining = list(order)
    done: set[str] = set()
    tiers: list[list[str]] = []
    while remaining:
        tier = [m for m in remaining if module_dependencies.get(m, set()).issubset(done)]
        if not tier:
            tier = [remaining[0]]
        tiers.append(sorted(tier))
        done.update(tier)
        remaining = [m for m in remaining if m not in done]
    return tiers


def generate(out_png: Path, out_svg: Path) -> None:
    project_root = out_png.resolve().parents[1]
    module_dependencies = _load_dependencies(project_root)
    order = _resolve_execution_order(
        module_dependencies,
        mandatory=["cellranger", "qc", "doublet_detection"],
        optional=[
            "clustering",
            "cell_cycle",
            "batch_correction",
            "differential_expression",
            "annotation",
            "trajectory",
            "pseudo_velocity",
            "rna_velocity",
            "cnv_inference",
            "pathway_analysis",
            "cell_communication",
            "gene_regulatory_network",
            "validate_cbioportal",
            "immune_phenotyping",
            "tumor_microenvironment",
            "gene_signature_scoring",
            "evolution",
            "pseudobulk_de",
            "cell_fate",
            "composition",
            "metacell",
        ],
    )
    tiers = _build_tiers(module_dependencies, order)

    x_gap = 2.9
    y_gap = 1.1
    node_pos: dict[str, tuple[float, float]] = {}
    for i, tier in enumerate(tiers):
        n = len(tier)
        y0 = (n - 1) * y_gap / 2
        for j, mod in enumerate(tier):
            node_pos[mod] = (i * x_gap, y0 - j * y_gap)

    fig_w = max(14, len(tiers) * 2.2)
    fig_h = max(8, max(len(t) for t in tiers) * 0.9 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off()

    for mod, deps in module_dependencies.items():
        x2, y2 = node_pos[mod]
        for dep in deps:
            x1, y1 = node_pos[dep]
            arrow = FancyArrowPatch(
                (x1 + 0.5, y1),
                (x2 - 0.5, y2),
                arrowstyle="-|>",
                mutation_scale=10,
                linewidth=1.2,
                color="#6b7280",
                alpha=0.9,
            )
            ax.add_patch(arrow)

    for mod, (x, y) in node_pos.items():
        if mod in {"cellranger", "qc", "doublet_detection"}:
            fc = "#d1fae5"
            ec = "#065f46"
        elif mod in {"batch_correction"}:
            fc = "#fee2e2"
            ec = "#991b1b"
        else:
            fc = "#dbeafe"
            ec = "#1e3a8a"
        ax.text(
            x,
            y,
            mod,
            ha="center",
            va="center",
            fontsize=9,
            bbox={
                "boxstyle": "round,pad=0.3,rounding_size=0.15",
                "facecolor": fc,
                "edgecolor": ec,
                "linewidth": 1.0,
            },
        )

    xs = [p[0] for p in node_pos.values()]
    ys = [p[1] for p in node_pos.values()]
    ax.set_xlim(min(xs) - 1.8, max(xs) + 1.8)
    ax.set_ylim(min(ys) - 1.5, max(ys) + 1.5)
    ax.set_title("singlecell_factory modular pipeline dependency flow", fontsize=13, pad=10)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_svg, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    project_root = Path(__file__).resolve().parents[1]
    generate(
        out_png=project_root / "docs" / "module_dependency_flow.png",
        out_svg=project_root / "docs" / "module_dependency_flow.svg",
    )
