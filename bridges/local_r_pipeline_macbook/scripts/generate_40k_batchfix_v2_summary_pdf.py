from __future__ import annotations

from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer


ROOT = Path("/Users/zerlinshen/Downloads/1. Antigravity/R pipeline for scRNA analysis")
PLOT_DIR = ROOT / "output" / "nsclc_40k_7donor_large_batchfix_v2_auto_bundle_rplots"
ASSET_DIR = ROOT / "output" / "nsclc_40k_7donor_large_batchfix_v2_auto_report_assets"
OUT_PDF = ROOT / "output" / "NSCLC_40K_7DONOR_LARGE_BATCHFIX_V2_summary.pdf"


def add_image(story, path: Path, width_cm: float = 15.5, max_height_cm: float = 8.5) -> None:
    if not path.exists():
        return
    img = Image(str(path))
    img.drawWidth = width_cm * cm
    img.drawHeight = img.imageHeight * img.drawWidth / img.imageWidth
    max_h = max_height_cm * cm
    if img.drawHeight > max_h:
        scale = max_h / img.drawHeight
        img.drawHeight *= scale
        img.drawWidth *= scale
    story.append(img)
    story.append(Spacer(1, 0.4 * cm))


def bulletize(lines: list[str]) -> str:
    return "<br/>".join(f"&bull; {line}" for line in lines)


def main() -> None:
    styles = getSampleStyleSheet()
    title = styles["Title"]
    h1 = styles["Heading1"]
    body = styles["BodyText"]
    body.leading = 15

    story = []
    story.append(Paragraph("NSCLC 40K 7-Donor Integrated Summary", title))
    story.append(Spacer(1, 0.4 * cm))
    story.append(
        Paragraph(
            "Remote run: <b>NSCLC_40K_7DONOR_LARGE_BATCHFIX_V2_AUTO</b><br/>"
            "This is the current clean multi-donor reference run after fixing Harmony CPU fallback and shape writeback.",
            body,
        )
    )
    story.append(Spacer(1, 0.4 * cm))

    story.append(Paragraph("Run Summary", h1))
    story.append(
        Paragraph(
            bulletize(
                [
                    "Raw cells: 26,027; cells after QC: 20,881; genes after QC: 28,706.",
                    "Seven donor batches were retained and integrated.",
                    "All requested modules completed successfully, including batch_correction.",
                    "Harmony now runs on CPU by design in normal auto mode on this workstation.",
                    "GPU remains useful in selected modules, but RAPIDS PCA and RAPIDS DE are currently unstable on the installed 5090/CUDA software stack.",
                ]
            ),
            body,
        )
    )
    story.append(Spacer(1, 0.4 * cm))

    story.append(Paragraph("Global Structure and Batch Handling", h1))
    add_image(story, PLOT_DIR / "umap_by_group.png")
    add_image(story, PLOT_DIR / "umap_by_cluster.png")
    add_image(story, ASSET_DIR / "umap_batch_before.png")
    add_image(story, ASSET_DIR / "umap_batch_after.png")

    story.append(PageBreak())
    story.append(Paragraph("QC, Composition, and Immune State", h1))
    add_image(story, PLOT_DIR / "qc_violin.png")
    add_image(story, PLOT_DIR / "qc_scatter.png")
    add_image(story, PLOT_DIR / "cell_fraction.png")
    add_image(story, ASSET_DIR / "umap_immune_subtype.png")

    story.append(PageBreak())
    story.append(Paragraph("Markers and Trajectory", h1))
    add_image(story, PLOT_DIR / "marker_feature.png")
    add_image(story, PLOT_DIR / "marker_dot.png")
    add_image(story, ASSET_DIR / "paga_trajectory.png")

    story.append(PageBreak())
    story.append(Paragraph("Interpretation", h1))
    story.append(
        Paragraph(
            bulletize(
                [
                    "This cohort is already suitable for a first-pass integrated NSCLC tumor microenvironment interpretation.",
                    "The run contains epithelial/tumor, immune, stromal, endothelial, plasma, mast, and myeloid structure.",
                    "Batch correction is now available, so donor-to-donor technical structure should be cleaner than in the earlier non-corrected run.",
                    "Trajectory outputs are useful as exploratory structure but should be refined on lineage-focused subsets before strong biological claims.",
                    "This 40k integrated object is a practical stable reference while the 900k pathway continues to be optimized.",
                ]
            ),
            body,
        )
    )
    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph("Suggested Next Directions", h1))
    story.append(
        Paragraph(
            bulletize(
                [
                    "Subset tumor epithelial cells for a focused tumor-state and CNV-oriented follow-up.",
                    "Subset T/NK/myeloid compartments for cleaner immune-state and exhaustion analysis.",
                    "Keep Harmony on CPU for now; use GPU selectively where validated stable.",
                    "Use this 40k result as the main remote-to-local reporting template for future cohort runs.",
                ]
            ),
            body,
        )
    )

    doc = SimpleDocTemplate(str(OUT_PDF), pagesize=A4, rightMargin=1.4 * cm, leftMargin=1.4 * cm, topMargin=1.4 * cm, bottomMargin=1.4 * cm)
    doc.build(story)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
