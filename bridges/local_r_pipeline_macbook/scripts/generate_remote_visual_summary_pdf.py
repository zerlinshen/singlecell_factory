from __future__ import annotations
from pathlib import Path
import json
import csv
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer


def add_image(story, path: Path, width_cm: float = 15.5, max_height_cm: float = 8.5) -> None:
    if not path.exists() or path.stat().st_size == 0:
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
    story.append(Spacer(1, 0.3 * cm))


def bulletize(lines: list[str]) -> str:
    return '<br/>'.join(f'&bull; {line}' for line in lines)


def main() -> None:
    import sys
    if len(sys.argv) < 4:
        raise SystemExit('Usage: generate_remote_visual_summary_pdf.py <asset_dir> <title> <out_pdf>')
    asset_dir = Path(sys.argv[1])
    title_text = sys.argv[2]
    out_pdf = Path(sys.argv[3])

    styles = getSampleStyleSheet()
    title = styles['Title']
    h1 = styles['Heading1']
    body = styles['BodyText']
    body.leading = 15

    manifest = {}
    if (asset_dir / 'run_manifest.json').exists():
        manifest = json.loads((asset_dir / 'run_manifest.json').read_text())
    module_rows = []
    if (asset_dir / 'module_status.csv').exists():
        with open(asset_dir / 'module_status.csv', newline='') as f:
            module_rows = list(csv.DictReader(f))

    metadata = manifest.get('metadata', {}) if isinstance(manifest, dict) else {}
    summary_lines = [
        f"Project: {manifest.get('project', 'unknown')}",
        f"Run dir: {manifest.get('run_dir', 'unknown')}",
    ]
    for key in ['raw_cells','cells_after_qc','cells_after_doublet_removal','n_clusters','sample_label_strategy','sample_label_nunique','css_status','css_n_reference_clusters','annotation_unknown_pct','composition_n_groups','composition_n_cell_types','tme_mean_cyt']:
        if key in metadata:
            summary_lines.append(f"{key}: {metadata[key]}")
    if module_rows:
        ok = [r['module'] for r in module_rows if r.get('status') == 'ok']
        summary_lines.append('Successful modules: ' + ', '.join(ok))

    story = []
    story.append(Paragraph(title_text, title))
    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph(bulletize(summary_lines), body))
    story.append(PageBreak())

    pages = [
        ('QC and clustering', [asset_dir/'qc'/'qc_violin_post_filter.png', asset_dir/'qc'/'qc_scatter_post_filter.png', asset_dir/'clustering'/'umap_leiden.png']),
        ('Annotation and composition', [asset_dir/'annotation'/'umap_cell_type.png', asset_dir/'annotation'/'umap_annotation_confidence.png', asset_dir/'annotation'/'cell_type_composition.png', asset_dir/'composition'/'composition_barplot.png', asset_dir/'composition'/'composition_boxplot.png']),
        ('Immune and TME views', [asset_dir/'immune_phenotyping'/'umap_immune_subtype.png', asset_dir/'immune_phenotyping'/'immune_signature_heatmap.png', asset_dir/'tumor_microenvironment'/'tme_cyt_umap.png', asset_dir/'tumor_microenvironment'/'tme_tis_umap.png', asset_dir/'tumor_microenvironment'/'tme_signature_heatmap.png']),
    ]
    for page_title, imgs in pages:
        story.append(Paragraph(page_title, h1))
        story.append(Spacer(1, 0.2 * cm))
        for img in imgs:
            add_image(story, img)
        story.append(PageBreak())

    doc = SimpleDocTemplate(str(out_pdf), pagesize=A4, rightMargin=1.4*cm, leftMargin=1.4*cm, topMargin=1.4*cm, bottomMargin=1.4*cm)
    doc.build(story)
    print(out_pdf)


if __name__ == '__main__':
    main()
