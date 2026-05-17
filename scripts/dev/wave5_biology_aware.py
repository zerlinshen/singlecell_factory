# STATUS: one-off, promote-on-reuse
#!/usr/bin/env python3
"""Wave-5 biology-aware summary composer.

Composes panel_recall output + strict-overlap value + summary templates into the
biology_validation_summary and trevino_s2f_caveat strings.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import string
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel-recall-json", required=True, type=Path,
                   help="output of wave5_panel_recall.py")
    p.add_argument("--templates-json", required=True, type=Path,
                   help="path to wave5_v4_summary_templates.json")
    p.add_argument("--strict-overlap", required=True, type=float,
                   help="v3 strict overlap value (e.g. 0.236)")
    p.add_argument("--trevino-artifact-examples", type=str,
                   default="MS4A12, FCRLA, SFTPC",
                   help="comma-separated example artifact genes")
    p.add_argument("--out-json", required=True, type=Path,
                   help="output JSON path")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    panel_recall = json.loads(args.panel_recall_json.read_text())
    templates = json.loads(args.templates_json.read_text())

    bio_template = templates["biology_validation_summary"]
    caveat_template = templates["trevino_s2f_caveat"]

    top10_genes = panel_recall.get("top10_genes", [])
    panel_hits = panel_recall.get("panel_hits", [])

    tokens = {
        "panel_recall": f"{panel_recall.get('panel_recall', 0.0):.4f}",
        "panel_hits_count": str(len(panel_hits)),
        "panel_total": str(panel_recall.get("panel_total", 0)),
        "top50_hit_rate": f"{panel_recall.get('top50_hit_rate', 0.0):.4f}",
        "top10_genes": ", ".join(top10_genes),
        "panel_sha": panel_recall.get("panel_sha", ""),
        "strict_overlap": f"{args.strict_overlap:.4f}",
        "trevino_artifact_examples": args.trevino_artifact_examples,
    }

    bio_str = string.Template(bio_template).safe_substitute(tokens)
    caveat_str = string.Template(caveat_template).safe_substitute(tokens)

    templates_sha = hashlib.sha256(args.templates_json.read_bytes()).hexdigest()

    out = {
        "biology_validation_summary": bio_str,
        "trevino_s2f_caveat": caveat_str,
        "_token_substitutions": tokens,
        "_templates_sha": templates_sha,
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(out, indent=2))
    print(f"[wave5_biology_aware] wrote {args.out_json}")


if __name__ == "__main__":
    main()
