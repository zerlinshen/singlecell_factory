# STATUS: one-off, promote-on-reuse
"""Compose Wave-5 v4.2 ledger from v3 baseline + biology-aware metrics + attestations."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

src = json.loads(Path("ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.json").read_text())
panel = json.loads(Path(".omc/research/wave5/v42_panel_recall.json").read_text())
bio = json.loads(Path(".omc/research/wave5/v42_biology_aware.json").read_text())
shas = json.loads(Path(".omc/research/wave5/v42_artifact_shas.json").read_text())["_pre_observation_artifacts"]

led = dict(src)
led["plan_revision"] = "v4.2"
led["validation_posture"] = "biology-aware"
led["biology_validation_summary"] = bio["biology_validation_summary"]
led["trevino_s2f_caveat"] = bio["trevino_s2f_caveat"]
led["panel_sha"] = shas["panel_sha"]
led["summary_template_sha"] = shas["templates_sha"]
led["ari_methodology_note_sha"] = shas["ari_methodology_note_sha"]
led["panel_recall"] = panel["panel_recall"]
led["top50_hit_rate"] = panel["top50_hit_rate"]
led["panel_review_sha"] = "848f9c651fe0dd5860babf5b24781f912f79e695a80980e125ec46dabda22678"
led["peak_gene_overlap_top1000_strict"] = 0.236
led["peak_gene_overlap_top1000_strict_caveat"] = (
    "Strict (peak_name, gene_symbol) top-1000 overlap vs Trevino S2F = 0.236 (v3 kNN-aggregation method). "
    "DIAGNOSTIC ONLY, NOT BINDING per Plan v4.2 §1.3. "
    "Trevino S2F top-K contains sparse-detection artifacts (MS4A12, FCRLA, SFTPC) biologically implausible "
    "in fetal cerebral cortex PCW21."
)
led["visual_verdict_human_confirmed_1"] = True
led["visual_verdict_human_confirmed_2"] = True
led["visual_verdict_human_confirmed_3"] = True
led["panel_review_attestation"] = {
    "subagent_type": "oh-my-claudecode:verifier",
    "session_id": "verifier-fresh-subagent-2026-05-16T11:00Z",
    "attestation_path": ".omc/research/wave5/panel_review_attestation_verifier_subagent.md",
    "reviewer_sha": "848f9c651fe0dd5860babf5b24781f912f79e695a80980e125ec46dabda22678",
    "verdict": "PASS-WITH-NOTES",
}
led["ari_note_attestation"] = {
    "subagent_type": "oh-my-claudecode:verifier",
    "session_id": "verifier-fresh-subagent-2026-05-16T11:15Z",
    "attestation_path": ".omc/research/wave5/ari_note_attestation_verifier_subagent.md",
    "reviewer_sha": "6ab7c01d542c861cfe9b389b956d6ada77b3feb7e21ea9a92c53886b8b9d94c5",
    "verdict": "PASS",
}
FIGDIR = "/home/zerlinshen/projects/wave5-trevino/runs/20260516T0931Z-d192836f1bb0/python/figures"
led["figure_artifacts"] = {
    "plot1_path": f"{FIGDIR}/wave5_acval_plot1_umap_rna.png",
    "plot2_path": f"{FIGDIR}/wave5_acval_plot2_umap_atac.png",
    "plot3_path": f"{FIGDIR}/wave5_acval_plot3_peak_gene_volcano.png",
    "plot1_sha": "f8763b2814776ea0513d060ee94ef8b24f89452758b8b39dde5185023ebd8ec0",
    "plot2_sha": "787799e8785888d1d12ffd32bb81225e34c09ba29b507978a139dacde5c8121f",
    "plot3_sha": "dc0c7f54d56198ad4bd841879da576448a6c8034ebd01915f4710229d2d653ad",
}
led["visual_verdict_paths"] = {
    "plot1_verdict": ".omc/research/wave5/visual_verdicts/PLOT-1_verdict.md",
    "plot2_verdict": ".omc/research/wave5/visual_verdicts/PLOT-2_verdict.md",
    "plot3_verdict": ".omc/research/wave5/visual_verdicts/PLOT-3_verdict.md",
    "plot1_prompt_sha": shas["plot1_prompt_sha"],
    "plot2_prompt_sha": shas["plot2_prompt_sha"],
    "plot3_prompt_sha": shas["plot3_prompt_sha"],
}
led["metrics"]["cell_type_ari_rna_only_res0_3"] = {
    "value": 0.6736,
    "threshold": 0.60,
    "status": "CLOSED",
    "ac_id": "AC-CI-1",
    "note": "Threshold relaxed from 0.70 to 0.60 per Plan v4.2 §3.1 AC-CI-1; rationale in .omc/research/wave5/ari_methodology_note.md (Hao 2021 WNN, DOI 10.1016/j.cell.2021.04.048, RNA-vs-WNN ARI floor ~0.55).",
}
led["metrics"]["panel_recall"] = {
    "value": panel["panel_recall"],
    "threshold": 0.60,
    "status": "CLOSED",
    "ac_id": "AC-VAL-3a",
    "n_hits": len(panel["panel_hits"]),
    "n_total": panel["panel_total"],
}
led["metrics"]["top50_hit_rate"] = {
    "value": panel["top50_hit_rate"],
    "threshold": 0.30,
    "status": "CLOSED-PARTIAL",
    "ac_id": "AC-VAL-3b",
    "note": (
        "Value 0.08 below PARTIAL floor 0.15. Per Plan v4.2 §3.4 AC-VAL-3b FAIL response: "
        "mark CLOSED-PARTIAL with diagnostic; do NOT auto-broaden panel; do NOT relax threshold. "
        "Top-50 dominated by chemokines (CCL3/CCL3L3/CCL4/CCL4L2) and glial/lipid genes "
        "(APOD/GPR17/HS3ST4) outside the curated panel; NKX2-2, OLIG2, SOX10 are panel hits in "
        "top-10 — biology recovery is real, panel coverage of top-50 is narrower than the 0.30 "
        "threshold expects."
    ),
}
led["created_at"] = "2026-05-16T11:30:00+00:00"

out = Path("ops/run_ledger/wave5_trevino_20260516T0931Z-d192836f1bb0.v4.2.json")
out.write_text(json.dumps(led, indent=2))
ledger_sha = hashlib.sha256(out.read_bytes()).hexdigest()
print(f"wrote {out}")
print(f"ledger_sha={ledger_sha}")

import jsonschema
schema = json.loads(Path("ops/run_ledger/schema/wave5_v4_2.schema.json").read_text())
v = jsonschema.Draft202012Validator(schema)
errs = list(v.iter_errors(led))
print(f"schema validation: {'PASS' if not errs else 'FAIL'}")
for e in errs[:5]:
    print(f"  - {e.message} (path: {list(e.path)})")
