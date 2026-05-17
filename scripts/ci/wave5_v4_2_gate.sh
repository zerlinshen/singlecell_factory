#!/usr/bin/env bash
set -euo pipefail
M="${1:?usage: $0 <ledger.json>}"

python - "$M" <<'PY'
import json
import os
import sys

try:
    import jsonschema
except ImportError:
    print("FAIL jsonschema package not available; activate sc_gpu_stable", file=sys.stderr)
    sys.exit(2)

ledger_path = sys.argv[1]
if not os.path.isfile(ledger_path):
    print(f"FAIL ledger not found: {ledger_path}", file=sys.stderr)
    sys.exit(2)

repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# When invoked as `python - <ledger>`, __file__ is "-"; fall back to CWD discovery.
schema_path = os.path.join(os.getcwd(), "ops/run_ledger/schema/wave5_v4_2.schema.json")
if not os.path.isfile(schema_path):
    for candidate in (
        "/home/zerlinshen/singlecell_factory/ops/run_ledger/schema/wave5_v4_2.schema.json",
    ):
        if os.path.isfile(candidate):
            schema_path = candidate
            break
if not os.path.isfile(schema_path):
    print(f"FAIL schema not found at {schema_path}", file=sys.stderr)
    sys.exit(2)

with open(ledger_path) as f:
    ledger = json.load(f)
with open(schema_path) as f:
    schema = json.load(f)

validator = jsonschema.Draft202012Validator(schema)
errors = list(validator.iter_errors(ledger))
if errors:
    msg = "; ".join(e.message for e in errors[:3])
    print(f"FAIL ledger schema invalid: {msg}", file=sys.stderr)
    sys.exit(2)

if ledger.get("plan_revision") != "v4.2" or ledger.get("validation_posture") != "biology-aware":
    print("FAIL plan_revision/validation_posture mismatch", file=sys.stderr)
    sys.exit(2)

verdicts = []  # list of (ac_id, status, detail)

def gate_threshold(ac_id, value, closed_thr, partial_thr, label):
    if value is None:
        verdicts.append((ac_id, "FAIL", f"{label} missing"))
        return
    if value >= closed_thr:
        verdicts.append((ac_id, "CLOSED", f"{label}={value:.3f} >= {closed_thr:.2f}"))
    elif value >= partial_thr:
        verdicts.append((ac_id, "PARTIAL", f"{label}={value:.3f} in [{partial_thr:.2f},{closed_thr:.2f})"))
    else:
        verdicts.append((ac_id, "FAIL", f"{label}={value:.3f} < {partial_thr:.2f}"))

# AC-VAL-3a panel_recall.
gate_threshold("AC-VAL-3a", ledger.get("panel_recall"), 0.60, 0.40, "panel_recall")

# AC-VAL-3b top50_hit_rate.
# Plan v4.2 §3.4 FAIL response: "If AC-VAL-3b fails (top-50 hit rate < 0.15): escalate to user with
# diagnostic dump. Do NOT retune. Do NOT swap to top-100. Mark CLOSED-PARTIAL." The AC matrix < 0.15
# threshold maps to PARTIAL (with diagnostic), not run-blocking FAIL — the AC is gating with a
# CLOSED-PARTIAL outcome rather than a hard FAIL.
_t50 = ledger.get("top50_hit_rate")
if _t50 is None:
    verdicts.append(("AC-VAL-3b", "FAIL", "top50_hit_rate missing"))
elif _t50 >= 0.30:
    verdicts.append(("AC-VAL-3b", "CLOSED", f"top50_hit_rate={_t50:.3f} >= 0.30"))
elif _t50 >= 0.15:
    verdicts.append(("AC-VAL-3b", "PARTIAL", f"top50_hit_rate={_t50:.3f} in [0.15,0.30) per AC matrix"))
else:
    verdicts.append(("AC-VAL-3b", "PARTIAL", f"top50_hit_rate={_t50:.3f} < 0.15; CLOSED-PARTIAL per §3.4 (escalate; no panel-broadening, no threshold relaxation)"))

# AC-CI-1 cell_type_ari.
# Plan v4.2 binds AC-CI-1 to RNA-only Leiden res=0.3 vs canonical labels (per ari_methodology_note.md).
# Prefer the v4.2 key `cell_type_ari_rna_only_res0_3`; fall back to legacy `cell_type_ari` only when the
# v4.2 key is absent (legacy v3 ledgers that wrote `cell_type_ari` as WNN-PRIMARY would otherwise be
# silently graded against the wrong methodology).
metrics = ledger.get("metrics") or {}
ari = metrics.get("cell_type_ari_rna_only_res0_3") or metrics.get("cell_type_ari")
if isinstance(ari, dict):
    ari_value = ari.get("value")
else:
    ari_value = ari
ari_note = ledger.get("ari_methodology_note_sha")
ari_att = ledger.get("ari_note_attestation") or {}
ari_subagent = ari_att.get("subagent_type", "")
if ari_note and ari_subagent:
    if ari_value is None:
        verdicts.append(("AC-CI-1", "FAIL", "cell_type_ari missing"))
    elif ari_value >= 0.60:
        verdicts.append(("AC-CI-1", "CLOSED", f"cell_type_ari={ari_value:.3f} >= 0.60 (biology-aware path)"))
    elif ari_value >= 0.70:
        # 0.70 is the fallback partial threshold; with biology-aware path < 0.60 is FAIL.
        verdicts.append(("AC-CI-1", "PARTIAL", f"cell_type_ari={ari_value:.3f}"))
    else:
        verdicts.append(("AC-CI-1", "FAIL", f"cell_type_ari={ari_value:.3f} < 0.60 (biology-aware)"))
else:
    if ari_value is None:
        verdicts.append(("AC-CI-1", "FAIL", "cell_type_ari missing; no ari_methodology_note_sha / attestation"))
    elif ari_value >= 0.70:
        verdicts.append(("AC-CI-1", "PARTIAL", f"cell_type_ari={ari_value:.3f} >= 0.70 (fallback partial; methodology note absent)"))
    else:
        verdicts.append(("AC-CI-1", "FAIL", f"cell_type_ari={ari_value:.3f} < 0.70 fallback"))

# AC-VAL-PLOT-1/2/3 binding visual verdicts + figure path existence on disk.
figs = ledger.get("figure_artifacts") or {}
for i in (1, 2, 3):
    ac_id = f"AC-VAL-PLOT-{i}"
    confirmed = ledger.get(f"visual_verdict_human_confirmed_{i}")
    plot_path = figs.get(f"plot{i}_path")
    reasons = []
    if confirmed is not True:
        reasons.append(f"visual_verdict_human_confirmed_{i}={confirmed!r}")
    if not plot_path:
        reasons.append(f"plot{i}_path missing")
    elif not os.path.isfile(plot_path):
        reasons.append(f"plot{i}_path not on disk: {plot_path}")
    if reasons:
        verdicts.append((ac_id, "FAIL", "; ".join(reasons)))
    else:
        verdicts.append((ac_id, "CLOSED", f"plot{i} confirmed and on disk"))

# AC-LEDGER-1 schema validation (already passed by step 2).
verdicts.append(("AC-LEDGER-1", "CLOSED", "schema validated"))

# AC-VAL-3c diagnostic strict overlap caveat.
strict = ledger.get("peak_gene_overlap_top1000_strict")
caveat = ledger.get("peak_gene_overlap_top1000_strict_caveat")
if strict is None:
    # Not reported; no obligation.
    verdicts.append(("AC-VAL-3c", "CLOSED", "strict diagnostic not reported"))
else:
    if caveat and len(caveat.strip()) >= 20:
        verdicts.append(("AC-VAL-3c", "CLOSED", f"strict={strict:.3f} with caveat"))
    else:
        verdicts.append(("AC-VAL-3c", "FAIL", f"strict={strict:.3f} reported without sufficient caveat"))

# Attestation enforcement (FIX-1 non-empty values).
panel_att = ledger.get("panel_review_attestation") or {}
panel_sub = panel_att.get("subagent_type", "")
panel_sess = panel_att.get("session_id", "")
if panel_sub and panel_sess:
    verdicts.append(("ATTESTATION", "CLOSED", "panel_review_attestation has subagent_type + session_id"))
else:
    verdicts.append(("ATTESTATION", "FAIL", "panel_review_attestation missing subagent_type and/or session_id"))

# Print per-AC verdict.
any_fail = False
any_partial = False
for ac_id, status, detail in verdicts:
    print(f"{status:8s} {ac_id}  {detail}")
    if status == "FAIL":
        any_fail = True
    elif status == "PARTIAL":
        any_partial = True

if any_fail:
    sys.exit(2)
if any_partial:
    sys.exit(1)
sys.exit(0)
PY
