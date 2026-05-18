"""CI verifier — ensure every public module has a non-empty __references__.

Walks workflow/modular/modules/*.py (excluding private modules starting with
underscore and the module_template). For each, asserts:
  1. Module has a top-level __references__ assignment (dict literal).
  2. Dict is non-empty.
  3. Each entry has a valid identifier (DOI, arXiv:NNNN.NNNNN, or PMID:NNNN+).
     Project-local utility entries (Tier-3 in the audit doc) are accepted via
     PMID:0 sentinel — see docs/SCIENTIFIC_AUDIT_2026-05-15.md.

Exit codes:
  0 — all modules pass.
  non-zero — drift report on stderr, list of modules with their failure reason.

Usage:
  python scripts/ci/verify_module_references.py
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MODULES_DIR = REPO_ROOT / "workflow" / "modular" / "modules"

# Public identifier formats:
#   DOI:           10.XXXX/anything-no-spaces
#   arXiv:         arXiv:YYMM.NNNNN[vN]
#   PMID:          PMID: NNN... (PMID: 0 sentinel = Tier-3 project-local utility)
#   URLs to canonical specs (https://...) are accepted for software/spec citations
_IDENT_RE = re.compile(r"^(10\.\d+/[^\s]+|arXiv:\d{4}\.\d{4,5}(v\d+)?|PMID:\s*\d+|https?://\S+)$")


def _extract_references(py_path: Path) -> tuple[dict | None, str | None]:
    """Return (references_dict, error_message) — exactly one is None."""
    try:
        src = py_path.read_text(encoding="utf-8")
        tree = ast.parse(src)
    except Exception as exc:
        return None, f"parse error: {exc}"
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == "__references__":
                    try:
                        return ast.literal_eval(node.value), None
                    except Exception as exc:
                        return None, f"__references__ is not a literal dict: {exc}"
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == "__references__":
            if node.value is not None:
                try:
                    return ast.literal_eval(node.value), None
                except Exception as exc:
                    return None, f"__references__ annotated assign is not a literal: {exc}"
    return None, "no top-level __references__ assignment"


def main() -> int:
    if not MODULES_DIR.exists():
        print(f"ERROR: modules dir not found: {MODULES_DIR}", file=sys.stderr)
        return 2

    failures: list[tuple[str, str]] = []
    total = 0
    passed = 0
    project_local = 0

    for py in sorted(MODULES_DIR.glob("*.py")):
        if py.name.startswith("_") or py.name == "module_template.py":
            continue
        total += 1
        refs, err = _extract_references(py)
        if err is not None:
            failures.append((py.name, err))
            continue
        assert refs is not None  # mypy hint
        if not isinstance(refs, dict) or not refs:
            failures.append((py.name, "__references__ is empty"))
            continue

        bad_ids: list[str] = []
        is_project_local = False
        for key, entry in refs.items():
            if not isinstance(entry, dict) or "doi" not in entry:
                bad_ids.append(f"{key}: missing 'doi' field")
                continue
            doi_val = entry["doi"].strip()
            if doi_val == "PMID: 0":
                is_project_local = True
                continue
            if not _IDENT_RE.match(doi_val):
                bad_ids.append(f"{key}: identifier '{doi_val}' does not match accepted patterns (DOI / arXiv:N.N / PMID:N / https URL)")

        if bad_ids:
            failures.append((py.name, "; ".join(bad_ids)))
        else:
            passed += 1
            if is_project_local:
                project_local += 1

    if failures:
        print("=== verify_module_references — FAIL ===", file=sys.stderr)
        for name, msg in failures:
            print(f"  {name}: {msg}", file=sys.stderr)
        print(f"\n{len(failures)} of {total} modules failed verification.", file=sys.stderr)
        return 1

    print(f"verify_module_references: OK ({passed}/{total} modules verified; {project_local} project-local)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
