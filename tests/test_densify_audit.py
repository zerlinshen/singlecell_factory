"""CI grep-ban: every .toarray() / .todense() in workflow/modular/modules/ must have a
# densify-allowed: <reason> comment on the same line or the immediately preceding line.

Fail fast so developers learn the policy before merging.
"""
import re
import textwrap
from pathlib import Path

MODULES_DIR = Path(__file__).parent.parent / "workflow" / "modular" / "modules"

# Matches the method call anywhere on a line (handles chained calls like .toarray().flatten())
_CALL_RE = re.compile(r"\.(toarray|todense)\s*\(")
_MARKER_RE = re.compile(r"#\s*densify-allowed\s*:", re.IGNORECASE)


def _check_file(path: Path) -> list[str]:
    violations: list[str] = []
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if not _CALL_RE.search(line):
            continue
        # Allow if this line carries the marker comment
        if _MARKER_RE.search(line):
            continue
        # Allow if the immediately preceding non-blank line carries the marker
        prev_idx = i - 1
        while prev_idx >= 0 and lines[prev_idx].strip() == "":
            prev_idx -= 1
        if prev_idx >= 0 and _MARKER_RE.search(lines[prev_idx]):
            continue
        violations.append(f"  {path.relative_to(MODULES_DIR.parent.parent)}:{i + 1}: {line.rstrip()}")
    return violations


def test_all_toarray_todense_have_densify_allowed_marker():
    """Every .toarray()/.todense() call in modules/ must carry a # densify-allowed comment."""
    py_files = sorted(MODULES_DIR.rglob("*.py"))
    assert py_files, f"No Python files found under {MODULES_DIR}"

    all_violations: list[str] = []
    for path in py_files:
        all_violations.extend(_check_file(path))

    if all_violations:
        report = "\n".join(all_violations)
        raise AssertionError(
            textwrap.dedent(f"""
            densify-allowed marker missing for {len(all_violations)} .toarray()/.todense() call(s).

            Add a comment like:
                # densify-allowed: <reason why this densification is bounded/safe>
            on the same line or the line immediately before the call.

            Violations:
            {report}
            """).strip()
        )


def test_modules_dir_exists():
    assert MODULES_DIR.is_dir(), f"modules directory not found: {MODULES_DIR}"
