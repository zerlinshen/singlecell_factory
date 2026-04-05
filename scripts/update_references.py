import ast
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")


def _fetch_crossref_meta(doi: str) -> dict:
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    req = urllib.request.Request(url, headers={"User-Agent": "singlecell_factory-reference-bot/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=8) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
        msg = payload.get("message", {})
        authors = msg.get("author", [])
        first_author = "Unknown"
        if authors:
            a0 = authors[0]
            first_author = a0.get("family") or a0.get("name") or "Unknown"
        journal = (msg.get("container-title") or ["Unknown"])[0]
        year = "Unknown"
        date_parts = msg.get("issued", {}).get("date-parts", [])
        if date_parts and date_parts[0]:
            year = str(date_parts[0][0])
        return {
            "authors": first_author,
            "journal": journal,
            "year": year,
            "title": msg.get("title", ["Unknown"])[0],
        }
    except (urllib.error.URLError, TimeoutError, ValueError, KeyError):
        return {"authors": "Unknown", "journal": "Unknown", "year": "Unknown", "title": "Unknown"}


def _extract_references_from_file(py_file: Path) -> list[dict]:
    refs: list[dict] = []
    content = py_file.read_text(encoding="utf-8")
    module_name = py_file.stem

    # 1) Structured references via __references__ dict.
    try:
        tree = ast.parse(content)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == "__references__" and isinstance(node.value, ast.Dict):
                        data = ast.literal_eval(node.value)
                        for _, ref in data.items():
                            if not isinstance(ref, dict):
                                continue
                            item = dict(ref)
                            item["module"] = module_name
                            refs.append(item)
    except Exception:
        pass

    # 2) DOI fallback from raw source text.
    found = set()
    for doi in DOI_RE.findall(content):
        d = doi.rstrip(").,;")
        if d in found:
            continue
        found.add(d)
        meta = _fetch_crossref_meta(d)
        refs.append(
            {
                "authors": meta["authors"],
                "journal": meta["journal"],
                "year": meta["year"],
                "title": meta["title"],
                "doi": d,
                "description": "auto-detected from module source",
                "module": module_name,
            }
        )
    return refs


def extract_references_from_modules(modules_dir: Path) -> list[dict]:
    refs = []
    for py_file in sorted(modules_dir.glob("*.py")):
        if py_file.name in {"__init__.py", "module_template.py"}:
            continue
        refs.extend(_extract_references_from_file(py_file))

    dedup = {}
    for r in refs:
        key = (r.get("doi", "").lower(), r.get("module", ""), r.get("description", ""))
        dedup[key] = r
    return list(dedup.values())


def _count_existing_rows(section_text: str) -> int:
    return sum(1 for line in section_text.splitlines() if line.strip().startswith("|") and line.count("|") >= 4) - 2


def update_readme_references(readme_path: Path, references: list[dict], force: bool = False) -> bool:
    if not readme_path.exists():
        raise FileNotFoundError(f"README not found: {readme_path}")

    content = readme_path.read_text(encoding="utf-8")
    references.sort(key=lambda x: (x.get("module", ""), x.get("year", ""), x.get("title", "")))

    lines = [
        "| # | Reference | DOI | Used by |",
        "|---|---|---|---|",
    ]
    for i, ref in enumerate(references, start=1):
        author = ref.get("authors", "Unknown")
        journal = ref.get("journal", "Unknown")
        year = ref.get("year", "Unknown")
        doi = ref.get("doi", "")
        doi_link = f"[{doi}](https://doi.org/{doi})" if doi else "N/A"
        module = f"`{ref.get('module', 'unknown')}`"
        desc = ref.get("description", "")
        if desc:
            module += f" ({desc})"
        lines.append(f"| {i} | {author} et al., *{journal}*, {year} | {doi_link} | {module} |")
    table = "\n".join(lines)

    pattern = r"(### Complete Citation List\n\n)([\s\S]*?)(\n---\n|\n## |\Z)"
    match = re.search(pattern, content)
    if not match:
        raise RuntimeError("Could not locate `### Complete Citation List` section in README.md")
    existing_rows = _count_existing_rows(match.group(2))
    if existing_rows > 0 and len(references) < int(existing_rows * 0.6) and not force:
        print(
            "Reference auto-update skipped: detected references are much fewer than existing table rows. "
            "Use --force to override."
        )
        return False
    new_content = re.sub(pattern, r"\1" + table + r"\3", content, count=1)
    readme_path.write_text(new_content, encoding="utf-8")
    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Force overwrite citation table even if detected refs are fewer.")
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    modules_dir = project_root / "workflow" / "modular" / "modules"
    readme_path = project_root / "README.md"
    refs = extract_references_from_modules(modules_dir)
    if refs:
        changed = update_readme_references(readme_path, refs, force=args.force)
        if changed:
            print(f"Updated references in {readme_path} with {len(refs)} module-level items.")
        else:
            print("README citation table kept unchanged.")
    else:
        print("No module-level references found; README was not modified.")
