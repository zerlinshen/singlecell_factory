import ast
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")

# Disk cache for Crossref responses. Without this, network jitter (intermittent
# DOI timeouts) makes README citation rows non-deterministic across consecutive
# runs, which breaks the pre-commit reference-manager gate.
_CACHE_PATH = Path(__file__).resolve().parents[1] / "data" / "external" / "crossref_cache.json"
_UNKNOWN_META = {"authors": "Unknown", "journal": "Unknown", "year": "Unknown", "title": "Unknown"}
_cache_loaded: bool = False
_cache: dict[str, dict] = {}
_cache_dirty: bool = False


def _load_cache() -> None:
    global _cache, _cache_loaded
    if _cache_loaded:
        return
    if _CACHE_PATH.exists():
        try:
            _cache = json.loads(_CACHE_PATH.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            _cache = {}
    _cache_loaded = True


def _save_cache() -> None:
    if not _cache_dirty:
        return
    _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    _CACHE_PATH.write_text(
        json.dumps(_cache, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _fetch_crossref_meta(doi: str) -> dict:
    global _cache_dirty
    _load_cache()
    if doi in _cache:
        return _cache[doi]
    if os.environ.get("UPDATE_REFERENCES_OFFLINE") == "1":
        # Offline mode: do not hit the network, return Unknown without caching.
        # Failures are not cached so the next online run can fill them in.
        return _UNKNOWN_META
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
        meta = {
            "authors": first_author,
            "journal": journal,
            "year": year,
            "title": msg.get("title", ["Unknown"])[0],
        }
        _cache[doi] = meta
        _cache_dirty = True
        return meta
    except (urllib.error.URLError, TimeoutError, ValueError, KeyError):
        # Don't cache transient failures — let the next online run retry.
        return dict(_UNKNOWN_META)


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
        key = (
            _ref_text(r.get("doi")).lower(),
            _ref_text(r.get("module")),
            _ref_text(r.get("description")),
        )
        dedup[key] = r
    return list(dedup.values())


def _ref_text(value, default: str = "") -> str:
    if value is None or value == "":
        return default
    return str(value)


def _count_existing_rows(section_text: str) -> int:
    return sum(1 for line in section_text.splitlines() if line.strip().startswith("|") and line.count("|") >= 4) - 2


def update_readme_references(readme_path: Path, references: list[dict], force: bool = False) -> bool:
    if not readme_path.exists():
        raise FileNotFoundError(f"README not found: {readme_path}")

    content = readme_path.read_text(encoding="utf-8")
    references.sort(
        key=lambda x: (
            _ref_text(x.get("module")),
            _ref_text(x.get("year")),
            _ref_text(x.get("title")),
        )
    )

    lines = [
        "| # | Reference | DOI | Used by |",
        "|---|---|---|---|",
    ]
    for i, ref in enumerate(references, start=1):
        author = _ref_text(ref.get("authors"), "Unknown")
        journal = _ref_text(ref.get("journal"), "Unknown")
        year = _ref_text(ref.get("year"), "Unknown")
        doi = _ref_text(ref.get("doi"))
        doi_link = f"[{doi}](https://doi.org/{doi})" if doi else "N/A"
        module = f"`{_ref_text(ref.get('module'), 'unknown')}`"
        desc = _ref_text(ref.get("description"))
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
    parser.add_argument("--refresh-cache", action="store_true", help="Discard the Crossref cache before running so every DOI is re-fetched.")
    args = parser.parse_args()

    if args.refresh_cache and _CACHE_PATH.exists():
        _CACHE_PATH.unlink()
        print(f"Refreshed Crossref cache (removed {_CACHE_PATH}).")

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
    _save_cache()
