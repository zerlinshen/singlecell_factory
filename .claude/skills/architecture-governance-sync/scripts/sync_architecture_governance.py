#!/usr/bin/env python3
"""Refresh README architecture governance text and a structure diagram PNG."""

from __future__ import annotations

import argparse
import html
import os
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

START = "<!-- ARCHITECTURE-GOVERNANCE:START -->"
END = "<!-- ARCHITECTURE-GOVERNANCE:END -->"
IGNORE_DIRS = {
    ".git",
    ".hg",
    ".svn",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    "node_modules",
    ".venv",
    "venv",
    "env",
}
IMPORTANT_FILES = {
    "README.md",
    "AGENTS.md",
    "CODEX_PROFILE.md",
    "SKILL_PLAYBOOK.md",
    "SKILL_PROTOCOL.md",
    "00_HUMAN_START_HERE.md",
}


@dataclass(frozen=True)
class Entry:
    name: str
    kind: str
    role: str


def classify(name: str, is_dir: bool) -> str:
    lower = name.lower()
    if name in IMPORTANT_FILES:
        return "governance / entrypoint"
    if lower in {"human_review", "reports", "figures", "figure_compare", "paper_assets"}:
        return "human-facing evidence"
    if lower in {"agent_runs", "before-every-run", "handoff", ".omx"}:
        return "agent-facing memory"
    if lower in {"scripts", "workflow", "workflows", "src", "bridges", "r", "tests"}:
        return "execution / code"
    if lower in {"docs", "documentation"}:
        return "documentation"
    if lower in {"results", "outputs", "data"}:
        return "data / generated artifacts"
    if is_dir:
        return "workspace area (inferred)"
    return "file artifact (inferred)"


def list_entries(root: Path, max_entries: int = 24) -> list[Entry]:
    entries: list[Entry] = []
    for child in sorted(root.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower())):
        if child.name in IGNORE_DIRS:
            continue
        if child.name.startswith(".") and child.name != ".omx":
            continue
        if child.name == "docs" and (child / "architecture").exists():
            role = "documentation; contains generated architecture governance"
        else:
            role = classify(child.name, child.is_dir())
        entries.append(Entry(child.name, "dir" if child.is_dir() else "file", role))
        if len(entries) >= max_entries:
            break
    return entries


def rel(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def current_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def existing_timestamp(text: str) -> str | None:
    for line in text.splitlines():
        if line.startswith("Last refreshed: `") and line.endswith("`"):
            return line.removeprefix("Last refreshed: `").removesuffix("`")
    return None


def markdown_snapshot(root: Path, diagram_path: Path, entries: list[Entry], generated_at: str) -> str:
    diagram_rel = rel(diagram_path, root)
    lines = [
        START,
        "## Architecture Governance",
        "",
        f"Last refreshed: `{generated_at}`",
        "",
        f"![Workspace structure]({diagram_rel})",
        "",
        "### Operating Map",
        "",
        "| Surface | Type | Governance role |",
        "| --- | --- | --- |",
    ]
    for entry in entries:
        lines.append(f"| `{entry.name}` | `{entry.kind}` | {entry.role} |")
    lines.extend(
        [
            "",
            "### Required Update Habit",
            "",
            "- Refresh this block and `docs/architecture/structure.png` whenever folders, workflow boundaries, run lanes, artifact locations, or remote/local contracts change.",
            "- Keep human-facing deliverables easy to find before agent-facing logs.",
            "- Keep run parameters, commands, source paths, outcomes, and residual risks in agent-facing run memory.",
            "- If remote behavior changes, update the remote README or before-every-run memory in the same workstream.",
            "",
            END,
            "",
        ]
    )
    return "\n".join(lines)


def replace_block(readme_text: str, block: str) -> str:
    if START in readme_text and END in readme_text:
        before = readme_text.split(START, 1)[0].rstrip()
        after = readme_text.split(END, 1)[1].lstrip()
        return f"{before}\n\n{block}{after}"
    insertion = "\n" if readme_text.endswith("\n") else "\n\n"
    return readme_text + insertion + block


def diagram_svg(root: Path, entries: list[Entry]) -> str:
    width = 1400
    row_h = 72
    height = 240 + max(1, len(entries)) * row_h
    esc_root = html.escape(root.name or str(root))
    rows = []
    y = 170
    colors = {
        "human-facing evidence": "#E8F3D6",
        "agent-facing memory": "#FCE7C8",
        "execution / code": "#DDEBFF",
        "documentation": "#E9DDFF",
        "data / generated artifacts": "#FFE1E1",
        "governance / entrypoint": "#D9F6F0",
    }
    for idx, entry in enumerate(entries):
        fill = colors.get(entry.role.split(";")[0], "#F3F4F6")
        safe_name = html.escape(entry.name)
        safe_role = html.escape(entry.role)
        rows.append(
            f'<rect x="90" y="{y}" width="1220" height="52" rx="14" fill="{fill}" stroke="#25313D" stroke-width="2"/>'
            f'<text x="120" y="{y + 33}" font-family="Avenir, Helvetica, Arial" font-size="24" fill="#17212B">{safe_name}</text>'
            f'<text x="610" y="{y + 33}" font-family="Avenir, Helvetica, Arial" font-size="21" fill="#43515F">{safe_role}</text>'
        )
        if idx < len(entries) - 1:
            rows.append(f'<line x1="700" y1="{y + 52}" x2="700" y2="{y + row_h}" stroke="#A7B0BA" stroke-width="2"/>')
        y += row_h
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#FAF7F0"/>
<text x="90" y="74" font-family="Avenir, Helvetica, Arial" font-size="42" font-weight="700" fill="#17212B">Architecture Governance Map</text>
<text x="90" y="118" font-family="Avenir, Helvetica, Arial" font-size="24" fill="#52616B">{esc_root}</text>
<rect x="90" y="140" width="1220" height="8" rx="4" fill="#25313D"/>
{''.join(rows)}
</svg>
"""


def write_png_with_pillow(png_path: Path, root: Path, entries: list[Entry]) -> bool:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return False
    width = 1400
    row_h = 78
    height = 250 + max(1, len(entries)) * row_h
    image = Image.new("RGB", (width, height), "#FAF7F0")
    draw = ImageDraw.Draw(image)
    try:
        title_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 42)
        body_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 24)
        small_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 21)
    except Exception:
        title_font = body_font = small_font = ImageFont.load_default()
    draw.text((90, 45), "Architecture Governance Map", fill="#17212B", font=title_font)
    draw.text((90, 100), root.name or str(root), fill="#52616B", font=body_font)
    draw.rounded_rectangle((90, 140, 1310, 148), radius=4, fill="#25313D")
    colors = {
        "human-facing evidence": "#E8F3D6",
        "agent-facing memory": "#FCE7C8",
        "execution / code": "#DDEBFF",
        "documentation": "#E9DDFF",
        "data / generated artifacts": "#FFE1E1",
        "governance / entrypoint": "#D9F6F0",
    }
    y = 175
    for entry in entries:
        fill = colors.get(entry.role.split(";")[0], "#F3F4F6")
        draw.rounded_rectangle((90, y, 1310, y + 56), radius=14, fill=fill, outline="#25313D", width=2)
        draw.text((120, y + 16), entry.name[:38], fill="#17212B", font=body_font)
        draw.text((610, y + 18), entry.role[:62], fill="#43515F", font=small_font)
        y += row_h
    image.save(png_path)
    return True


def convert_svg_to_png(svg_path: Path, png_path: Path) -> bool:
    commands = [
        ["rsvg-convert", str(svg_path), "-o", str(png_path)],
        ["magick", str(svg_path), str(png_path)],
        ["convert", str(svg_path), str(png_path)],
        ["qlmanage", "-t", "-s", "1400", "-o", str(png_path.parent), str(svg_path)],
    ]
    for command in commands:
        if not shutil.which(command[0]):
            continue
        try:
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            ql_out = png_path.parent / f"{svg_path.name}.png"
            if ql_out.exists() and ql_out != png_path:
                ql_out.replace(png_path)
            return png_path.exists() and png_path.stat().st_size > 0
        except Exception:
            continue
    return False


def sync(root: Path, check: bool) -> int:
    root = root.resolve()
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"Workspace root not found or not a directory: {root}")
    arch_dir = root / "docs" / "architecture"
    readme = root / "README.md"
    arch_dir.mkdir(parents=True, exist_ok=True)
    entries = list_entries(root)
    svg_path = arch_dir / "structure.svg"
    png_path = arch_dir / "structure.png"
    snapshot_path = arch_dir / "ARCHITECTURE_GOVERNANCE.md"
    readme_text = readme.read_text(encoding="utf-8") if readme.exists() else f"# {root.name}\n"
    generated_at = existing_timestamp(readme_text) if check else None
    block = markdown_snapshot(root, png_path, entries, generated_at or current_timestamp())
    new_readme = replace_block(readme_text, block)
    svg = diagram_svg(root, entries)
    snapshot = "# Architecture Governance Snapshot\n\n" + block

    drift = []
    if readme_text != new_readme:
        drift.append(str(readme))
    if not svg_path.exists() or svg_path.read_text(encoding="utf-8") != svg:
        drift.append(str(svg_path))
    if not snapshot_path.exists() or snapshot_path.read_text(encoding="utf-8") != snapshot:
        drift.append(str(snapshot_path))
    if not png_path.exists() or png_path.stat().st_size == 0:
        drift.append(str(png_path))

    if check:
        if drift:
            print("Architecture governance drift detected:")
            for item in drift:
                print(f"- {item}")
            return 1
        print("Architecture governance is up to date.")
        return 0

    readme.write_text(new_readme, encoding="utf-8")
    svg_path.write_text(svg, encoding="utf-8")
    snapshot_path.write_text(snapshot, encoding="utf-8")
    if not write_png_with_pillow(png_path, root, entries):
        if not convert_svg_to_png(svg_path, png_path):
            raise SystemExit("Could not generate structure.png; install Pillow, rsvg-convert, ImageMagick, or qlmanage.")
    print(f"Updated README: {readme}")
    print(f"Updated snapshot: {snapshot_path}")
    print(f"Updated diagram SVG: {svg_path}")
    print(f"Updated diagram PNG: {png_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Refresh README architecture governance and structure diagram artifacts."
    )
    parser.add_argument("workspace_root", help="Workspace root containing or needing README.md")
    parser.add_argument("--check", action="store_true", help="Fail if generated artifacts are stale")
    args = parser.parse_args()
    return sync(Path(args.workspace_root), args.check)


if __name__ == "__main__":
    raise SystemExit(main())
