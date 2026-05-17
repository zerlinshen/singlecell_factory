# STATUS: one-off, promote-on-reuse
#!/usr/bin/env python3
"""Wave-5 visual verdict scaffolder.

Authors a verdict MD with SHA-pinned prompt + plot, ready for human or
second-agent fill-in. Does NOT call an LLM.
"""
from __future__ import annotations

import argparse
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_checklist_items(prompt_md: str) -> list[str]:
    """Extract biological-substance checklist items from the prompt MD.

    Heuristic: collects bullet/numbered lines under a "biological substance"
    or "checklist" heading. Falls back to 4 generic items if nothing found.
    """
    lines = prompt_md.splitlines()
    items: list[str] = []
    capturing = False
    heading_re = re.compile(r"^(#+)\s+(.*)$")
    bullet_re = re.compile(r"^\s*(?:[-*+]|\d+\.)\s+(.+?)\s*$")

    for line in lines:
        m = heading_re.match(line)
        if m:
            heading_text = m.group(2).lower()
            if "biolog" in heading_text or "checklist" in heading_text or "check" in heading_text:
                capturing = True
                continue
            if capturing and m.group(1):
                # New section heading ends capture.
                capturing = False
                continue
        if capturing:
            b = bullet_re.match(line)
            if b:
                items.append(b.group(1).strip())

    if not items:
        items = [
            "biological substance check 1",
            "biological substance check 2",
            "biological substance check 3",
            "biological substance check 4",
        ]
    # Always cap at 4 per the brief.
    return items[:4] if len(items) >= 4 else items + [
        f"biological substance check {i}" for i in range(len(items) + 1, 5)
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prompt-md", required=True, type=Path,
                   help="path to PLOT-<i>.md prompt")
    p.add_argument("--plot-png", required=True, type=Path,
                   help="path to PNG plot")
    p.add_argument("--plot-id", required=True, type=str,
                   help="plot id, e.g. PLOT-1")
    p.add_argument("--out-md", required=True, type=Path,
                   help="output verdict MD path")
    p.add_argument("--auto-verdict", choices=["pass", "partial", "fail", "pending"],
                   default="pending",
                   help="if not pending, prefill verdict")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    prompt_text = args.prompt_md.read_text()
    prompt_sha = sha256_file(args.prompt_md)
    plot_sha = sha256_file(args.plot_png)
    items = extract_checklist_items(prompt_text)

    authored_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

    frontmatter = (
        "---\n"
        f"plot_id: {args.plot_id}\n"
        f"prompt_md: {args.prompt_md}\n"
        f"prompt_sha: {prompt_sha}\n"
        f"plot_png: {args.plot_png}\n"
        f"plot_sha: {plot_sha}\n"
        f"verdict: {args.auto_verdict}\n"
        f"authored_at: {authored_at}\n"
        "---\n"
    )

    body_lines = ["", f"# Visual Verdict: {args.plot_id}", ""]
    body_lines.append("## Biological-substance checklist")
    body_lines.append("")
    auto = args.auto_verdict != "pending"
    for item in items:
        mark = "x" if auto else " "
        if auto:
            reason = f"{args.auto_verdict}: auto-verdict via CLI; reviewer should manually confirm"
        else:
            reason = "pending"
        body_lines.append(f"- [{mark}] {item}: {reason}")
    body_lines.append("")
    body_lines.append("## Reviewer notes")
    body_lines.append("")
    body_lines.append("_(fill in reasoning, citations, screenshots references)_")
    body_lines.append("")

    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text(frontmatter + "\n".join(body_lines))
    print(f"[wave5_visual_verdict] wrote {args.out_md}")
    print(f"  prompt_sha={prompt_sha}  plot_sha={plot_sha}  verdict={args.auto_verdict}")


if __name__ == "__main__":
    main()
