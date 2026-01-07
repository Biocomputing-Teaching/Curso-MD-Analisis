#!/usr/bin/env python3
"""Update TOC blocks in docs/episodes/*.md."""

from __future__ import annotations

from pathlib import Path
import re
import unicodedata

ROOT_DIR = Path(__file__).resolve().parents[1]
EPISODES_DIR = ROOT_DIR / "docs" / "episodes"

TOC_START = "<!-- toc:start -->"
TOC_END = "<!-- toc:end -->"
HEADING_RE = re.compile(r"^(##)\s+(.*)")
CODE_FENCE_RE = re.compile(r"^```")


def slugify(text: str) -> str:
    normalized = unicodedata.normalize("NFKD", text)
    ascii_text = normalized.encode("ascii", "ignore").decode("ascii")
    cleaned = re.sub(r"[^a-zA-Z0-9\s-]", "", ascii_text)
    collapsed = re.sub(r"[\s_-]+", "-", cleaned.strip().lower())
    return collapsed.strip("-")


def build_toc(lines: list[str]) -> list[str]:
    toc_lines = [TOC_START, "## Tabla de contenidos"]
    in_code = False
    for line in lines:
        if CODE_FENCE_RE.match(line.strip()):
            in_code = not in_code
            continue
        if in_code:
            continue
        match = HEADING_RE.match(line)
        if not match:
            continue
        title = match.group(2).strip()
        if title.lower() == "tabla de contenidos":
            continue
        anchor = slugify(title)
        if not anchor:
            continue
        toc_lines.append(f"- [{title}](#{anchor})")
    toc_lines.append(TOC_END)
    return toc_lines


def update_file(path: Path) -> None:
    content = path.read_text(encoding="utf-8")
    lines = content.splitlines()

    toc_block = build_toc(lines)
    if TOC_START in content and TOC_END in content:
        pre, rest = content.split(TOC_START, 1)
        _, post = rest.split(TOC_END, 1)
        new_content = pre.rstrip("\n") + "\n\n" + "\n".join(toc_block) + post
        path.write_text(new_content.rstrip("\n") + "\n", encoding="utf-8")
        return

    insert_idx = None
    for i, line in enumerate(lines):
        if line.startswith("## "):
            insert_idx = i
            break
    if insert_idx is None:
        return

    new_lines = lines[:insert_idx] + [""] + toc_block + [""] + lines[insert_idx:]
    path.write_text("\n".join(new_lines).rstrip("\n") + "\n", encoding="utf-8")


def main() -> int:
    for path in sorted(EPISODES_DIR.glob("*.md")):
        update_file(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
