#!/usr/bin/env python3
"""
Standardize notebook structure:
- Ensure a single H1 in the intro cell.
- Insert a Table of Contents based on H2 headings.
- Add H2 step headings before each code cell.
"""

from __future__ import annotations

from pathlib import Path
import json
import re
import unicodedata
import uuid

ROOT_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT_DIR / "docs"
NOTEBOOKS_GLOB = "**/notebooks/*.ipynb"

HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)")


def slugify(text: str) -> str:
    normalized = unicodedata.normalize("NFKD", text)
    ascii_text = normalized.encode("ascii", "ignore").decode("ascii")
    cleaned = re.sub(r"[^a-zA-Z0-9\s-]", "", ascii_text)
    collapsed = re.sub(r"[\s_-]+", "-", cleaned.strip().lower())
    return collapsed.strip("-")


def normalize_headings(source: list[str], allow_h1: bool) -> list[str]:
    new_lines = []
    seen_h1 = False
    for raw in source:
        line = raw.rstrip("\n")
        match = HEADING_RE.match(line.strip())
        if match and match.group(1) == "#":
            if allow_h1 and not seen_h1:
                seen_h1 = True
                new_lines.append(line)
            else:
                new_lines.append("## " + match.group(2))
        else:
            new_lines.append(line)
    return [line + "\n" for line in new_lines]


def extract_intro_cell(cells: list[dict]) -> dict:
    for cell in cells:
        if cell.get("cell_type") == "markdown":
            cell = dict(cell)
            cell["source"] = normalize_headings(cell.get("source", []), allow_h1=True)
            return cell
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": ["# Notebook\n", "\n", "## Parte\n", "\n", "Descripcion pendiente.\n"],
    }


def find_intro_h2(source: list[str]) -> list[str]:
    headings = []
    for line in "".join(source).splitlines():
        match = HEADING_RE.match(line.strip())
        if match and match.group(1) == "##":
            headings.append(match.group(2).strip())
    return headings


def step_title_from_code(cell: dict, index: int) -> str:
    source = cell.get("source", [])
    title = None
    for raw in source:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#!"):
            continue
        if line.startswith("# -*-"):
            continue
        if line.startswith("#"):
            text = line.lstrip("#").strip()
            if text.startswith("%%"):
                text = text.lstrip("%").strip()
            if text:
                title = text
                break
            continue
        break
    if title:
        return f"Paso {index}: {title}"
    return f"Paso {index}"


def build_markdown_cell(text_lines: list[str], metadata: dict | None = None) -> dict:
    return {
        "cell_type": "markdown",
        "id": uuid.uuid4().hex,
        "metadata": metadata or {},
        "source": [line + "\n" for line in text_lines],
    }


def build_toc(headings: list[str]) -> dict:
    lines = ["## Tabla de contenidos", ""]
    for heading in headings:
        anchor = slugify(heading)
        if not anchor:
            continue
        lines.append(f"- [{heading}](#{anchor})")
    return build_markdown_cell(lines, {"codex": {"auto_toc": True}})


def is_auto_toc(cell: dict) -> bool:
    meta = cell.get("metadata", {})
    return bool(meta.get("codex", {}).get("auto_toc"))


def is_auto_step(cell: dict) -> bool:
    meta = cell.get("metadata", {})
    return bool(meta.get("codex", {}).get("auto_step"))


def has_h1(cell: dict) -> bool:
    if cell.get("cell_type") != "markdown":
        return False
    for line in "".join(cell.get("source", [])).splitlines():
        match = HEADING_RE.match(line.strip())
        if match and match.group(1) == "#":
            return True
    return False


def markdown_has_h2(cells: list[dict]) -> bool:
    for cell in cells:
        if cell.get("cell_type") != "markdown":
            continue
        for line in "".join(cell.get("source", [])).splitlines():
            match = HEADING_RE.match(line.strip())
            if match and match.group(1) == "##":
                if match.group(2).strip().lower() == "tabla de contenidos":
                    continue
                return True
    return False


def collect_h2_headings(cells: list[dict]) -> list[str]:
    headings = []
    for cell in cells:
        if cell.get("cell_type") != "markdown":
            continue
        if is_auto_toc(cell):
            continue
        for line in "".join(cell.get("source", [])).splitlines():
            match = HEADING_RE.match(line.strip())
            if match and match.group(1) == "##":
                title = match.group(2).strip()
                if title.lower() == "tabla de contenidos":
                    continue
                headings.append(title)
    return headings


def standardize_notebook(path: Path) -> None:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    cells = notebook.get("cells", [])
    for cell in cells:
        if not cell.get("id"):
            cell["id"] = uuid.uuid4().hex

    intro_cell = extract_intro_cell(cells)

    pending_markdown: list[dict] = []
    step_index = 0
    new_cells: list[dict] = [intro_cell]

    toc_placeholder = build_toc([])
    new_cells.append(toc_placeholder)

    found_intro = False
    for cell in cells:
        if cell.get("cell_type") == "markdown":
            if not found_intro:
                found_intro = True
                continue
        if not found_intro:
            continue
        if is_auto_toc(cell):
            continue
        if cell.get("cell_type") == "markdown":
            normalized = dict(cell)
            normalized["source"] = normalize_headings(normalized.get("source", []), allow_h1=False)
            pending_markdown.append(normalized)
            continue
        if cell.get("cell_type") == "code":
            step_index += 1
            if not markdown_has_h2(pending_markdown):
                title = step_title_from_code(cell, step_index)
                pending_markdown.insert(0, build_markdown_cell([f"## {title}"], {"codex": {"auto_step": True}}))
            new_cells.extend(pending_markdown)
            pending_markdown = []
            new_cells.append(cell)
            continue
        new_cells.append(cell)

    if pending_markdown:
        new_cells.extend(pending_markdown)

    toc_headings = find_intro_h2(intro_cell.get("source", [])) + collect_h2_headings(new_cells)
    toc_cell = build_toc(toc_headings)
    new_cells[1] = toc_cell

    notebook["cells"] = new_cells
    path.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    for path in sorted(DOCS_DIR.glob(NOTEBOOKS_GLOB)):
        standardize_notebook(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
