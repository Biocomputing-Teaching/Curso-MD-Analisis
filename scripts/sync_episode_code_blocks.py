#!/usr/bin/env python3
"""
Sync embedded Python code blocks in markdown files with source scripts/notebooks.
"""

from __future__ import annotations

from pathlib import Path
import json
import re


ROOT_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT_DIR / "docs"
CODE_BLOCK_PATTERN = re.compile(r"```python\n.*?\n```", re.DOTALL)
SYNC_FROM_PATTERN = re.compile(r"<!--\s*sync-from:\s*(.*?)\s*-->")
LINK_LINE_PATTERN = re.compile(r"\n\nFuente del script: .*")
BASE_WEB = "https://biocomputing-teaching.github.io/Curso-MD-Analisis"
BASE_GH = "https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main"

SCRIPTS_DIR_GLOB = "**/scripts/*.py"
NOTEBOOKS_DIR_GLOB = "**/notebooks/*.ipynb"


def load_notebook_code(notebook_path: Path) -> str:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    code_cells = [
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "code"
    ]
    return "\n\n".join(cell.rstrip() for cell in code_cells).rstrip()


def find_source_for_markdown(markdown_path: Path, markdown_text: str, sources: dict[str, Path]) -> Path | None:
    match = SYNC_FROM_PATTERN.search(markdown_text)
    if match:
        return (ROOT_DIR / match.group(1)).resolve()

    return sources.get(markdown_path.stem)


def load_source_text(source_path: Path) -> str:
    if source_path.suffix == ".ipynb":
        return load_notebook_code(source_path)
    return source_path.read_text(encoding="utf-8").rstrip()


def build_source_link(source_path: Path) -> str:
    if DOCS_DIR in source_path.parents:
        rel = source_path.relative_to(DOCS_DIR).as_posix()
        href = f"{BASE_WEB}/{rel}"
    else:
        rel = source_path.relative_to(ROOT_DIR).as_posix()
        href = f"{BASE_GH}/{rel}"
    return f'\n\nFuente del script: <a href="{href}">{source_path.name}</a>'


def sync_markdown(markdown_path: Path, sources: dict[str, Path]) -> bool:
    markdown_text = markdown_path.read_text(encoding="utf-8")
    match = CODE_BLOCK_PATTERN.search(markdown_text)
    if not match:
        return False

    source_path = find_source_for_markdown(markdown_path, markdown_text, sources)
    if source_path is None:
        print(f"Skipping {markdown_path}: no matching source")
        return False
    if not source_path.exists():
        print(f"Skipping {markdown_path}: missing {source_path}")
        return False

    source_text = load_source_text(source_path)

    new_block = f"```python\n{source_text}\n```"
    updated_text = CODE_BLOCK_PATTERN.sub(new_block, markdown_text, count=1)
    updated_text = LINK_LINE_PATTERN.sub("", updated_text, count=1)
    updated_text = updated_text.replace(new_block, new_block + build_source_link(source_path), 1)

    if updated_text == markdown_text:
        return False

    markdown_path.write_text(updated_text, encoding="utf-8")
    return True


def main() -> int:
    if not DOCS_DIR.exists():
        print(f"Missing {DOCS_DIR}")
        return 1

    scripts = list(DOCS_DIR.glob(SCRIPTS_DIR_GLOB))
    notebooks = list(DOCS_DIR.glob(NOTEBOOKS_DIR_GLOB))
    sources: dict[str, Path] = {path.stem: path for path in notebooks}
    sources.update({path.stem: path for path in scripts})

    updated = 0
    updated_files: list[Path] = []
    for markdown_path in sorted(DOCS_DIR.rglob("*.md")):
        if sync_markdown(markdown_path, sources):
            updated += 1
            updated_files.append(markdown_path)

    print(f"Updated {updated} file(s)")
    if updated_files:
        for path in updated_files:
            print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
