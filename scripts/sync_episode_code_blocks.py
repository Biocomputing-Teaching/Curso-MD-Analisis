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
PY_BLOCK_PATTERN = re.compile(r"```python\n.*?\n```", re.DOTALL)
EMBED_PATTERN = re.compile(r'<div class="notebook-embed">.*?</div>(?:\s*</div>)+', re.DOTALL)
SYNC_BLOCK_PATTERN = re.compile(r'(?:```python\n.*?\n```|<div class="notebook-embed">.*?</div>(?:\s*</div>)+)', re.DOTALL)
SYNC_FROM_PATTERN = re.compile(r"<!--\s*sync-from:\s*(.*?)\s*-->")
LINK_LINE_PATTERN = re.compile(r"\n*Fuente del script: .*")
BASE_WEB = "{{ site.baseurl }}"
BASE_GH = "https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main"
NOTEBOOK_HTML_DIR = DOCS_DIR / "episodes" / "notebooks" / "rendered"
NBVIEWER_BASE = "https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis"
RESOURCES_PATH = DOCS_DIR / "resources.md"
RESOURCES_SCRIPTS_START = "<!-- resources:scripts-start -->"
RESOURCES_SCRIPTS_END = "<!-- resources:scripts-end -->"
RESOURCES_NOTEBOOKS_START = "<!-- resources:notebooks-start -->"
RESOURCES_NOTEBOOKS_END = "<!-- resources:notebooks-end -->"

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


def notebook_url(notebook_path: Path) -> str:
    html_path = NOTEBOOK_HTML_DIR / f"{notebook_path.stem}.html"
    rel = html_path.relative_to(DOCS_DIR).as_posix()
    return f"{BASE_WEB}/{rel}"


def nbviewer_url(notebook_path: Path) -> str:
    rel = notebook_path.relative_to(DOCS_DIR).as_posix()
    return f"{NBVIEWER_BASE}/{rel}"


def notebook_download_url(notebook_path: Path) -> str:
    rel = notebook_path.relative_to(DOCS_DIR).as_posix()
    return f"{BASE_WEB}/{rel}"


def script_download_url(script_path: Path) -> str:
    rel = script_path.relative_to(DOCS_DIR).as_posix()
    return f"{BASE_WEB}/{rel}"


def build_notebook_embed(notebook_path: Path, script_path: Path | None) -> str:
    view_url = notebook_url(notebook_path)
    download_nb = notebook_download_url(notebook_path)
    download_py = script_download_url(script_path) if script_path else None
    links = [f'<a href="{download_nb}" download>Descargar notebook</a>']
    if download_py:
        links.append(f'<a href="{download_py}" download>Descargar script (.py)</a>')
    links_html = " | ".join(links)
    return (
        f'<div class="notebook-embed">'
        f'<iframe src="{view_url}" loading="lazy"></iframe>'
        f'<div class="notebook-links">{links_html}</div>'
        f"</div>"
    )


def extract_title(markdown_text: str) -> str | None:
    if not markdown_text.startswith("---"):
        return None
    lines = markdown_text.splitlines()
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            break
        if lines[i].startswith("title:"):
            return lines[i].split(":", 1)[1].strip().strip('"').strip("'")
    return None


def page_title_for(markdown_path: Path) -> str:
    text = markdown_path.read_text(encoding="utf-8")
    title = extract_title(text)
    if title:
        return title
    return markdown_path.stem.replace("-", " ").title()


def page_url_for(markdown_path: Path) -> str:
    rel = markdown_path.relative_to(DOCS_DIR).as_posix()
    if rel == "index.md":
        return f"{BASE_WEB}/"
    if rel.endswith("/index.md"):
        return f"{BASE_WEB}/{rel[:-len('index.md')]}"
    return f"{BASE_WEB}/{rel[:-3]}/"


def replace_block_lines(text: str, start_marker: str, end_marker: str, new_lines: list[str]) -> str:
    lines = text.splitlines()
    start_idx = next((i for i, line in enumerate(lines) if start_marker in line), None)
    end_idx = next((i for i, line in enumerate(lines) if end_marker in line), None)
    if start_idx is None or end_idx is None or end_idx < start_idx:
        raise ValueError(f"Markers not found for block {start_marker} ... {end_marker}")
    indent = lines[start_idx].split(start_marker)[0]
    replaced = [f"{indent}{line}" if line else "" for line in new_lines]
    result = lines[:start_idx] + replaced + lines[end_idx + 1 :]
    return "\n".join(result) + ("\n" if text.endswith("\n") else "")


def build_resources_table(
    items: list[Path],
    source_pages: dict[Path, set[Path]],
    include_view: bool,
) -> list[str]:
    rows: list[str] = []
    for item in sorted(items, key=lambda path: path.name):
        pages = source_pages.get(item.resolve(), set()) or source_pages.get(item, set())
        if pages:
            page_links = [
                f'<a href="{page_url_for(page)}">{page_title_for(page)}</a>'
                for page in sorted(pages, key=lambda path: path.as_posix())
            ]
            pages_html = "<br>".join(page_links)
        else:
            pages_html = "No asociado"
        if include_view:
            view_link = f'<a href="{nbviewer_url(item)}">Ver (nbviewer)</a>'
            html_link = f'<a href="{notebook_url(item)}">Ver (HTML)</a>'
            download_link = f'<a href="{notebook_download_url(item)}" download>Descargar</a>'
            rows.append(
                f"<tr><td><code>{item.name}</code></td><td>{view_link}<br>{html_link}</td>"
                f"<td>{download_link}</td><td>{pages_html}</td></tr>"
            )
        else:
            download_link = f'<a href="{script_download_url(item)}" download>Descargar</a>'
            rows.append(
                f"<tr><td><code>{item.name}</code></td><td>{download_link}</td><td>{pages_html}</td></tr>"
            )
    return rows


def update_resources_index(
    markdown_paths: list[Path],
    scripts: list[Path],
    notebooks: list[Path],
    sources: dict[str, Path],
    notebooks_by_stem: dict[str, Path],
) -> bool:
    if not RESOURCES_PATH.exists():
        return False

    scripts = [path.resolve() for path in scripts]
    notebooks = [path.resolve() for path in notebooks]

    source_pages: dict[Path, set[Path]] = {}
    for markdown_path in markdown_paths:
        if markdown_path == RESOURCES_PATH:
            continue
        markdown_text = markdown_path.read_text(encoding="utf-8")
        markers = list(SYNC_FROM_PATTERN.finditer(markdown_text))
        for marker in markers:
            source_path = (ROOT_DIR / marker.group(1)).resolve()
            if source_path.exists():
                source_pages.setdefault(source_path.resolve(), set()).add(markdown_path)
        if not markers and SYNC_BLOCK_PATTERN.search(markdown_text):
            default_source = sources.get(markdown_path.stem)
            if default_source is not None:
                source_pages.setdefault(default_source.resolve(), set()).add(markdown_path)
        for match in re.findall(r'href="([^"]+\\.(ipynb|py))"|\\(([^)]+\\.(ipynb|py))\\)', markdown_text):
            href = next((item for item in match if item and item.endswith((".ipynb", ".py"))), "")
            if not href or href.startswith("http"):
                continue
            cleaned = href.replace("{{ site.baseurl }}/", "").lstrip("/")
            linked_path = DOCS_DIR / cleaned
            if linked_path.exists():
                source_pages.setdefault(linked_path.resolve(), set()).add(markdown_path)
        for script_path in scripts:
            if script_path.name in markdown_text:
                source_pages.setdefault(script_path.resolve(), set()).add(markdown_path)
                notebook_match = notebooks_by_stem.get(script_path.stem)
                if notebook_match is not None:
                    source_pages.setdefault(notebook_match.resolve(), set()).add(markdown_path)

    notebooks_index = DOCS_DIR / "episodes" / "notebooks" / "index.md"
    if notebooks_index.exists():
        for notebook_path in notebooks:
            source_pages.setdefault(notebook_path.resolve(), set()).add(notebooks_index)

    resources_text = RESOURCES_PATH.read_text(encoding="utf-8")
    script_rows = build_resources_table(scripts, source_pages, include_view=False)
    notebook_rows = build_resources_table(notebooks, source_pages, include_view=True)
    scripts_block = [
        RESOURCES_SCRIPTS_START,
        "<table class=\"resources-table\">",
        "<thead><tr><th>Script</th><th>Descargar</th><th>Página</th></tr></thead>",
        "<tbody>",
        *script_rows,
        "</tbody>",
        "</table>",
        RESOURCES_SCRIPTS_END,
    ]
    notebooks_block = [
        RESOURCES_NOTEBOOKS_START,
        "<table class=\"resources-table\">",
        "<thead><tr><th>Notebook</th><th>Ver</th><th>Descargar</th><th>Página</th></tr></thead>",
        "<tbody>",
        *notebook_rows,
        "</tbody>",
        "</table>",
        RESOURCES_NOTEBOOKS_END,
    ]
    updated = replace_block_lines(resources_text, RESOURCES_SCRIPTS_START, RESOURCES_SCRIPTS_END, scripts_block)
    updated = replace_block_lines(updated, RESOURCES_NOTEBOOKS_START, RESOURCES_NOTEBOOKS_END, notebooks_block)
    changed = updated != resources_text
    RESOURCES_PATH.write_text(updated, encoding="utf-8")
    return changed


def build_source_link(source_path: Path) -> str:
    if DOCS_DIR in source_path.parents:
        rel = source_path.relative_to(DOCS_DIR).as_posix()
        href = f"{BASE_WEB}/{rel}"
    else:
        rel = source_path.relative_to(ROOT_DIR).as_posix()
        href = f"{BASE_GH}/{rel}"
    return f'\n\nFuente del script: <a href="{href}">{source_path.name}</a>'


def sync_markdown(markdown_path: Path, sources: dict[str, Path], notebooks_by_stem: dict[str, Path]) -> bool:
    markdown_text = markdown_path.read_text(encoding="utf-8")
    if not SYNC_BLOCK_PATTERN.search(markdown_text):
        return False

    updated_text = markdown_text
    updated_text = LINK_LINE_PATTERN.sub("", updated_text)

    default_source = find_source_for_markdown(markdown_path, markdown_text, sources)
    used_default = False

    def replace_block(match: re.Match[str]) -> str:
        nonlocal used_default
        block = match.group(0)
        prefix = markdown_text[: match.start()]
        source_path = None
        for marker in SYNC_FROM_PATTERN.finditer(prefix):
            source_path = (ROOT_DIR / marker.group(1)).resolve()
        if source_path is None and not used_default and default_source is not None:
            source_path = default_source
            used_default = True

        if source_path is None:
            return block
        if not source_path.exists():
            print(f"Skipping {markdown_path}: missing {source_path}")
            return block

        if source_path.suffix == ".py":
            notebook_match = notebooks_by_stem.get(source_path.stem)
            if notebook_match is not None:
                source_path = notebook_match

        if source_path.suffix == ".ipynb":
            script_candidate = DOCS_DIR / "episodes" / "scripts" / f"{source_path.stem}.py"
            script_path = script_candidate if script_candidate.exists() else None
            return build_notebook_embed(source_path, script_path)

        source_text = load_source_text(source_path)
        new_block = f"```python\n{source_text}\n```"
        return new_block + build_source_link(source_path)

    updated_text = SYNC_BLOCK_PATTERN.sub(replace_block, updated_text)
    updated_text = updated_text.replace("https://biocomputing-teaching.github.io/Curso-MD-Analisis", BASE_WEB)

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
    sources: dict[str, Path] = {path.stem: path for path in scripts}
    sources.update({path.stem: path for path in notebooks})
    notebooks_by_stem = {path.stem: path for path in notebooks}

    updated = 0
    updated_files: list[Path] = []
    markdown_paths = sorted(DOCS_DIR.rglob("*.md"))
    for markdown_path in markdown_paths:
        if sync_markdown(markdown_path, sources, notebooks_by_stem):
            updated += 1
            updated_files.append(markdown_path)

    if update_resources_index(markdown_paths, scripts, notebooks, sources, notebooks_by_stem):
        updated += 1
        updated_files.append(RESOURCES_PATH)

    print(f"Updated {updated} file(s)")
    if updated_files:
        for path in updated_files:
            print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
