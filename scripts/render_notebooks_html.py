#!/usr/bin/env python3
"""
Render notebooks to static HTML files under docs/episodes/notebooks/rendered.
"""

from __future__ import annotations

import argparse
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT_DIR / "docs"
NOTEBOOKS_DIR = DOCS_DIR / "episodes" / "notebooks"
OUTPUT_DIR = NOTEBOOKS_DIR / "rendered"


def render_notebook(notebook_path: Path, output_path: Path) -> None:
    from nbconvert import HTMLExporter
    import nbformat

    notebook = nbformat.read(str(notebook_path), as_version=4)
    exporter = HTMLExporter()
    body, _ = exporter.from_notebook_node(notebook)
    output_path.write_text(body, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render notebooks to HTML.")
    parser.add_argument("--force", action="store_true", help="Re-render all notebooks.")
    args = parser.parse_args()

    if not NOTEBOOKS_DIR.exists():
        raise SystemExit(f"Missing notebooks directory: {NOTEBOOKS_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    notebooks = sorted(NOTEBOOKS_DIR.glob("*.ipynb"))
    if not notebooks:
        print("No notebooks found.")
        return 0

    updated = 0
    for notebook in notebooks:
        output_path = OUTPUT_DIR / f"{notebook.stem}.html"
        if not args.force and output_path.exists():
            if output_path.stat().st_mtime >= notebook.stat().st_mtime:
                continue
        render_notebook(notebook, output_path)
        updated += 1
        print(f"Rendered {notebook.name} -> {output_path}")

    print(f"Updated {updated} notebook(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
