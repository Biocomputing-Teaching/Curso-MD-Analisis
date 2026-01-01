#!/usr/bin/env python3
"""
Update notebooks under docs/**/notebooks/ from matching scripts in docs/**/scripts/.
"""

from __future__ import annotations

from pathlib import Path
import json
import re


ROOT_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT_DIR / "docs"
SCRIPTS_GLOB = "**/scripts/*.py"
NOTEBOOKS_GLOB = "**/notebooks/*.ipynb"
CELL_SPLIT_PATTERN = re.compile(r"^#\\s*%%.*$", re.MULTILINE)


def split_script_into_cells(script_text: str) -> list[str]:
    if CELL_SPLIT_PATTERN.search(script_text):
        parts = CELL_SPLIT_PATTERN.split(script_text)
        return [part.strip("\n") for part in parts if part.strip()]
    return [script_text.strip("\n")]


def build_code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [line if line.endswith("\n") else f"{line}\n" for line in source.splitlines()] or [""],
    }


def sync_notebook(notebook_path: Path, script_path: Path) -> bool:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    script_text = script_path.read_text(encoding="utf-8")
    code_chunks = split_script_into_cells(script_text)
    new_code_cells = [build_code_cell(chunk) for chunk in code_chunks]

    cells = notebook.get("cells", [])
    first_code_idx = next((i for i, cell in enumerate(cells) if cell.get("cell_type") == "code"), None)
    if first_code_idx is None:
        cells = cells + new_code_cells
    else:
        cells = [cell for cell in cells if cell.get("cell_type") != "code"]
        cells[first_code_idx:first_code_idx] = new_code_cells

    notebook["cells"] = cells
    notebook_path.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return True


def main() -> int:
    scripts = {path.stem: path for path in DOCS_DIR.glob(SCRIPTS_GLOB)}
    updated = 0
    updated_files: list[Path] = []
    missing_scripts: list[Path] = []
    failed_syncs: list[Path] = []
    missing_notebooks: list[Path] = []

    for notebook_path in sorted(DOCS_DIR.glob(NOTEBOOKS_GLOB)):
        script_path = scripts.get(notebook_path.stem)
        if not script_path:
            missing_scripts.append(notebook_path)
            continue
        try:
            if sync_notebook(notebook_path, script_path):
                updated += 1
                updated_files.append(notebook_path)
        except (OSError, json.JSONDecodeError, UnicodeError) as exc:
            print(f"Failed {notebook_path}: {exc}")
            failed_syncs.append(notebook_path)

    print(f"Updated {updated} notebook(s)")
    if updated_files:
        for path in updated_files:
            print(f"- {path}")
    if missing_scripts:
        print(f"Missing scripts for {len(missing_scripts)} notebook(s)")
        for path in missing_scripts:
            print(f"- {path}")
    if failed_syncs:
        print(f"Failed to sync {len(failed_syncs)} notebook(s)")
        for path in failed_syncs:
            print(f"- {path}")
    notebooks = {path.stem: path for path in DOCS_DIR.glob(NOTEBOOKS_GLOB)}
    for stem, script_path in scripts.items():
        if stem not in notebooks:
            missing_notebooks.append(script_path)
    if missing_notebooks:
        print(f"Missing notebooks for {len(missing_notebooks)} script(s)")
        for path in sorted(missing_notebooks):
            print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
