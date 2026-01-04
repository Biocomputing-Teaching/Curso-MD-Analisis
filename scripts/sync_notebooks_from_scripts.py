#!/usr/bin/env python3
"""
Update notebooks under docs/**/notebooks/ from matching scripts in docs/**/scripts/.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json
import re


ROOT_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT_DIR / "docs"
SCRIPTS_GLOB = "**/scripts/*.py"
NOTEBOOKS_GLOB = "**/notebooks/*.ipynb"
CELL_SPLIT_PATTERN = re.compile(r"^#\\s*%%.*$", re.MULTILINE)
IGNORE_SCRIPT_SUFFIXES = ()
IGNORE_NOTEBOOKS = {"06-pyemma.ipynb", "07-deeptime.ipynb"}
IGNORE_SCRIPT_NAMES = {"course_paths.py"}


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


def script_from_notebook(notebook_path: Path) -> str:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    code_cells = [
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "code"
    ]
    if not code_cells:
        return ""
    return "\n\n# %%\n\n".join(cell.rstrip() for cell in code_cells).rstrip() + "\n"


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


def sync_script_from_notebook(script_path: Path, notebook_path: Path) -> bool:
    script_text = script_from_notebook(notebook_path)
    if not script_text:
        return False
    if script_path.exists() and script_path.read_text(encoding="utf-8") == script_text:
        return False
    script_path.write_text(script_text, encoding="utf-8")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description="Sync scripts and notebooks using timestamps.")
    parser.add_argument("--check", action="store_true", help="Report pending sync without writing files.")
    args = parser.parse_args()

    scripts = {
        path.stem: path
        for path in DOCS_DIR.glob(SCRIPTS_GLOB)
        if not any(path.name.endswith(suffix) for suffix in IGNORE_SCRIPT_SUFFIXES)
        and path.name not in IGNORE_SCRIPT_NAMES
    }
    notebooks = {
        path.stem: path
        for path in DOCS_DIR.glob(NOTEBOOKS_GLOB)
        if path.name not in IGNORE_NOTEBOOKS
    }
    updated = 0
    updated_files: list[Path] = []
    pending_updates: list[Path] = []
    missing_scripts: list[Path] = []
    failed_syncs: list[Path] = []
    missing_notebooks: list[Path] = []

    # Sync scripts from notebooks (notebooks are source of truth).
    stems = sorted(set(scripts.keys()) | set(notebooks.keys()))
    for stem in stems:
        script_path = scripts.get(stem)
        notebook_path = notebooks.get(stem)
        if script_path is None and notebook_path is None:
            continue
        if script_path is None and notebook_path is not None:
            missing_scripts.append(notebook_path)
            continue
        if notebook_path is None and script_path is not None:
            missing_notebooks.append(script_path)
            continue

        try:
            if args.check:
                notebook_script_text = script_from_notebook(notebook_path)
                script_text = script_path.read_text(encoding="utf-8") if script_path.exists() else ""
                if notebook_script_text and notebook_script_text != script_text:
                    pending_updates.append(script_path)
            elif sync_script_from_notebook(script_path, notebook_path):
                updated += 1
                updated_files.append(script_path)
        except (OSError, json.JSONDecodeError, UnicodeError) as exc:
            print(f"Failed {notebook_path}: {exc}")
            failed_syncs.append(notebook_path)

    if args.check:
        if pending_updates:
            print("Sync needed for:")
            for path in pending_updates:
                print(f"- {path}")
            return 1
        print("No sync needed.")
        return 0

    print(f"Updated {updated} script(s)")
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
    if missing_notebooks:
        print(f"Missing notebooks for {len(missing_notebooks)} script(s)")
        for path in sorted(missing_notebooks):
            print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
