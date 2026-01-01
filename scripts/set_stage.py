#!/usr/bin/env python3
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LAYOUT_PATH = ROOT / "docs" / "_layouts" / "default.html"

NAV_START = "<!-- stage:nav-start -->"
NAV_END = "<!-- stage:nav-end -->"

NAV_ITEMS = {
    "pre": [
        ("Inicio", "/"),
        ("Preparación", "/setup/"),
        ("Referencia", "/reference/"),
    ],
    "course": [
        ("Inicio", "/"),
        ("Preparación", "/setup/"),
        ("Episodios", "/episodes/"),
        ("Referencia", "/reference/"),
        ("Datos", "/data/"),
    ],
}


def replace_block(text: str, start_marker: str, end_marker: str, new_block: str) -> str:
    start_index = text.find(start_marker)
    end_index = text.find(end_marker)
    if start_index == -1 or end_index == -1 or end_index < start_index:
        raise ValueError(f"Markers not found for block {start_marker} ... {end_marker}")
    end_index += len(end_marker)
    return text[:start_index] + new_block + text[end_index:]


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


def update_nav(stage: str) -> None:
    layout = LAYOUT_PATH.read_text()
    items = NAV_ITEMS[stage]
    lines = [NAV_START]
    for label, href in items:
        lines.append(f'<a href="{{{{ site.baseurl }}}}{href}">{label}</a>')
    lines.append(NAV_END)
    updated = replace_block_lines(layout, NAV_START, NAV_END, lines)
    LAYOUT_PATH.write_text(updated)




def main() -> None:
    parser = argparse.ArgumentParser(description="Set site stage for navigation and setup content.")
    parser.add_argument("stage", choices=sorted(NAV_ITEMS.keys()))
    args = parser.parse_args()

    update_nav(args.stage)
    print(f"Stage updated to: {args.stage}")


if __name__ == "__main__":
    main()
