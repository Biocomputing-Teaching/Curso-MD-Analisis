#!/usr/bin/env python3
"""Audit figure references inside docs and point out missing or unused assets."""

from collections import defaultdict
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
DOC_ROOT = ROOT / "docs"
FIGURE_DIR = DOC_ROOT / "figures"

REFERENCE_PATTERN = re.compile(r"figures/([^\s)\"']+\.[A-Za-z0-9_-]+)")
TEXT_EXTENSIONS = {".md", ".markdown", ".html", ".htm", ".yml", ".yaml", ".css", ".js", ".txt"}
LATEX_DIR = ROOT / "LaTeX"
INCLUDEGRAPHICS_PATTERN = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")


def normalize_latex_reference(raw_path: str) -> str:
    """Normalize LaTeX \includegraphics paths into docs/figures relative paths."""
    path = raw_path.replace("\\DIRFigBook", "").replace("\\DIRFig", "")
    path = path.strip().replace("\\", "/")
    if not path:
        return ""
    segments = [segment for segment in path.split("/") if segment and segment != "."]
    while segments and segments[0] in {"..", "docs", "figures", "episodes"}:
        segments.pop(0)
    return "/".join(segments)


def collect_latex_references():
    """Gather figure references from LaTeX sources so they are treated as used."""
    references = set()
    ref_locations = defaultdict(list)

    if not LATEX_DIR.exists():
        return references, ref_locations

    for path in LATEX_DIR.rglob("*.tex"):
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue

        for match in INCLUDEGRAPHICS_PATTERN.finditer(text):
            normalized = normalize_latex_reference(match.group(1))
            if not normalized or normalized == ".":
                continue
            resolved = resolve_figure_reference(normalized)
            references.add(resolved)
            ref_locations[resolved].append(path.relative_to(ROOT))

    return references, ref_locations


def collect_references():
    """Scan document sources for figure references and remember the files that cite them."""
    references = set()
    ref_locations = defaultdict(list)

    for path in DOC_ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in TEXT_EXTENSIONS:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue

        for match in REFERENCE_PATTERN.finditer(text):
            relative_path = match.group(1)
            references.add(relative_path)
            ref_locations[relative_path].append(path.relative_to(ROOT))

    return references, ref_locations


def collect_actual_figures():
    """Return every file path that lives under docs/figures, relative to that directory."""
    actual = set()
    if not FIGURE_DIR.exists():
        return actual

    for path in FIGURE_DIR.rglob("*"):
        if path.is_file():
            actual.add(path.relative_to(FIGURE_DIR).as_posix())
    return actual


def resolve_figure_reference(path_str: str) -> str:
    """Resolve references without an extension by looking under docs/figures/."""
    candidate = FIGURE_DIR / path_str
    if candidate.exists():
        return Path(path_str).as_posix()
    parsed = Path(path_str)
    if parsed.suffix:
        return parsed.as_posix()
    if not FIGURE_DIR.exists():
        return path_str
    matches = sorted(FIGURE_DIR.rglob(f"{parsed.name}.*"))
    if not matches:
        return path_str
    return matches[0].relative_to(FIGURE_DIR).as_posix()


def main():
    references, reference_locations = collect_references()
    latex_references, latex_locations = collect_latex_references()
    references.update(latex_references)
    for name, locations in latex_locations.items():
        reference_locations[name].extend(locations)
    actual_figures = collect_actual_figures()

    missing = sorted(references - actual_figures)
    unused = sorted(actual_figures - references)

    if missing:
        print("Missing figure assets referenced in docs:")
        for ref in missing:
            locations = reference_locations.get(ref, [])
            print(f"  • {ref}")
            for location in locations:
                print(f"     → referenced from {location}")
            print("     → suggestion: add the file under docs/figures/ or update the reference")
        print()

    if unused:
        print("Figure files under docs/figures/ that are never referenced:")
        for path in unused:
            print(f"  • {path}")
        print("  → suggestion: either reference them (for new logos/assets) or remove them to avoid stale content\n")

    if not missing and not unused:
        print("All /figures/ references have corresponding files in docs/figures/.")

    if missing:
        sys.exit(1)


if __name__ == "__main__":
    main()
