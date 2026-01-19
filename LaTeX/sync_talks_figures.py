#!/usr/bin/env python3
"""Sync LaTeX figures referenced in the beamer to docs/figures/talks."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LATEX_DIR = ROOT / "LaTeX"
FIG_SRC = ROOT / "docs" / "figures"
DOCS_FIG_DIR = ROOT / "docs" / "figures" / "talks"

EXTENSIONS = [".png", ".jpg", ".jpeg", ".pdf", ".gif", ".tif", ".tiff", ".eps", ".svg"]
INCLUDE_RE = re.compile(r"\\includegraphics(?:\\[(.*?)\\])?\\{([^}]+)\\}")


def normalize_path(raw_path: str) -> Path:
    path = raw_path.replace("\\DIRFigBook", "").replace("\\DIRFig", "")
    path = path.lstrip("/")
    if path.startswith("figures/"):
        path = path[len("figures/") :]
    if path.startswith("./"):
        path = path[2:]
    return Path(path)


def resolve_path(path: Path) -> Path | None:
    if path.suffix:
        candidate = FIG_SRC / path
        return candidate if candidate.exists() else None
    for ext in EXTENSIONS:
        candidate = FIG_SRC / f"{path}{ext}"
        if candidate.exists():
            return candidate
    search_dir = FIG_SRC / path.parent if path.parent != Path(".") else FIG_SRC
    matches = list(search_dir.glob(f"{path.name}.*"))
    return matches[0] if matches else None


def main() -> int:
    parser = argparse.ArgumentParser(description="Copy LaTeX figures used by the beamer to docs.")
    parser.add_argument("--tex", default="curso-md-analisis.tex", help="LaTeX source file")
    args = parser.parse_args()

    tex_path = LATEX_DIR / args.tex
    if not tex_path.exists():
        print(f"[sync_talks_figures] ERROR: No existe el archivo LaTeX: {tex_path}")
        return 1

    print(f"[sync_talks_figures] Leyendo LaTeX: {tex_path}")
    latex_text = tex_path.read_text(encoding="utf-8", errors="ignore")
    raw_paths = INCLUDE_RE.findall(latex_text)
    if not raw_paths:
        print("[sync_talks_figures] No se encontraron figuras.")
        return 0

    DOCS_FIG_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[sync_talks_figures] Copiando figuras a: {DOCS_FIG_DIR}")
    copied = 0
    missing = []
    seen = set()

    for raw in raw_paths:
        _, path = raw
        normalized = normalize_path(path)
        source = resolve_path(normalized)
        if not source:
            missing.append(raw)
            continue
        if source in seen:
            continue
        seen.add(source)
        try:
            rel = source.relative_to(FIG_SRC)
        except ValueError:
            rel = source.name
        dest = DOCS_FIG_DIR / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        copied += 1

    if missing:
        print(f"[sync_talks_figures] Aviso: faltan {len(missing)} figuras.")
    if copied:
        print(f"[sync_talks_figures] OK: sincronizadas {copied} figuras.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
