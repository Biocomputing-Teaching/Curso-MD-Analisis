#!/usr/bin/env python3
"""Generate markdown talks from the LaTeX beamer source using pandoc."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LATEX_DIR = ROOT / "LaTeX"
DOCS_TALKS_DIR = ROOT / "docs" / "talks"

SECTIONS = [
    ("Episode 1: Getting Started (environment and data)", "01-getting-started"),
    ("Episode 2: Running Simulations", "02-running-simulations"),
    ("Episode 3: Model Building & Editing", "03-model-building"),
    ("Episode 4: Advanced Simulation Examples", "04-advanced-examples"),
    ("Episode 5: Trajectory analysis", "05-trajectory-analysis"),
    ("Episode 6: MSMs with PyEMMA", "06-pyemma"),
    ("Episode 7: MSMs with deeptime", "07-deeptime"),
]

SECTION_RE = re.compile(r"\\section\[(.*?)\]\{([^}]*)\}")


def extract_sections(latex_text: str) -> list[tuple[str, str]]:
    sections = []
    matches = list(SECTION_RE.finditer(latex_text))
    print(f"[render_talks_md] Detected sections: {len(matches)}")
    for idx, match in enumerate(matches):
        title = match.group(2)
        start = match.end()
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(latex_text)
        body = latex_text[start:end].strip()
        sections.append((title, body))
    return sections


def sanitize_fragment(latex_fragment: str) -> str:
    cleaned_lines = []
    for line in latex_fragment.splitlines():
        stripped = line.strip()
        if stripped in ("}", "}%"):
            continue
        cleaned_lines.append(line)
    return "\n".join(cleaned_lines)


def fallback_markdown(latex_fragment: str) -> str:
    text = latex_fragment
    text = text.replace("\\DIRFigBook", "../docs/figures/talks").replace("\\DIRFig", "../docs/figures/talks")
    text = re.sub(r"\\frametitle\\{([^}]*)\\}", r"## \\1", text)
    text = re.sub(r"\\framesubtitle\\{([^}]*)\\}", r"### \\1", text)
    text = re.sub(r"\\subsection\\[(.*?)\\]\\{([^}]*)\\}", r"### \\2", text)
    text = re.sub(r"\\subsection\\{([^}]*)\\}", r"### \\1", text)
    text = re.sub(r"\\subsubsection\\[(.*?)\\]\\{([^}]*)\\}", r"#### \\2", text)
    text = re.sub(r"\\subsubsection\\{([^}]*)\\}", r"#### \\1", text)
    text = re.sub(r"\\begin\\{itemize\\}", "", text)
    text = re.sub(r"\\end\\{itemize\\}", "", text)
    text = re.sub(r"^\\s*\\\\item\\s*", "- ", text, flags=re.M)
    text = re.sub(r"\\begin\\{frame\\}", "", text)
    text = re.sub(r"\\end\\{frame\\}", "", text)
    text = re.sub(r"\\begin\\{[^}]+\\}", "", text)
    text = re.sub(r"\\end\\{[^}]+\\}", "", text)
    text = re.sub(r"\\includegraphics(?:\\[(.*?)\\])?\\{[^}]+\\}", "", text)
    text = re.sub(r"\\only<[^>]+>\\{([^}]*)\\}", r"\\1", text)
    text = re.sub(r"\\texttt\\{([^}]*)\\}", r"`\\1`", text)
    text = re.sub(r"\\emph\\{([^}]*)\\}", r"*\\1*", text)
    text = re.sub(r"\\cite\\{[^}]+\\}", "", text)
    text = re.sub(r"\\cita\\{[^}]+\\}", "", text)
    text = re.sub(r"\\parencite\\{[^}]+\\}", "", text)
    text = re.sub(r"\\textcite\\{[^}]+\\}", "", text)
    text = re.sub(r"\\begin\\{center\\}|\\end\\{center\\}", "", text)
    text = re.sub(r"\\begin\\{columns\\}.*?\\end\\{columns\\}", "", text, flags=re.S)
    text = re.sub(r"\\begin\\{figure\\}.*?\\end\\{figure\\}", "", text, flags=re.S)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip() + "\n"


def sanitize_title(title: str) -> str:
    cleaned = title.replace("\\&", "&")
    cleaned = re.sub(r"\\[\"'`^~]\\{?([A-Za-z])\\}?", r"\\1", cleaned)
    cleaned = re.sub(r"\\c\\{([A-Za-z])\\}", r"\\1", cleaned)
    cleaned = cleaned.replace("{", "").replace("}", "")
    cleaned = cleaned.replace("\\", "")
    return cleaned


def render_markdown(latex_fragment: str) -> str:
    normalized = latex_fragment.replace("\\DIRFigBook", "../docs/figures/talks").replace(
        "\\DIRFig", "../docs/figures/talks"
    )
    normalized = sanitize_fragment(normalized)
    return fallback_markdown(normalized)


def write_markdown(title: str, slug: str, body: str) -> None:
    safe_title = sanitize_title(title)
    front_matter = (
        "---\n"
        "layout: default\n"
        f"title: \"{safe_title}\"\n"
        "type: talk\n"
        f"permalink: /talks/{slug}/\n"
        "---\n\n"
    )
    header = f"# {safe_title}\n\n"
    wrapped_body = "{% raw %}\n" + body + "\n{% endraw %}\n"
    output_path = DOCS_TALKS_DIR / f"{slug}.md"
    output_path.write_text(front_matter + header + wrapped_body, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render talks markdown from LaTeX using pandoc.")
    parser.add_argument("--tex", default="curso-md-analisis.tex", help="LaTeX source file")
    args = parser.parse_args()

    tex_path = LATEX_DIR / args.tex
    if not tex_path.exists():
        print(f"[render_talks_md] ERROR: LaTeX file does not exist: {tex_path}")
        return 1

    print(f"[render_talks_md] Reading LaTeX: {tex_path}")
    latex_text = tex_path.read_text(encoding="utf-8")
    sections = extract_sections(latex_text)
    if not sections:
        print("[render_talks_md] ERROR: No sections found in the LaTeX.")
        return 1

    slug_by_title = {sanitize_title(title): slug for title, slug in SECTIONS}
    DOCS_TALKS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[render_talks_md] Output in: {DOCS_TALKS_DIR}")

    written = 0
    for title, body in sections:
        slug = slug_by_title.get(sanitize_title(title))
        if not slug:
            print(f"[render_talks_md] Warning: Section without slug: {title}")
            continue
        print(f"[render_talks_md] Processing: {title} -> {slug}.md")
        md_body = render_markdown(body)
        if not md_body:
            md_body = "Content pending generation.\n"
        write_markdown(title, slug, md_body)
        written += 1

    print(f"[render_talks_md] OK: wrote {written} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
