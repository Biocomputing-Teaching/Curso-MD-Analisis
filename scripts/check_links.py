import pathlib
import re
import sys

DOCS_ROOT = pathlib.Path(__file__).resolve().parents[1] / "docs"

LINK_RE = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
HTML_LINK_RE = re.compile(r"<a\\s+[^>]*href=[\"']([^\"']+)[\"'][^>]*>")
BASE_WEB = "https://biocomputing-teaching.github.io/Curso-MD-Analisis"


def normalize_link(link):
    link = link.strip()
    if link.startswith("http://") or link.startswith("https://"):
        return None
    if link.startswith("mailto:") or link.startswith("#"):
        return None
    link = link.replace("{{ site.baseurl }}", "").replace("{{site.baseurl}}", "")
    return link


def resolve_path(md_path, link):
    if link.startswith("/"):
        rel = link.lstrip("/")
        base = DOCS_ROOT / rel
    elif link.startswith("docs/"):
        base = DOCS_ROOT / link[len("docs/"):]
    else:
        base = md_path.parent / link

    if base.is_dir():
        index_md = base / "index.md"
        if index_md.exists():
            return index_md
        alt_md = base.with_suffix(".md")
        if alt_md.exists():
            return alt_md
    if base.exists():
        return base

    if str(base).endswith("/"):
        alt_md = pathlib.Path(str(base).rstrip("/") + ".md")
        if alt_md.exists():
            return alt_md

    if not base.suffix and (base.with_suffix(".md")).exists():
        return base.with_suffix(".md")

    return None


def main():
    missing = []
    externalized = []
    for md_path in DOCS_ROOT.rglob("*.md"):
        text = md_path.read_text(encoding="utf-8")
        links = LINK_RE.findall(text) + HTML_LINK_RE.findall(text)
        for match in links:
            if match.startswith(BASE_WEB):
                externalized.append((md_path, match))
                continue
            link = normalize_link(match)
            if not link:
                continue
            resolved = resolve_path(md_path, link)
            if resolved is None:
                missing.append((md_path, match))

    if externalized:
        print("Links should be relative or use {{ site.baseurl }}:")
        for md_path, link in externalized:
            print(f"- {md_path.relative_to(DOCS_ROOT)} -> {link}")
        sys.exit(1)

    if missing:
        print("Missing links:")
        for md_path, link in missing:
            print(f"- {md_path.relative_to(DOCS_ROOT)} -> {link}")
        sys.exit(1)

    print("All local links resolved.")


if __name__ == "__main__":
    main()
