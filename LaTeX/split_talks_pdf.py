#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys
import tempfile


EPISODES = [
    (1, "getting-started"),
    (2, "running-simulations"),
    (3, "model-building"),
    (4, "advanced-examples"),
    (5, "trajectory-analysis"),
    (6, "pyemma"),
    (7, "deeptime"),
]


def read_pages(pdf_path):
    text = subprocess.check_output(["pdftotext", "-layout", pdf_path, "-"])
    pages = text.decode("utf-8", errors="ignore").split("\f")
    while pages and not pages[-1].strip():
        pages.pop()
    return pages


def find_episode_starts(pages):
    starts = {}
    for episode, _ in EPISODES:
        pattern = re.compile(rf"Episode\s+{episode}[^\n]*Summary")
        for i, page in enumerate(pages, start=1):
            if pattern.search(page):
                starts[episode] = i
                break
    return starts


def extract_range(pdf_path, start, end, output_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        pattern = os.path.join(tmpdir, "page-%03d.pdf")
        subprocess.check_call(
            ["pdfseparate", "-f", str(start), "-l", str(end), pdf_path, pattern]
        )
        pages = [pattern % i for i in range(start, end + 1)]
        subprocess.check_call(["pdfunite", *pages, output_path])


def main():
    parser = argparse.ArgumentParser(
        description="Split the main course PDF into per-episode PDFs."
    )
    parser.add_argument("--pdf", required=True, help="Input PDF path.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    args = parser.parse_args()

    pages = read_pages(args.pdf)
    total_pages = len(pages)
    starts = find_episode_starts(pages)

    missing = [ep for ep, _ in EPISODES if ep not in starts]
    if missing:
        missing_str = ", ".join(str(ep) for ep in missing)
        print(
            f"Missing Summary page for episodes: {missing_str}. "
            "Check the LaTeX source for missing sections.",
            file=sys.stderr,
        )
        return 1

    for idx, (episode, slug) in enumerate(EPISODES):
        start = starts[episode]
        if idx + 1 < len(EPISODES):
            end = starts[EPISODES[idx + 1][0]] - 1
        else:
            end = total_pages

        out_name = f"curso-md-analisis-ep{episode:02d}-{slug}.pdf"
        out_path = os.path.join(args.outdir, out_name)
        extract_range(args.pdf, start, end, out_path)
        print(f"Wrote {out_path} ({start}-{end}).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
