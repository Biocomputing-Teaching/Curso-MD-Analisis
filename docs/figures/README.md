# Figure assets

`docs/figures/` is now the single source of truth for every image used by the Jekyll site and the LaTeX slides.

- The Beamer drivers set `\graphicspath{{../docs/figures/}}`, so a bare `\includegraphics{logo-cozyme-160x34.png}` loads `docs/figures/logo-cozyme-160x34.png` without duplicating files.
- When a slide references a figure, you can keep the asset inside `docs/figures/` (or `docs/figures/talks/` for talk-specific imagery) and the site will pick it up automatically.
- If you regenerate the talk markdown, `python LaTeX/sync_talks_figures.py` copies every figure used by the LaTeX beamer into `docs/figures/talks/` so the web talk pages have the same images.
