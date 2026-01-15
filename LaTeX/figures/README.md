# LaTeX figures

This directory now hosts all of the figures used by the LaTeX slides and PDFs.

- Web-only assets live under `docs/figures/`, so they do not clutter this tree.
- There is no longer a nested `LaTeX/figures/episodes/` folder; every slide simply pulls from `LaTeX/figures/`.
- The Beamer drivers already set `\graphicspath{{./figures/}{../docs/figures/}}` so each `\includegraphics` call resolves whether it points at a local slide figure or one published through the docs site.
