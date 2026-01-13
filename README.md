# MD Course - Repository maintenance

The public course content lives in
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/docs/index.md">docs/index.md</a>
and is published as a website. This README collects the technical maintenance
notes for the repository.

## Structure

- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a> website (Jekyll).
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/docs/index.md">docs/index.md</a> website index (the only content visible on the homepage).
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes">docs/episodes/</a> course episodes.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/scripts">docs/episodes/scripts/</a> episode source scripts.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/notebooks">docs/episodes/notebooks/</a> renderable notebooks.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/data">docs/data/</a> sample data.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/figures">docs/figures/</a> figures.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/scripts">scripts/</a> maintenance utilities.

## Maintenance workflows

### Makefile (targets)

Episode maintenance tasks can be run with `make` (or `make all`). All targets
operate on the files under
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes">docs/episodes/</a>
or their derived outputs.

- `make` / `make all`: run the full chain (`sync-code`, `sync-notebooks`,
  `update-episode-toc`, `render-notebooks`, `check-links`).
- `make sync-code`: sync the episode code blocks with the source scripts or
  notebooks.
- `make sync-notebooks`: sync scripts from notebooks and standardize notebook
  cell structure.
- `make standardize-notebooks`: only standardize notebook cell structure.
- `make update-episode-toc`: update the episode table of contents.
- `make render-notebooks`: render notebooks as HTML
  (`docs/episodes/notebooks/rendered`).
- `make check-links`: validate internal site links.
- `make check-sync`: verify notebooks match scripts.
- `make check`: run `sync-code`, `sync-notebooks`, `render-notebooks`,
  `check-links`, and `check-sync`.

### Menu state (pre/course)

```bash
python scripts/set_stage.py pre
python scripts/set_stage.py course
```

Updates the block between `<!-- stage:nav-start -->` and `<!-- stage:nav-end -->`
in <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/docs/_layouts/default.html">docs/_layouts/default.html</a>.

### Sync embedded code blocks in episodes

```bash
python scripts/sync_episode_code_blocks.py
```

Updates the ```python blocks in markdown files under
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a> and replaces them with
whichever script or notebook matches the same stem.
If you need to force the source, add to the markdown:

```html
<!-- sync-from: path/to/script.py -->
```

The script adds a link to the source file under the code block.

### Sync scripts from notebooks

```bash
python scripts/sync_notebooks_from_scripts.py
```

Uses the notebooks in
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/notebooks">docs/episodes/notebooks/</a>
as the source of truth and updates matching scripts in
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/scripts">docs/episodes/scripts/</a>.
Respects cells separated with `# %%`.

### Validate internal links

```bash
python scripts/check_links.py
```

Validates relative links inside
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a> and fails if any do not exist.

### Generate a sample DCD

```bash
python scripts/generate_example_dcd.py
```

Generates a sample DCD under
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/data">docs/data/</a> and requires
an environment with OpenMM.

## View the site locally (Jekyll)

This repository does not pin Jekyll versions, so the simplest option is
using a globally installed Jekyll. If you prefer pinned versions, create a
Gemfile and use `bundle exec` (out of scope for this README).

### Install (Linux, Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y ruby-full build-essential zlib1g-dev
gem install jekyll bundler webrick
```

### Install (macOS, Homebrew)

```bash
brew install ruby
echo 'export PATH="/opt/homebrew/opt/ruby/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
gem install jekyll bundler webrick
```

### Run locally

```bash
jekyll serve --source docs --livereload --baseurl /Curso-MD-Analisis
```

The site will be available at
<a href="http://127.0.0.1:4000/Curso-MD-Analisis/">http://127.0.0.1:4000/Curso-MD-Analisis/</a>.

### Common issues

- **`webrick` not found (Ruby 3+)**: run `gem install webrick`.
- **Port in use**: use `--port 4001` (or whatever port you want).
- **Gem permission errors**: install Ruby via `rbenv` or `brew` and avoid `sudo gem`.

## Best practices

- Edit public content under
  <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a>, not
  <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/README.md">README.md</a>.
- Avoid manually editing embedded code blocks in episodes: use the sync
  scripts to keep consistency.
- Keep links with `{{ site.baseurl }}` when pointing to resources under
  <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a>.
- Log maintenance changes in this README as the repository evolves with new
  scripts or workflows.
