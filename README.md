# Curso MD - Manteniment del repositori

El contingut públic del curs és a `docs/index.md` i es publica com a web.
Aquest README recull els detalls tècnics de manteniment del repositori.

## Estructura

- `docs/` lloc web (Jekyll).
- `docs/index.md` índex de la web (únic contingut visible a la portada).
- `docs/episodes/` episodis del curs.
- `docs/episodes/scripts/` scripts font dels episodis.
- `docs/episodes/notebooks/` notebooks renderitzables.
- `docs/data/` dades d'exemple.
- `docs/figures/` figures.
- `scripts/` utilitats de manteniment.

## Fluxos de manteniment

### Estat del menú (pre/curs)

```bash
python scripts/set_stage.py pre
python scripts/set_stage.py course
```

Actualitza el bloc entre `<!-- stage:nav-start -->` i `<!-- stage:nav-end -->`
a `docs/_layouts/default.html`.

### Sincronitzar codi incrustat als episodis

```bash
python scripts/sync_episode_code_blocks.py
```

Busca blocs ```python a `docs/**/*.md` i els reemplaça pel script o notebook
de referència que tingui el mateix nom (per `stem`).
Si cal forçar la font, afegeix al markdown:

```html
<!-- sync-from: ruta/al/script.py -->
```

L'script afegeix un enllaç al fitxer font sota el bloc de codi.

### Regenerar notebooks des de scripts

```bash
python scripts/sync_notebooks_from_scripts.py
```

Converteix cada `docs/**/scripts/*.py` en el notebook homònim a
`docs/**/notebooks/*.ipynb`. Respecta les cel·les separades amb `# %%`.

### Verificar enllaços interns

```bash
python scripts/check_links.py
```

Valida enllaços relatius a `docs/` i falla si algun no existeix.

### Generar un DCD d'exemple

```bash
python scripts/generate_example_dcd.py
```

Genera `docs/data/alanine-dipeptide.dcd`. Requereix un entorn amb OpenMM.

## Visualitzar la web localment (Jekyll)

Aquest repositori no inclou `Gemfile`, així que l'opció més simple és usar
Jekyll instal·lat com a gem global. Si prefereixes versions fixades, crea un
`Gemfile` i usa `bundle exec` (fora de l'abast d'aquest README).

### Instal·lació (Linux, Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y ruby-full build-essential zlib1g-dev
gem install jekyll bundler webrick
```

### Instal·lació (macOS, Homebrew)

```bash
brew install ruby
echo 'export PATH="/opt/homebrew/opt/ruby/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
gem install jekyll bundler webrick
```

### Execució local

```bash
jekyll serve --source docs --livereload --baseurl /Curso-MD-Analisis
```

La web queda disponible a [http://127.0.0.1:4000/Curso-MD-Analisis/](http://127.0.0.1:4000/Curso-MD-Analisis/).

### Problemes típics

- **`webrick` no trobat (Ruby 3+)**: executa `gem install webrick`.
- **Port ocupat**: usa `--port 4001` (o el port que vulguis).
- **Errors de permisos en gem**: instal·la Ruby via `rbenv` o `brew` i evita `sudo gem`.

## Bones pràctiques

- Editar el contingut públic a `docs/`, no a `README.md`.
- Evitar tocar a mà els blocs de codi incrustats als episodis: usa la
  sincronització per mantenir consistència.
- Mantenir enllaços amb `{{ site.baseurl }}` quan apuntin a recursos de `docs/`.
- Guardar els canvis de manteniment en aquest README quan el repositori
  evolucioni amb nous scripts o fluxos.
