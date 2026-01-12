# Curso MD - Manteniment del repositori

El contingut públic del curs és a
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/docs/index.md">docs/index.md</a>
i es publica com a web. Aquest README recull els detalls tècnics de manteniment
del repositori.

## Estructura

- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a> lloc web (Jekyll).
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/docs/index.md">docs/index.md</a> índex de la web (únic contingut visible a la portada).
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes">docs/episodes/</a> episodis del curs.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/scripts">docs/episodes/scripts/</a> scripts font dels episodis.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/notebooks">docs/episodes/notebooks/</a> notebooks renderitzables.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/data">docs/data/</a> dades d'exemple.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/figures">docs/figures/</a> figures.
- <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/scripts">scripts/</a> utilitats de manteniment.

## Fluxos de manteniment

### Makefile (targets)

Les tasques del manteniment dels episodis es poden llançar amb `make` (o
`make all`). Tots els objectius operen sobre els fitxers de
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes">docs/episodes/</a>
o les seves sortides derivades.

- `make` / `make all`: executa la cadena completa (`sync-code`, `sync-notebooks`,
  `update-episode-toc`, `render-notebooks`, `check-links`).
- `make sync-code`: sincronitza els blocs de codi dels episodis amb els scripts
  o notebooks font.
- `make sync-notebooks`: regenera els notebooks a partir dels scripts i
  estandarditza l'estructura de cel·les.
- `make standardize-notebooks`: només estandarditza l'estructura de cel·les
  dels notebooks.
- `make update-episode-toc`: actualitza la taula de continguts dels episodis.
- `make render-notebooks`: renderitza els notebooks com a HTML
  (carpeta `docs/episodes/notebooks/rendered`).
- `make check-links`: valida els enllaços interns del lloc web.
- `make check-sync`: comprova que els notebooks coincideixen amb els scripts.
- `make check`: executa `sync-code`, `sync-notebooks`, `render-notebooks`,
  `check-links` i `check-sync`.

### Estat del menú (pre/curs)

```bash
python scripts/set_stage.py pre
python scripts/set_stage.py course
```

Actualitza el bloc entre `<!-- stage:nav-start -->` i `<!-- stage:nav-end -->`
a <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/docs/_layouts/default.html">docs/_layouts/default.html</a>.

### Sincronitzar codi incrustat als episodis

```bash
python scripts/sync_episode_code_blocks.py
```

Actualitza els blocs ```python dels fitxers markdown dins de
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a> i els reemplaça pel
script o notebook de referència que tingui el mateix nom (per `stem`).
Si cal forçar la font, afegeix al markdown:

```html
<!-- sync-from: ruta/al/script.py -->
```

L'script afegeix un enllaç al fitxer font sota el bloc de codi.

### Regenerar notebooks des de scripts

```bash
python scripts/sync_notebooks_from_scripts.py
```

Converteix els scripts de
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/scripts">docs/episodes/scripts/</a>
en els notebooks homònims a
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/episodes/notebooks">docs/episodes/notebooks/</a>.
Respecta les cel·les separades amb `# %%`.

### Verificar enllaços interns

```bash
python scripts/check_links.py
```

Valida enllaços relatius dins de
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a> i falla si algun no existeix.

### Generar un DCD d'exemple

```bash
python scripts/generate_example_dcd.py
```

Genera un DCD d'exemple a
<a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs/data">docs/data/</a> i requereix
un entorn amb OpenMM.

## Visualitzar la web localment (Jekyll)

Aquest repositori no fixa versions de Jekyll, així que l'opció més simple és
usar Jekyll instal·lat com a gem global. Si prefereixes versions fixades,
crea un fitxer de dependències i usa `bundle exec` (fora de l'abast d'aquest README).

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

La web queda disponible a
<a href="http://127.0.0.1:4000/Curso-MD-Analisis/">http://127.0.0.1:4000/Curso-MD-Analisis/</a>.

### Problemes típics

- **`webrick` no trobat (Ruby 3+)**: executa `gem install webrick`.
- **Port ocupat**: usa `--port 4001` (o el port que vulguis).
- **Errors de permisos en gem**: instal·la Ruby via `rbenv` o `brew` i evita `sudo gem`.

## Bones pràctiques

- Editar el contingut públic a
  <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a>, no a
  <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/README.md">README.md</a>.
- Evitar tocar a mà els blocs de codi incrustats als episodis: usa la
  sincronització per mantenir consistència.
- Mantenir enllaços amb `{{ site.baseurl }}` quan apuntin a recursos de
  <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/tree/main/docs">docs/</a>.
- Guardar els canvis de manteniment en aquest README quan el repositori
  evolucioni amb nous scripts o fluxos.
