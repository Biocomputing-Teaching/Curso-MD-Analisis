---
layout: default
title: Episodio 2 - Preparación de sistemas en OpenMM
permalink: /episodes/02-preparacion-sistema/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Siguiente</a>
</div>

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Objetivos

- Preparar un sistema simple (alanina) y uno complejo (proteína-ligando).
- Generar estructuras listas para simulación.
- Guardar archivos de salida reproducibles.

## Contenido

- Preparación simple con Modeller y solvatación ligera.
- PDBFixer para reparar proteínas complejas.
- Minimización con OpenMM y SystemGenerator.

## Preparar la carpeta de trabajo

Antes de ejecutar los scripts, crea la carpeta madre donde se guardarán todos los resultados (independientemente de dónde lances los cálculos) con <a href="{{ site.baseurl }}/episodes/scripts/course_paths.py">course_paths.py</a>:

```bash
python docs/episodes/scripts/course_paths.py
```

Esto genera la carpeta principal (por defecto `~/Concepcion26`), crea `data/` y `results/`, y te muestra el `export COURSE_DIR="..."` que debes ejecutar en la terminal para que todos los scripts usen ese mismo destino.

Después, copia los datos dentro de `data/`:

```bash
cp -R Curso-MD-Data-main/* "$COURSE_DIR/data"
```

## Demo guiada

### Parte simple: preparación de alanina

<!-- sync-from: docs/episodes/scripts/02-preparacion-sistema_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema_simple.py" download>Descargar script (.py)</a></div></div>

### Parte compleja: preparación de proteína

<!-- sync-from: docs/episodes/scripts/02-preparacion-sistema.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

## Ejercicio

- Ejecutar el script simple y revisar `"$COURSE_DIR/results/02-preparacion-sistema/simple/alanine_solvated.pdb"`.
- Ejecutar el script complejo con `--output-dir outputs`.
- Probar un `--output-base` diferente.

## Puntos clave

- La preparación depende del sistema objetivo.
- Guardar versiones preparadas facilita la trazabilidad.

## Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema_simple.ipynb">02-preparacion-sistema_simple.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb">02-preparacion-sistema.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema_simple.py">02-preparacion-sistema_simple.py</a>
- <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py">02-preparacion-sistema.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Siguiente</a>
</div>
