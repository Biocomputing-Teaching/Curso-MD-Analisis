---
layout: default
title: Episodio 3 - Preparación y edición del sistema en OpenMM
permalink: /episodes/02-preparacion-sistema/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Siguiente</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Preparar la carpeta de trabajo](#preparar-la-carpeta-de-trabajo)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
- [Puntos clave](#puntos-clave)
- [Modeller (guía OpenMM)](#modeller-guia-openmm)
<!-- toc:end -->

## Duración

- **Sesión:** 60 min
- **Ejercicios:** 45 min

## Objetivos

- Preparar un sistema simple (alanina) y uno complejo (proteína-ligando).
- Replicar las operaciones de Modeller de la guía OpenMM.
- Generar estructuras listas para simulación y guardarlas de forma reproducible.

## Contenido

- Preparación simple con Modeller y solvatación ligera.
- PDBFixer para reparar proteínas complejas.
- Modeller: hidrógenos, solvente, membranas, partículas extra, eliminar agua.
- Ejemplo de guardado de resultados (Saving the Results).

## Preparar la carpeta de trabajo

Antes de ejecutar los scripts, crea la carpeta madre donde se guardarán todos los resultados (independientemente de dónde lances los cálculos) con <a href="{{ site.baseurl }}/episodes/scripts/course_paths.py">course_paths.py</a>:

```bash
python docs/episodes/scripts/course_paths.py
```

Esto genera la carpeta principal (por defecto `~/Concepcion26`), crea `data/` y `results/`, y te muestra el `export COURSE_DIR="..."` que debes ejecutar en la terminal para que todos los scripts usen ese mismo destino.

Después, copia los datos dentro de `data/`:

```bash
cp -R Course-MD-Data-main/* "$COURSE_DIR/data"
```

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/02-preparacion-sistema_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema_simple.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el script simple y revisar `"$COURSE_DIR/results/02-preparacion-sistema/simple/alanine_solvated.pdb"`.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema_simple.ipynb">02-preparacion-sistema_simple.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema_simple.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/02-preparacion-sistema_simple.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema_simple.py">02-preparacion-sistema_simple.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/02-preparacion-sistema.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el script complejo con `--output-dir outputs`.
- Probar un `--output-base` diferente.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb">02-preparacion-sistema.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/02-preparacion-sistema.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py">02-preparacion-sistema.py</a>

## Puntos clave

- La preparación depende del sistema objetivo.
- Guardar versiones preparadas facilita la trazabilidad.

## Modeller (guía OpenMM)

### Añadir hidrógenos

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

### Añadir solvente

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

### Añadir membrana

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

### Añadir partículas extra

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

### Eliminar agua

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

### Guardar resultados (Saving the Results)

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/02-preparacion-sistema.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py" download>Descargar script (.py)</a></div></div>

Script del curso: <a href="{{ site.baseurl }}/episodes/scripts/openmm_modeller_save.py">openmm_modeller_save.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Siguiente</a>
</div>
