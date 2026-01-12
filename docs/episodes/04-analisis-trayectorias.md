---
layout: default
title: Episodio 5 - Análisis de trayectorias con MDTraj
permalink: /episodes/04-analisis-trayectorias/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Siguiente</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
<!-- toc:end -->

## Duración

- **Sesión:** 60 min
- **Ejercicios:** 45 min

## Objetivos

- Analizar trayectorias simples y complejas.
- Reimaginar y realinear trayectorias.
- Generar RMSD con MDTraj.

## Contenido

- Lectura de DCD con MDTraj.
- Reimaging y superposición al frame inicial.
- Gráficas para comparar sistemas.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el análisis simple con <a href="{{ site.baseurl }}/data/">traj.dcd</a>.

### Puntos clave

- Las métricas son comparables entre sistemas.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias_simple.ipynb">04-analisis-trayectorias_simple.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias_simple.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/04-analisis-trayectorias_simple.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el análisis del complejo con <a href="{{ site.baseurl }}/reference/">output_traj.dcd</a>.
- Guardar las figuras resultantes.

### Puntos clave

- MDTraj permite análisis reproducible con pocas líneas.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb">04-analisis-trayectorias.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/04-analisis-trayectorias.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Siguiente</a>
</div>
