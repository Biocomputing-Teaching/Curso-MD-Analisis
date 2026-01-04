---
layout: default
title: Episodio 4 - Análisis de trayectorias con MDTraj
permalink: /episodes/04-analisis-trayectorias/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Siguiente</a>
</div>

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Objetivos

- Analizar trayectorias simples y complejas.
- Reimaginar y realinear trayectorias.
- Generar RMSD con MDTraj.

## Contenido

- Lectura de DCD con MDTraj.
- Reimaging y superposición al frame inicial.
- Gráficas para comparar sistemas.

## Demo guiada

### Parte simple: RMSD de alanina

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py" download>Descargar script (.py)</a></div></div>

### Parte compleja: RMSD del complejo

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py" download>Descargar script (.py)</a></div></div>

## Ejercicio

- Ejecutar el análisis simple con <a href="{{ site.baseurl }}/data/">traj.dcd</a>.
- Ejecutar el análisis del complejo con <a href="{{ site.baseurl }}/reference/">output_traj.dcd</a>.
- Guardar las figuras resultantes.

## Puntos clave

- Las métricas son comparables entre sistemas.
- MDTraj permite análisis reproducible con pocas líneas.

## Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias_simple.ipynb">04-analisis-trayectorias_simple.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb">04-analisis-trayectorias.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">04-analisis-trayectorias_simple.py</a>
- <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">04-analisis-trayectorias.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Siguiente</a>
</div>
