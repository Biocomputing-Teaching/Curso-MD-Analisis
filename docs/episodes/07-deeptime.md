---
layout: default
title: Episodio 7 - MSM con deeptime
permalink: /episodes/07-deeptime/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
</div>

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Objetivos

- Construir MSM con deeptime para un sistema simple y uno complejo.
- Aplicar TICA, clustering y MSM en un flujo reproducible.
- Comparar resultados con PyEMMA.

## Contenido

- Features basados en distancias.
- Reducción de dimensionalidad con TICA.
- Estimación de MSM con MaximumLikelihoodMSM.
- Flujo inspirado en el ejemplo de <a href="https://github.com/joanmp-uoc/TFM">Joan Mora</a> (PCA + clustering + MSM).

## Demo guiada

### Parte simple: alanina dipeptido

<!-- sync-from: docs/episodes/scripts/07-deeptime-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-alanine.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-alanine.py" download>Descargar script (.py)</a></div></div>

### Parte compleja: proteína-ligando

<!-- sync-from: docs/episodes/scripts/07-deeptime-complex.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-complex.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-complex.py" download>Descargar script (.py)</a></div></div>

## Ejercicio

- Ejecutar deeptime en alanina con <a href="{{ site.baseurl }}/data/">traj.dcd</a> y <a href="{{ site.baseurl }}/data/alanine-dipeptide.pdb">alanine-dipeptide.pdb</a>.
- Ejecutar deeptime en el complejo con <a href="{{ site.baseurl }}/reference/">output_traj.dcd</a> y <a href="{{ site.baseurl }}/reference/">output_minimised.pdb</a>.
- Comparar las escalas de tiempo con PyEMMA.

## Puntos clave

- deeptime ofrece un pipeline MSM ligero y modular.
- La comparación con PyEMMA ayuda a validar resultados.

## Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-alanine.ipynb">07-deeptime-alanine.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-complex.ipynb">07-deeptime-complex.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-alanine.py">07-deeptime-alanine.py</a>
- <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-complex.py">07-deeptime-complex.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
</div>
