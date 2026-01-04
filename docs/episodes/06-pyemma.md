---
layout: default
title: Episodio 6 - MSM con PyEMMA
permalink: /episodes/06-pyemma/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Siguiente</a>
</div>

## Overview

- **Teaching:** 70 min
- **Exercises:** 50 min

## Objetivos

- Construir MSM con PyEMMA para un sistema simple y uno complejo.
- Seguir los tutoriales oficiales de PyEMMA en ambos sistemas.
- Comparar escalas de tiempo y validaciones básicas.

## Contenido

- Featurización, TICA y clustering.
- Estimación y validación del MSM.
- Análisis, PCCA/TPT y HMM.

## Demo guiada

### Parte simple: alanina dipeptido

<!-- sync-from: docs/episodes/scripts/06-pyemma-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py" download>Descargar script (.py)</a></div></div>

### Parte compleja: proteína-ligando

<!-- sync-from: docs/episodes/scripts/06-pyemma-complex.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py" download>Descargar script (.py)</a></div></div>

## Tutoriales PyEMMA cubiertos (ambos sistemas)

- 00 - Showcase pentapeptide: flujo completo de preparación.
- 01 - Data IO & featurization: lectura y features.
- 02 - Dimension reduction & discretization: TICA + clustering.
- 03 - MSM estimation & validation: ITS y tiempos de relajación.
- 04 - MSM analysis: escalas de tiempo y estados activos.
- 05 - PCCA & TPT: macroestados y flujos.
- 06 - Expectations & observables: observables simples.
- 07 - Hidden Markov models: HMM bayesiano.
- 08 - Common problems: fracción activa y conectividad.

## Ejercicio

- Ejecutar PyEMMA en alanina con <a href="{{ site.baseurl }}/data/">traj.dcd</a> y <a href="{{ site.baseurl }}/data/alanine-dipeptide.pdb">alanine-dipeptide.pdb</a>.
- Ejecutar PyEMMA en el complejo con <a href="{{ site.baseurl }}/reference/">output_traj.dcd</a> y <a href="{{ site.baseurl }}/reference/">output_minimised.pdb</a>.
- Comparar las escalas de tiempo e ITS.

## Puntos clave

- PyEMMA permite un pipeline completo y reproducible.
- El sistema simple valida el flujo antes del complejo.

## Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb">06-pyemma-alanine.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb">06-pyemma-complex.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py">06-pyemma-alanine.py</a>
- <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py">06-pyemma-complex.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Siguiente</a>
</div>
