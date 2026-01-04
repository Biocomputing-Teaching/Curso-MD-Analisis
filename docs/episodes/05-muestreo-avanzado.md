---
layout: default
title: Episodio 5 - Muestreo avanzado en OpenMM
permalink: /episodes/05-muestreo-avanzado/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Siguiente</a>
</div>

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Objetivos

- Aplicar una restricción simple en alanina.
- Simular un complejo con solvente explícito.
- Usar barostato para sistemas periódicos.

## Contenido

- Restricciones armónicas en sistemas pequeños.
- Solvatación con Modeller para sistemas grandes.
- Caja periódica y condiciones de contorno.

## Demo guiada

### Parte simple: restricción en alanina

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py" download>Descargar script (.py)</a></div></div>

### Parte compleja: solvatación y barostato

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py" download>Descargar script (.py)</a></div></div>

## Ejercicio

- Ejecutar la restricción y revisar la estabilidad.
- Ejecutar con `--solvate` y `--padding 12`.
- Probar otro `--water-model` y comparar el número de átomos.

## Puntos clave

- Las restricciones permiten explorar ventanas simples.
- La solvatación convierte el sistema en periódico.

## Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb">05-muestreo-avanzado_simple.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb">05-muestreo-avanzado.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>
- <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Siguiente</a>
</div>
