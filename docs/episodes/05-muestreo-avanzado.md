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

- Aplicar una restricción simple en alanina.
- Simular un complejo con solvente explícito.
- Usar barostato para sistemas periódicos.

## Contenido

- Restricciones armónicas en sistemas pequeños.
- Solvatación con Modeller para sistemas grandes.
- Caja periódica y condiciones de contorno.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar la restricción y revisar la estabilidad.

### Puntos clave

- Las restricciones permiten explorar ventanas simples.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb">05-muestreo-avanzado_simple.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/05-muestreo-avanzado_simple.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar con `--solvate` y `--padding 12`.
- Probar otro `--water-model` y comparar el número de átomos.

### Puntos clave

- La solvatación convierte el sistema en periódico.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb">05-muestreo-avanzado.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/05-muestreo-avanzado.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Siguiente</a>
</div>
