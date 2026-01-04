---
layout: default
title: Episodio 3 - Simulaciones clásicas y control de calidad
permalink: /episodes/03-simulaciones-clasicas/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Siguiente</a>
</div>

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Objetivos

- Simular alanina y un complejo proteína-ligando.
- Configurar integrador de Langevin y reportes.
- Generar PDB minimizado y trayectorias DCD.

## Contenido

- Simulación simple con alanina.
- SystemGenerator con FF14SB y GAFF para el complejo.
- Reporters para energía y trayectoria.

## Demo guiada

### Parte simple: simulación de alanina

<!-- sync-from: docs/episodes/scripts/03-simulaciones-clasicas_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py" download>Descargar script (.py)</a></div></div>

### Parte compleja: simulación del complejo

<!-- sync-from: docs/episodes/scripts/03-simulaciones-clasicas.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py" download>Descargar script (.py)</a></div></div>

## Ejercicio

- Cambiar `--steps` y `--interval` en ambos sistemas.
- Verificar los archivos de salida de alanina y del complejo.
- Probar con otro `--output` base en el complejo.

## Puntos clave

- Los conceptos básicos se prueban con alanina.
- El complejo introduce más átomos y mayor coste.

## Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb">03-simulaciones-clasicas_simple.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb">03-simulaciones-clasicas.ipynb</a>
- <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py">03-simulaciones-clasicas_simple.py</a>
- <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py">03-simulaciones-clasicas.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Siguiente</a>
</div>
