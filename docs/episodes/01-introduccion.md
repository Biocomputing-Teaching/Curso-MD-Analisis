---
layout: default
title: Episodio 1 - Introducción (entorno y datos)
permalink: /episodes/01-introduccion/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Siguiente</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
- [Puntos clave](#puntos-clave)
<!-- toc:end -->

## Duración

- **Sesión:** 45 min
- **Ejercicios:** 30 min



## Objetivos

- Identificar los archivos de entrada para alanina y complejo proteína-ligando.
- Verificar que OpenMM y OpenFF Toolkit funcionen en el entorno.
- Cargar datos base para ambos sistemas siguiendo el flujo de primeros pasos.

## Contenido

- PDB y MOL/SDF como entradas.
- Lectura de alanina con OpenMM.
- Lectura de proteína y ligando con OpenFF Toolkit.
- Ruta de trabajo `COURSE_DIR` y validación del entorno.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/01-introduccion_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion_simple.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el script simple y confirmar los conteos de átomos.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion_simple.ipynb">01-introduccion_simple.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion_simple.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/01-introduccion_simple.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion_simple.py">01-introduccion_simple.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/01-introduccion.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el script complejo y cambiar el ligando por <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.sdf">ligand1.sdf</a>.
- Guardar la salida en un archivo de texto.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb">01-introduccion.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/01-introduccion.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion.py">01-introduccion.py</a>



## Puntos clave

- El curso trabaja con un sistema simple y uno complejo en paralelo.
- OpenMM y OpenFF Toolkit cubren la lectura de PDB y MOL/SDF.
  
<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Siguiente</a>
</div>
