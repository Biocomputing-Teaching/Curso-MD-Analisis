---
layout: default
title: Episodio 1 - Introducción y entorno de trabajo
permalink: /episodes/01-introduccion/
---

## Overview

- **Teaching:** 45 min
- **Exercises:** 30 min

## Notebook

- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/01-introduccion.ipynb">01-introduccion.ipynb</a>
- <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/01-introduccion.ipynb">Ver en nbviewer</a>

## Scripts

<div class="episode-nav">
  
  <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/">Todos los episodios</a>
  <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/02-preparacion-sistema/">Siguiente</a>
</div>


- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/scripts/01-introduccion.py">01-introduccion.py</a>

## Objetivos

- Identificar los componentes de una simulación MD.
- Configurar un entorno de trabajo reproducible.
- Ejecutar un primer script de OpenMM.

## Contenido

- Conceptos básicos: potencial, integrador, condiciones de borde.
- Unidades y convenciones en OpenMM.
- Estructura mínima de un script en OpenMM.

## Demo guiada

### Crear un sistema mínimo en OpenMM

```python
from openmm import app
import openmm as mm

pdb = app.PDBFile('../../data/alanine-dipeptide.pdb')
print('OpenMM', mm.__version__)
print('Atoms:', pdb.topology.getNumAtoms())
```

Fuente del script: <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/scripts/01-introduccion.py">01-introduccion.py</a>

Fuente del script: <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/scripts/01-introduccion.py">01-introduccion.py</a>

## Ejercicio

- Ejecutar el script en CPU y GPU (si está disponible).
- Registrar energía potencial y temperatura inicial.
- Guardar el entorno `conda list` en un archivo de texto.

## Puntos clave

- OpenMM permite flujos de trabajo reproducibles en Python.
- Es vital documentar unidades y parámetros desde el inicio.
