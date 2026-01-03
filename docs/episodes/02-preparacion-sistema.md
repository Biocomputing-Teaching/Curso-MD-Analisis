---
layout: default
title: Episodio 2 - Preparación de sistemas en OpenMM
permalink: /episodes/02-preparacion-sistema/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/02-preparacion-sistema.ipynb">02-preparacion-sistema.ipynb</a>
- <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/02-preparacion-sistema.ipynb">Ver en nbviewer</a>

## Scripts

<div class="episode-nav">
  <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/01-introduccion/">Anterior</a>
  <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/">Todos los episodios</a>
  <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/03-simulaciones-clasicas/">Siguiente</a>
</div>


- <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/scripts/02-preparacion-sistema.py">02-preparacion-sistema.py</a>

## Objetivos

- Construir un sistema solvado con OpenMM.
- Seleccionar campos de fuerza y parámetros básicos.
- Generar archivos iniciales para simulación.

## Contenido

- Topologías y archivos PDB.
- ForceField y Modeller en OpenMM.
- Solvatación, iones y caja periódica.

## Demo guiada

### Solvatar un sistema con Modeller

```python
from openmm import app, unit

pdb = app.PDBFile('../../data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0 * unit.nanometer)

print('Atoms after solvation:', modeller.topology.getNumAtoms())
```

Fuente del script: <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/scripts/02-preparacion-sistema.py">02-preparacion-sistema.py</a>

Fuente del script: <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/scripts/02-preparacion-sistema.py">02-preparacion-sistema.py</a>

## Ejercicio

- Crear una caja con padding de 1.2 nm y 0.1 M.
- Comparar el número de átomos entre dos condiciones.
- Guardar el sistema en <a href="https://biocomputing-teaching.github.io/Curso-MD-Analisis/data/">solvated.pdb</a>.

## Puntos clave

- La preparación determina la calidad de la simulación.
- Mantener registros de cada paso facilita la reproducibilidad.
