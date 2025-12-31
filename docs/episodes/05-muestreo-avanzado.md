---
layout: default
title: Episodio 5 - Muestreo avanzado en OpenMM
permalink: /episodes/05-muestreo-avanzado/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [05-muestreo-avanzado.ipynb]({{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb)
- [Ver en nbviewer](https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/05-muestreo-avanzado.ipynb)

## Scripts

- [05-muestreo-avanzado.py]({{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py)

## Objetivos

- Implementar SMD y umbrella sampling en OpenMM.
- Diseñar un protocolo de metadynamics.
- Evaluar convergencia basica.

## Contenido

- Potenciales de restriccion y CVs.
- Ventanas y combinacion de perfiles.
- Buenas practicas de muestreo.

## Demo guiada

### Umbrella sampling con una restriccion armonica

```python
import openmm as mm
from openmm import app, unit

# Sistema base ya preparado
# Se crea una restriccion armonica sobre una distancia

force = mm.CustomBondForce('0.5*k*(r-r0)^2')
force.addPerBondParameter('k')
force.addPerBondParameter('r0')

# Atomos i, j (indices de ejemplo)
force.addBond(10, 120, [1000.0 * unit.kilojoule_per_mole / unit.nanometer**2, 0.35 * unit.nanometer])

system.addForce(force)
```

## Ejercicio

- Definir 5 ventanas con diferentes `r0`.
- Ejecutar 0.5 ns por ventana.
- Comparar perfiles de energia potencial.

## Puntos clave

- El muestreo avanzado requiere diagnosticos de convergencia.
- Documentar cada ventana y parametro es esencial.
