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
- Evaluar convergencia básica.

## Contenido

- Potenciales de restricción y CVs.
- Ventanas y combinación de perfiles.
- Buenas prácticas de muestreo.

## Demo guiada

### Umbrella sampling con una restricción armónica

```python
import openmm as mm
from openmm import app, unit

pdb = app.PDBFile('../../data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)

force = mm.CustomBondForce('0.5*k*(r-r0)^2')
force.addPerBondParameter('k')
force.addPerBondParameter('r0')
force.addBond(1, 4, [1000.0 * unit.kilojoule_per_mole / unit.nanometer**2, 0.35 * unit.nanometer])
system.addForce(force)

print('Custom forces:', system.getNumForces())
```

Fuente del script: [05-muestreo-avanzado.py]({{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py)

## Ejercicio

- Definir 5 ventanas con diferentes `r0`.
- Ejecutar 0.5 ns por ventana.
- Comparar perfiles de energía potencial.

## Puntos clave

- El muestreo avanzado requiere diagnósticos de convergencia.
- Documentar cada ventana y parámetro es esencial.
