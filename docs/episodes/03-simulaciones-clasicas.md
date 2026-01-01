---
layout: default
title: Episodio 3 - Simulaciones clásicas y control de calidad
permalink: /episodes/03-simulaciones-clasicas/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [03-simulaciones-clasicas.ipynb]({{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb)
- [Ver en nbviewer](https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/03-simulaciones-clasicas.ipynb)

## Scripts

- [03-simulaciones-clasicas.py]({{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py)

## Objetivos

- Ejecutar minimización, calentamiento y producción.
- Configurar integradores y termostatos.
- Evaluar estabilidad y calidad de la trayectoria.

## Contenido

- Integradores (Langevin, Verlet).
- Control de temperatura y presión.
- Reporters para energía, trayectoria y log.

## Demo guiada

### Minimizar y correr una simulación corta

```python
from openmm import app, unit
import openmm as mm

pdb = app.PDBFile('../../data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)

simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy(maxIterations=100)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
simulation.step(200)

state = simulation.context.getState(getEnergy=True)
print('Potential energy:', state.getPotentialEnergy())
```

Fuente del script: [03-simulaciones-clasicas.py]({{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py)

## Ejercicio

- Correr 10,000 pasos y guardar `log.csv`.
- Graficar energía potencial vs tiempo.
- Identificar si hay deriva de energía.

## Puntos clave

- El control de calidad es parte del flujo de trabajo.
- Reporters bien configurados ahorran tiempo de depuración.
