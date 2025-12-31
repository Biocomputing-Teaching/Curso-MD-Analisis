---
layout: default
title: Episodio 3 - Simulaciones clasicas y control de calidad
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [03-simulaciones-clasicas.ipynb](../notebooks/03-simulaciones-clasicas.ipynb)

## Scripts

- [03-simulaciones-clasicas.py](../scripts/03-simulaciones-clasicas.py)

## Objetivos

- Ejecutar minimizacion, calentamiento y produccion.
- Configurar integradores y termostatos.
- Evaluar estabilidad y calidad de la trayectoria.

## Contenido

- Integradores (Langevin, Verlet).
- Control de temperatura y presion.
- Reporters para energia, trayectoria y log.

## Demo guiada

### Minimizar y correr una simulacion corta

```python
from openmm import app
import openmm as mm
from openmm import unit

pdb = app.PDBFile('solvated.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds
)

integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)

simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

print('Minimizing...')
simulation.minimizeEnergy(maxIterations=500)

simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

simulation.reporters.append(app.StateDataReporter('log.csv', 1000, step=True, temperature=True, potentialEnergy=True, density=True))
simulation.reporters.append(app.DCDReporter('traj.dcd', 1000))

simulation.step(5000)
```

## Ejercicio

- Correr 10,000 pasos y guardar `log.csv`.
- Graficar energia potencial vs tiempo.
- Identificar si hay deriva de energia.

## Puntos clave

- El control de calidad es parte del flujo de trabajo.
- Reporters bien configurados ahorran tiempo de depuracion.
