---
layout: default
title: Episodio 1 - Introduccion y entorno de trabajo
permalink: /episodes/01-introduccion/
---

## Overview

- **Teaching:** 45 min
- **Exercises:** 30 min

## Notebook

- [01-introduccion.ipynb]({{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb)

## Scripts

- [01-introduccion.py]({{ site.baseurl }}/episodes/scripts/01-introduccion.py)

## Objetivos

- Identificar los componentes de una simulacion MD.
- Configurar un entorno de trabajo reproducible.
- Ejecutar un primer script de OpenMM.

## Contenido

- Conceptos basicos: potencial, integrador, condiciones de borde.
- Unidades y convenciones en OpenMM.
- Estructura minima de un script en OpenMM.

## Demo guiada

### Crear un sistema minimo en OpenMM

```python
from openmm import app
import openmm as mm
from openmm import unit

pdb = app.PDBFile('alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds
)

integrator = mm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds
)

platform = mm.Platform.getPlatformByName('CUDA')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
state = simulation.context.getState(getEnergy=True)
print('Potential energy:', state.getPotentialEnergy())
```

## Ejercicio

- Ejecutar el script en CPU y GPU (si esta disponible).
- Registrar energia potencial y temperatura inicial.
- Guardar el entorno `conda list` en un archivo de texto.

## Puntos clave

- OpenMM permite flujos de trabajo reproducibles en Python.
- Es vital documentar unidades y parametros desde el inicio.
