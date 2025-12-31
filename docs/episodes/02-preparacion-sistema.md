---
layout: default
title: Episodio 2 - Preparacion de sistemas en OpenMM
permalink: /episodes/02-preparacion-sistema/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [02-preparacion-sistema.ipynb]({{ site.baseurl }}/episodes/notebooks/02-preparacion-sistema.ipynb)

## Scripts

- [02-preparacion-sistema.py]({{ site.baseurl }}/episodes/scripts/02-preparacion-sistema.py)

## Objetivos

- Construir un sistema solvado con OpenMM.
- Seleccionar campos de fuerza y parametros basicos.
- Generar archivos iniciales para simulacion.

## Contenido

- Topologias y archivos PDB.
- ForceField y Modeller en OpenMM.
- Solvatacion, iones y caja periodica.

## Demo guiada

### Solvatar un sistema con Modeller

```python
from openmm import app
from openmm import unit

pdb = app.PDBFile('protein.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0 * unit.nanometer,
    ionicStrength=0.15 * unit.molar
)

app.PDBFile.writeFile(modeller.topology, modeller.positions, open('solvated.pdb', 'w'))
print('Atoms after solvation:', modeller.topology.getNumAtoms())
```

## Ejercicio

- Crear una caja con padding de 1.2 nm y 0.1 M.
- Comparar el numero de atomos entre dos condiciones.
- Guardar el sistema en `solvated.pdb`.

## Puntos clave

- La preparacion determina la calidad de la simulacion.
- Mantener registros de cada paso facilita la reproducibilidad.
