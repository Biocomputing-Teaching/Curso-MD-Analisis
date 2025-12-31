---
layout: default
title: Episodio 4 - Analisis de trayectorias con MDAnalysis
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [04-analisis-trayectorias.ipynb](../notebooks/04-analisis-trayectorias.ipynb)

## Scripts

- [04-analisis-trayectorias.py](../scripts/04-analisis-trayectorias.py)

## Objetivos

- Cargar trayectorias con MDAnalysis.
- Calcular RMSD, RMSF y distancias internas.
- Generar graficas reproducibles.

## Contenido

- Universes, selections y atom groups.
- Analisis basico y estadisticas.
- Exportacion de resultados.

## Demo guiada

### RMSD con MDAnalysis

```python
import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

pdb_path = '../../data/alanine-dipeptide.pdb'
dcd_path = '../../data/alanine-dipeptide.dcd'

if os.path.exists(dcd_path):
    u = mda.Universe(pdb_path, dcd_path)
    output_name = 'rmsd_dcd.png'
else:
    u = mda.Universe('../../data/alanine-dipeptide-multi.pdb')
    output_name = 'rmsd_multimodel.png'

atoms = u.select_atoms('all')

R = rms.RMSD(atoms, atoms)
R.run()

plt.plot(R.results.rmsd[:, 1], R.results.rmsd[:, 2])
plt.xlabel('Frame')
plt.ylabel('RMSD (A)')
plt.title('RMSD')
plt.savefig(output_name, dpi=150)
```

## Ejercicio

- Calcular RMSF por residuo.
- Comparar RMSD entre dos replicas.
- Guardar figuras en `docs/figures/`.

## Puntos clave

- La seleccion correcta de atomos es critica.
- Automatizar analisis mejora la trazabilidad.
