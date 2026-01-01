---
layout: default
title: Episodio 4 - Análisis de trayectorias con MDAnalysis
permalink: /episodes/04-analisis-trayectorias/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [04-analisis-trayectorias.ipynb]({{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb)
- [Ver en nbviewer](https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/04-analisis-trayectorias.ipynb)

## Scripts

- [04-analisis-trayectorias.py]({{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py)

## Objetivos

- Cargar trayectorias con MDAnalysis.
- Calcular RMSD, RMSF y distancias internas.
- Generar gráficas reproducibles.

## Contenido

- Universes, selections y atom groups.
- Análisis básico y estadísticas.
- Exportación de resultados.

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

Fuente del script: [04-analisis-trayectorias.py]({{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py)

## Ejercicio

- Calcular RMSF por residuo.
- Comparar RMSD entre dos réplicas.
- Guardar figuras en `docs/figures/`.

## Puntos clave

- La selección correcta de átomos es crítica.
- Automatizar análisis mejora la trazabilidad.
