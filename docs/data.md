---
layout: default
title: Datos
permalink: /data/
---

## Datos del curso

Los datos de ejemplo están en `docs/data/` y se usan en los ejercicios y notebooks.

### Archivos principales

- `alanine-dipeptide.pdb`: estructura mínima para pruebas con OpenMM.
- `alanine-dipeptide-multi.pdb`: PDB multi-model para análisis de trayectorias sin DCD.
- `example_rmsd.csv`: serie de RMSD de ejemplo para gráficas.

### Generar trayectoria DCD

Para crear un DCD corto con OpenMM:

```bash
python scripts/generate_example_dcd.py
```

Esto genera `docs/data/alanine-dipeptide.dcd`.
