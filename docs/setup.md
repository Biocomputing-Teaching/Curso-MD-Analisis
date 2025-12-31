---
layout: default
title: Preparacion
---

## Requisitos de software

- Python 3.10 o superior
- OpenMM (con soporte CUDA si hay GPU)
- MDAnalysis
- PyEMMA
- JupyterLab o Jupyter Notebook
- Git

## Instalacion recomendada (conda)

```bash
conda create -n md-openmm python=3.10
conda activate md-openmm
conda install -c conda-forge openmm mdanalysis pyemma jupyterlab numpy scipy matplotlib
```

## Verificacion rapida

```bash
python - <<'PY'
import openmm
import mdanalysis as mda
import pyemma
print('openmm', openmm.__version__)
print('mdanalysis', mda.__version__)
print('pyemma', pyemma.__version__)
PY
```

## Datos y materiales

Los ejemplos del curso se almacenaran en `docs/data/` y `docs/figures/`.
Los notebooks estan en `docs/episodes/notebooks/`. Si se usan datasets grandes se compartiran por un enlace externo durante el curso.
Los scripts equivalentes estan en `docs/episodes/scripts/`.
Para generar una trayectoria corta en DCD, ejecutar `python scripts/generate_example_dcd.py`.
