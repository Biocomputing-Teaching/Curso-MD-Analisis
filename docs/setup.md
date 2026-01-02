---
layout: default
title: Preparación
permalink: /setup/
---

## Requisitos de software

- Python 3.10 o superior
- OpenMM (con soporte CUDA si hay GPU)
- MDAnalysis
- PyEMMA
- JupyterLab o Jupyter Notebook
- Git

## Instalación recomendada (conda)

```bash
conda create -n md-openmm python=3.10
conda activate md-openmm
conda install -c conda-forge openmm mdanalysis pyemma jupyterlab numpy scipy matplotlib
```

## Verificación rápida

```bash
python - <<'PY'
import openmm
import MDAnalysis as mda
import pyemma
print('openmm', openmm.__version__)
print('MDAnalysis', mda.__version__)
print('pyemma', pyemma.__version__)
PY
```

## Guía para participantes

- Traer computador con acceso al clúster o GPU local.
- Instalar el entorno indicado en la página de preparación.
- Confirmar acceso a datos y notebooks antes del inicio.

## Requisitos previos

- Manejo básico de Linux.
- Python básico.
