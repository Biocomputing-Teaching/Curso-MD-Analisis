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
conda install -c conda-forge openmm mdanalysis deeptime jupyterlab numpy scipy matplotlib
```

## Verificación rápida

```bash
python - <<'PY'
import openmm
import MDAnalysis as mda
import deeptime
print('openmm', openmm.__version__)
print('MDAnalysis', mda.__version__)
print('deeptime', deeptime.__version__)
PY
```

## Guía para participantes

- Traer computador con acceso al clúster o GPU local.
- Instalar el entorno indicado en la página de preparación.
- Confirmar acceso a datos y notebooks antes del inicio.

## Requisitos previos

- Manejo básico de Linux.
- Se asume un manejo básico de python. Explorar, por ejemplo, este enlace: [Biomolecular simulation: python basics](https://emleddin.github.io/2020-06-05-py-tutorial/01-introduction/index.html)
