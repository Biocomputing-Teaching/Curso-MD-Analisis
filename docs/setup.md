---
layout: default
title: Preparación
permalink: /setup/
---

## Instalación rápida (OpenMM)

OpenMM se instala con `conda` o `pip`. La guía oficial recomienda Miniconda para un entorno aislado y reproducible.

### Opción recomendada (conda)

```bash
conda create -n md-openmm python=3.10
conda activate md-openmm
conda install -c conda-forge openmm
```

Si tienes GPU NVIDIA y quieres CUDA:

```bash
conda install -c conda-forge openmm cuda-version=12
```

### Alternativa (pip)

```bash
pip install openmm
```

Si usas GPU NVIDIA con CUDA 12:

```bash
pip install "openmm[cuda12]"
```

Si usas GPU AMD (HIP):

```bash
pip install "openmm[hip6]"
```

## OpenMM-Setup (GUI)

La aplicación `openmm-setup` genera scripts listos para correr y corrige problemas comunes en estructuras.

```bash
conda install -c conda-forge openmm-setup
```

Para abrirla:

```bash
openmm-setup
```

## Paquetes del curso

Dentro del entorno `md-openmm`, instala los paquetes que usamos en los episodios:

```bash
conda install -c conda-forge jupyterlab mdanalysis mdtraj deeptime openff-toolkit openmmforcefields pdbfixer
```

PyEMMA es opcional (el proyecto está congelado). Si lo necesitas:

```bash
pip install pyemma
```

## Ejecutar simulaciones (scripts de OpenMM)

La guía de OpenMM incluye scripts en `python-examples`:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py) (PDB directo)
- [`simulateAmber.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateAmber.py) (prmtop + inpcrd)
- [`simulateCharmm.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateCharmm.py) (psf + crd)
- [`simulateGromacs.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateGromacs.py) (top + gro)
- [`simulateTinker.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateTinker.py) (AMOEBA)

Ejemplo genérico:

```bash
python simulatePdb.py entrada.pdb
```

## Datos del curso

Usamos un directorio externo controlado por `COURSE_DIR`:

```bash
export COURSE_DIR=~/Concepcion26
```

Estructura esperada:

- `data/` archivos de entrada
- `results/` salidas y análisis

Descarga de datos:

```bash
mkdir -p "$COURSE_DIR"
cd "$COURSE_DIR"
curl -L -o Course-MD-Data.zip https://github.com/Biocomputing-Teaching/Course-MD-Data/archive/refs/heads/main.zip
unzip -q Course-MD-Data.zip
rm -rf "$COURSE_DIR/data"
mkdir -p "$COURSE_DIR/data"
cp -R Course-MD-Data-main/* "$COURSE_DIR/data"
rm -rf Course-MD-Data-main Course-MD-Data.zip
```

## Verificación rápida

```bash
python - <<'PY'
import openmm
print('openmm', openmm.__version__)
PY
```
