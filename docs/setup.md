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
- OpenFF Toolkit
- OpenMMForceFields
- PDBFixer
- MDTraj
- Plotly (con Kaleido para exportar SVG)
- Biopython (lectura de PDB en notebooks)
- nglview (visualización interactiva)

## Instalación recomendada (conda)

```bash
conda create -n md-openmm python=3.10
conda activate md-openmm
conda install -c conda-forge openmm mdanalysis pyemma deeptime jupyterlab numpy scipy matplotlib openff-toolkit openmmforcefields pdbfixer mdtraj plotly python-kaleido biopython nglview ipywidgets jupyterlab_widgets
```

Nota: el desarrollo de pyemma se ha abandonado, con lo que su instalación es más problemática, y ahora los esfuerzos se destinan al desarollo de deeptime con una visión más amplia del análisis de datos en simulación.

## Notas

- Si `nglview` falla por un conflicto entre `bleach` y `tinycss2`, prueba:

```bash
conda install -c conda-forge "tinycss2<1.5" bleach
```

## Directorio de trabajo del curso

Todos los archivos generados (trayectorias, gráficos, salidas) se guardan fuera del repositorio en un directorio externo
definido por la variable de entorno `COURSE_DIR`. Valor por defecto: `~/Concepcion26`.

```bash
export COURSE_DIR=~/Concepcion26
```

Es recomendable que esta línea esté incluída en `$HOME/.zshrc` (Mac OS) o en `$HOME/.bashrc` (linux).

Estructura esperada dentro de `COURSE_DIR`:

- `data/` con todos los archivos de entrada.
- `results/` para todos los resultados generados por los episodios.

Descarga el repositorio <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data">Course-MD-Data</a> (los archivos están en la raíz) y colócalo en `COURSE_DIR/data`:

```bash
mkdir -p $COURSE_DIR
cd $COURSE_DIR
curl -L -o Course-MD-Data.zip https://github.com/Biocomputing-Teaching/Course-MD-Data/archive/refs/heads/main.zip
unzip -q Course-MD-Data.zip
rm -rf "$COURSE_DIR/data"
mkdir -p "$COURSE_DIR/data"
cp -R Course-MD-Data-main/* "$COURSE_DIR/data"
rm -rf Course-MD-Data-main Course-MD-Data.zip
```

Si lo prefieres, puedes clonar `Course-MD-Data` con `git` y copiar su contenido a `COURSE_DIR/data`.

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
- Se asume un manejo básico de python. Explorar, por ejemplo, este enlace: <a href="https://emleddin.github.io/2020-06-05-py-tutorial/01-introduction/index.html">Biomolecular simulation: python basics</a>
