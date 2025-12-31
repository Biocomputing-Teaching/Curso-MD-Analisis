---
layout: default
title: Episodio 6 - Modelos de Markov con PyEMMA
permalink: /episodes/06-msm-pyemma/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [06-msm-pyemma.ipynb](../notebooks/06-msm-pyemma.ipynb)

## Scripts

- [06-msm-pyemma.py](../scripts/06-msm-pyemma.py)

## Objetivos

- Extraer features para MSM.
- Construir y validar un modelo de Markov.
- Interpretar estados metaestables.

## Contenido

- Seleccion de features y preprocesamiento.
- Clustering y matriz de transicion.
- Tiempos de relajacion y estados.

## Demo guiada

### Pipeline minimo en PyEMMA

```python
import pyemma

# Cargar features (ejemplo: distancias)
feat = pyemma.coordinates.featurizer('solvated.pdb')
feat.add_distances(feat.select_Backbone(), periodic=False)

traj_files = ['traj_1.dcd', 'traj_2.dcd']
X = pyemma.coordinates.load(traj_files, features=feat)

# TICA
tica = pyemma.coordinates.tica(X, lag=10)
Y = tica.get_output()

# Clustering
cl = pyemma.coordinates.cluster_kmeans(Y, k=50, max_iter=50)

# MSM
msm = pyemma.msm.estimate_markov_model(cl.dtrajs, lag=10)
print('Timescales:', msm.timescales_[:5])
```

## Ejercicio

- Probar k=30 y k=80 y comparar timescales.
- Identificar estados metaestables con PCCA.

## Puntos clave

- La eleccion de features domina la calidad del MSM.
- Validar con implied timescales evita sobreajuste.
