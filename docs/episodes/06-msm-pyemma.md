---
layout: default
title: Episodio 6 - Modelos de Markov con PyEMMA
permalink: /episodes/06-msm-pyemma/
---

## Overview

- **Teaching:** 60 min
- **Exercises:** 45 min

## Notebook

- [06-msm-pyemma.ipynb]({{ site.baseurl }}/episodes/notebooks/06-msm-pyemma.ipynb)
- [Ver en nbviewer](https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/06-msm-pyemma.ipynb)

## Scripts

- [06-msm-pyemma.py]({{ site.baseurl }}/episodes/scripts/06-msm-pyemma.py)

## Objetivos

- Extraer features para MSM.
- Construir y validar un modelo de Markov.
- Interpretar estados metaestables.

## Contenido

- Selección de features y preprocesamiento.
- Clustering y matriz de transición.
- Tiempos de relajación y estados.

## Demo guiada

### Pipeline mínimo en PyEMMA

```python
import numpy as np
import pyemma

x1 = np.random.normal(-1.0, 0.2, size=(1000, 1))
x2 = np.random.normal(1.0, 0.2, size=(1000, 1))
X = [np.vstack([x1, x2])]

tica = pyemma.coordinates.tica(X, lag=10)
Y = tica.get_output()

cl = pyemma.coordinates.cluster_kmeans(Y, k=20, max_iter=50)
msm = pyemma.msm.estimate_markov_model(cl.dtrajs, lag=10)
print('Timescales:', msm.timescales_[:5])
```

Fuente del script: [06-msm-pyemma.py]({{ site.baseurl }}/episodes/scripts/06-msm-pyemma.py)

## Ejercicio

- Probar k=30 y k=80 y comparar timescales.
- Identificar estados metaestables con PCCA.

## Puntos clave

- La elección de features domina la calidad del MSM.
- Validar con implied timescales evita sobreajuste.
