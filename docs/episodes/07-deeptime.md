---
layout: default
title: Episodio 7 - MSM con deeptime
permalink: /episodes/07-deeptime/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Scripts de la guía para preparar datos](#scripts-de-la-guia-para-preparar-datos)
- [Fundamentos operatoriales](#fundamentos-operatoriales)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
<!-- toc:end -->

## Duración

- **Sesión:** 60 min
- **Ejercicios:** 45 min

## Objetivos

- Explorar MSM y análisis espectral utilizando el paquete `deeptime`.
- Introducir la interpretación del operador de Koopman y su relación con las dinámicas visibles.
- Comparar los resultados con PyEMMA y reforzar la elección de features.

## Contenido

- Features basadas en distancias y torsiones.
- Reducción con TICA y transformaciones lineales.
- Estimación de MSM con `MaximumLikelihoodMSM` y validación.
- Estudios espectrales con `KoopmanModel` y análisis de macrostados.

## Scripts de la guía para preparar datos

En este episodio reutilizamos los scripts del Application Layer de OpenMM para construir datos homogéneos:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py): trayectorias rápidas para validación de Koopman.
- [`simulateAmber.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateAmber.py): series largas para estimar espectros y validación CK.
- [`simulateTinker.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateTinker.py): ejemplos AMOEBA para comparar con modelos clásicos.

Estas salidas se procesan con `deeptime` para mantener consistencia con PyEMMA.

## Fundamentos operatoriales

El operador de Koopman actúa sobre una observación $g(\mathbf{x})$ propagándola en el tiempo

$$
(\mathcal{K}_\tau g)(\mathbf{x}) = \mathbb{E}[g(\mathbf{x}_{t+\tau}) | \mathbf{x}_t = \mathbf{x}],
$$

que en la práctica se aproxima con matrices $K_{ij}$ que relacionan microestados. La diagonalización de $K$ proporciona valores propios $\lambda_k$ y funciones propias $\psi_k$, de modo que se puede reconstruir el estado estacionario

$$
\rho(\mathbf{x}) \approx \sum_k \lambda_k \psi_k(\mathbf{x}).
$$

La simulación de `deeptime` complementa esta visión al permitir estimar la matriz de transferencia y calcular proyecciones espectrales para distinguir macrostados.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/07-deeptime-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-alanine.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-alanine.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el canal simple en alanina y ajustar la descomposición SVD para TICA.
- Comparar los primeros tres valores propios con el modelo PyEMMA: ¿coinciden los tiempos de relajación?

### Puntos clave

- La normalización de features garantiza que las funciones propias sean significativas.
- El uso de `KoopmanModel` permite extraer observables que capturan la dinámica lenta.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-alanine.ipynb">07-deeptime-alanine.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-alanine.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/07-deeptime-alanine.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-alanine.py">07-deeptime-alanine.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/07-deeptime-complex.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-complex.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-complex.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Usar el complejo proteína-ligando y construir una MSM con `MaximumLikelihoodMSM`.
- Realizar análisis de `ChapmanKolmogorov` y comparar con los resultados de PyEMMA, enfocándose en la conservación de probabilidades.

### Puntos clave

- El operador de Koopman describe cómo la observación se propaga con $\tau$ y se compara con la matriz de transición clásica.
- Los espacios latentes permiten distinguir macrostados y estudiar flujos de transición.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-complex.ipynb">07-deeptime-complex.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-complex.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/07-deeptime-complex.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-complex.py">07-deeptime-complex.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
</div>
