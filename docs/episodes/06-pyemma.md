---
layout: default
title: Episodio 6 - MSM con PyEMMA
permalink: /episodes/06-pyemma/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Siguiente</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Scripts de la guía para preparar datos](#scripts-de-la-guia-para-preparar-datos)
- [Fundamentos matemáticos](#fundamentos-matematicos)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
<!-- toc:end -->

## Duración

- **Sesión:** 70 min
- **Ejercicios:** 50 min

## Objetivos

- Formular el proceso completo de MSM (featurización, TICA, clustering y estimación).
- Reforzar cómo las transiciones entre microestados se traducen en escalas de tiempo internacionales.
- Comparar los resultados entre el sistema de alanina y el complejo proteína-ligando.

## Contenido

- Lectura y normalización de datos (DCD + PDB).
- Extracción de features y reducción de dimensión con TICA.
- Estimación de la matriz de transición y cálculo de tiempos de relajación.
- Validación (ITS, Chapman-Kolmogorov) y análisis de macrostados (PCCA/TPT).

## Scripts de la guía para preparar datos

Los scripts oficiales de OpenMM generan las trayectorias que luego featurizamos con PyEMMA:

- `simulatePdb.py`: entrada mínima para validar el pipeline con alanina.
- `simulateAmber.py`: producción con prmtop/inpcrd para el complejo proteína-ligando.
- `simulateGromacs.py` y `simulateCharmm.py`: soporte multi-formato cuando el sistema viene de otros paquetes.

Todas las trayectorias se convierten a DCD y se acompañan de reportes de energía para verificar estabilidad antes de discretizar.

## Fundamentos matemáticos

Denotamos con $C_{ij}$ el conteo de transiciones del microestado $i$ al $j$ en lag time $\tau$. La estimación de la matriz de transición normalizada

$$
T_{ij} = \frac{C_{ij}}{\sum_j C_{ij}}
$$

es el núcleo del modelo. Sus valores propios $\lambda_k$ dan escalas de tiempo mediante

$$
\tau_k = -\frac{\tau}{\ln \lambda_k},
$$

y permiten interpretar las dinámicas lentas como procesos casi invariantes. La densidad de flujo se analiza sumando probabilidades estacionarias $\pi_i$ y transiciones de PCCA, mientras que los observables $(O)$ se evalúan como

$$
\langle O \rangle = \sum_i \pi_i O_i.
$$

La verificación de Chapman-Kolmogorov asegura que $T(\tau)^n \approx T(n\tau)$ dentro del margen de error estadístico.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/06-pyemma-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el flujo completo en alanina (featurización, TICA y clustering) y calcular los tres tiempos de relajación más lentos.
- Validar el modelo mostrando la correspondencia de $\tau_k$ con las pendientes de ITS.

### Puntos clave

- La elección de las distancias y los ángulos influye en la matriz $C$ y, por tanto, en las escalas de tiempo.
- La normalización de las features y la regularización del clustering evitan sobreestima de $\lambda_k$.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb">06-pyemma-alanine.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-alanine.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/06-pyemma-alanine.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py">06-pyemma-alanine.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/06-pyemma-complex.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Aplicar MSM completo con el complejo proteico: featurización, TICA, clustering y estimación.
- Calcular los tiempos de relajación y comparar con el sistema simple. ¿Qué estados dominan $\lambda_2$ y $\lambda_3$?

### Puntos clave

- Las trayectorias de salida se validan con las TPT y la matriz de transición extendida.
- Usar relevantes observables (distancias, energías locales) para reconstruir la densidad de estados.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb">06-pyemma-complex.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-complex.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/06-pyemma-complex.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py">06-pyemma-complex.py</a>
<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Siguiente</a>
</div>
