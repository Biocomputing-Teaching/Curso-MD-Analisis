---
layout: default
title: Episodio 5 - Muestreo Avanzado y fundamentos matemáticos
permalink: /episodes/05-muestreo-avanzado/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Siguiente</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Fundamentos termodinámicos](#fundamentos-termodinamicos)
- [Conexión con la guía oficial de OpenMM](#conexion-con-la-guia-oficial-de-openmm)
- [Scripts de la guía (OpenMM Application Layer)](#scripts-de-la-guia-openmm-application-layer)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
<!-- toc:end -->

## Duración

- **Sesión:** 60 min
- **Ejercicios:** 45 min

## Objetivos

- Establecer la relación entre la estadística de Boltzmann y los algoritmos de muestreo avanzado.
- Mostrar cómo cada script del curso representa una estrategia matemática (recocidos, fuerzas externas y reportes cuantitativos).
- Validar la solvatación periódica en las trayectorias de alanina y del complejo proteína-ligando.

## Contenido

- Formulación de la distribución canónica y energía libre.
- Control de la temperatura efectiva mediante integradores adaptativos.
- Reportes de energías y fuerzas para estimar derivadas del potencial.
- Estrategias de solvatación y caja periódica con OpenMM.

## Fundamentos termodinámicos

La distribución canónica

$$
\rho(\mathbf{x}) = \frac{1}{Z} \exp\left(-\beta U(\mathbf{x})\right),
$$

donde $\beta = 1/(k_B T)$ fija la relación entre energía y entropía. Las estrategias de recocido simulado modifican \(T\) gradualmente para recorrer regiones de alto y bajo potencial sin saltos bruscos, y permiten estimar diferencias de energía libre mediante

$$
\Delta G = -k_B T \ln \frac{Z_{\text{solv}}}{Z_{\text{no\ soluto}}}.
$$

Para estimar la estabilidad de las restricciones también vigila la norma del gradiente del potencial

$$
\mathbf{g} = \nabla U(\mathbf{x}), \qquad \langle \|\mathbf{g}\|^2 \rangle = \frac{1}{N} \sum_{i=1}^N \|\mathbf{g}_i\|^2,
$$

que se reporta periódicamente para identificar drift numérico.

![Embudo energético y estrategias de muestreo]({{ site.baseurl }}/figures/function_funnel.png)

## Conexión con la guía oficial de OpenMM

Los capítulos de [OpenMM Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/loading_and_reporting.html) y sus secciones relacionadas (building systems, simulation parameters, alchemical free energy) describen pasos coherentes con nuestros scripts: definir el sistema, personalizar fuerzas externas y registrar eventos. De allí adoptamos tres pilares:

- Recocido simulado ajustando la temperatura en los integradores para favorecer saltos entre pozos.
- Fuerzas externas `CustomExternalForce` para acotar regiones y estudiar respuestas mecánicas.
- Reportes basados en `State` para recuperar energías potenciales y fuerzas y calibrar la convergencia.

Todos los scripts enlazados utilizan esta arquitectura y sirven como ejemplos guiados para cada sistema.

## Scripts de la guía (OpenMM Application Layer)

Para alimentar el análisis, apoyamos los scripts oficiales descritos en la guía de OpenMM:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py): simulaciones rápidas desde PDB, útil para pruebas de recocido.
- [`simulateAmber.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateAmber.py): entrada AMBER (prmtop/inpcrd) para comparar energías con y sin solvente.
- [`simulateCharmm.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateCharmm.py): topologías PSF/CRD para evaluar restricciones y reportes.
- [`simulateGromacs.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateGromacs.py): entrada GROMACS (top/gro), ideal para validar preprocesado.
- [`simulateTinker.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateTinker.py): AMOEBA, útil para contrastar efectos de polarización.

En cada caso, las salidas (DCD, CSV y reportes de energía) se usan en los ejercicios de este episodio.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar el flujo simple con las restricciones esféricas y registrar el valor de $\langle \|\mathbf{g}\|^2 \rangle$ antes y después del recocido.
- Comparar la energía potencial promedio con y sin CustomExternalForce.

### Puntos clave

- Mantener los reportes cada 10 pasos para monitorizar la deriva del integrador.
- Ajustar $\beta$ para conservar la estabilidad numérica en sistemas pequeños.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb">05-muestreo-avanzado_simple.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/05-muestreo-avanzado_simple.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Ejecutar con `--solvate` y `--padding 12` registrando la densidad del solvente y la energía de Coulomb.
- Probar otro modelo de agua y analizar la variación en $\Delta G$ estimado mediante la diferencia de particiones.

### Puntos clave

- La solvatación periódica estabiliza $\Delta G$ en escalas largas y mejora la convergencia del recocido.
- Comparar distribuciones de energía para validar el muestreo extendido.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb">05-muestreo-avanzado.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/05-muestreo-avanzado.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>
