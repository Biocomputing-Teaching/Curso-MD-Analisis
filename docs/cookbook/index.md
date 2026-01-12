---
layout: default
title: Libro de Cocina
permalink: /cookbook/
---

Este capítulo complementa los episodios clásicos con recetas cortas que explican los tutoriales más recientes del OpenMM Cookbook. Cada sección resume el contenido esencial y enlaza al notebook original para que el alumno pueda profundizar desde un punto de partida en castellano.

## 1. Cargar archivos y reportar resultados
El tutorial **Carga de archivos y reporte de resultados** ("Loading Input Files and Reporting Results") explica cómo aprovechar `OpenMM` para leer entradas generadas con AmberTools (prmtop+inpcrd) y, de forma analógica, con archivos de CHARMM, GROMACS o Tinker. Se muestra cómo registrar fácilmente información de energía, presión y coordenadas usando los reporteros construidos en la capa de aplicación, y también cómo definir reporteros personalizados que se disparen en cada paso de integración. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/loading_and_reporting.html)

## 2. Construir sistemas desde cero
En **Construcción de sistemas desde cero** ("Building Systems from Scratch") se abandona OpenMM Setup para trabajar directamente con la capa de librería. Se crea un sistema binario de Lennard-Jones, se definen las reglas de mezcla (Lorentz-Berthelot) y se elaboran manualmente las fuerzas base (`Nonbonded`, `HarmonicBond`, etc.). Este enfoque es ideal para experimentar con modelos sencillos o inspeccionar los objetos `System` antes de usarlos en simulaciones biomoleculares. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/building_systems.html)

## 3. Parámetros de simulación
El tutorial **Selección de parámetros de simulación** ("Selecting Values for Simulation Parameters") profundiza en las compensaciones entre velocidad y precisión. Detalla cómo elegir el tamaño del paso de integración según el modo rápido más importante, cómo jugar con las restricciones (`HBonds`, `AllBonds`, `HAngles`), cómo aplicar repartición de masa al hidrógeno y qué configuraciones de termostato/barostato son adecuadas para distintos ensambles. Es una hoja de ruta para ajustar la exactitud numérica incluso cuando trabajas con scripts generados previamente. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/simulation_parameters.html)

## 4. Energías libres alquímicas
**Cálculo de energías libres alquímicas** ("Alchemical Free Energy Calculations") muestra cómo usar `CustomForce` y parámetros globales para representar estados intermedios entre dos sistemas, y cómo combinarlo con `OpenMMTools` y `PyMBAR` para estimar diferencias de energía libre. Incluye recomendaciones para instalar dependencias (`openmmtools`, `pymbar`) y enlaza al código de referencia `YANK`. Es ideal para extender los episodios centrados en sistemas proteína-ligando. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/Alchemical_free_energy_calculations.html)

## 5. Polímeros de grano grueso
El capítulo **Polímeros de grano grueso** ("Implementing a Coarse-Grained Polymer Force Field") enseña cómo construir topologías de polímero con beads de Lennard-Jones, asignar masa y parámetros físicos realistas, y definir la fuerza tanto desde API como desde archivos XML. Combina técnicas de modelado redox y parametrización manual para explicar cómo se pueden adaptar simulaciones de proteínas a marcos de grano grueso. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/coarse_grained_polymer.html)

## 6. Intercambio de réplicas con temperado del soluto (REST)
**Simulación REST** ("Running a REST simulation") detalla la creación de un sistema REST: se copia un `System` "vanilla", se definen los átomos REST, se agregan `CustomBondForce`, `CustomAngleForce` y `CustomTorsionForce` con factores globales, y se escala el `NonbondedForce` para las interacciones modificadas. Explica la ejecución de réplicas y la coordinación con `OpenMMTools`. Este material encaja perfectamente después del episodio de MSM/fenómenos lentos. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/Running_a_REST_simulation.html)

## 7. Muestreo umbrella
Finalmente, **Muestreo umbrella** ("Umbrella Sampling") resume los fundamentos teóricos (CVs, potenciales de sesgo `bias` y perfil de energía libre), y explica cómo aplicar sesgos gaussianos para explorar regiones de alta energía. Incluye la relación \(F(x)=-k_B T \log(p(x))\) y cómo reconstruir el perfil combinando histogramas. [Abrir en el Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/umbrella_sampling.html)
