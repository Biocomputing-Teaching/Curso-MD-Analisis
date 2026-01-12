---
showheader: true
layout: default
permalink: /
---

# Curso práctico de Dinámica Molecular y Análisis de Trayectorias

## Organización

- **Fechas:** 15, 16 y 19 de enero de 2026
- **Lugar:** Universidad Andrés Bello, Concepción, Chile
- **Modalidad:** Presencial, práctico (mañana y tarde)
- **Idioma:** Español
- **Profesor:** <a href="mailto:jordi.villa@uvic.cat">Jordi Villà i Freixa</a>, <a href="https://mon.uvic.cat/cbbl">Computational Biochemistry and Biophysics Lab</a>, <a href="https://www.uvic.cat">Universitat de Vic - Universitat Central de Catalunya</a>, <a href="https://iris-cc.cat">IRIS-CC</a>
- **Cupo limitado:** 15-20 participantes
- **Inscripción previa obligatoria**
- **Contacto:** Verónica Andrea Jiménez Curihual <a href="mailto:veronica.jimenez@unab.cl">veronica.jimenez@unab.cl</a>

## Requisitos previos

- Conocimientos básicos de simulaciones moleculares o química computacional
- Familiaridad con la línea de comandos de Linux
- Conocimientos básicos de Python
- Acceso a un clúster o workstation con GPU, idealmente con OpenMM

## Programa del curso

### Día 1: Jueves 15 de enero
**Tema central:** Fundamentos de simulaciones clásicas con OpenMM

**Mañana (09:00 - 13:00):**
- Bienvenida y presentación del curso
- Teoría mínima de dinámica molecular
- Preparación de sistemas en OpenMM (topología, fuerza, solvente)
- Parametrización con OpenMM y OpenFF
- Ejecución de simulaciones cortas en GPU

**Tarde (14:30 - 18:00):**
- Análisis básico de trayectorias con MDAnalysis
- Ejercicios: RMSD, RMSF, radio de giro, distancias internas
- Buenas prácticas de control de calidad

### Día 2: Viernes 16 de enero
**Tema central:** Técnicas avanzadas y modelos de Markov

**Mañana (09:00 - 13:00):**
- Técnicas avanzadas en OpenMM:
  - **Steered MD (SMD)**
  - **Umbrella Sampling**
  - **Metadynamics**
- Configuración y ejecución de protocolos avanzados
- Consideraciones de eficiencia y escalabilidad

**Tarde (14:30 - 18:00):**
- Modelos de Markov con PyEMMA
  - Introducción a MSM
  - Extracción de features con MDAnalysis + PyEMMA
  - Clustering y construcción del modelo
  - Tiempos de relajación y estados metaestables
- Caso aplicado: transiciones conformacionales en una proteína

### Día 3: Lunes 19 de enero (solo mañana)
**Tema central:** Proyectos individuales y consultoría

**Mañana (09:30 - 13:00):**
- Discusión de proyectos de los participantes
- Consultoría técnica personalizada
- Recomendaciones de buenas prácticas en simulaciones y análisis

## Charla del curso

Material ampliado basado en el beamer del curso, con un capítulo HTML por sesión:

- <a href="{{ site.baseurl }}/talks/">Acceder a la charla y a los capítulos</a>


## Colaboradores

Con la colaboración de [![COZYME COST ACTION: Pan-European Network on Computational Redesign of Enzymes]({{ site.baseurl }}/figures/logo-cozyme-160x34.png)](https://cozyme.eu/). Pan-European Network on Computational Redesign of Enzymes.

![COST-Actions]({{ site.baseurl }}/figures/COST_LOGO_white_tranparent_small-300x93.png)

Con el soporte de Fondecyt Fondo Nacional de Desarrollo Científico y Tecnológico

![Conicyt]({{ site.baseurl }}/figures/logo-conicyt-chile-en-marcha.png)
