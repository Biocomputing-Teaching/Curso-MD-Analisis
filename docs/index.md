---
showheader: true
layout: default
title: Home
---
# Curso práctico de Dinámica Molecular y Análisis de Trayectorias

**Fechas:** 15, 16 y 19 de enero de 2026  
**Lugar:** Universidad Andrés Bello, Concepción, Chile  
**Modalidad:** Presencial, práctico (mañana y tarde)  
**Idioma:** Español

## Requisitos previos
- Conocimientos básicos de bioinformática o química computacional
- Familiaridad con la línea de comandos de Linux
- Conocimientos básicos de Python
- Acceso a un clúster de cálculo habilitado con AMBER, Python, PyEMMA, etc.

---

## Programa del curso

### Día 1: Jueves 15 de enero
**Tema central:** Introducción práctica a la dinámica molecular clásica

**Mañana (09:00 - 13:00):**
- Bienvenida y presentación del curso
- Fundamentos de dinámica molecular (MD): teoría mínima imprescindible
- Preparación de sistemas con `tleap` (AMBER)
- Configuración de simulaciones clásicas (input files, topología, parámetros)
- Ejecución de simulaciones cortas en el clúster (`pmemd/cuda` o `pmemd.MPI`)

**Tarde (14:30 - 18:00):**
- Visualización y análisis básico de trayectorias con `cpptraj`
- Introducción al análisis con `MDAnalysis` (Python)
- Ejercicios prácticos: cálculo de RMSD, RMSF, radio de giro, distancias internas

---

### Día 2: Viernes 16 de enero
**Tema central:** Simulaciones avanzadas y modelos de Markov

**Mañana (09:00 - 13:00):**
- Técnicas avanzadas de simulación:
  - **Accelerated MD (aMD)**
  - **Steered MD (SMD)**
- Preparación de simulaciones aMD y SMD con AMBER
- Consideraciones de eficiencia y escalabilidad al clúster

**Tarde (14:30 - 18:00):**
- Construcción de modelos de Markov con `PyEMMA`:
  - Introducción a los **Markov State Models (MSM)**
  - Extracción de features con `MDAnalysis` + `PyEMMA`
  - Clustering y construcción del modelo
  - Análisis de tiempos de relajación, estados metaestables
- Ejemplo aplicado: detección de transiciones conformacionales en una proteína

---

### Día 3: Lunes 19 de enero (solo mañana)
**Tema central:** Proyectos individuales y consultoría

**Mañana (09:30 - 13:00):**
- Espacio para discusión de proyectos de los participantes
- Consultoría técnica personalizada:
  - Diseño experimental computacional
  - Análisis de trayectorias
  - Preparación de simulaciones en el clúster
  - Revisión de scripts y código Python
- Recomendaciones de buenas prácticas en simulaciones y análisis

---

## Herramientas utilizadas
- **AMBER** (tleap, pmemd, cpptraj)
- **MDAnalysis**
- **PyEMMA**
- **Python 3**, Jupyter Notebooks
- Clúster de cálculo con sistema de gestión de trabajos (ej. SLURM)

---

## Inscripción y contacto
- **Cupo limitado:** 15-20 participantes
- **Inscripción previa obligatoria**
- **Contacto:** Nombre del coordinador/a – `correo@ejemplo.cl`

---

## Observaciones
- Los participantes deben traer su propio computador portátil con acceso al clúster habilitado
- Se facilitarán notebooks y scripts de análisis
- El curso está orientado a estudiantes de posgrado, investigadores jóvenes y profesionales en bioquímica, biofísica o biología estructural

