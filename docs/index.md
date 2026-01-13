---
showheader: true
layout: default
permalink: /
---

# Practical Course on Molecular Dynamics and Trajectory Analysis

## Organization

- **Dates:** January 15, 16, and 19, 2026
- **Location:** Universidad Andres Bello, Concepcion, Chile
- **Format:** In-person, hands-on (morning and afternoon)
- **Language:** Spanish
- **Instructor:** <a href="mailto:jordi.villa@uvic.cat">Jordi Villa i Freixa</a>, <a href="https://mon.uvic.cat/cbbl">Computational Biochemistry and Biophysics Lab</a>, <a href="https://www.uvic.cat">Universitat de Vic - Universitat Central de Catalunya</a>, <a href="https://iris-cc.cat">IRIS-CC</a>
- **Limited spots:** 15-20 participants
- **Registration required in advance**
- **Contact:** Veronica Andrea Jimenez Curihual <a href="mailto:veronica.jimenez@unab.cl">veronica.jimenez@unab.cl</a>

## Prerequisites

- Basic knowledge of molecular simulations or computational chemistry
- Familiarity with the Linux command line
- Basic Python knowledge
- Access to a GPU cluster or workstation, ideally with OpenMM

## Course program

### Day 1: Thursday, January 15
**Main theme:** Fundamentals of classical simulations with OpenMM

**Morning (09:00 - 13:00):**
- Welcome and course overview
- Minimal theory of molecular dynamics
- System preparation in OpenMM (topology, force field, solvent)
- Parameterization with OpenMM and OpenFF
- Running short GPU simulations

**Afternoon (14:30 - 18:00):**
- Basic trajectory analysis with MDAnalysis
- Exercises: RMSD, RMSF, radius of gyration, internal distances
- Quality control best practices

### Day 2: Friday, January 16
**Main theme:** Advanced techniques and Markov models

**Morning (09:00 - 13:00):**
- Advanced techniques in OpenMM:
  - **Steered MD (SMD)**
  - **Umbrella Sampling**
  - **Metadynamics**
- Setup and execution of advanced protocols
- Efficiency and scalability considerations

**Afternoon (14:30 - 18:00):**
- Markov models with PyEMMA
  - MSM introduction
  - Feature extraction with MDAnalysis + PyEMMA
  - Clustering and model construction
  - Relaxation times and metastable states
- Applied case: conformational transitions in a protein

### Day 3: Monday, January 19 (morning only)
**Main theme:** Individual projects and consultation

**Morning (09:30 - 13:00):**
- Participant project discussions
- Personalized technical consultation
- Best-practice recommendations for simulations and analysis

## Course talk

Expanded material based on the course beamer, with one HTML chapter per session:

- <a href="{{ site.baseurl }}/talks/">Access the talk and chapters</a>


## Collaborators

With collaboration from [![COZYME COST ACTION: Pan-European Network on Computational Redesign of Enzymes]({{ site.baseurl }}/figures/logo-cozyme-160x34.png)](https://cozyme.eu/). Pan-European Network on Computational Redesign of Enzymes.

![COST-Actions]({{ site.baseurl }}/figures/COST_LOGO_white_tranparent_small-300x93.png)

Supported by Fondecyt Fondo Nacional de Desarrollo Cientifico y Tecnologico

![Conicyt]({{ site.baseurl }}/figures/logo-conicyt-chile-en-marcha.png)
