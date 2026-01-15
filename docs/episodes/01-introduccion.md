---
layout: default
title: Episode 1 - Introduction (environment and data)
permalink: /episodes/01-introduccion/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Table of contents](#table-of-contents)
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Preliminary sequence and ligand checks](#preliminary-sequence-and-ligand-checks)
- [Alanine dipeptide](#alanine-dipeptide)
  - [Guided demo](#guided-demo)
  - [Exercise](#exercise)
  - [Notebooks and scripts](#notebooks-and-scripts)
- [Protein-ligand complex](#protein-ligand-complex)
  - [Guided demo](#guided-demo-1)
  - [Exercise](#exercise-1)
  - [Notebooks and scripts](#notebooks-and-scripts-1)
- [Key points](#key-points)
<!-- toc:end -->

## Duration

- **Session:** 45 min
- **Exercises:** 30 min

## Objectives

- Identify input files for alanine and the protein-ligand complex.
- Verify OpenMM and OpenFF Toolkit work in the environment.
- Load base data for both systems following the first-steps workflow.

## Content

- PDB and MOL/SDF as inputs.
- Reading alanine with OpenMM.
- Reading protein and ligand with OpenFF Toolkit.
- Working path `COURSE_DIR` and environment validation.

## Preliminary sequence and ligand checks

Before we launch the OpenMM demos, review two focused notebooks to explore the inputs.

- Sequence exploration: `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduction-sequence-check.ipynb">01-introduction-sequence-check.ipynb</a>` walks through parsing the PDB sequence, running UniProt and PDB BLAST queries, and collecting database summaries without touching OpenMM.
- Ligand exploration: `<a href="{{ site.baseurl }}/episodes/notebooks/01-introduction-ligand-check.ipynb">01-introduction-ligand-check.ipynb</a>` inspects the `$COURSE_DIR/data/complex/ligand1.sdf` molecule with RDKit to report formula, properties, and a quick sketch before any molecular dynamics setup.


## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/01-introduccion_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion_simple.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the simple script and confirm atom counts.

### Notebooks and scripts

- Recorrer esta versión guiada del episodio muestra la estructura del directorio `COURSE_DIR`, explica cómo leer un PDB con Biopython y prepara el sistema de alanina para el flujo de trabajo básico. (<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/01-introduccion.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/01-introduccion.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the complex script and swap the ligand for <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.sdf">ligand1.sdf</a>.
- Save the output to a text file.

### Notebooks and scripts

- Esta versión completa repasa la preparación del complejo proteína-ligando: secuencia rápida con Biopython, visualización con nglview y escritura de archivos Amber para la simulación. (<a href="{{ site.baseurl }}/episodes/notebooks/01-introduccion.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/01-introduccion.py">script</a>)

## Key points

- The course works with a simple and a complex system in parallel.
- OpenMM and OpenFF Toolkit cover reading PDB and MOL/SDF.

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/03-simulaciones-clasicas/">Next</a>
</div>
