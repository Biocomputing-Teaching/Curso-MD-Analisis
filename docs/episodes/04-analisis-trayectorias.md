---
layout: default
title: Episode 5 - Trajectory analysis with MDTraj
permalink: /episodes/04-analisis-trayectorias/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
<!-- toc:end -->

## Duration

- **Session:** 60 min
- **Exercises:** 45 min

## Objectives

- Analyze simple and complex trajectories.
- Reimage and align trajectories.
- Generate RMSD with MDTraj.

## Content

- Read DCD files with MDTraj.
- Reimage and align to the initial frame.
- Plots to compare systems.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the simple analysis with <a href="{{ site.baseurl }}/data/">traj.dcd</a>.

### Key points

- The metrics are comparable across systems.

### Notebooks and scripts

- This notebook reads the alanine trajectory with MDTraj, reimages frames, and compares RMSD metrics for the simple system. (<a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/04-analisis-trayectorias.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-analisis-trayectorias.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the complex analysis with <a href="{{ site.baseurl }}/reference/">output_traj.dcd</a>.
- Save the resulting figures.

### Key points

- MDTraj enables reproducible analysis with a few lines.

### Notebooks and scripts

- This notebook performs MDTraj analyses on the protein-ligand trajectory, re-imaging the complex and exporting RMSD plots to highlight differences. (<a href="{{ site.baseurl }}/episodes/notebooks/04-analisis-trayectorias.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-analisis-trayectorias.py">script</a>)

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>
