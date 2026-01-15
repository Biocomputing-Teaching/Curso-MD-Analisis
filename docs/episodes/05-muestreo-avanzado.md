---
layout: default
title: Episode 5 - Advanced sampling and mathematical foundations
permalink: /episodes/05-muestreo-avanzado/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>

> **System focus** The advanced sampling episode keeps using the `protein2-ligand2.pdb` pair. Run `bash docs/episodes/scripts/01-split-protein2-ligand2.sh` after refreshing `$COURSE_DIR/data/complex` to regenerate the canonical `protein2.pdb` and `ligand2.pdb` before launching any sampling notebooks.

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Thermodynamic foundations](#thermodynamic-foundations)
- [Connection to the official OpenMM guide](#connection-to-the-official-openmm-guide)
- [Guide scripts (OpenMM Application Layer)](#guide-scripts-openmm-application-layer)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
<!-- toc:end -->

## Duration

- **Session:** 60 min
- **Exercises:** 45 min

## Objectives

- Connect Boltzmann statistics with advanced sampling algorithms.
- Show how each course script represents a mathematical strategy (annealing, external forces, and quantitative reporting).
- Validate periodic solvation in alanine and protein-ligand complex trajectories.

## Content

- Canonical distribution and free energy formulation.
- Effective temperature control via adaptive integrators.
- Energy and force reports to estimate potential derivatives.
- Solvation and periodic box strategies with OpenMM.

## Thermodynamic foundations

The canonical distribution

$$
\rho(\mathbf{x}) = \frac{1}{Z} \exp\left(-\beta U(\mathbf{x})\right),
$$

where $\beta = 1/(k_B T)$ sets the relation between energy and entropy. Simulated annealing strategies modify \(T\) gradually to explore high- and low-potential regions without abrupt jumps, and allow estimating free energy differences via

$$
\Delta G = -k_B T \ln \frac{Z_{\text{solv}}}{Z_{\text{no\ soluto}}}.
$$

To estimate restraint stability, also monitor the norm of the potential gradient

$$
\mathbf{g} = \nabla U(\mathbf{x}), \qquad \langle \|\mathbf{g}\|^2 \rangle = \frac{1}{N} \sum_{i=1}^N \|\mathbf{g}_i\|^2,
$$

which is reported periodically to identify numerical drift.

![Energy funnel and sampling strategies]({{ site.baseurl }}/figures/function_funnel.png)

## Connection to the official OpenMM guide

The [OpenMM Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/loading_and_reporting.html) chapters and related sections (building systems, simulation parameters, alchemical free energy) describe steps consistent with our scripts: defining the system, customizing external forces, and logging events. From there we adopt three pillars:

- Simulated annealing by adjusting temperature in integrators to favor jumps between wells.
- `CustomExternalForce` external forces to constrain regions and study mechanical responses.
- `State`-based reporting to retrieve potential energies and forces and calibrate convergence.

All linked scripts use this architecture and serve as guided examples for each system.

## Guide scripts (OpenMM Application Layer)

To feed the analysis, we rely on the official scripts described in the OpenMM guide:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py): fast simulations from PDB, useful for annealing tests.
- [`simulateAmber.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateAmber.py): AMBER input (prmtop/inpcrd) to compare energies with and without solvent.
- [`simulateCharmm.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateCharmm.py): PSF/CRD topologies to evaluate restraints and reporting.
- [`simulateGromacs.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateGromacs.py): GROMACS input (top/gro), ideal for validating preprocessing.
- [`simulateTinker.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateTinker.py): AMOEBA, useful for contrasting polarization effects.

In each case, the outputs (DCD, CSV, and energy reports) are used in this episode's exercises.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the simple workflow with spherical restraints and record the value of $\langle \|\mathbf{g}\|^2 \rangle$ before and after annealing.
- Compare the average potential energy with and without CustomExternalForce.

### Key points

- Keep reports every 10 steps to monitor integrator drift.
- Adjust $\beta$ to preserve numerical stability in small systems.

### Notebooks and scripts

- This notebook executes the advanced sampling routine for alanine (restraints, custom forces, reporting energies) and tracks gradient norms. (<a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py" download>Download script (.py)</a></div></div>

### Exercise

- Run with `--solvate` and `--padding 12`, recording solvent density and Coulomb energy.
- Try another water model and analyze the variation in $\Delta G$ estimated from partition differences.

### Key points

- Periodic solvation stabilizes $\Delta G$ on long time scales and improves annealing convergence.
- Compare energy distributions to validate extended sampling.

### Notebooks and scripts

- This notebook covers the complex system with solvation options, barostat control, and energy/force logging to validate enhanced sampling workflows. (<a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">script</a>)

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>
