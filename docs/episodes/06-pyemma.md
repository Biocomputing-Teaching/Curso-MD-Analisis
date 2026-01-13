---
layout: default
title: Episode 6 - MSMs with PyEMMA
permalink: /episodes/06-pyemma/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Guide scripts to prepare data](#guide-scripts-to-prepare-data)
- [Mathematical foundations](#mathematical-foundations)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
<!-- toc:end -->

## Duration

- **Session:** 70 min
- **Exercises:** 50 min

## Objectives

- Formulate the full MSM workflow (featurization, TICA, clustering, and estimation).
- Reinforce how transitions between microstates map to implied timescales.
- Compare results between the alanine system and the protein-ligand complex.

## Content

- Read and normalize data (DCD + PDB).
- Feature extraction and dimensionality reduction with TICA.
- Estimate the transition matrix and compute relaxation times.
- Validation (ITS, Chapman-Kolmogorov) and macrostate analysis (PCCA/TPT).

## Guide scripts to prepare data

The official OpenMM scripts generate the trajectories that we featurize with PyEMMA:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py): minimal input to validate the pipeline with alanine.
- [`simulateAmber.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateAmber.py): production with prmtop/inpcrd for the protein-ligand complex.
- [`simulateGromacs.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateGromacs.py) and [`simulateCharmm.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateCharmm.py): multi-format support when the system comes from other packages.

All trajectories are converted to DCD and paired with energy reports to verify stability before discretization.

## Mathematical foundations

We denote by $C_{ij}$ the count of transitions from microstate $i$ to $j$ at lag time $\tau$. The normalized transition matrix estimate

$$
T_{ij} = \frac{C_{ij}}{\sum_j C_{ij}}
$$

is the core of the model. Its eigenvalues $\lambda_k$ give timescales via

$$
\tau_k = -\frac{\tau}{\ln \lambda_k},
$$

and allow interpreting slow dynamics as near-invariant processes. The flux density is analyzed by combining stationary probabilities $\pi_i$ and PCCA transitions, while observables $(O)$ are evaluated as

$$
\langle O \rangle = \sum_i \pi_i O_i.
$$

The Chapman-Kolmogorov test ensures that $T(\tau)^n \approx T(n\tau)$ within statistical error.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/06-pyemma-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the full workflow on alanine (featurization, TICA, and clustering) and compute the three slowest relaxation times.
- Validate the model by showing the correspondence of $\tau_k$ with the ITS slopes.

### Key points

- The choice of distances and angles affects matrix $C$ and therefore the timescales.
- Feature normalization and clustering regularization prevent overestimation of $\lambda_k$.

### Notebooks and scripts

- This notebook builds MSMs with PyEMMA for alanine, covering featurization, TICA, clustering, and implied timescales. (<a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/06-pyemma-complex.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py" download>Download script (.py)</a></div></div>

### Exercise

- Apply the full MSM workflow to the protein complex: featurization, TICA, clustering, and estimation.
- Compute relaxation times and compare with the simple system. Which states dominate $\lambda_2$ and $\lambda_3$?

### Key points

- Output trajectories are validated with TPT and the extended transition matrix.
- Use relevant observables (distances, local energies) to reconstruct the state density.

### Notebooks and scripts

- This notebook repeats the MSM workflow on the protein-ligand trajectories, highlighting longer timescales and macrostate analysis. (<a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py">script</a>)

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Next</a>
</div>
