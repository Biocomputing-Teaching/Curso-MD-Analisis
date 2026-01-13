---
layout: default
title: Episode 7 - MSMs with deeptime
permalink: /episodes/07-deeptime/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Guide scripts to prepare data](#guide-scripts-to-prepare-data)
- [Operator foundations](#operator-foundations)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
<!-- toc:end -->

## Duration

- **Session:** 60 min
- **Exercises:** 45 min

## Objectives

- Explore MSMs and spectral analysis using the `deeptime` package.
- Introduce the Koopman operator interpretation and its relation to observable dynamics.
- Compare results with PyEMMA and reinforce feature selection.

## Content

- Features based on distances and torsions.
- TICA reduction and linear transformations.
- MSM estimation with `MaximumLikelihoodMSM` and validation.
- Spectral studies with `KoopmanModel` and macrostate analysis.

## Guide scripts to prepare data

In this episode we reuse OpenMM Application Layer scripts to build consistent data:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py): fast trajectories for Koopman validation.
- [`simulateAmber.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateAmber.py): long series for spectrum estimation and CK validation.
- [`simulateTinker.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulateTinker.py): AMOEBA examples to compare with classical models.

These outputs are processed with `deeptime` to stay consistent with PyEMMA.

## Operator foundations

The Koopman operator acts on an observable $g(\mathbf{x})$ by propagating it in time

$$
(\mathcal{K}_\tau g)(\mathbf{x}) = \mathbb{E}[g(\mathbf{x}_{t+\tau}) | \mathbf{x}_t = \mathbf{x}],
$$

which in practice is approximated with matrices $K_{ij}$ relating microstates. Diagonalizing $K$ yields eigenvalues $\lambda_k$ and eigenfunctions $\psi_k$, allowing reconstruction of the stationary state

$$
\rho(\mathbf{x}) \approx \sum_k \lambda_k \psi_k(\mathbf{x}).
$$

`deeptime` complements this view by estimating the transfer matrix and computing spectral projections to distinguish macrostates.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/07-deeptime-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-alanine.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-alanine.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the simple pipeline on alanine and tune the SVD decomposition for TICA.
- Compare the first three eigenvalues with the PyEMMA model: do the relaxation times match?

### Key points

- Feature normalization ensures the eigenfunctions are meaningful.
- Using `KoopmanModel` extracts observables that capture slow dynamics.

### Notebooks and scripts

- This notebook applies Deeptime workflows to the alanine dataset, tuning the SVD/TICA decomposition and comparing eigenvalues with the PyEMMA model. (<a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-alanine.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-alanine.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/07-deeptime-complex.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-complex.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-complex.py" download>Download script (.py)</a></div></div>

### Exercise

- Use the protein-ligand complex and build an MSM with `MaximumLikelihoodMSM`.
- Run `ChapmanKolmogorov` analysis and compare with PyEMMA, focusing on probability conservation.

### Key points

- The Koopman operator describes how the observable propagates with $\tau$ and compares with the classical transition matrix.
- Latent spaces help distinguish macrostates and study transition fluxes.

### Notebooks and scripts

- This notebook performs Deeptime spectral analysis on the protein-ligand complex, including Koopman models and Chapman-Kolmogorov validation. (<a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-complex.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/07-deeptime-complex.py">script</a>)

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
</div>
