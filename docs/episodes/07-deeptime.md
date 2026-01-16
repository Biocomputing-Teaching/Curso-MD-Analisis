---
layout: default
title: Episode 7 - MSMs with deeptime
permalink: /episodes/07-deeptime/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
</div>

> **System focus** Deeptime sessions now follow the Alanine dipeptide example from `deeptime` (see https://deeptime-ml.github.io/latest/notebooks/examples/ala2-example.html); all required hydra data is fetched from `mdshare`, so no additional preprocessing is needed.

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Guide scripts to prepare data](#guide-scripts-to-prepare-data)
- [Operator foundations](#operator-foundations)
- [Deeptime showcase](#deeptime-showcase)
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

## Deeptime showcase

### Guided demo

This episode now follows the Deeptime-managed Alanine dipeptide example (`Ala2 Example <https://deeptime-ml.github.io/latest/notebooks/examples/ala2-example.html>`__). We imported that notebook into the repository (see the iframe below) so you can open it locally in the same environment we use for the course. It fetches the heavy-atom positions and backbone dihedrals from `mdshare`, runs the same Torch-let Koopman/TICA flow, and highlights how the spectral models compare with PyEMMA’s timelines.

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/07-deeptime-ala2.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/07-deeptime-ala2.ipynb" download>Download notebook</a></div></div>

### Exercise

- Re-run the deeptime example and experiment with different projector sizes (SVD truncation) in the TICA step.  
- Compare the first three Koopman eigenvalues with the PyEMMA model discussed in Episode 6—are the slow-time scales in agreement?

### Key points

- `mdshare` provides the alanine datasets (`alanine-dipeptide-3x250ns-*`) used by the official example, so no extra downloads are needed.
- The notebook blends Torch-based Koopman models with deeptime’s `KoopmanModel`/`MaximumLikelihoodMSM` estimators for the same relaxation timescales you saw in the PyEMMA workflow.

### Notebooks and scripts

- The episode now ships a single Deeptime notebook that mirrors the upstream `Ala2 Example` (see the iframe and download link above).

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
</div>
