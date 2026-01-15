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

> **System focus** The PyEMMA episode works with the files derived from `protein2-ligand2.pdb`. Whenever `$COURSE_DIR/data/complex` is refreshed, re-run `bash docs/episodes/scripts/01-split-protein2-ligand2.sh` so `protein2.pdb` and `ligand2.pdb` match the source used in all MSM builders.

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Guide scripts to prepare data](#guide-scripts-to-prepare-data)
- [Mathematical foundations](#mathematical-foundations)
- [PyEMMA configuration and API primers](#pyemma-configuration-and-api-primers)
- [Tutorial pillars from the official PyEMMA guide](#tutorial-pillars-from-the-official-pyemma-guide)
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

## PyEMMA configuration and API primers

### Adjust runtime values step by step

PyEMMA exposes a `config` object that controls logging, progress bars, trajectory caching, parallelism, and even automatic version checks. Use it interactively (see [runtime configuration](http://www.emma-project.org/latest/Configuration.html)) to peek at the defaults and tweak them:

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py" download>Download script (.py)</a></div></div>

Every attribute listed under `pyemma.util._config.Config` in the [API index](http://www.emma-project.org/latest/api/index.html) can be inspected (e.g., `.mute`, `.traj_info_max_entries`, `.use_trajectory_lengths_cache`) and saved for reproduction.

### Persisting tweaks and logging

Save your tuned configuration with `config.save("/path/to/pyemma.cfg")`, then point `PYEMMA_CFG_DIR` at the directory that contains the file so every import reads the same settings. The configuration loader searches:

1. `./pyemma.cfg`
2. `$HOME/.pyemma/pyemma.cfg`
3. `~/pyemma.cfg`
4. `$PYTHONPATH/pyemma/pyemma.cfg` (always the default)

This lets you version-control different logging settings (the `.logging_config` attribute points to the `logging.yml` file) and reuse them across the simple and complex notebooks.

### Parallelism and resources

PyEMMA honors `PYEMMA_NJOBS` and `OMP_NUM_THREADS` (along with `SLURM_CPUS_ON_NODE` when available) to cap cluster usage. The tutorials for MSM estimation and clustering spawn both threads and processes, so set `PYEMMA_NJOBS` before launching Jupyter if you need to limit the total worker count.

### API highlights

- `pyemma.coordinates`: handles trajectory loading, caching, and transformation pipelines (see [`pyemma.coordinates.api`](http://www.emma-project.org/latest/api/coordinates/api.html)). Always feed it the same normalized data used in the notebooks so the coordinate cache is coherent.
- `pyemma.msm`: powers the transition matrix, MSM, and PCCA/TPT routines. The methods there (`estimate_markov_model`, `.timescales`, `.pcca`) are referenced in both scripts, so compare their docstrings with the API page to understand additional arguments (lag time, reversible, etc.).
- `pyemma.analysis` and `pyemma.plots`: supply helpers for implied timescale plotting, ITS, and Chapman-Kolmogorov validation; follow the API links when you extend the example workflows.

## Tutorial pillars from the official PyEMMA guide

### Tutorial 1 – simple MSM pipeline

The [official tutorial](http://www.emma-project.org/latest/tutorial.html#simple_msm) demonstrates the same featurization → tICA → clustering → MSM steps that the alanine notebook covers. Pay attention to the explanation of `pyemma.coordinates.source` (for loading DCD/PDB pairs) and the bootstrap/ITS sections so you can mirror the same plots without errors.

### Tutorial 2 – protein complex workflow

The advanced tutorial ([“Protein-ligand transitions”](http://www.emma-project.org/latest/tutorial.html#protein_ligand_complex)) walks through building MSMs on larger datasets with `pyemma.msm.markov_model` and `pyemma.plots.plot_cktest`. Use it to validate the complex notebook’s clustering parameters and to create the same diagnostic figures with the stored transition matrix.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/06-pyemma-alanine.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-alanine.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-alanine.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-alanine.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the full workflow on alanine (featurization, TICA, and clustering) and compute the three slowest relaxation times.
- Validate the model by showing the correspondence of $\tau_k$ with the ITS slopes.
### Tutorial check
Follow [Tutorial 1: simple MSM workflow](http://www.emma-project.org/latest/tutorial.html#simple_msm) to confirm you reuse the same coordinate pipelines and clustering diagnostics before moving on to the protein complex.

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
### Tutorial check
Use [Tutorial 2: protein-ligand complex](http://www.emma-project.org/latest/tutorial.html#protein_ligand_complex) to ensure the extended transition matrix and CK test match the notebook examples.

### Notebooks and scripts

- This notebook repeats the MSM workflow on the protein-ligand trajectories, highlighting longer timescales and macrostate analysis. (<a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-complex.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-complex.py">script</a>)

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Next</a>
</div>
