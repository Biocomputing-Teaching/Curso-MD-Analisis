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
- [PyEMMA configuration and API primers](#pyemma-configuration-and-api-primers)
- [Tutorial pillars from the official PyEMMA guide](#tutorial-pillars-from-the-official-pyemma-guide)
- [Pentapeptide showcase](#pentapeptide-showcase)
- [Project](#project)
- [References](#references)
<!-- toc:end -->

## Duration

- **Session:** 70 min
- **Exercises:** 50 min

## Objectives

- Formulate the full MSM workflow (featurization, TICA, clustering, and estimation).
- Reinforce how transitions between microstates map to implied timescales.
- Compare results between the pentapeptide showcase and the protein-ligand complex.

## Content

- Read and normalize data (DCD + PDB).
- Feature extraction and dimensionality reduction with TICA.
- Estimate the transition matrix and compute relaxation times.
- Validation (ITS, Chapman-Kolmogorov) and macrostate analysis (PCCA/TPT).

## Guide scripts to prepare data

The official OpenMM scripts generate the trajectories that we featurize with PyEMMA:

- [`simulatePdb.py`](https://github.com/openmm/openmm/blob/master/examples/python-examples/simulatePdb.py): minimal input to validate the pipeline with a pentapeptide demo.
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

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/06-pyemma-pentapeptide.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/06-pyemma-pentapeptide.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/06-pyemma-pentapeptide.py" download>Download script (.py)</a></div></div>

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

The [official tutorial](http://www.emma-project.org/latest/tutorial.html#simple_msm) demonstrates the same featurization → tICA → clustering → MSM steps that the pentapeptide example covers. Pay attention to the explanation of `pyemma.coordinates.source` (for loading DCD/PDB pairs) and the bootstrap/ITS sections so you can mirror the same plots without errors.

For supplementary recipes and worked examples on MSM construction and visualization, see the `msmhelper tutorials <https://moldyn.github.io/msmhelper/tutorials/>`__ provided by MolDyn; they revisit the same ideas (featurization, dimension reduction, clustering, observables) and occasionally feature alternative observables you can compare against the PyEMMA output.

### Tutorial 2 – protein complex workflow

The advanced tutorial ([“Protein-ligand transitions”](http://www.emma-project.org/latest/tutorial.html#protein_ligand_complex)) walks through building MSMs on larger datasets with `pyemma.msm.markov_model` and `pyemma.plots.plot_cktest`. Use it to validate the complex notebook’s clustering parameters and to create the same diagnostic figures with the stored transition matrix.

## Pentapeptide showcase

### Guided demo

We now follow the official PyEMMA **00 – Showcase pentapeptide: a PyEMMA walkthrough** (PyEMMA 2.5.5+9.ge9257b9.dirty) as our guided notebook. The walkthrough loads the pentapeptide dataset from `mdshare`, adds backbone torsions/positions/distances, scores the feature sets with VAMP2, applies TICA, clusters the projected data, builds a Bayesian MSM, validates it, and finally inspects PCCA/TPT, observables, and manuscript-quality figures.

#### Data input and featurization

PyEMMA internally inspects `mdtraj.version.version`. When you install MDTraj 1.9+ the `mdtraj` module does not expose a `version` attribute, which triggers `AttributeError`. The first cell in the notebook (and our helper script) patches the module so PyEMMA can access the already-installed `version` module:

```python
import importlib
import mdtraj

if not hasattr(mdtraj, 'version'):
    mdtraj.version = importlib.import_module('mdtraj.version')
```

```python
import mdshare
import pyemma

pdb = mdshare.fetch('pentapeptide-impl-solv.pdb', working_directory='data')
files = mdshare.fetch('pentapeptide-*-500ns-impl-solv.xtc', working_directory='data')
```

Three feature groups are built with `pyemma.coordinates.featurizer`: backbone torsions (with cossin expansion), backbone heavy atom positions, and backbone-heavy atom distances, all computed with `periodic=False` after a global alignment that removes the simulation box.

#### Feature selection

Each feature list is ranked with a cross-validated VAMP2 score that splits the 25 trajectories into training/validation, fits `pyemma.coordinates.vamp`, and scores the held-out half. Across lags of 0.1–0.5 ns the backbone torsions consistently return the highest kinetic variance, so the notebook continues with that featurization.

#### Coordinate transform and discretization

Time-lagged independent component analysis compresses the torsion data into a few kinetic components, and k-means discretizes that low-dimensional projection into microstates.

##### TICA

Time-lagged Independent Component Analysis (TICA) pulls out the kinetically slow degrees of freedom by searching for linear combinations of features with maximum time-lagged autocorrelation; the resulting ICs (``IC1``, ``IC2``…) form the low-dimensional space where MSM transitions are best resolved.

The torsion-based data are projected with `pyemma.coordinates.tica(..., lag=5)` (0.5 ns) using the default kinetic map scaling and enough components to capture 95% of the kinetic variance. Histograms and densities on the first four ICs show well-separated metastable basins, and one trajectory is plotted over time to illustrate discrete jumps.

##### Discretization

k-means clustering is applied to the TICA coordinates. A VAMP2 scan over `[5, 10, 30, 75, 200, 450]` cluster centers reveals saturation near `k=75`, which is therefore fixed for the remainder of the tutorial (`fixed_seed=1` for reproducibility). The centers and discrete trajectories are shown against the first two TICA dimensions.

#### MSM estimation and validation

Implied timescales are computed with `pyemma.msm.its(cluster.dtrajs, lags=50, nits=10, errors='bayes')` and plotted in ns. The timescales converge above 0.5 ns, which justifies a lag time of 5 steps for `pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=5, dt_traj='0.1 ns')`. Active-state and active-count fractions are printed, and a Chapman–Kolmogorov test with `mlags=6` over five metastable sets confirms Markovian consistency.

The Chapman–Kolmogorov (CK) test checks whether the transition probabilities propagate consistently across multiples of the lag time. By comparing estimates obtained at `kτ` directly to the k-step propagation of the original lag-τ transition matrix, we verify the Markov assumption: a passing CK test indicates the kinetics remain invariant under coarse-graining, lending confidence to the MSM's lag choice.

#### MSM spectral analysis

The sample mean/std of the first 15 timescales are plotted together with their separations, showing the gap between the 4th and 5th process. Stationary distributions, free energies, and the first four right eigenvectors are visualized as contour plots on the first two TICA coordinates, revealing how slow processes connect the dense clusters.

#### PCCA & TPT

`msm.pcca(5)` produces fuzzy memberships that are contoured across IC1/IC2, and the argument-wise `argmax` provides crisp metastable assignments. Sample structures are saved via `pyemma.coordinates.save_traj` and visualized with `nglview`, while the tutorial also prints state free energies, MFPT tables, and direction-specific MFPTs, showing metastable state 1’s short lifetime. Transition path theory between states 2 and 4 uses `pyemma.msm.tpt`, and the committor is plotted as a contour map.

Transition Path Theory (TPT) is the mathematical framework that uses these metastable sets to define reactive pathways, generate committor probabilities, and compute net fluxes between states. The committor tells us the chance that a trajectory starting in a configuration reaches the designated product set before returning to the reactant set, while the reactive flux highlights which microstate transitions actually carry probability current under equilibrium assumptions[^noe_tpt]. Together they transform the MSM into a pathway-level picture of the slow dynamics and guide which states to sample for experimental observables.

The committor gives the probability that a configuration will reach a specified final macrostate before returning to the initial macrostate, so contour maps of the committor reveal the transition pathways’ location, while the reactive fluxes computed by TPT count the net probability current between states. Examining both the committor surface and the fluxes ensures we understand the preferred transition routes in addition to the equilibrium populations[^noe_tpt].

#### Expectations and observables

Samples of 20 frames per Markov state are converted to MDTraj objects for SASA (`shrake_rupley`) and radius of gyration (`compute_rg`). The Bayesian MSM evaluates ensemble expectations, standard deviations, and confidence intervals for the radius of gyration, while Trp-1 SASA is projected on the TICA plane to produce contour plots, autocorrelation functions, and relaxation curves for ML and Bayesian MSMs. Comparing metastable state expectations against the global average highlights state 1 as the most distinct ensemble.

#### Hidden Markov models

The notebook notes that hidden Markov models, initialized from the PCCA macrostates, further reduce discretization error and provide a compact kinetic description; readers are directed to Notebook 07 for hands-on HMM exercises.

#### Assembling manuscript figures & wrapping up

The final chapters reconfigure Matplotlib defaults and reproduce Figures 2–6 from the Showcase notebook, including the system overview, ITS + CK panels, free-energy and state maps, flux sketches, and Trp-1 autocorrelation/relaxation plots. Saved PDFs live in `data/figure_[2-6].pdf`. A Wrapping up section gestures to the remaining notebooks for deeper theory, TPT, HMMs, observables, and debugging tips.

### Exercise

- Reproduce the full PyEMMA 00 – Showcase pentapeptide workflow: featurization, VAMP2 scoring, TICA, clustering, MSM estimation, ITS/CK validation, PCCA/TPT, and observables.
- Confirm that the three slowest implied timescales are stable across MSM re-estimations and that backbone torsions remain the best feature set.

### Tutorial check

Use [PyEMMA Tutorial 1: simple MSM workflow](http://www.emma-project.org/latest/tutorial.html#simple_msm) to cross-check the coordinate pipelines and diagnostics before moving on to the protein complex.

### Key points

- Cross-validated VAMP2 scores guard against premature feature selection.
- TICA lag, clustering size, and MSM lag should be tuned by ITS convergence and Chapman–Kolmogorov tests.
- Observables, PCCA/TPT, and HMM post-processing translate MSMs into experimentally accessible quantities.

### Notebooks and scripts

- This episode now follows the PyEMMA 00 – Showcase pentapeptide notebook (PyEMMA 2.5.5+9.ge9257b9.dirty). View it at <https://pyemma.org/latest/tutorial.html#showcase> and rerun it locally using `mdshare`’s pentapeptide dataset.

## Project

Before the session ends, launch a project that covers the entire PyEMMA workflow and submit a ZIP archive containing the script/notebook you ran, any required trajectory/topology files (or download instructions), and the key plots that demonstrate your MSM analysis. Choose one of the following:

1. **Small custom system** – pick a peptide or fragment that runs within a few hours on your machine, generate trajectories (aggregated ≥50 ns), featurize with torsions/positions, then build and validate an MSM with at least three timescales plotted.
2. **Stanford Alanine decapeptide** – recreate the published Ala2 example, compute PCCA memberships, and display how the slowest implied timescales compare with the PyEMMA showcase.
3. **Your existing large simulation** – take any dataset you already have, subsample it as needed, apply PyEMMA/TICA/clustering, and show that the implied timescales/VAMP2 scores converge with the lag and cluster choices you pick.

Deliverables:

- A zipped folder containing your runnable script (or notebook), supporting data (or commands to pull it), and the MSM plots (timescales, state maps, etc.).
- A short README describing which path you chose, simulation length, MSM lag/cluster parameters, and the takeaways from the plots.
- If the data is hosted online, include the exact download steps so we can reproduce your simulation.

The goal is to get comfortable with the full simulation‑to‑MSM pipeline and produce evidence (scripts + plots) you can share by this afternoon.

## References

- [^noe_tpt]: Noé, F.; Schütte, C.; Vanden-Eijnden, E.; Reich, L.; Weikl, T. R. (2009). *Constructing the equilibrium ensemble of folding pathways from short off-equilibrium simulations*. Proc. Natl. Acad. Sci. U.S.A., 106, 19011–19016.

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-muestreo-avanzado/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/07-deeptime/">Next</a>
</div>
