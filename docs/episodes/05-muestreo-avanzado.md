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
- [Simulated annealing](#simulated-annealing)
- [Umbrella sampling](#umbrella-sampling)
- [SBMOpenMM resources](#sbmopenmm-resources)
- [Guide scripts (OpenMM Application Layer)](#guide-scripts-openmm-application-layer)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
- [References](#references)
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

## Simulated annealing

Refer to Example 5-1 in the “Simulated annealing” block of the [OpenMM User Guide: Advanced Simulation Examples](https://docs.openmm.org/latest/userguide/application/04_advanced_sim_examples.html#simulated-annealing) page. The sample is a minimal annealing loop that lowers the integrator temperature from 300 K to 0 K in 100 steps while running 1,000 time steps at each temperature:

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/05-muestreo-avanzado_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/05-muestreo-avanzado_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py" download>Download script (.py)</a></div></div>

The accompanying explanation notes that the loop simply updates the `LangevinMiddleIntegrator` temperature before each batch of 1,000 steps, which mirrors the ramping schedules we describe here. Keep that block as a quick reference when tuning integrator parameters or defining temperature sequences for your exercises so the practical code stays aligned with the theory in this episode.

### Replica exchange (REMD) context

The same temperature-ramping intuition underpins replica-exchange (REMD) methods. Instead of moving a single trajectory between high and low temperatures, REMD maintains multiple replicas at different thermodynamic states and swaps them through Metropolis trials to jump across barriers. `OpenMMTools` exposes `ReplicaExchangeSampler` (see https://openmmtools.readthedocs.io/en/stable/multistate.html#replicaexchangesampler-replica-exchange-among-thermodynamic-states)
to manage state definitions, collect swap statistics, and maintain detailed balance for arbitrary Hamiltonians.

For a broader explanation of why REMD is statistically valid and how the exchanges let you preserve canonical sampling, see the review `PMC6484850 <https://pmc.ncbi.nlm.nih.gov/articles/PMC6484850/>`__ [^remd_review], which discusses accepted protocols, ensemble averages, and implementations. Pairing replica exchange with the annealing schedule above gives two complementary views on negotiating rough energy landscapes: direct temperature ramps and ensemble-wide swaps that both respect the Boltzmann distribution.

### Gaussian accelerated MD (GaMD) with OpenMM

Gaussian accelerated molecular dynamics (GaMD) adds a boost potential that smooths the energy surface and reduces barriers without requiring pre-selected reaction coordinates. The `gamd-openmm` repository (https://github.com/MiaoLab20/gamd-openmm/tree/main) wraps OpenMM integrators to apply GaMD on top of existing routines, making it easy to run the same systems we already prepare for this lecture. Running GaMD alongside the simulated annealing loop lets you sample the same conformations while still estimating unbiased observables via the provided reweighting formulas [^gamd_paper].

### Martini coarse-graining with OpenMM

The MacCallum/Tieleman group has implemented both Martini 2 and Martini 3 directly in OpenMM (https://github.com/maccallumlab/martini_openmm/tree/martini3), taking advantage of OpenMM’s flexible force field layer to encode all Martini interactions, custom scripts to translate GROMACS topology files, and the GPU performance that makes OpenMM competitive with other engines. OpenMM’s extensibility—custom forces, fields, and integrators—lets these implementations reproduce the full Martini potential while still running inside the same infrastructure we use for simulated annealing and GaMD. See the Biophysical Journal paper by MacCallum et al. (2023) for the complete description of this OpenMM-friendly Martini pipeline [^martini_openmm].

## Umbrella sampling

The [Umbrella Sampling tutorial](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/umbrella_sampling.html) from the OpenMM Cookbook summarizes how to add Gaussian bias potentials along a collective variable, combine per-window histograms, and reconstruct \(F(x) = -k_B T \ln p(x)\). We adopt that setup by biasing the distance between two backbone atoms in alanine dipeptide, recording the sampled distances, and plotting a coarse free-energy-like curve from the histograms collected across the windows.

The notebook <a href="{{ site.baseurl }}/episodes/notebooks/05-umbrella-sampling.ipynb">05-umbrella-sampling.ipynb</a> reproduces those windows inline, gathers the sampled distances, stores the collected points in `$COURSE_DIR/data/umbrella_samples.csv`, and visualizes the emergent profile so you can poke at bin edges and the \(F(x) = -k_B T \ln p(x)\) conversion without leaving the browser. Its cells call the same systems described in the cookbook block, which demonstrates how each window’s bias shifts the sampled coordinate and why the umbrella reconstruction formula holds when you combine overlapping histograms.

The cookbook’s Bash section explains how to compile and run the [WHAM](https://github.com/choderalab/wham) binary to combine the windowed histograms and recover the free energy via the `<rmin> <rmax> <dx> <temperature> <histogram-files>` interface. Follow those steps locally with commands such as:

```bash
git clone https://github.com/choderalab/wham.git
cd wham
cmake .
make
./wham 0.18 0.35 0.005 300 $COURSE_DIR/data/umbrella_samples.csv > umbrella_free_energy.dat
```

These commands mirror the cookbook block: compile the OpenMM-provided WHAM executable, then pass your collected `umbrella_samples.csv` (or other histogram files you generated) to produce the final `umbrella_free_energy.dat` curve, matching the plot described on the page.

## SBMOpenMM resources

This subsection highlights how `sbmopenmm` (repository: https://github.com/CompBiochBiophLab/sbm-openmm) and its tutorials help explore advanced sampling through structure-based models. The [official guides](https://compbiochbiophlab.github.io/sbm-openmm/build/html/index.html) show how to build a Hamiltonian from the protein topology and apply external forces, a workflow that aligns with the methodology of this episode.

- [Quickstart and single-basin setup](https://github.com/CompBiochBiophLab/sbm-openmm#quick-start) – the top-level README section explains how to import a PDB, generate contact lists, and run simulations with custom forces (perfect for experimenting with temperature control and energy barriers).
- [Repository tutorials](https://github.com/CompBiochBiophLab/sbm-openmm/tree/master/tutorials) – a directory of folding and reaction-path scripts showing how to add restraint forces, log energies and gradients, and combine them with enhanced observation techniques to estimate free energy landscapes.

These references can serve as a baseline for seeing how simplified SBM setups follow the same sampling and reporting patterns we implement with a full OpenMM stack.

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

## References

- [^remd_review]: Sugita, Y. and Okamoto, Y., "Replica-exchange molecular dynamics method for protein folding", *Chem. Phys. Lett.* **2000**, 314, 141–151, and review `PMC6484850 <https://pmc.ncbi.nlm.nih.gov/articles/PMC6484850/>`__.
- [^gamd_paper]: Miao, Y.; Feher, V. A.; McCammon, J. A., "Gaussian accelerated molecular dynamics: Unconstrained enhanced sampling and free energy calculation", *J. Chem. Theory Comput.* **2015**, 11, 3584–3595.
- [^martini_openmm]: MacCallum, J. L.; Tieleman, D. P. et al., "Martini implementation in OpenMM," *Biophys. J.* **2023**, and the open-source repository at https://github.com/maccallumlab/martini_openmm/tree/martini3; the OpenMM adaptation covers both Martini 2 and 3 force fields and leverages OpenMM’s extensible API to process GROMACS topologies with the same flexibility as the native package.

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/04-analisis-trayectorias/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>
