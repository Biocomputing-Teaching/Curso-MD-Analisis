---
layout: default
title: Episode 4 - Advanced sampling and mathematical foundations
permalink: /episodes/04-muestreo-avanzado/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/05-analisis-trayectorias/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [Thermodynamic foundations](#thermodynamic-foundations)
- [Official OpenMM tutorials](#official-openmm-tutorials)
- [Simulated annealing](#simulated-annealing)
- [Gaussian accelerated MD (GaMD) with OpenMM](#gaussian-accelerated-md-gamd-with-openmm)
- [Martini coarse-graining with OpenMM](#martini-coarse-graining-with-openmm)
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

## Official OpenMM tutorials

The [OpenMM Cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/loading_and_reporting.html) chapters and related sections (building systems, simulation parameters, alchemical free energy) describe steps consistent with our scripts: defining the system, customizing external forces, and logging events. From there we adopt three pillars:

- Simulated annealing by adjusting temperature in integrators to favor jumps between wells.
- `CustomExternalForce` external forces to constrain regions and study mechanical responses.
- `State`-based reporting to retrieve potential energies and forces and calibrate convergence.

All linked scripts use this architecture and serve as guided examples for each system.

## Simulated annealing

Refer to Example 5-1 in the “Simulated annealing” block of the [OpenMM User Guide: Advanced Simulation Examples](https://docs.openmm.org/latest/userguide/application/04_advanced_sim_examples.html#simulated-annealing) page. The sample here ramps the integrator temperature from 300 K toward lower values in configurable increments while running a fixed number of steps per temperature, mirroring the temperature schedule in that guide.

Script source: <a href="{{ site.baseurl }}/episodes/scripts/openmm_advanced_annealing.py">openmm_advanced_annealing.py</a>

### Guided demo

<!-- sync-from: docs/episodes/notebooks/04-simulated-annealing.ipynb -->

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-simulated-annealing.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-simulated-annealing.ipynb" download>Download notebook</a></div></div>

### Exercise

- Run the notebook from within your `COURSE_DIR` so the annealing loop reads `data/alanine-dipeptide.pdb` and dumps snapshots into `results/05-muestreo-avanzado/openmm-annealing`; tweak `steps` or `increments` to explore more aggressive ramps.


### Notebooks

- Simulated annealing with OpenMM and configurable temperature ramps. (<a href="{{ site.baseurl }}/episodes/notebooks/04-simulated-annealing.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/openmm_advanced_annealing.py">script derived from this notebook</a>)

The accompanying explanation notes that the loop simply updates the `LangevinMiddleIntegrator` temperature before each batch of 1,000 steps, which mirrors the ramping schedules we describe here. Keep that block as a quick reference when tuning integrator parameters or defining temperature sequences for your exercises so the practical code stays aligned with the theory in this episode.

### Replica exchange (REMD) context

The same temperature-ramping intuition underpins replica-exchange (REMD) methods. Instead of moving a single trajectory between high and low temperatures, REMD maintains multiple replicas at different thermodynamic states and swaps them through Metropolis trials to jump across barriers. `OpenMMTools` exposes `ReplicaExchangeSampler` (see https://openmmtools.readthedocs.io/en/stable/multistate.html#replicaexchangesampler-replica-exchange-among-thermodynamic-states)
to manage state definitions, collect swap statistics, and maintain detailed balance for arbitrary Hamiltonians.

For a broader explanation of why REMD is statistically valid and how the exchanges let you preserve canonical sampling, see the review `PMC6484850 <https://pmc.ncbi.nlm.nih.gov/articles/PMC6484850/>`__ [^remd_review], which discusses accepted protocols, ensemble averages, and implementations. Pairing replica exchange with the annealing schedule above gives two complementary views on negotiating rough energy landscapes: direct temperature ramps and ensemble-wide swaps that both respect the Boltzmann distribution.

### REMD for alanine dipeptide

This guided block runs REMD on the alanine dipeptide system already prepared for the previous exercises. Four Langevin trajectories at staggered temperatures attempt swaps between nearest neighbors, so the algorithm crosses barriers without needing a single long trajectory at the highest temperature. The notebook below is the canonical source; the Python script shown beside it is generated from the same cells and can be used directly in batch pipelines.

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-remd.py">04-remd.py</a>

<!-- sync-from: docs/episodes/notebooks/04-remd.ipynb -->

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-remd.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-remd.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-remd.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the notebook from your `COURSE_DIR`, let it read `data/alanine-dipeptide.pdb`, and track the swaps logged at `results/05-muestreo-avanzado/remd/remd_swaps.csv`. Increase the number of iterations and explorer `steps-per-iteration` to see how swap acceptance changes with longer integration windows.


### Notebooks

- Replica exchange with multiple replicas (alanine dipeptide). (<a href="{{ site.baseurl }}/episodes/notebooks/04-remd.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-remd.py">script derived from this notebook</a>)

## Gaussian accelerated MD (GaMD) with OpenMM

Gaussian accelerated molecular dynamics (GaMD) adds a boost potential that smooths the energy surface and reduces barriers without requiring pre-selected reaction coordinates. The `gamd-openmm` repository (https://github.com/MiaoLab20/gamd-openmm/tree/main) wraps OpenMM integrators to apply GaMD on top of existing routines, making it easy to run the same systems we already prepare for this lecture. Running GaMD alongside the simulated annealing loop lets you sample the same conformations while still estimating unbiased observables via the provided reweighting formulas [^gamd_paper].

### Guided demo: GaMD on alanine dipeptide

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-gamd.py">04-gamd.py</a>

<!-- sync-from: docs/episodes/notebooks/04-gamd.ipynb -->

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-gamd.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-gamd.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-gamd.py" download>Download script (.py)</a></div></div>

### Exercise

- Launch the notebook from inside your `COURSE_DIR`, keep `data/alanine-dipeptide.pdb` as the input, and inspect `results/05-muestreo-avanzado/gamd/gamd_log.csv` to see how the boost amplitude adapts across iterations.

### Notebooks

- GaMD with OpenMM and adaptive boost logging. (<a href="{{ site.baseurl }}/episodes/notebooks/04-gamd.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-gamd.py">script derived from this notebook</a>)

## Martini coarse-graining with OpenMM

The MacCallum/Tieleman group has implemented both Martini 2 and Martini 3 directly in OpenMM (https://github.com/maccallumlab/martini_openmm/tree/martini3), taking advantage of OpenMM’s flexible force field layer to encode all Martini interactions, custom scripts to translate GROMACS topology files, and the GPU performance that makes OpenMM competitive with other engines. OpenMM’s extensibility—custom forces, fields, and integrators—lets these implementations reproduce the full Martini potential while still running inside the same infrastructure we use for simulated annealing and GaMD. See the Biophysical Journal paper by MacCallum et al. (2023) for the complete description of this OpenMM-friendly Martini pipeline [^martini_openmm].

### Guided demo: Martini coarse-graining for the complex system

Setting up the complex protein-ligand system for Martini teaches how coarse-grained beads map back to the all-atom coordinates and how a Martini topology feeds into the rest of the pipeline. The notebook below is canonical; the helper script records the Martini command line plus a JSON summary so you can reproduce conversions later.

Script source: <a href="{{ site.baseurl }}/episodes/scripts/04-martini-complex.py">04-martini-complex.py</a>

<!-- sync-from: docs/episodes/notebooks/04-martini-complex.ipynb -->

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/04-martini-complex.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/04-martini-complex.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-martini-complex.py" download>Download script (.py)</a></div></div>

### Exercise

- Run the notebook from your `COURSE_DIR`, ensure `martini_openmm` is installed, and compare the bead count/output topology with the original `protein-ligand` complex. Check that the `.json` summary under `results/05-muestreo-avanzado/martini/` documents the executed command.


### Notebooks

- Martini coarse-graining for the complex protein-ligand system. (<a href="{{ site.baseurl }}/episodes/notebooks/04-martini-complex.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/04-martini-complex.py">script derived from this notebook</a>)

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
```python
#!/usr/bin/env python3
import os
from pathlib import Path

from openmm import unit, app
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
PDB_IN = DATA_DIR / "alanine-dipeptide.pdb"
OUT_DIR = COURSE_DIR / "results" / "05-muestreo-avanzado" / "simple"
OUT_DIR.mkdir(parents=True, exist_ok=True)

pdb = PDBFile(str(PDB_IN))
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
)

# Constrain the distance between two atoms as a simple example.
force = mm.CustomBondForce("0.5*k*(r-r0)^2")
force.addPerBondParameter("k")
force.addPerBondParameter("r0")
force.addBond(0, 1, [500.0 * unit.kilojoule_per_mole / unit.nanometer**2, 0.25 * unit.nanometer])
system.addForce(force)

integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy(maxIterations=200)
simulation.step(2000)

print("Simulation with restraint finished. Output dir:", OUT_DIR)
```

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">05-muestreo-avanzado_simple.py</a>

### Exercise

- Run the simple workflow with spherical restraints and record the value of $\langle \|\mathbf{g}\|^2 \rangle$ before and after annealing.
- Compare the average potential energy with and without CustomExternalForce.

### Key points

- Keep reports every 10 steps to monitor integrator drift.
- Adjust $\beta$ to preserve numerical stability in small systems.

### Notebooks and scripts

- This notebook executes the advanced sampling routine for alanine (restraints, custom forces, reporting energies) and tracks gradient norms. (<a href="{{ site.baseurl }}/episodes/notebooks/04-muestreo-avanzado_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/05-muestreo-avanzado.py -->
```python
#!/usr/bin/env python3
import argparse
import os
import sys
import time
from pathlib import Path

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import app, unit, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller, StateDataReporter, DCDReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "complex"
DEFAULT_PROTEIN = DATA_DIR / "protein.pdb"
DEFAULT_LIGAND = DATA_DIR / "ligand1.mol"


def get_platform():
    speed = 0
    platform = None
    for i in range(openmm.Platform.getNumPlatforms()):
        candidate = openmm.Platform.getPlatform(i)
        if candidate.getSpeed() > speed:
            platform = candidate
            speed = candidate.getSpeed()
    print("Using platform", platform.getName())
    if platform.getName() in {"CUDA", "OpenCL"}:
        platform.setPropertyDefaultValue("Precision", "mixed")
        print("Set precision for platform", platform.getName(), "to mixed")
    return platform


def main() -> None:
    parser = argparse.ArgumentParser(description="Simulate protein-ligand complex with optional solvation")
    parser.add_argument("-p", "--protein", default=str(DEFAULT_PROTEIN), help="Protein PDB file")
    parser.add_argument("-l", "--ligand", default=str(DEFAULT_LIGAND), help="Ligand MOL file")
    parser.add_argument("-o", "--output", default="solvated", help="Base name for output files")
    parser.add_argument("-s", "--steps", type=int, default=5000, help="Number of steps")
    parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps)")
    parser.add_argument("-f", "--friction-coeff", type=float, default=1.0, help="Friction coefficient (ps)")
    parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
    parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
    parser.add_argument("--solvate", action="store_true", help="Add solvent box")
    parser.add_argument("--padding", type=float, default=10.0, help="Padding for solvent box (A)")
    parser.add_argument("--water-model", default="tip3p", choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"], help="Water model")
    parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
    parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
    parser.add_argument("--ionic-strength", type=float, default=0.0, help="Ionic strength (M)")
    parser.add_argument("--no-neutralize", action="store_true", help="Don't neutralize")
    parser.add_argument("-e", "--equilibration-steps", type=int, default=200, help="Equilibration steps")
    args = parser.parse_args()

    t0 = time.time()
    out_dir = COURSE_DIR / "results" / "05-muestreo-avanzado" / "complex"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_base = str(out_dir / args.output)
    output_complex = output_base + "_complex.pdb"
    output_min = output_base + "_minimised.pdb"
    output_traj = output_base + "_traj.dcd"

    print("Reading ligand")
    ligand_mol = Molecule.from_file(args.ligand)

    print("Preparing system")
    forcefield_kwargs = {
        "constraints": app.HBonds,
        "rigidWater": True,
        "removeCMMotion": False,
        "hydrogenMass": 4 * unit.amu,
    }
    system_generator = SystemGenerator(
        forcefields=["amber/ff14SB.xml", "amber/tip3p_standard.xml"],
        small_molecule_forcefield="gaff-2.11",
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs,
    )

    print("Reading protein")
    protein_pdb = PDBFile(args.protein)

    print("Preparing complex")
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
    lig_top = ligand_mol.to_topology()
    modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())

    if args.solvate:
        print("Adding solvent")
        modeller.addSolvent(
            system_generator.forcefield,
            model=args.water_model,
            padding=args.padding * unit.angstroms,
            positiveIon=args.positive_ion,
            negativeIon=args.negative_ion,
            ionicStrength=args.ionic_strength * unit.molar,
            neutralize=not args.no_neutralize,
        )

    with open(output_complex, "w") as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
    step_size = args.step_size * unit.picoseconds
    friction = args.friction_coeff / unit.picosecond
    temperature = args.temperature * unit.kelvin
    duration = (step_size * args.steps).value_in_unit(unit.nanoseconds)

    if system.usesPeriodicBoundaryConditions():
        system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))

    integrator = LangevinIntegrator(temperature, friction, step_size)
    platform = get_platform()

    simulation = Simulation(modeller.topology, system, integrator, platform=platform)
    simulation.context.setPositions(modeller.positions)

    print("Minimising ...")
    simulation.minimizeEnergy()

    with open(output_min, "w") as outfile:
        PDBFile.writeFile(
            modeller.topology,
            simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
            file=outfile,
            keepIds=True,
        )

    simulation.context.setVelocitiesToTemperature(temperature)
    print("Equilibrating ...")
    simulation.step(args.equilibration_steps)

    simulation.reporters.append(DCDReporter(output_traj, args.interval, enforcePeriodicBox=True))
    simulation.reporters.append(StateDataReporter(sys.stdout, args.interval * 5, step=True, potentialEnergy=True, temperature=True))

    print("Starting simulation with", args.steps, "steps ...")
    t1 = time.time()
    simulation.step(args.steps)
    t2 = time.time()
    print("Simulation complete in", round((t2 - t1) / 60, 3), "mins")
    print("Simulation time was", round(duration, 3), "ns")
    print("Total wall clock time was", round((t2 - t0) / 60, 3), "mins")


if __name__ == "__main__":
    main()
```

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

Script source: <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">05-muestreo-avanzado.py</a>

### Exercise

- Run with `--solvate` and `--padding 12`, recording solvent density and Coulomb energy.
- Try another water model and analyze the variation in $\Delta G$ estimated from partition differences.

### Key points

- Periodic solvation stabilizes $\Delta G$ on long time scales and improves annealing convergence.
- Compare energy distributions to validate extended sampling.

### Notebooks and scripts

- This notebook covers the complex system with solvation options, barostat control, and energy/force logging to validate enhanced sampling workflows. (<a href="{{ site.baseurl }}/episodes/notebooks/04-muestreo-avanzado.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/05-muestreo-avanzado.py">script</a>)

## References

- [^remd_review]: Sugita, Y. and Okamoto, Y., "Replica-exchange molecular dynamics method for protein folding", *Chem. Phys. Lett.* **2000**, 314, 141–151, and review `PMC6484850 <https://pmc.ncbi.nlm.nih.gov/articles/PMC6484850/>`__.
- [^gamd_paper]: Miao, Y.; Feher, V. A.; McCammon, J. A., "Gaussian accelerated molecular dynamics: Unconstrained enhanced sampling and free energy calculation", *J. Chem. Theory Comput.* **2015**, 11, 3584–3595.
- [^martini_openmm]: MacCallum, J. L.; Tieleman, D. P. et al., "Martini implementation in OpenMM," *Biophys. J.* **2023**, and the open-source repository at https://github.com/maccallumlab/martini_openmm/tree/martini3; the OpenMM adaptation covers both Martini 2 and 3 force fields and leverages OpenMM’s extensible API to process GROMACS topologies with the same flexibility as the native package.

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/05-analisis-trayectorias/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/06-pyemma/">Next</a>
</div>
