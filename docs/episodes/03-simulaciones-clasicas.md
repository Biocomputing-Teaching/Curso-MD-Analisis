---
layout: default
title: Episode 2 - Running simulations
permalink: /episodes/03-simulaciones-clasicas/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Next</a>
</div>

<!-- toc:start -->
## Table of contents
- [Duration](#duration)
- [Objectives](#objectives)
- [Content](#content)
- [OpenMM examples (official guide)](#openmm-examples-official-guide)
- [Alanine dipeptide](#alanine-dipeptide)
- [Protein-ligand complex](#protein-ligand-complex)
<!-- toc:end -->

## Duration

- **Session:** 60 min
- **Exercises:** 45 min

## Objectives

- Reproduce the base OpenMM workflow to run simulations.
- Run examples with PDB, AMBER, Gromacs, CHARMM, and Tinker files.
- Simulate alanine and a protein-ligand complex with the course pipeline.
- Generate DCD trajectories and energy reports.

## Content

- OpenMM base example with PDB (simulatePdb).
- AMBER, Gromacs, CHARMM, and Tinker input examples.
- Integrators, barostats, and basic reporting.
- Course examples (alanine and protein-ligand complex).

## OpenMM examples (official guide)

Before running these scripts, download the OpenMM example files to `COURSE_DIR/data/openmm-examples` (see <a href="{{ site.baseurl }}/data/">Data</a>).

### A First Example (PDB)

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py" download>Download script (.py)</a></div></div>

Course script: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_pdb.py">openmm_running_pdb.py</a>
 (defaults to `alanine-dipeptide.pdb`).

### Using AMBER Files

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

inpcrd = AmberInpcrdFile('input.inpcrd')
prmtop = AmberPrmtopFile('input.prmtop', periodicBoxVectors=inpcrd.boxVectors)
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
```

Course script: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_amber.py">openmm_running_amber.py</a>

### Using Gromacs Files

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile('input.gro')
top = GromacsTopFile('input.top', periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir='/usr/local/gromacs/share/gromacs/top')
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
```

Course script: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_gromacs.py">openmm_running_gromacs.py</a>

### Using CHARMM Files

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

psf = CharmmPsfFile('input.psf')
pdb = PDBFile('input.pdb')
params = CharmmParameterSet('charmm22.rtf', 'charmm22.prm')
system = psf.createSystem(params, nonbondedMethod=NoCutoff,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
```

Course script: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_charmm.py">openmm_running_charmm.py</a>

### Using Tinker Files

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

tinker = TinkerFiles('amoeba_solvated_phenol.xyz', ['amoeba_phenol.prm', 'amoebabio18.prm'])
system = tinker.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.7*nanometer, vdwCutoff=0.9*nanometer)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(tinker.topology, system, integrator)
simulation.context.setPositions(tinker.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
```

Course script: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_tinker.py">openmm_running_tinker.py</a>

### OpenMM-Setup (script generator)

```bash
conda install -c conda-forge openmm-setup
openmm-setup
```

This assistant generates scripts equivalent to the examples above and lets you validate options without editing code.

## Alanine dipeptide

### Guided demo

<!-- sync-from: docs/episodes/scripts/03-simulaciones-clasicas_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py" download>Download script (.py)</a></div></div>

### Exercise

- Change `--steps` and `--interval` in the simple system.
- Verify alanine output files.

### Key points

- Basic concepts are tested on alanine.

### Notebooks and scripts

- This notebook runs the simple alanine simulation (OpenMM PDB example) and shows how changing `--steps`/`--interval` affects the output. (<a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py">script</a>)

## Protein-ligand complex

### Guided demo

<!-- sync-from: docs/episodes/scripts/03-simulaciones-clasicas.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb" download>Download notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py" download>Download script (.py)</a></div></div>

### Exercise

- Change `--steps` and `--interval` in the complex system.
- Verify complex output files.
- Try a different `--output` base for the complex.

### Key points

- The complex introduces more atoms and higher cost.

### Notebooks and scripts

- This notebook executes the complex protein-ligand simulation (AMBER/Gromacs/CHARMM/Tinker inputs) and emphasizes file selection and reporting. (<a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb">notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py">script</a>)

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Previous</a>
  <a href="{{ site.baseurl }}/episodes/">All episodes</a>
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Next</a>
</div>
