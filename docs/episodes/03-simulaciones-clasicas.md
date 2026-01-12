---
layout: default
title: Episodio 2 - Ejecución de simulaciones
permalink: /episodes/03-simulaciones-clasicas/
---

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Siguiente</a>
</div>

<!-- toc:start -->
## Tabla de contenidos
- [Duración](#duracion)
- [Objetivos](#objetivos)
- [Contenido](#contenido)
- [Ejemplos OpenMM (guía oficial)](#ejemplos-openmm-guia-oficial)
- [Parte simple](#parte-simple)
- [Parte compleja](#parte-compleja)
<!-- toc:end -->

## Duración

- **Sesión:** 60 min
- **Ejercicios:** 45 min

## Objetivos

- Replicar el flujo base de OpenMM para ejecutar simulaciones.
- Ejecutar ejemplos con archivos PDB, AMBER, Gromacs, CHARMM y Tinker.
- Simular alanina y un complejo proteína-ligando con el pipeline del curso.
- Generar trayectorias DCD y reportes de energía.

## Contenido

- Ejemplo base de OpenMM con PDB (simulatePdb).
- Ejemplos de entrada AMBER, Gromacs, CHARMM y Tinker.
- Integradores, barostatos y reportes básicos.
- Ejemplos del curso (alanina y complejo proteína-ligando).

## Ejemplos OpenMM (guía oficial)

Antes de ejecutar estos scripts, descarga los archivos de ejemplo de OpenMM a `COURSE_DIR/data/openmm-examples` (ver <a href="{{ site.baseurl }}/data/">Datos</a>).

### A First Example (PDB)

<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py" download>Descargar script (.py)</a></div></div>

Script del curso: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_pdb.py">openmm_running_pdb.py</a>
 (por defecto usa `alanine-dipeptide.pdb`).

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

Script del curso: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_amber.py">openmm_running_amber.py</a>

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

Script del curso: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_gromacs.py">openmm_running_gromacs.py</a>

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

Script del curso: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_charmm.py">openmm_running_charmm.py</a>

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

Script del curso: <a href="{{ site.baseurl }}/episodes/scripts/openmm_running_tinker.py">openmm_running_tinker.py</a>

### OpenMM-Setup (generador de scripts)

```bash
conda install -c conda-forge openmm-setup
openmm-setup
```

Este asistente genera scripts equivalentes a los ejemplos anteriores y permite validar opciones sin editar código.

## Parte simple

### Demo guiada

<!-- sync-from: docs/episodes/scripts/03-simulaciones-clasicas_simple.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas_simple.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Cambiar `--steps` y `--interval` en el sistema simple.
- Verificar los archivos de salida de alanina.

### Puntos clave

- Los conceptos básicos se prueban con alanina.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb">03-simulaciones-clasicas_simple.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas_simple.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/03-simulaciones-clasicas_simple.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas_simple.py">03-simulaciones-clasicas_simple.py</a>

## Parte compleja

### Demo guiada

<!-- sync-from: docs/episodes/scripts/03-simulaciones-clasicas.py -->
<div class="notebook-embed"><iframe src="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas.html" loading="lazy"></iframe><div class="notebook-links"><a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb" download>Descargar notebook</a> | <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py" download>Descargar script (.py)</a></div></div>

### Ejercicio

- Cambiar `--steps` y `--interval` en el sistema complejo.
- Verificar los archivos de salida del complejo.
- Probar con otro `--output` base en el complejo.

### Puntos clave

- El complejo introduce más átomos y mayor coste.

### Notebooks y scripts

- <a href="{{ site.baseurl }}/episodes/notebooks/03-simulaciones-clasicas.ipynb">03-simulaciones-clasicas.ipynb</a> (<a href="{{ site.baseurl }}/episodes/notebooks/rendered/03-simulaciones-clasicas.html">HTML</a> | <a href="https://nbviewer.org/url/biocomputing-teaching.github.io/Curso-MD-Analisis/episodes/notebooks/03-simulaciones-clasicas.ipynb">nbviewer</a>)
- <a href="{{ site.baseurl }}/episodes/scripts/03-simulaciones-clasicas.py">03-simulaciones-clasicas.py</a>

<div class="episode-nav">
  <a href="{{ site.baseurl }}/episodes/01-introduccion/">Anterior</a>
  <a href="{{ site.baseurl }}/episodes/">Todos los episodios</a>
  <a href="{{ site.baseurl }}/episodes/02-preparacion-sistema/">Siguiente</a>
</div>
