---
layout: default
title: Datos
permalink: /data/
---

### Descargar materiales

- [Descargar `alanine-dipeptide.pdb`]({{ site.baseurl }}/data/alanine-dipeptide.pdb): estructura mínima para pruebas con OpenMM.
- [Descargar `alanine-dipeptide-multi.pdb`]({{ site.baseurl }}/data/alanine-dipeptide-multi.pdb): PDB multi-model para análisis de trayectorias sin DCD.
- [Descargar `example_rmsd.csv`]({{ site.baseurl }}/data/example_rmsd.csv): serie de RMSD de ejemplo para gráficas.
Para obtener el DCD, genera la trayectoria localmente (ver [Referencia]({{ site.baseurl }}/reference/)).

### Código completo del generador DCD

```python
from openmm import app, unit
import openmm as mm

pdb = app.PDBFile('docs/data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds
)

integrator = mm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds
)

simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy(maxIterations=100)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

simulation.reporters.append(app.DCDReporter('docs/data/alanine-dipeptide.dcd', 10))
simulation.step(200)

print('Wrote docs/data/alanine-dipeptide.dcd')
```

Nota: el comando para generar el DCD está en [Referencia]({{ site.baseurl }}/reference/).
