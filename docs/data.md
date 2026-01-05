---
layout: default
title: Datos
permalink: /data/
---
<!-- sync-from: scripts/generate_example_dcd.py -->

### Descargar materiales

Los datos del curso viven en <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data">Course-MD-Data</a>. Descárgalos a `COURSE_DIR/data` siguiendo las instrucciones de <a href="{{ site.baseurl }}/setup/">Preparación</a>.

- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/alanine-dipeptide.pdb">Descargar alanine-dipeptide.pdb</a>: estructura mínima para pruebas con OpenMM.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/alanine-dipeptide-multi.pdb">Descargar alanine-dipeptide-multi.pdb</a>: PDB multi-model para análisis de trayectorias sin DCD.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/example_rmsd.csv">Descargar example_rmsd.csv</a>: serie de RMSD de ejemplo para gráficas.
Para obtener el DCD, genera la trayectoria localmente (ver <a href="{{ site.baseurl }}/reference/">Referencia</a>).

### Datos para complejo proteína-ligando (OpenMM)

- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/protein_orig.pdb">Descargar protein_orig.pdb</a>: proteína con ligando (para preparación).
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/protein.pdb">Descargar protein.pdb</a>: proteína preparada para simulación.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.mol">Descargar ligand1.mol</a>: ligando en formato MOL.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.sdf">Descargar ligand1.sdf</a>: ligando en formato SDF.

### Código completo del generador DCD

```python
from openmm import app, unit
import openmm as mm

pdb = app.PDBFile('docs/data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds
)

integrator = mm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds
)

simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy(maxIterations=100)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

simulation.reporters.append(app.DCDReporter('docs/data/alanine-dipeptide.dcd', 10))
simulation.step(200)

print('Wrote docs/data/alanine-dipeptide.dcd')
```

Fuente del script: <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/scripts/generate_example_dcd.py">generate_example_dcd.py</a>

Nota: el comando para generar el DCD está en <a href="{{ site.baseurl }}/reference/">Referencia</a>.
