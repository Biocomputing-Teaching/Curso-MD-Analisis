---
layout: default
title: Data
permalink: /data/
---
<!-- sync-from: scripts/generate_example_dcd.py -->

### Download materials

Course data lives in <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data">Course-MD-Data</a>. Download it to `COURSE_DIR/data` by following the instructions in <a href="{{ site.baseurl }}/setup/">Setup</a>.

- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/alanine-dipeptide.pdb">Download alanine-dipeptide.pdb</a>: minimal structure for OpenMM testing.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/alanine-dipeptide-multi.pdb">Download alanine-dipeptide-multi.pdb</a>: multi-model PDB for trajectory analysis without DCD.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/example_rmsd.csv">Download example_rmsd.csv</a>: sample RMSD series for plots.
To get the DCD, generate the trajectory locally (see <a href="{{ site.baseurl }}/reference/">Reference</a>).

### Protein-ligand complex data (OpenMM)

- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/protein_orig.pdb">Download protein_orig.pdb</a>: protein with ligand (for preparation).
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/protein.pdb">Download protein.pdb</a>: protein prepared for simulation.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.mol">Download ligand1.mol</a>: ligand in MOL format.
- <a href="https://github.com/Biocomputing-Teaching/Course-MD-Data/raw/main/complex/ligand1.sdf">Download ligand1.sdf</a>: ligand in SDF format.

### OpenMM data (python-examples)

The official OpenMM examples use input files you can download from the OpenMM repo.
Store them in `COURSE_DIR/data/openmm-examples` for the Episode 2 scripts (Running simulations).

```bash
mkdir -p "$COURSE_DIR/data/openmm-examples"
base="https://raw.githubusercontent.com/openmm/openmm/master/examples/python-examples"
for file in \
  input.pdb input.inpcrd input.prmtop input.gro input.top \
  ala_ala_ala.pdb ala_ala_ala.psf charmm22.rtf charmm22.par \
  amoeba_solvated_phenol.xyz amoeba_phenol.prm amoebabio18.prm; do
  curl -L -o "$COURSE_DIR/data/openmm-examples/$file" "$base/$file"
done
```

Note: the Gromacs example requires access to force field files (`includeDir`). If you do not have them installed, use only the PDB/AMBER/CHARMM/Tinker examples.

### Full DCD generator code

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

Script source: <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/scripts/generate_example_dcd.py">generate_example_dcd.py</a>

Script source: <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/scripts/generate_example_dcd.py">generate_example_dcd.py</a>

Script source: <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/scripts/generate_example_dcd.py">generate_example_dcd.py</a>

Script source: <a href="https://github.com/Biocomputing-Teaching/Curso-MD-Analisis/blob/main/scripts/generate_example_dcd.py">generate_example_dcd.py</a>

Note: the command to generate the DCD is in <a href="{{ site.baseurl }}/reference/">Reference</a>.
