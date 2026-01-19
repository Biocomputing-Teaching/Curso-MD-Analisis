from pathlib import Path
import os

COURSE_DIR = Path(os.environ.get('COURSE_DIR', str(Path.home() / 'Concepcion26'))).expanduser()
COMPLEX_DIR = COURSE_DIR / 'data' / 'complex'
MARTINI_DIR = COURSE_DIR / 'results' / '05-muestreo-avanzado' / 'martini'
MARTINI_DIR.mkdir(parents=True, exist_ok=True)
MARTINI_TOP = MARTINI_DIR / 'complex-martini.top'
MARTINI_GRO = MARTINI_DIR / 'complex-martini.gro'

print('Course directory:', COURSE_DIR)
print('Complex directory contents:')
for path in sorted(COMPLEX_DIR.glob('*')):
    print(' ', path.name)
print('Martini outputs:', MARTINI_TOP.name, MARTINI_GRO.name)

# %%

command = [
    'martini_openmm',
    '--protein', str(COMPLEX_DIR / 'protein.pdb'),
    '--ligand', str(COMPLEX_DIR / 'ligand1.mol'),
    '--output-top', str(MARTINI_TOP),
    '--output-gro', str(MARTINI_GRO),
]
print('Run this command to build the Martini topology:')
print(' ', ' '.join(command))

# %%

import shutil
import subprocess

if shutil.which('martini_openmm'):
    print('martini_openmm available; running conversion now...')
    subprocess.run(command, check=True)
else:
    print('martini_openmm is not installed. Follow https://github.com/maccallumlab/martini_openmm to get it.')

# %%

from openmm import app, unit
from openmm.app import GromacsTopFile, GromacsGroFile, Simulation, DCDReporter
from openmm import LangevinIntegrator

if MARTINI_TOP.exists() and MARTINI_GRO.exists():
    top = GromacsTopFile(str(MARTINI_TOP))
    gro = GromacsGroFile(str(MARTINI_GRO))
    system = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 20 * unit.femtoseconds)
    simulation = Simulation(top.topology, system, integrator)
    simulation.context.setPositions(gro.getPositions())
    simulation.minimizeEnergy()
    simulation.reporters.append(DCDReporter(str(MARTINI_DIR / 'martini_bg.dcd'), 100, enforcePeriodicBox=True))
    simulation.step(500)
    print('Martini simulation finished; DCD saved to', MARTINI_DIR / 'martini_bg.dcd')
else:
    print('Missing Martini files. Run the conversion first, then rerun this cell.')
