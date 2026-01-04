#!/usr/bin/env python3
import os
from pathlib import Path

from openmm import unit, app
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation, StateDataReporter, DCDReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
PDB_IN = DATA_DIR / "alanine-dipeptide.pdb"
OUT_DIR = COURSE_DIR / "results" / "03-simulaciones-clasicas" / "simple"
OUT_DIR.mkdir(parents=True, exist_ok=True)

pdb = PDBFile(str(PDB_IN))
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
)

integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy(maxIterations=200)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

simulation.reporters.append(StateDataReporter(str(OUT_DIR / "log.csv"), 1000, step=True, temperature=True, potentialEnergy=True))
simulation.reporters.append(DCDReporter(str(OUT_DIR / "traj.dcd"), 1000))

simulation.step(5000)
print("Written", OUT_DIR / "log.csv", "and", OUT_DIR / "traj.dcd")
