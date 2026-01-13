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
