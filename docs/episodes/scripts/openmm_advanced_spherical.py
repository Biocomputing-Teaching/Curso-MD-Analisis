#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
OUT_DIR = COURSE_DIR / "results" / "05-muestreo-avanzado" / "openmm-spherical"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: spherical container example")
    parser.add_argument(
        "-p",
        "--pdb",
        default=str(DATA_DIR / "alanine-dipeptide.pdb"),
        help="Input PDB file",
    )
    parser.add_argument("--steps", type=int, default=1000, help="Number of steps")
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    ensure_file(pdb_path, "PDB input")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=None,
    )
    force = mm.CustomExternalForce("100*max(0, r-2)^2; r=sqrt(x*x+y*y+z*z)")
    system.addForce(force)
    for i in range(system.getNumParticles()):
        force.addParticle(i, [])

    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 91 / unit.picosecond, 0.004 * unit.picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.step(args.steps)

    print("Spherical container run complete")


if __name__ == "__main__":
    main()
