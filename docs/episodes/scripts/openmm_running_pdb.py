#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from sys import stdout

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
OUT_DIR = COURSE_DIR / "results" / "03-simulaciones-clasicas" / "openmm-pdb"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: first PDB example")
    parser.add_argument(
        "-p",
        "--pdb",
        default=str(DATA_DIR / "alanine-dipeptide.pdb"),
        help="Input PDB file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output.dcd",
        help="Output DCD filename",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=10000,
        help="Number of simulation steps",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=1000,
        help="Reporting interval (steps)",
    )
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    ensure_file(pdb_path, "PDB input")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber19-all.xml", "amber19/tip3pfb.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=app.HBonds,
    )
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
    )

    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()

    dcd_path = OUT_DIR / args.output
    simulation.reporters.append(DCDReporter(str(dcd_path), args.interval))
    simulation.reporters.append(
        StateDataReporter(stdout, args.interval, step=True, potentialEnergy=True, temperature=True)
    )
    simulation.step(args.steps)

    print("Written", dcd_path)


if __name__ == "__main__":
    main()
