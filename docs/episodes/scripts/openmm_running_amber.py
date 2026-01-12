#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from sys import stdout

from openmm import app, unit
import openmm as mm
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, Simulation, DCDReporter, StateDataReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "openmm-examples"
OUT_DIR = COURSE_DIR / "results" / "03-simulaciones-clasicas" / "openmm-amber"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: AMBER input example")
    parser.add_argument(
        "-c",
        "--inpcrd",
        default=str(DATA_DIR / "input.inpcrd"),
        help="AMBER inpcrd file",
    )
    parser.add_argument(
        "-p",
        "--prmtop",
        default=str(DATA_DIR / "input.prmtop"),
        help="AMBER prmtop file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output.dcd",
        help="Output DCD filename",
    )
    parser.add_argument("--steps", type=int, default=10000, help="Number of steps")
    parser.add_argument("--interval", type=int, default=1000, help="Reporting interval")
    args = parser.parse_args()

    inpcrd_path = Path(args.inpcrd)
    prmtop_path = Path(args.prmtop)
    ensure_file(inpcrd_path, "AMBER inpcrd")
    ensure_file(prmtop_path, "AMBER prmtop")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    inpcrd = AmberInpcrdFile(str(inpcrd_path))
    prmtop = AmberPrmtopFile(str(prmtop_path), periodicBoxVectors=inpcrd.boxVectors)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=app.HBonds,
    )
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
    )

    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
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
