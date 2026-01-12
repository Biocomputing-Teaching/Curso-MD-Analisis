#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from sys import stdout

from openmm import app, unit
import openmm as mm
from openmm.app import GromacsGroFile, GromacsTopFile, Simulation, DCDReporter, StateDataReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "openmm-examples"
OUT_DIR = COURSE_DIR / "results" / "03-simulaciones-clasicas" / "openmm-gromacs"
DEFAULT_GROMACS_TOP = os.environ.get("GROMACS_TOP", "/usr/local/gromacs/share/gromacs/top")


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: Gromacs input example")
    parser.add_argument(
        "-g",
        "--gro",
        default=str(DATA_DIR / "input.gro"),
        help="Gromacs GRO file",
    )
    parser.add_argument(
        "-t",
        "--top",
        default=str(DATA_DIR / "input.top"),
        help="Gromacs TOP file",
    )
    parser.add_argument(
        "--include-dir",
        default=DEFAULT_GROMACS_TOP,
        help="Directory with Gromacs force field files",
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

    gro_path = Path(args.gro)
    top_path = Path(args.top)
    include_dir = Path(args.include_dir)

    ensure_file(gro_path, "GRO input")
    ensure_file(top_path, "TOP input")
    if not include_dir.exists():
        raise SystemExit(f"Missing Gromacs includeDir: {include_dir}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    gro = GromacsGroFile(str(gro_path))
    top = GromacsTopFile(
        str(top_path),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(include_dir),
    )
    system = top.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=app.HBonds,
    )
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
    )

    simulation = Simulation(top.topology, system, integrator)
    simulation.context.setPositions(gro.positions)
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
