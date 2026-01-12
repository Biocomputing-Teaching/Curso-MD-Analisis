#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from sys import stdout

from openmm import app, unit
import openmm as mm
from openmm.app import TinkerFiles, Simulation, DCDReporter, StateDataReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "openmm-examples"
OUT_DIR = COURSE_DIR / "results" / "03-simulaciones-clasicas" / "openmm-tinker"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: Tinker input example")
    parser.add_argument(
        "--xyz",
        default=str(DATA_DIR / "amoeba_solvated_phenol.xyz"),
        help="Tinker XYZ file",
    )
    parser.add_argument(
        "--prm",
        default=str(DATA_DIR / "amoeba_phenol.prm"),
        help="Tinker PRM file",
    )
    parser.add_argument(
        "--prm-bio",
        default=str(DATA_DIR / "amoebabio18.prm"),
        help="Tinker PRM (bio) file",
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

    xyz_path = Path(args.xyz)
    prm_path = Path(args.prm)
    prm_bio_path = Path(args.prm_bio)
    ensure_file(xyz_path, "XYZ input")
    ensure_file(prm_path, "PRM input")
    ensure_file(prm_bio_path, "PRM bio input")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tinker = TinkerFiles(str(xyz_path), [str(prm_path), str(prm_bio_path)])
    system = tinker.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.7 * unit.nanometer,
        vdwCutoff=0.9 * unit.nanometer,
    )
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.001 * unit.picoseconds
    )

    simulation = Simulation(tinker.topology, system, integrator)
    simulation.context.setPositions(tinker.positions)
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
