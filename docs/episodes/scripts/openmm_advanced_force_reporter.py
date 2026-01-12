#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
OUT_DIR = COURSE_DIR / "results" / "05-muestreo-avanzado" / "openmm-force-reporter"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


class ForceReporter:
    def __init__(self, file, report_interval):
        self._out = open(file, "w")
        self._report_interval = report_interval

    def __del__(self):
        try:
            self._out.close()
        except Exception:
            pass

    def describeNextReport(self, simulation):
        steps = self._report_interval - simulation.currentStep % self._report_interval
        return {"steps": steps, "periodic": None, "include": ["forces"]}

    def report(self, simulation, state):
        forces = state.getForces().value_in_unit(unit.kilojoules / unit.mole / unit.nanometer)
        for f in forces:
            self._out.write(f"{f[0]} {f[1]} {f[2]}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: custom force reporter")
    parser.add_argument(
        "-p",
        "--pdb",
        default=str(DATA_DIR / "alanine-dipeptide.pdb"),
        help="Input PDB file",
    )
    parser.add_argument("--steps", type=int, default=100, help="Number of steps")
    parser.add_argument("--interval", type=int, default=10, help="Reporting interval")
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    ensure_file(pdb_path, "PDB input")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds)

    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()

    forces_path = OUT_DIR / "forces.txt"
    simulation.reporters.append(ForceReporter(str(forces_path), args.interval))
    simulation.step(args.steps)

    print("Written", forces_path)


if __name__ == "__main__":
    main()
