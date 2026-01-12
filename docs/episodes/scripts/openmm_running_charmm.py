#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from sys import stdout

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, CharmmPsfFile, CharmmParameterSet, Simulation, DCDReporter, StateDataReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "openmm-examples"
OUT_DIR = COURSE_DIR / "results" / "03-simulaciones-clasicas" / "openmm-charmm"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: CHARMM input example")
    parser.add_argument(
        "--psf",
        default=str(DATA_DIR / "ala_ala_ala.psf"),
        help="CHARMM PSF file",
    )
    parser.add_argument(
        "--pdb",
        default=str(DATA_DIR / "ala_ala_ala.pdb"),
        help="CHARMM PDB file",
    )
    parser.add_argument(
        "--rtf",
        default=str(DATA_DIR / "charmm22.rtf"),
        help="CHARMM RTF file",
    )
    parser.add_argument(
        "--prm",
        default=str(DATA_DIR / "charmm22.par"),
        help="CHARMM PAR file",
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

    psf_path = Path(args.psf)
    pdb_path = Path(args.pdb)
    rtf_path = Path(args.rtf)
    prm_path = Path(args.prm)
    ensure_file(psf_path, "PSF input")
    ensure_file(pdb_path, "PDB input")
    ensure_file(rtf_path, "RTF input")
    ensure_file(prm_path, "PAR input")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    psf = CharmmPsfFile(str(psf_path))
    pdb = PDBFile(str(pdb_path))
    params = CharmmParameterSet(str(rtf_path), str(prm_path))
    system = psf.createSystem(
        params,
        nonbondedMethod=app.NoCutoff,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=app.HBonds,
    )
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
    )

    simulation = Simulation(psf.topology, system, integrator)
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
