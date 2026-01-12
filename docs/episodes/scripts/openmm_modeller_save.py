#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, ForceField, Modeller, Simulation

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
OUT_DIR = COURSE_DIR / "results" / "02-preparacion-sistema" / "openmm-modeller"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: Modeller save results")
    parser.add_argument(
        "-p",
        "--pdb",
        default=str(DATA_DIR / "alanine-dipeptide.pdb"),
        help="Input PDB file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output.pdb",
        help="Output PDB filename",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=1.0,
        help="Solvent padding (nm)",
    )
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    ensure_file(pdb_path, "PDB input")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading...")
    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber99sb.xml", "tip3p.xml")
    modeller = Modeller(pdb.topology, pdb.positions)

    print("Adding hydrogens...")
    modeller.addHydrogens(forcefield)
    print("Adding solvent...")
    modeller.addSolvent(forcefield, model="tip3p", padding=args.padding * unit.nanometer)

    print("Minimizing...")
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=100)

    output_path = OUT_DIR / args.output
    print("Saving...", output_path)
    positions = simulation.context.getState(positions=True).getPositions()
    with open(output_path, "w") as handle:
        PDBFile.writeFile(simulation.topology, positions, handle)

    print("Done")


if __name__ == "__main__":
    main()
