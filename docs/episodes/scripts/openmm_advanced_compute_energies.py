#!/usr/bin/env python3
import argparse
import os
import random
from pathlib import Path

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
OUT_DIR = COURSE_DIR / "results" / "05-muestreo-avanzado" / "openmm-energies"


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def generate_structures(base_pdb: Path, output_dir: Path, count: int) -> list[Path]:
    pdb = PDBFile(str(base_pdb))
    positions = pdb.positions
    output_dir.mkdir(parents=True, exist_ok=True)

    random.seed(42)
    paths = []
    for i in range(count):
        jitter = [
            (pos[0] + random.uniform(-0.005, 0.005) * unit.nanometer,
             pos[1] + random.uniform(-0.005, 0.005) * unit.nanometer,
             pos[2] + random.uniform(-0.005, 0.005) * unit.nanometer)
            for pos in positions
        ]
        out_path = output_dir / f"structure_{i+1:02d}.pdb"
        with open(out_path, "w") as handle:
            PDBFile.writeFile(pdb.topology, jitter, handle)
        paths.append(out_path)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="OpenMM guide: compute energies for multiple structures")
    parser.add_argument(
        "-p",
        "--pdb",
        default=str(DATA_DIR / "alanine-dipeptide.pdb"),
        help="Base PDB file",
    )
    parser.add_argument(
        "-d",
        "--structures-dir",
        default=str(DATA_DIR / "energy-structures"),
        help="Directory with PDB structures",
    )
    parser.add_argument("--count", type=int, default=5, help="Structures to generate if directory is empty")
    args = parser.parse_args()

    base_pdb = Path(args.pdb)
    structures_dir = Path(args.structures_dir)
    ensure_file(base_pdb, "Base PDB input")

    if structures_dir.exists():
        pdb_files = sorted(structures_dir.glob("*.pdb"))
    else:
        pdb_files = []

    if not pdb_files:
        pdb_files = generate_structures(base_pdb, structures_dir, args.count)
        print(f"Generated {len(pdb_files)} structures in {structures_dir}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        PDBFile(str(base_pdb)).topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )
    integrator = mm.VerletIntegrator(0.002 * unit.picoseconds)
    simulation = Simulation(PDBFile(str(base_pdb)).topology, system, integrator)

    for pdb_path in pdb_files:
        pdb = PDBFile(str(pdb_path))
        simulation.context.setPositions(pdb.positions)
        state = simulation.context.getState(energy=True)
        print(pdb_path.name, state.getPotentialEnergy())


if __name__ == "__main__":
    main()
