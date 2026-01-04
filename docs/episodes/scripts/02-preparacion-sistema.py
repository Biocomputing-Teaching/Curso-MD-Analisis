#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

from pdbfixer import PDBFixer
from openmmforcefields.generators import SystemGenerator
from openmm import unit
from openmm.app import PDBFile, Simulation
from openmm import LangevinIntegrator

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "complex"
DEFAULT_INPUT = DATA_DIR / "protein_orig.pdb"


def prepare_protein(pdb_in: Path, output_dir: Path, output_base: str) -> None:
    print("Processing", pdb_in)
    fixer = PDBFixer(filename=str(pdb_in))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    print("Residues:", fixer.missingResidues)
    print("Atoms:", fixer.missingAtoms)
    print("Terminals:", fixer.missingTerminals)
    print("Non-standard:", fixer.nonstandardResidues)

    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    fixer.removeHeterogens(False)

    output_dir.mkdir(parents=True, exist_ok=True)
    fixed_path = output_dir / f"{output_base}_fixed.pdb"
    min_path = output_dir / f"{output_base}_minimised.pdb"

    with fixed_path.open("w") as outfile:
        PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile, keepIds=True)

    system_generator = SystemGenerator(forcefields=["amber/ff14SB.xml"])
    system = system_generator.create_system(fixer.topology)
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = Simulation(fixer.topology, system, integrator)
    simulation.context.setPositions(fixer.positions)
    print("Minimising")
    simulation.minimizeEnergy()

    with min_path.open("w") as outfile:
        PDBFile.writeFile(
            fixer.topology,
            simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(),
            file=outfile,
            keepIds=True,
        )

    print("Written", fixed_path)
    print("Written", min_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare protein with PDBFixer and OpenMM")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input PDB file")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--output-base", default="protein_prepared", help="Output base name")
    args = parser.parse_args()

    out_dir = Path(args.output_dir) if args.output_dir else (COURSE_DIR / "results" / "02-preparacion-sistema" / "complex")
    out_dir.mkdir(parents=True, exist_ok=True)
    prepare_protein(Path(args.input), out_dir, args.output_base)


if __name__ == "__main__":
    main()
