#!/usr/bin/env python3
import os
from pathlib import Path

from openmm import unit
from openmm.app import PDBFile, Modeller, ForceField

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
PDB_IN = DATA_DIR / "alanine-dipeptide.pdb"
OUT_DIR = COURSE_DIR / "results" / "02-preparacion-sistema" / "simple"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT = OUT_DIR / "alanine_solvated.pdb"

pdb = PDBFile(str(PDB_IN))
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
modeller.addSolvent(
    forcefield,
    model="tip3p",
    padding=1.0 * unit.nanometer,
)

with OUTPUT.open("w") as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

print("Written", OUTPUT)

# %%
