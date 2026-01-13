
# Quick residue scan before visualization (BioPython)
import os
from pathlib import Path

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
PDB_PATH = COURSE_DIR / "data" / "alanine-dipeptide.pdb"

try:
    from Bio.PDB import PDBParser
except ImportError as exc:
    raise SystemExit(
        "Biopython is required for this step. Install with:"
        "  conda install -c conda-forge biopython"
    ) from exc

parser = PDBParser(QUIET=True)
structure = parser.get_structure("alanine", str(PDB_PATH))

residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            residues.append((chain.id, residue.id[1], residue.resname))
    break

print("PDB:", PDB_PATH)
print("Residues:", " - ".join(f"{r[2]}({r[0]}{r[1]})" for r in residues))

# %%


# Interactive visualization (requires nglview)
import os
from pathlib import Path

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
PDB_PATH = COURSE_DIR / "data" / "alanine-dipeptide.pdb"

try:
    import nglview as nv
except ImportError as exc:
    raise SystemExit(
        "nglview is required for visualization. Install with:"
        "  conda install -c conda-forge nglview"
    ) from exc

view = nv.show_file(str(PDB_PATH))
view

# %%


import os
from collections import Counter
from pathlib import Path

import openmm as mm
import openmm.app as app
from openmm.app import ForceField, Modeller, PDBFile
from openmmforcefields.generators import SystemGenerator

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
PDB_PATH = DATA_DIR / "alanine-dipeptide.pdb"
OUT_DIR = COURSE_DIR / "results" / "01-introduccion" / "simple"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WATER_NAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}


def summarize_topology(label, topology):
    residues = list(topology.residues())
    chains = list(topology.chains())
    residue_counts = Counter(res.name for res in residues)
    water_count = sum(1 for res in residues if res.name in WATER_NAMES)

    print(f"{label}:")
    print("  Atoms:", topology.getNumAtoms())
    print("  Residues:", len(residues))
    print("  Chains:", len(chains))
    if water_count:
        print("  Water residues:", water_count)
    for name, count in residue_counts.most_common(10):
        print(f"  Residue {name}: {count}")


def resolve_ff19sb():
    candidates = [
        "amber/ff19SB.xml",
        "amber14/ff19SB.xml",
        "amber14/protein.ff19SB.xml",
    ]
    for ff in candidates:
        try:
            ForceField(ff)
        except Exception:
            continue
        print("Using protein force field:", ff)
        return ff

    data_dir = Path(app.__file__).resolve().parent / "data"
    matches = sorted(p for p in data_dir.rglob("*.xml") if "ff19" in p.name.lower())
    if matches:
        ff_path = str(matches[0])
        print("Using protein force field:", ff_path)
        return ff_path

    raise SystemExit(
        "Could not locate ff19SB. Tried: " + ", ".join(candidates) + "\n"
        "Install with: conda install -c conda-forge openmmforcefields ambertools"
    )


def write_amber_files(topology, system, positions, out_dir, base_name):
    try:
        import parmed as pmd
    except ImportError as exc:
        raise SystemExit(
            "Parmed is required to write Amber files. Install with:\n"
            "  conda install -c conda-forge parmed"
        ) from exc

    structure = pmd.openmm.load_topology(topology, system, positions)
    prmtop_path = out_dir / f"{base_name}.prmtop"
    inpcrd_path = out_dir / f"{base_name}.inpcrd"
    structure.save(str(prmtop_path), overwrite=True)
    structure.save(str(inpcrd_path), overwrite=True)
    print("Written:", prmtop_path)
    print("Written:", inpcrd_path)


print("OpenMM", mm.__version__)
print("PDB:", PDB_PATH)
print("Output dir:", OUT_DIR)

pdb = PDBFile(str(PDB_PATH))
summarize_topology("Alanine dipeptide (input)", pdb.topology)

protein_ff = resolve_ff19sb()
protein_forcefield = ForceField(protein_ff)

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(protein_forcefield)
summarize_topology("Alanine dipeptide (with H)", modeller.topology)

system_generator = SystemGenerator(forcefields=[protein_ff])
system = system_generator.create_system(modeller.topology)

write_amber_files(modeller.topology, system, modeller.positions, OUT_DIR, "alanine-dipeptide")

# %%
