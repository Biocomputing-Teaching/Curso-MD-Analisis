
    # Visualizacion interactiva (requiere nglview)
    import os
    from pathlib import Path

    COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
    PROTEIN_PDB = COURSE_DIR / "data" / "complex" / "protein.pdb"
    LIGAND_SDF = COURSE_DIR / "data" / "complex" / "ligand1.sdf"

    try:
        import nglview as nv
    except ImportError as exc:
        raise SystemExit(
            "nglview is required for visualization. Install with:
"
            "  conda install -c conda-forge nglview"
        ) from exc

    view = nv.show_file(str(PROTEIN_PDB))
    try:
        view.add_component(str(LIGAND_SDF))
    except Exception:
        pass
    view

# %%


import os
from collections import Counter
from pathlib import Path

import openmm as mm
import openmm.app as app
from openff.toolkit import Molecule
from openmm.app import ForceField, Modeller, PDBFile
from openmmforcefields.generators import SystemGenerator

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "complex"
PROTEIN_PDB = DATA_DIR / "protein.pdb"
LIGAND_SDF = DATA_DIR / "ligand1.sdf"
OUT_DIR = COURSE_DIR / "results" / "01-introduccion" / "complex"
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
print("Protein PDB:", PROTEIN_PDB)
print("Ligand SDF:", LIGAND_SDF)
print("Output dir:", OUT_DIR)

protein = PDBFile(str(PROTEIN_PDB))
ligand = Molecule.from_file(str(LIGAND_SDF))

summarize_topology("Protein (input PDB)", protein.topology)
print("Ligand atoms:", ligand.n_atoms)

protein_ff = resolve_ff19sb()
protein_forcefield = ForceField(protein_ff)

protein_modeller = Modeller(protein.topology, protein.positions)
protein_modeller.addHydrogens(protein_forcefield)
summarize_topology("Protein (with H)", protein_modeller.topology)

ligand_top = ligand.to_topology()
modeller = Modeller(protein_modeller.topology, protein_modeller.positions)
modeller.add(ligand_top.to_openmm(), ligand_top.get_positions().to_openmm())

summarize_topology("Protein + ligand (combined)", modeller.topology)

system_generator = SystemGenerator(
    forcefields=[protein_ff],
    small_molecule_forcefield="gaff-2.11",
    molecules=[ligand],
)
system = system_generator.create_system(modeller.topology, molecules=ligand)

write_amber_files(modeller.topology, system, modeller.positions, OUT_DIR, "protein_ligand")
