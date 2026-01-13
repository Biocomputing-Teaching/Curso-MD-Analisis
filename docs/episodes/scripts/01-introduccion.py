
# Quick sequence scan before visualization (BioPython)
import os
from pathlib import Path

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
PROTEIN_PDB = COURSE_DIR / "data" / "complex" / "protein.pdb"

try:
    from Bio.PDB import PDBParser
except ImportError as exc:
    raise SystemExit(
        "Biopython is required for this step. Install with:"
        "  conda install -c conda-forge biopython"
    ) from exc

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", str(PROTEIN_PDB))

three_to_one = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H",
}

sequences = {}
for model in structure:
    for chain in model:
        seq = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            seq.append(three_to_one.get(residue.resname, "X"))
        if seq:
            sequences[chain.id] = "".join(seq)
    break

print("Protein PDB:", PROTEIN_PDB)
for chain, seq in sequences.items():
    print(f"Chain {chain}: {len(seq)} residues")
    print(seq)

# %%


# Interactive visualization (requires nglview)
import os
from pathlib import Path

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
PROTEIN_PDB = COURSE_DIR / "data" / "complex" / "protein.pdb"
LIGAND_SDF = COURSE_DIR / "data" / "complex" / "ligand1.sdf"

try:
    import nglview as nv
except ImportError as exc:
    raise SystemExit(
        "nglview is required for visualization. Install with:"
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

# %%


import os
from pathlib import Path
import requests
from Bio.Data import IUPACData
from Bio.PDB import PDBParser
from rdkit import Chem

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
PROTEIN_PDB = COURSE_DIR / "data" / "complex" / "protein.pdb"
LIGAND_SDF = COURSE_DIR / "data" / "complex" / "ligand1.sdf"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", str(PROTEIN_PDB))
chain_sequences = {}
for model in structure:
    for chain in model:
        seq = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            seq.append(IUPACData.protein_letters_3to1.get(residue.resname, "X"))
        chain_sequences[chain.id] = "".join(seq)

print("Chains found:", list(chain_sequences.keys()))

SEQ = next(iter(chain_sequences.values()), "")
if not SEQ:
    print("No sequence available for the protein chains.")
else:
    def query_uniprot(sequence):
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": f"sequence:{sequence}", "format": "json", "size": 1}
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        hits = resp.json().get("results", [])
        return hits[0] if hits else None

    def search_pdb(sequence):
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        query = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "value": sequence,
                    "identityCutoff": 0.9,
                    "evalueCutoff": 1e-5,
                    "sequenceType": "protein",
                    "matrix": "BLOSUM62",
                    "gapsAllowed": True,
                    "target": "pdb_protein_sequence"
                }
            },
            "return_type": "entry"
        }
        resp = requests.post(url, json=query, timeout=10)
        resp.raise_for_status()
        result_set = resp.json().get("result_set", [])
        return result_set[0]["identifier"] if result_set else None

    uniprot_hit = None
    pdb_hit = None
    try:
        uniprot_hit = query_uniprot(SEQ)
        if uniprot_hit:
            name = uniprot_hit.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
            acc = uniprot_hit.get("primaryAccession")
            print("UniProt match:", acc, name)
        else:
            print("No UniProt match for the sequence.")
    except Exception as exc:
        print("UniProt query failed:", exc)
    try:
        pdb_hit = search_pdb(SEQ)
        print("Closest PDB:", pdb_hit or "none")
    except Exception as exc:
        print("PDB lookup failed:", exc)

supplier = Chem.SDMolSupplier(str(LIGAND_SDF), removeHs=False)
ligand = next((mol for mol in supplier if mol), None)
if ligand:
    smiles = Chem.MolToSmiles(ligand)
    name = ligand.GetProp("_Name") if ligand.HasProp("_Name") else "(unnamed)"
    print("Ligand:", name)
    print("SMILES:", smiles)
else:
    print("Could not load ligand file.")

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
