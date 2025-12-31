from openmm import app
import openmm as mm

pdb = app.PDBFile('../../data/alanine-dipeptide.pdb')
print('OpenMM', mm.__version__)
print('Atoms:', pdb.topology.getNumAtoms())
