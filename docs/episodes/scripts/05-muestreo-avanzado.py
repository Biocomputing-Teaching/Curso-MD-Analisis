import openmm as mm
from openmm import app, unit

pdb = app.PDBFile('../../data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)

force = mm.CustomBondForce('0.5*k*(r-r0)^2')
force.addPerBondParameter('k')
force.addPerBondParameter('r0')
force.addBond(1, 4, [1000.0 * unit.kilojoule_per_mole / unit.nanometer**2, 0.35 * unit.nanometer])
system.addForce(force)

print('Custom forces:', system.getNumForces())
