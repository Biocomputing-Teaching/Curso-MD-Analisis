from openmm import app, unit
import openmm as mm

pdb = app.PDBFile('docs/data/alanine-dipeptide.pdb')
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds
)

integrator = mm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds
)

simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy(maxIterations=100)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

simulation.reporters.append(app.DCDReporter('docs/data/alanine-dipeptide.dcd', 10))
simulation.step(200)

print('Wrote docs/data/alanine-dipeptide.dcd')
