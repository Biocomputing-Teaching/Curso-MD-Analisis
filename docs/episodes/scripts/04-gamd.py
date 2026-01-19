import os
from pathlib import Path
from openmm import app, unit, mm
from openmm.app import PDBFile, ForceField, Simulation

COURSE_DIR = Path(os.environ.get('COURSE_DIR', str(Path.home() / 'Concepcion26'))).expanduser()
DATA_DIR = COURSE_DIR / 'data'
OUT_DIR = COURSE_DIR / 'results' / '05-muestreo-avanzado' / 'gamd'
OUT_DIR.mkdir(parents=True, exist_ok=True)

PDB_PATH = DATA_DIR / 'alanine-dipeptide.pdb'
TEMPERATURE = 300
STEPS = 1000
ITERATIONS = 12
BOOST_FACTOR = 0.3
print('Running GaMD for', PDB_PATH)
print('Outputs in', OUT_DIR)

# %%

pdb = PDBFile(str(PDB_PATH))
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
)

boost_force = mm.CustomExternalForce('boost_amplitude')
boost_force.addGlobalParameter('boost_amplitude', 0.0)
system.addForce(boost_force)

integrator = mm.LangevinIntegrator(TEMPERATURE * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
print('Setup complete')

# %%

def compute_boost(potential_energy, history, boost_factor):
    energies = history + [potential_energy]
    v_min = min(energies)
    v_max = max(energies)
    width = v_max - v_min
    if width.magnitude == 0:
        return 0 * unit.kilojoule_per_mole
    boost = boost_factor * max(0, ((potential_energy - v_min) / width)) * unit.kilojoule_per_mole
    return boost

energy_history = []

# %%

for iteration in range(ITERATIONS):
    simulation.step(STEPS)
    state = simulation.context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy()
    boost = compute_boost(potential_energy, energy_history, BOOST_FACTOR)
    energy_history.append(potential_energy)
    simulation.context.setParameter('boost_amplitude', boost)
    print(
        f'Iteration {iteration + 1}/{ITERATIONS}: V={potential_energy:.2f}, boost={boost:.2f}'
    )

log_path = OUT_DIR / 'gamd_log.csv'
with open(log_path, 'w') as handle:
    handle.write('iteration,potential_energy,kilojoule_per_mole,boost
')
    for idx, energy in enumerate(energy_history):
        boost_level = compute_boost(energy, energy_history[: idx], BOOST_FACTOR)
        handle.write(f"{idx},{energy.value_in_unit(unit.kilojoule_per_mole):.4f},{boost_level.value_in_unit(unit.kilojoule_per_mole):.4f}
")

print('Log written to', log_path)
