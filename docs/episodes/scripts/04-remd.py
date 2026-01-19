import os
from pathlib import Path
import math

from openmm import app, unit, mm
from openmm.app import PDBFile, ForceField, Simulation
import numpy as np

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data"
OUT_DIR = COURSE_DIR / "results" / "05-muestreo-avanzado" / "remd"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PDB_PATH = DATA_DIR / "alanine-dipeptide.pdb"
TEMPERATURES_K = [300, 320, 340, 360]
STEPS_PER_ITERATION = 500
N_ITERATIONS = 24

print("REMD will run for temperatures:", TEMPERATURES_K)
print("Outputs will land in", OUT_DIR)

# %%

pdb = PDBFile(str(PDB_PATH))
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds,
)

kB = unit.BOLTZMANN_CONSTANT_kB


def build_simulations(system, pdb, temperatures):
    simulations = []
    for temp in temperatures:
        integrator = mm.LangevinIntegrator(temp * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulations.append(simulation)
    return simulations


def beta(temp_kelvin):
    return 1.0 / (kB * temp_kelvin * unit.kelvin)


def get_potential_energy(simulation):
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    return state.getPotentialEnergy(), state.getPositions(asNumpy=True)

# %%

rng = np.random.default_rng(42)
simulations = build_simulations(system, pdb, TEMPERATURES_K)
swap_history = []

for iteration in range(N_ITERATIONS):
    for sim in simulations:
        sim.step(STEPS_PER_ITERATION)

    for idx in range(len(simulations) - 1):
        energy_i, positions_i = get_potential_energy(simulations[idx])
        energy_j, positions_j = get_potential_energy(simulations[idx + 1])
        beta_i = beta(TEMPERATURES_K[idx])
        beta_j = beta(TEMPERATURES_K[idx + 1])
        delta = (beta_i - beta_j) * (energy_j - energy_i)
        acceptance = min(1.0, math.exp(delta.value_in_unit(unit.dimensionless)))
        accepted = rng.random() < acceptance
        swap_history.append((iteration, idx, float(delta), accepted))
        if accepted:
            simulations[idx].context.setPositions(positions_j)
            simulations[idx + 1].context.setPositions(positions_i)

print("REMD iterations complete. Swap attempts:", len(swap_history))
print("Accepted swaps:", sum(1 for entry in swap_history if entry[3]))

log_path = OUT_DIR / "remd_swaps.csv"
with open(log_path, "w") as handle:
    handle.write("iteration,pair,delta,accepted
")
    for iteration, pair_idx, delta_val, accepted in swap_history:
        handle.write(f"{iteration},{pair_idx},{delta_val:.4f},{int(accepted)}
")

print("Swap details logged to", log_path)
