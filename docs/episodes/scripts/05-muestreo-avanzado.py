#!/usr/bin/env python3
import argparse
import os
import sys
import time
from pathlib import Path

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import app, unit, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller, StateDataReporter, DCDReporter

COURSE_DIR = Path(os.environ.get("COURSE_DIR", str(Path.home() / "Concepcion26"))).expanduser()
DATA_DIR = COURSE_DIR / "data" / "complex"
DEFAULT_PROTEIN = DATA_DIR / "protein.pdb"
DEFAULT_LIGAND = DATA_DIR / "ligand1.mol"


def get_platform():
    speed = 0
    platform = None
    for i in range(openmm.Platform.getNumPlatforms()):
        candidate = openmm.Platform.getPlatform(i)
        if candidate.getSpeed() > speed:
            platform = candidate
            speed = candidate.getSpeed()
    print("Using platform", platform.getName())
    if platform.getName() in {"CUDA", "OpenCL"}:
        platform.setPropertyDefaultValue("Precision", "mixed")
        print("Set precision for platform", platform.getName(), "to mixed")
    return platform


def main() -> None:
    parser = argparse.ArgumentParser(description="Simulate protein-ligand complex with optional solvation")
    parser.add_argument("-p", "--protein", default=str(DEFAULT_PROTEIN), help="Protein PDB file")
    parser.add_argument("-l", "--ligand", default=str(DEFAULT_LIGAND), help="Ligand MOL file")
    parser.add_argument("-o", "--output", default="solvated", help="Base name for output files")
    parser.add_argument("-s", "--steps", type=int, default=5000, help="Number of steps")
    parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps)")
    parser.add_argument("-f", "--friction-coeff", type=float, default=1.0, help="Friction coefficient (ps)")
    parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
    parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
    parser.add_argument("--solvate", action="store_true", help="Add solvent box")
    parser.add_argument("--padding", type=float, default=10.0, help="Padding for solvent box (A)")
    parser.add_argument("--water-model", default="tip3p", choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"], help="Water model")
    parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
    parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
    parser.add_argument("--ionic-strength", type=float, default=0.0, help="Ionic strength (M)")
    parser.add_argument("--no-neutralize", action="store_true", help="Don't neutralize")
    parser.add_argument("-e", "--equilibration-steps", type=int, default=200, help="Equilibration steps")
    args = parser.parse_args()

    t0 = time.time()
    out_dir = COURSE_DIR / "results" / "05-muestreo-avanzado" / "complex"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_base = str(out_dir / args.output)
    output_complex = output_base + "_complex.pdb"
    output_min = output_base + "_minimised.pdb"
    output_traj = output_base + "_traj.dcd"

    print("Reading ligand")
    ligand_mol = Molecule.from_file(args.ligand)

    print("Preparing system")
    forcefield_kwargs = {
        "constraints": app.HBonds,
        "rigidWater": True,
        "removeCMMotion": False,
        "hydrogenMass": 4 * unit.amu,
    }
    system_generator = SystemGenerator(
        forcefields=["amber/ff14SB.xml", "amber/tip3p_standard.xml"],
        small_molecule_forcefield="gaff-2.11",
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs,
    )

    print("Reading protein")
    protein_pdb = PDBFile(args.protein)

    print("Preparing complex")
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
    lig_top = ligand_mol.to_topology()
    modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())

    if args.solvate:
        print("Adding solvent")
        modeller.addSolvent(
            system_generator.forcefield,
            model=args.water_model,
            padding=args.padding * unit.angstroms,
            positiveIon=args.positive_ion,
            negativeIon=args.negative_ion,
            ionicStrength=args.ionic_strength * unit.molar,
            neutralize=not args.no_neutralize,
        )

    with open(output_complex, "w") as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
    step_size = args.step_size * unit.picoseconds
    friction = args.friction_coeff / unit.picosecond
    temperature = args.temperature * unit.kelvin
    duration = (step_size * args.steps).value_in_unit(unit.nanoseconds)

    if system.usesPeriodicBoundaryConditions():
        system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))

    integrator = LangevinIntegrator(temperature, friction, step_size)
    platform = get_platform()

    simulation = Simulation(modeller.topology, system, integrator, platform=platform)
    simulation.context.setPositions(modeller.positions)

    print("Minimising ...")
    simulation.minimizeEnergy()

    with open(output_min, "w") as outfile:
        PDBFile.writeFile(
            modeller.topology,
            simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
            file=outfile,
            keepIds=True,
        )

    simulation.context.setVelocitiesToTemperature(temperature)
    print("Equilibrating ...")
    simulation.step(args.equilibration_steps)

    simulation.reporters.append(DCDReporter(output_traj, args.interval, enforcePeriodicBox=True))
    simulation.reporters.append(StateDataReporter(sys.stdout, args.interval * 5, step=True, potentialEnergy=True, temperature=True))

    print("Starting simulation with", args.steps, "steps ...")
    t1 = time.time()
    simulation.step(args.steps)
    t2 = time.time()
    print("Simulation complete in", round((t2 - t1) / 60, 3), "mins")
    print("Simulation time was", round(duration, 3), "ns")
    print("Total wall clock time was", round((t2 - t0) / 60, 3), "mins")


if __name__ == "__main__":
    main()
